#
# Copyright (C) 2020 Sylvain Marsat.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
# This code is based on the example code in lisabeta:
# lisabeta/lisabeta/inference/ptemcee_smbh.py
# Make sure you check regularly for updates etc as this copy is outside of the
# git for lisabeta.


"""
    Python functions for ptemcee simulated inference of MBHBs with LISA.
"""


import os
import h5py
import itertools
import copy
import numpy as np
import json

from tqdm import tqdm as tqdm

from astropy.cosmology import Planck15 as cosmo

import lisabeta
import lisabeta.pyconstants as pyconstants
import lisabeta.tools.pytools as pytools
import lisabeta.tools.pyspline as pyspline
import lisabeta.tools.pyoverlap as pyoverlap
import lisabeta.lisa.pyresponse as pyresponse
import lisabeta.lisa.snrtools as snrtools
import lisabeta.lisa.pyLISAnoise as pyLISAnoise
import lisabeta.lisa.lisatools as lisatools
import lisabeta.lisa.lisa as lisa
import lisabeta.lisa.lisa_fisher as lisa_fisher
import lisabeta.utils.plotutils as plotutils
import lisabeta.inference.inference as inference

import ptemcee
import ptemcee.mpi_pool as mpi_pool
import time

try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None



# Physical signal parameters
params = {
    # Total *redshifted* mass M=m1+m2, solar masses
    "M": 4e6,
    # Mass ratio q=m1/m2
    "q": 3.0,
    # Dimensionless spin component 1 along orbital momentum
    "chi1": 0.5,
    # Dimensionless spin component 2 along orbital momentum
    "chi2": 0.2,
    # Time shift of coalescence, s -- coalescence is at t0*yr + Deltat*s, t0 in waveform_params
    "Deltat": 0.0,
    # Luminosity distance, Mpc
    "dist": 3.65943e+04,
    # Inclination, observer's colatitude in source-frame
    "inc": 1.0471975511965976,
    # Phase, observer's longitude in source-frame
    "phi": 1.2,
    # Longitude in the sky
    "lambda": 0.8,
    # Latitude in the sky
    "beta": 0.3,
    # Polarization angle
    "psi": 1.7,
    # Flag indicating whether angles and Deltat pertain to the L-frame or SSB-frame
    "Lframe": True
  }

# Parameters for the waveform generation and other options
waveform_params = {
    # Frequency range
    "minf": 1e-5,
    "maxf": 0.5,
    # Reference epoch of coalescence, yr -- coalescence is at t0*yr + Deltat*s, Deltat in params
    "t0": 0.0,
    # Always cut signals timetomerger_max*yr before merger -- to avoid needlessly long signals using minf
    "timetomerger_max": 1.0,
    # Option to cut the signal pre-merger -- must be in L-frame
    "DeltatL_cut": None,
    # Further options to cut signals
    "fend": None,
    "tmin": None,
    "tmax": None,
    # Options for the time and phase alignment -- development/testing
    "phiref": 0.0,
    "fref_for_phiref": 0.0,
    "tref": 0.0,
    "fref_for_tref": 0.0,
    "force_phiref_fref": True,
    "toffset": 0.0,
    # TDI channels to generate
    "TDI": "TDIAET",
    # Internal accuracy params
    "acc": 1e-4,
    "order_fresnel_stencil": 0,
    # Waveform approximant and set of harmonics to use
    "approximant": "IMRPhenomHM",
    "modes": None,
    # LISA response options
    "LISAconst": "Proposal",
    "responseapprox": "full",
    "frozenLISA": False,
    "TDIrescaled": True,
    # Noise options -- can also be given as a numpy array for interpolation
    "LISAnoise": {
        "InstrumentalNoise": "SciRDv1",
        "WDbackground": False,
        "WDduration" : 0.0,
        "lowf_add_pm_noise_f0": 0.0,
        "lowf_add_pm_noise_alpha": 2.0
    }
  }

run_params_musthave = ['out_dir', 'out_name']
list_params = [
    "M",
    "q",
    "chi1",
    "chi2",
    "Deltat",
    "dist",
    "inc",
    "phi",
    "lambda",
    "beta",
    "psi"]


# Default pymultinest params
run_params_default = {
    "sampler": "ptemcee",
    "sample_Lframe": True,
    "multimodal": False,
    "multimodal_pattern": "8modes",
    "p_jump": 0.5,
    "zerolike": False,
    "n_walkers": 16,
    "n_iter": 1000,
    "burn_in": 300,
    "autocor_method": "autocor_new", # "acor", #
    "thin_samples": False,
    "upsample": 1,
    "seed": None,
    "print_info": False,
    "n_iter_info": 10,
    "output": True,
    "output_raw": True,
    "simple_likelihood": False,
    "skip_fisher": False,
    "init_method": "fisher", # "prior", #
    "init_scale_cov": 1.,
    "n_temps": 5, # up to 10
    "temp_max": None,
    "adaptive_temp": False,
    "likelihood_method": "residuals",
    "likelihood_residuals_ngrid": 128,
}

# Default parameter range for angular parameters
params_range_default = {
    "chi1": [-1., 1.],
    "chi2": [-1., 1.],
    "inc": [0.000001, np.pi],
    "phi": [-np.pi, np.pi],
    "lambda": [-np.pi, np.pi],
    "beta": [-np.pi/2., np.pi/2.],
    "psi": [0., np.pi]
}


################################################################################
# Main

def RunPTEMCEE(input_file):

    """
    Helper function to run ptemcee within SYNEX.

    parameters
    ----------
    input_file : str
        json input file with all parameters saved.
    """

    # Using MPI or not
    # We allow for different cases: MPI used or not, master or worker
    use_mpi=False
    is_master = True
    mapper = map
    if MPI is not None:
        MPI_size = MPI.COMM_WORLD.Get_size()
        MPI_rank = MPI.COMM_WORLD.Get_rank()
        use_mpi = (MPI_size > 1)
        if use_mpi:
            print("MPI rank/size: %d / %d" % (MPI_rank, MPI_size), flush=True)
            pool = ptemcee.mpi_pool.MPIPool(debug=False)
            is_master = pool.is_master()
            mapper = pool.map
        else:
            print("No MPI", flush=True)
            is_master = True
            mapper = map

    # Define globals if using MPI to avoid cluster problems
    if use_mpi: global lnprior, prior, lnlikelihood

    # Load json file with all parameters
    with open(input_file, 'r') as input_file:
        input_params = json.load(input_file)
    # Source params and prior params must be given
    source_params = input_params['source_params']
    source_params = pytools.complete_mass_params(source_params)
    prior_params = input_params['prior_params']
    # Source params can be specified in the Lframe
    source_params_are_Lframe = source_params["Lframe"] # source_params.pop('Lframe', False)
    # Waveform params and run params have predefined default values
    # Need to cast from json to python, list of modes list[list] -> list[tuple]
    # waveform_params = waveform_params_default.copy()
    input_params['waveform_params'] = inference.waveform_params_json2py(
                                                input_params['waveform_params'])
    # waveform_params.update(input_params['waveform_params'])
    waveform_params = input_params['waveform_params'].copy()
    run_params = run_params_default.copy()
    run_params.update(input_params['run_params'])
    zerolike = run_params.pop('zerolike', False)
    sample_Lframe = run_params.pop('sample_Lframe', False)
    simple_likelihood = run_params.pop('simple_likelihood', False)

    # Check run params are meant to be for the ptemcee sampler
    if not (run_params['sampler'] == 'ptemcee'):
        raise ValueError('run_params sampler flag is %s instead of ptemcee.'\
                          % (run_params['sampler']))
    # run params must contain output prefix for files -- no default here
    for musthave in run_params_musthave:
        if not (musthave in run_params):
            raise ValueError('run_params must contain %s.' % (musthave))
    # Check that the output directory exists
    if (run_params['output'] or run_params['output_raw']):
        if not os.path.isdir(run_params['out_dir']):
            raise ValueError('Output dir %s does not exist.' % run_params['out_dir'])

    # If source params are given in the Lframe, convert to SSB
    if source_params_are_Lframe:
        source_params_SSBframe = lisatools.convert_Lframe_to_SSBframe(
                                    source_params,
                                    t0=waveform_params['t0'],
                                    frozenLISA=waveform_params['frozenLISA'])
        source_params_Lframe = source_params.copy()
    else:
        source_params_SSBframe = source_params.copy()
        source_params_Lframe = lisatools.convert_SSBframe_to_Lframe(
                                    source_params,
                                    t0=waveform_params['t0'],
                                    frozenLISA=waveform_params['frozenLISA'])

    # Use default param ranges to complete input ranges -- check completeness
    prior_type = prior_params['prior_type']
    params_range = prior_params['params_range']
    infer_params = prior_params['infer_params']
    for i, param in enumerate(infer_params):
        if params_range[i] == []:
            params_range[i] = params_range_default[param]
    if [] in params_range:
        raise ValueError('Some parameter ranges are neither given \
                      nor covered by the defaults.')

    # Parameters inferred and parameters fixed
    infer_params = prior_params['infer_params']
    fixed_params = [p for p in list_params if p not in infer_params]
    n_dim = len(infer_params)
    n_fixed = len(fixed_params)
    # Vector x of infer_params for the injection
    # Also store values of fixed params in dict for later use
    if sample_Lframe:
        fixed_params_dict = dict([(p, source_params_Lframe[p]) for p in fixed_params])
        x_inj = np.array([source_params_Lframe[p] for p in infer_params])
    else:
        fixed_params_dict = dict([(p, source_params_SSBframe[p]) for p in fixed_params])
        x_inj = np.array([source_params_SSBframe[p] for p in infer_params])

    # TODO: implement wrapped parameters
    # wrapped_params = [0] * n_dim
    # for i in range(n_dim):
    #     wrapped_params[i] = int(prior_params['wrap_params'][i])
    # run_params['wrapped_params'] = wrapped_params

    # Multimodal jumps
    # TODO: support both 2-mode and 8-mode jumps
    # TODO: be more flexible when some parameters are pinned to their value
    # For now we require all ['inc', 'lambda', 'beta', 'psi'] to be
    # inferred params
    # TODO: improve the back-and-forth between vector and dict
    multimodal = run_params['multimodal']
    multimodal_pattern = run_params['multimodal_pattern']
    extra_proposal_prob = run_params['p_jump']
    extra_proposal_jump = None
    if multimodal:
        if not sample_Lframe:
            raise ValueError('The multimodal option is only supported \
                                together with sample_Lframe.')

        list_angle_params = ['inc', 'lambda', 'beta', 'psi']
        index_map_angles = {}
        for param in list_angle_params:
            if not param in infer_params:
                raise ValueError('For multimodality, infer_params must contain \
                                        [inc, phi, lambda, beta, psi]')
            index_map_angles[param] = infer_params.index(param)
        # x is (n, ndim) array with physical parameters (for n walkers)
        # returns x_jump
        def jump_extra(x, random_state=np.random):
            angle_params = {}
            for param in list_angle_params:
                angle_params[param] = x[:,index_map_angles[param]]
            skymodes = inference.proposal_skymode(size=len(x), pattern=multimodal_pattern, random_state=random_state)
            angle_params_map = inference.map_skymode(angle_params, skymodes[:,0], skymodes[:,1])
            x_jump = x.copy()
            for param in list_angle_params:
                x_jump[:,index_map_angles[param]] = angle_params_map[param]
            return x_jump
        extra_proposal_jump = jump_extra

    # Prior function (assumes independent multiplicative prior in all params)
    def prior(x):
        p = 1.
        for i in range(n_dim):
            p *= inference.compute_prior(prior_type[i], params_range[i], x[i])
        return p
    # Set up prior bounds for x vector of parameters, dimension (n_dim, 2)
    prior_bounds = np.array(params_range)
    # Check the injection is within the prior bounds
    if not inference.prior_check(x_inj, prior_bounds):
        print(x_inj, prior_bounds)
        raise ValueError('Source parameters not within the prior bounds!')

    # Set up likelihood class
    if simple_likelihood:
        likelihoodClass = lisa.LikelihoodLISASMBH_LinearResiduals(source_params_Lframe, ngrid=128, **waveform_params)
    else:
        # Historic 0-noise likelihood
        likelihoodClass = lisa.LikelihoodLISASMBH(source_params_Lframe, **waveform_params)

    # Waveform params to be used for template waveforms
    # May differ from the waveform params for injection, given as an update
    waveform_params_template = copy.deepcopy(waveform_params)
    if 'waveform_params_update_for_template' in waveform_params:
        waveform_params_template.update(
                         waveform_params['waveform_params_update_for_template'])

    # Define ln-likelihood and ln-prior function, wrapped for use in ptemcee
    # TODO: add m1>m2 to prior range checks
    def lnlikelihood(x):
        if not inference.prior_check(x, prior_bounds):
            return -1e99
        template_params = {}
        for i in range(n_dim):
            template_params[infer_params[i]] = x[i]
        for p in fixed_params:
            template_params[p] = fixed_params_dict[p]
        template_params = pytools.complete_mass_params(template_params)
        template_params = pytools.complete_spin_params(template_params)
        template_params['Lframe'] = sample_Lframe
        lnL = likelihoodClass.lnL(template_params, **waveform_params_template)
        return lnL
    def lnprior(x):
        if not inference.prior_check(x, prior_bounds):
            return -1e99
        p = prior(x)
        if p==0.: return -1e99
        return np.log(p)

    # Read run_params
    # p_extra = run_params['p_jump']
    init_method = run_params['init_method']
    init_scale_cov = run_params['init_scale_cov']
    n_temps = run_params['n_temps']
    temp_max = run_params['temp_max']
    adaptive_temp = run_params['adaptive_temp']
    n_walkers = run_params['n_walkers']
    n_iter = run_params['n_iter']
    burn_in = run_params['burn_in']
    autocor_method = run_params['autocor_method']
    upsample = run_params['upsample']
    seed = run_params['seed']
    # Only master will print info
    print_info = run_params['print_info'] and is_master
    n_iter_info = run_params['n_iter_info']
    n_dim = len(infer_params)

    if print_info:
        print('Source parameters SSBframe:')
        print(source_params_SSBframe)
        print('Source parameters Lframe:')
        print(source_params_Lframe)

    # Compute Fisher matrix - not available for all waveform models
    if not run_params['skip_fisher']:
        # NOTE: might be slow if timetomerger_max is set at generic value like 1yr
        # NOTE: numerical problems in Fisher when timetomerger_max is too large
        # For Fisher only, we set timetomerger_max based on SNR=1 threshold
        waveform_params_22 = waveform_params.copy()
        waveform_params_22['modes'] = [(2,2)]
        t_SNR1_s = snrtools.lisa_tofSNR(1., source_params_SSBframe, **waveform_params_22)
        t_SNR1_yr = t_SNR1_s / pyconstants.YRSID_SI
        # Ensure at least 1d of signal for fisher regardless of SNR threshold
        t_threshold_yr = np.fmax(t_SNR1_yr, 1./365.25)
        waveform_params_fisher = waveform_params.copy()
        waveform_params_fisher['timetomerger_max'] = t_threshold_yr
        # TODO: improve handling of different choices of infer_params
        if print_info:
            print('Computing Fisher matrix...')
            print('Threshold max(t(SNR=1), 1d) used for Fisher: %g yr' % (t_threshold_yr))
            # print('Fisher params. \n Source params ', source_params_SSBframe, '\n steps: ', lisa_fisher.default_steps,
            #     '\n freqs: ', ['log', None], '\n list_params: ', infer_params, '\n list_fixed_params: ', fixed_params,
            #     '\n Lframe: ', sample_Lframe, '\n waveform_params_fisher: ', waveform_params_fisher)
        fishercov = lisa_fisher.fisher_covariance_smbh(source_params_SSBframe,
                            steps=lisa_fisher.default_steps, freqs=['log', None],
                            list_params=infer_params, list_fixed_params=fixed_params,
                            Lframe=sample_Lframe,
                            **waveform_params_fisher)
        if print_info:
            print('Fisher covariance matrix:')
            print(fishercov['cov'])

    # Initialize walkers, drawing from the prior or from Fisher matrix
    x_ini = np.zeros((n_temps, n_walkers, n_dim))
    if init_method=='prior':
        # Drawing from the prior -- assumes independent priors in all params
        for k in range(n_temps):
            for j in range(n_dim):
                x_ini[k,:,j] = inference.draw_prior(prior_type[j], prior_bounds[j], n_walkers)
    elif init_method=='fisher':
        if run_params['skip_fisher']:
            raise ValueError('Cannot have skip_fisher and init_method=fisher.')
        # Drawing from the Fisher matrix - same for each temperature
        if sample_Lframe:
            mean = np.array([source_params_Lframe[p] for p in infer_params])
        else:
            mean = np.array([source_params_SSBframe[p] for p in infer_params])
        # We allow for the possibility of scaling the covariance for initialization
        cov = fishercov['cov'] * init_scale_cov
        try:
            for k in range(n_temps):
                for i in range(n_walkers):
                    x_ini[k,i,:] = inference.draw_multivariate_normal_constrained(
                                          mean, cov, params_range, infer_params)
        except:
            print('Fisher initialization failed, initializing from prior.')
            # Drawing from the prior -- assumes independent priors in all params
            for k in range(n_temps):
                for j in range(n_dim):
                    x_ini[k,:,j] = inference.draw_prior(prior_type[j], prior_bounds[j], n_walkers)
    else:
        raise ValueError('Option init_method = %s not recognized.' % init_method)

    # Initial points all start on the main mode (1,0)
    # if multimodal:
    #     x_ini[:,n_dim] = 1
    #     x_ini[:,n_dim+1] = 0

    # Enforce range of phase parameters in case Fisher uncertainties are large
    # TODO: we should have proper rejection based on the prior ranges
    if 'inc' in infer_params:
        iinc = infer_params.index('inc')
        x_ini[:,:,iinc] = pytools.modpi(x_ini[:,:,iinc])
    if 'phi' in infer_params:
        iphi = infer_params.index('phi')
        x_ini[:,:,iphi] = pytools.mod2pi(x_ini[:,:,iphi])
    if 'lambda' in infer_params:
        ilambda = infer_params.index('lambda')
        x_ini[:,:,ilambda] = pytools.mod2pi(x_ini[:,:,ilambda])
    if 'beta' in infer_params:
        ibeta = infer_params.index('beta')
        x_ini[:,:,ibeta] = pytools.modpi2(x_ini[:,:,ibeta])
    if 'psi' in infer_params:
        ipsi = infer_params.index('psi')
        x_ini[:,:,ipsi] = pytools.mod2pi(x_ini[:,:,ipsi])

    # Temperature ladder
    betas = ptemcee.sampler.make_ladder(n_dim, ntemps=n_temps, Tmax=temp_max)

    # If using mpi, workers are simply told to wait
    # They will be sollicited by master for likelihood computations via pool.map
    if use_mpi and (not is_master):
        pool.wait()
        # Is that exit necessary ? Was in the emcee/schwimmbad example
        # Maybe there to ensure exits if pool.close() is not called
        sys.exit(0)

    # Whether using mpi and being master or simply not using mpi
    # Sampler set-up, running, post-processing and output all done by master
    if is_master:
        # Set up and initialize sampler
        sampler = ptemcee.sampler.Sampler(n_walkers, n_dim, lnlikelihood, lnprior,
                                          betas=betas, mapper=mapper,
                                          adaptive=adaptive_temp,
                                          extra_proposal_prob=extra_proposal_prob,
                                          extra_proposal_jump=extra_proposal_jump)
        chain = sampler.chain(x_ini)

        # A la grace de dieu
        if print_info:
            print('Running ptemcee...')
        print("outdir checks 4:",run_params["out_dir"], run_params["out_name"])
        print("outdir checks 5:", run_params["output"], run_params["output_raw"])
        t1 = time.time()
        chain.run(n_iter)
        t2 = time.time()
        print("Time to execute ptemcee: ", round((t2-t1)*10.)/10., "s")

        # TODO: output meta information, acceptance rate and likelihood speed

        # Post-processing: unfold if multimodal was used
        if print_info:
            print('Post-processing...')

        # Only keep 0-temperature chain for post-processing
        # TODO: allow to output other temperatures, and temp. adjustments, for debug
        # ptemcee format: (niter, ntemp, nwalker, ndim)
        # We translate to multiemcee format: (nwalker, niter, ndim) for temp=0
        chain_unfold = np.transpose(chain.x[:,0,:,:], axes=(1,0,2))
        # Same format without ndim for like and post values
        lnlikevals = np.transpose(chain.logl[:,0,:], axes=(1,0))
        lnpostvals = np.transpose(chain.logP[:,0,:], axes=(1,0))
        # Format: (niter, ntemp)
        betavals = chain.betas
        # ptemcee does not seem to store individual acceptance, only their count
        acceptance_ratios = chain.jump_acceptance_ratio[0,:]

        # Post-processing: burn-in, autocor. length and thinning, merge walkers
        # Also complete with params that were fixed to their values
        if burn_in>=n_iter:
            print('WARNING: burn_in is larger than n_iter, ignoring burn_in.')
            burn_in = 0
        if run_params['thin_samples']:
            autocor_len = {}
            for i,param in enumerate(infer_params):
                # Note: autocorr_new computes a mean across walkers
                autocor_len[param] = inference.get_autocor(chain_unfold[:,burn_in:,i], method=autocor_method)
                if print_info:
                    print('Autocorrelation length for %s: %g' % (param,autocor_len[param]))
            # The thinning length is the worse autocor length for all params
            thin_len = int(np.ceil(np.max([autocor_len[param] for param in infer_params])))
            thin_len = int(np.ceil(thin_len / float(upsample)))
        else:
            thin_len = 1
        if print_info:
            print('Autocorrelation thinning: %d' % (thin_len,))
        # Processed output: thinned and merged chain (completed with fixed params)
        n_out = n_walkers * ((n_iter-burn_in-1) // thin_len + 1)
        chain_processed = np.zeros((n_out, n_dim + n_fixed), dtype=float)
        # Here ravel with order='F' flattens by reading along columns first
        for i in range(n_dim):
            i_out = list_params.index(infer_params[i])
            chain_processed[:,i_out] = np.ravel(chain_unfold[:,burn_in::thin_len,i],
                                                order='F')
        for i in range(n_fixed):
            i_out = list_params.index(fixed_params[i])
            if sample_Lframe:
                chain_processed[:,i_out] = source_params_Lframe[fixed_params[i]]
            else:
                chain_processed[:,i_out] = source_params_SSBframe[fixed_params[i]]
        # Thin and merge likelihood, posterior, acceptance vals
        lnlikevals_processed = np.ravel(lnlikevals[:,burn_in::thin_len], order='F')
        lnpostvals_processed = np.ravel(lnpostvals[:,burn_in::thin_len], order='F')
        # acceptvals_processed = np.ravel(acceptvals[:,burn_in::thin_len], order='F')

        # Output
        basename = run_params['out_dir'] + run_params['out_name']
        # Main output: processed, thinned chains and Fisher matrix
        # TODO: hierarchy, create a level (group ?) samples
        if run_params['output']:
            with h5py.File(basename + '.h5', 'w') as f:
                # Source params
                source_params_Lframe_gr = f.create_group('source_params_Lframe')
                pytools.write_h5py_dict(source_params_Lframe_gr, source_params_Lframe)
                source_params_SSBframe_gr = f.create_group('source_params_SSBframe')
                pytools.write_h5py_dict(source_params_SSBframe_gr, source_params_SSBframe)
                for i,param in enumerate(list_params):
                    f.create_dataset(param, data=chain_processed[:,i])
                f.create_dataset('lnlike', data=lnlikevals_processed)
                f.create_dataset('lnpost', data=lnpostvals_processed)
                f.create_dataset('acceptance_ratios', data=acceptance_ratios)
                f.create_dataset('betas', data=betavals)
                # f.create_dataset('accept', data=acceptvals_processed)
                if not run_params['skip_fisher']:
                    fishercov_gr = f.create_group('fishercov')
                    lisa_fisher.write_h5py_fishercov(fishercov_gr, fishercov)
        # Extra output: full chains, walker-separated, for testing and diagnostics
        if run_params['output_raw']:
            with h5py.File(basename + '_raw.h5', 'w') as f:
                # Source params
                source_params_Lframe_gr = f.create_group('source_params_Lframe')
                pytools.write_h5py_dict(source_params_Lframe_gr, source_params_Lframe)
                source_params_SSBframe_gr = f.create_group('source_params_SSBframe')
                pytools.write_h5py_dict(source_params_SSBframe_gr, source_params_SSBframe)
                for i,param in enumerate(infer_params):
                    f.create_dataset(param, data=chain_unfold[:,:,i])
                f.create_dataset('lnlike', data=lnlikevals)
                f.create_dataset('lnpost', data=lnpostvals)
                # f.create_dataset('accept', data=acceptvals)

        # If using mpi, close the pool
        if use_mpi:
            pool.close()
