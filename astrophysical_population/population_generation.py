#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 14:02:11 2023

@author: duverne
"""

import numpy as np
import argparse, os
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.table import Table

def generate_uniform_quantity(N_sim, inf=0., sup=1., plot=False):
    """
    Generate an array of floats evenly genrated in [inf, sup].

    Parameters
    ----------
    N_sim : int
        Number of generated float.
    inf : float, optional
        Inferior limit of the interval to generate numbers. The default is 0.
    sup : int, optional
        Inferior limit of the interval to generate numbers. The default is 1.
    plot : bool, optional
        Plot the generated distribution histogram. The default is False.

    Returns
    -------
    quantity : numpy array size N_sim
        Array of the generated numbers.

    """
    quantity = (sup - inf) * np.random.random_sample(N_sim) + inf
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 12.5))
        ax.hist(quantity, bins=100)

    return quantity

def generate_normal_quantity(N_sim, mu=0., std_dev=1., plot=False):
    """
    Generate an array of floats with a normal distribution with mean=mu\
    and std=std_dev.

    Parameters
    ----------
    N_sim : int
        Number of generated float.
    mu : float, optional
        Mean of the normal distribution. The default is 0..
    sup : int, optional
        Standard deviation of the normal distribution. The default is 1..
    plot : bool, optional
        Plot the generated distribution histogram. The default is False.

    Returns
    -------
    quantity : numpy array size N_sim
        Array of the generated numbers.

    """
    quantity = np.random.normal(mu, std_dev, N_sim)
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 12.5))
        ax.hist(quantity, bins=100)

    return quantity

def generate_uniform_position(min_radius, max_radius):
    """
    Generate a position in (d, rad, dec) evenly choosen in volume.

    Parameters
    ----------
    min_radius : Float
        Radius of the sphere in which the position is excluded.
    max_radius : Float
        Radius of the sphere in which the position is generated.

    Returns
    -------
    d : Float
        Distance.
    src_lambda : Float
        Right ascension.
    src_beta : Float
        Declination.

    """
    x = generate_uniform_quantity(1, -max_radius, max_radius)
    y = generate_uniform_quantity(1, -max_radius, max_radius)
    z = generate_uniform_quantity(1, -max_radius, max_radius)

    d = np.sqrt(x**2 + y**2 + z**2)
    while  d > max_radius or d < min_radius:
        x = generate_uniform_quantity(1, -max_radius, max_radius)
        y = generate_uniform_quantity(1, -max_radius, max_radius)
        z = generate_uniform_quantity(1, -max_radius, max_radius)
        d = np.sqrt(x**2 + y**2 + z**2)
    src_lambda =  np.arctan2(x, y)
    src_beta = np.arcsin(z / d)

    return d, src_lambda, src_beta

def generate_mass(min_mass, max_mass, mass_ratio_cut = 20.):

    mass_1 = generate_uniform_quantity(1, min_mass, max_mass)
    mass_2 = generate_uniform_quantity(1, min_mass, max_mass)
    q = max(mass_1, mass_2) / min(mass_1, mass_2)

    while q > mass_ratio_cut:
        mass_1 = generate_uniform_quantity(1, min_mass, max_mass)
        mass_2 = generate_uniform_quantity(1, min_mass, max_mass)
        q = np.max(mass_1, mass_2) / np.min(mass_1, mass_2)

    return mass_1, mass_2

def generate_cbc(N_cbc, save=False, path2save='./parameters_cbc.dat'):

    # Polarisation uniform in [0, pi] (phi param in lisabeta)
    polarization = generate_uniform_quantity(N_cbc, 0., np.pi)
    # Coa phase uniform in[0, 2 pi] (psi param in lisabeta)
    coa_phase = generate_uniform_quantity(N_cbc, -np.pi, np.pi)
    # cos(inclination) uniform in [-1, 1]
    inclination = np.arccos(generate_uniform_quantity(N_cbc, -1., 1.))

    # Instanciate the arrays for the distance and angles parameters
    distance = np.zeros(N_cbc)
    src_lambda = np.zeros(N_cbc)
    src_beta = np.zeros(N_cbc)

    # Masses are chosen for given mass ratio and total mass.
    # No need to generate them for our purpose then.
    # mass1, mass2 = np.zeros(N_cbc), np.zeros(N_cbc)

    # For spinsz : Uniform distribution between -1, 1
    spin1z = generate_uniform_quantity(N_cbc, -1., 1.)
    spin2z = generate_uniform_quantity(N_cbc, -1., 1.)

    # Generate the starting time t0 of the merger
    t0 = generate_uniform_quantity(N_cbc, 0., 1.)

    # Generate the positions of the binaries
    for i in tqdm(range(N_cbc), desc="Generating the CBC's positions"):
        distance[i], src_lambda[i], src_beta[i] = generate_uniform_position(min_radius=1000.,
                                                                            max_radius=20000.)

    parameters_cbc = Table([inclination,
                            coa_phase,
                            polarization,
                            src_lambda,
                            src_beta,
                            spin1z,
                            spin2z,
                            t0,
                            distance],
                           names=['inclination',
                                  'coa_phase',
                                  'polarization',
                                  'lambda',
                                  'beta',
                                  'spin1z',
                                  'spin2z',
                                  't0',
                                  'distance'])
    if save:
        parameters_cbc.write(path2save,
                             format='ascii.commented_header',
                             overwrite=True)
    return parameters_cbc


def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate a SMBBH population.')
    parser.add_argument('--N_BBH',
                        type=int,
                        default=10,
                        help='Number of CBC in the population.')

    parser.add_argument('--path',
                        type=str,
                        default='./',
                        help='Path to save the injection file and the plots.')

    args = parser.parse_args()

    filename = os.path.join(args.path, 'parameters_cbc.dat')
    parameters_cbc = generate_cbc(args.N_BBH,
                                  save=True,
                                  path2save=filename)

    # Create the directory for saving the plots if it does not exists
    path = os.path.join(args.path, 'plots')
    if not os.path.exists(path):
        os.mkdir(path)

    # Plotting interresting quantities
    label_size = 25
    fontsize = 30
    ########### PLOT RaDec ###############
    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9,
                              wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)


    main_axes.hist(parameters_cbc['lambda'], alpha = 0.5,
                   label=r'$\lambda$', bins=100)

    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'$\lambda$ [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Lambda', fontsize=fontsize)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)
    fig.savefig(path+'/BBH_lambda_distribution.png')
    plt.close(fig)


    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9,
                              wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)


    main_axes.hist(parameters_cbc['beta'], alpha = 0.5,
                   label=r'$\beta$', bins=100)
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'$\beta$ [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Beta', fontsize=fontsize)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)
    fig.savefig(path+'/BBH_beta_distribution.png')
    plt.close(fig)

    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9,
                              wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.scatter(parameters_cbc['lambda'], parameters_cbc['beta'],
                      color='blue', marker='o', s=60)
    main_axes.set_xlabel(r'$\lambda$ [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'$\beta$ [rad]]', fontsize=fontsize)
    fig.savefig(path+'/BBH_lambda_vs_beta.png')
    plt.close(fig)

    ########### PLOT distance ###############

    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9, wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.hist(parameters_cbc['distance'], alpha = 0.5,
                   label='Distance', bins=100)
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'Distance [Mpc]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Distance', fontsize=fontsize)
    fig.savefig(path+'/BBH_distance_distribution.png')
    plt.close(fig)


    # ########### PLOT inclination ###############
    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9, wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.hist(parameters_cbc['inclination'], alpha = 0.5,
                   label='Inclination' ,bins=100)
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'$\iota$ [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Inclinations', fontsize=fontsize)
    fig.savefig(path+'/BBH_inclination_distribution.png')
    plt.close(fig)

    # ########### PLOT coa_phase ###############
    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9,
                              wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.hist(parameters_cbc['coa_phase'], alpha = 0.5,
                   label='phase at Coalescence', bins=100)
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'Coalescence Phase [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Phase at Coalescence', fontsize=fontsize)
    fig.savefig(path+'/BBH_coa_phase_distribution.png')
    plt.close(fig)


    # ########### PLOT polarization ###############
    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9, wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.hist(parameters_cbc['polarization'], alpha = 0.5,
                   label='polarization', bins=100)
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'Polarization [rad]', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Polarization', fontsize=fontsize)
    fig.savefig(path+'/BBH_polarization_distribution.png')
    plt.close(fig)

    # ########### PLOT masses ###############
    # fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    # plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
    #                           right = 0.9, top = 0.9,
    #                           wspace = 2.5, hspace = 0.5)

    # main_axes.grid(visible=True, which='major',
    #                color='#666666', linestyle='-')
    # main_axes.minorticks_on()
    # main_axes.grid(visible=True, which='minor',
    #                color='#999999', linestyle='-', alpha=0.5)
    # main_axes.xaxis.set_tick_params(labelsize=label_size)
    # main_axes.yaxis.set_tick_params(labelsize=label_size)

    # main_axes.hist(parameters_cbc['mass1'], alpha = 0.5,
    #                 label=r'Mass 1', bins=100, color='blue')
    # main_axes.hist(parameters_cbc['mass2'], alpha = 0.5,
    #                 label=r'Mass 2', bins=100, color='red')
    # main_axes.legend(fontsize=fontsize)
    # main_axes.set_xlabel(r'Mass [M$_{\odot}$]', fontsize=fontsize)
    # main_axes.set_title(r'Distribution of Masses', fontsize=fontsize)
    # main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    # fig.savefig(path+'/BBH_masses_distribution.png')
    # plt.close(fig)

    # fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    # plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
    #                           right = 0.9, top = 0.9,
    #                           wspace = 2.5, hspace = 0.5)

    # main_axes.grid(visible=True, which='major',
    #                color='#666666', linestyle='-')
    # main_axes.minorticks_on()
    # main_axes.grid(visible=True, which='minor',
    #                color='#999999', linestyle='-', alpha=0.5)
    # main_axes.xaxis.set_tick_params(labelsize=label_size)
    # main_axes.yaxis.set_tick_params(labelsize=label_size)

    # main_axes.scatter(parameters_cbc['mass1'], parameters_cbc['mass2'],
    #                   color='blue', marker='o', s=60)
    # main_axes.set_xlabel(r'Mass 1 [M$_{\odot}$]', fontsize=fontsize)
    # main_axes.set_ylabel(r'Mass 2 [M$_{\odot}$]', fontsize=fontsize)
    # fig.savefig(path+'/BBH_M1_vs_M2.png')
    # plt.close(fig)

    # ########### PLOT spinz ###############
    fig, main_axes = plt.subplots(1,1, figsize=(15,15))
    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                              right = 0.9, top = 0.9,
                              wspace = 2.5, hspace = 0.5)

    main_axes.grid(visible=True, which='major',
                   color='#666666', linestyle='-')
    main_axes.minorticks_on()
    main_axes.grid(visible=True, which='minor',
                   color='#999999', linestyle='-', alpha=0.5)
    main_axes.xaxis.set_tick_params(labelsize=label_size)
    main_axes.yaxis.set_tick_params(labelsize=label_size)

    main_axes.hist(parameters_cbc['spin1z'], alpha = 0.5,
                    label=r'Spinz 1',bins=100, color='blue')
    main_axes.hist(parameters_cbc['spin2z'], alpha = 0.5,
                    label=r'Spinz 2',bins=100, color='red')
    main_axes.legend(fontsize=fontsize)
    main_axes.set_xlabel(r'Spinz', fontsize=fontsize)
    main_axes.set_title(r'Distribution of Spinz', fontsize=fontsize)
    main_axes.set_ylabel(r'Number of BBH', fontsize=fontsize)
    fig.savefig(path+'/BBH_spinz_distribution.png')
    plt.close(fig)


if __name__ == "__main__":
    main()