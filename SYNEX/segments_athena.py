"""
Function to replace gwemopt segments when considering a non-geostationary telescope
(e.g. Athena).
Functions are kept as close to original gwemopt.segments.py script so that it can
be proposed for integration into gwemopt later as either new script or added
to original.

J. Baird 31/01/2022
"""
import os, sys
import numpy as np
import healpy as hp
import copy

import astropy.coordinates
import astropy.constants as consts
from astropy.time import Time
import astropy.units as u
from joblib import Parallel, delayed
import ephem
import pickle

import ligo.segments as segments
import gwemopt.utils

def get_telescope_segments(params):
    """
    No other changes needed here. gwemopt.utils.get_exposures returns the times
    per tile after reducing total time available by overheads for slew, filter changes
    etc.
    """

    for telescope in params["telescopes"]:

        params["config"][telescope]["segmentlist"] = get_segments(params, params["config"][telescope])

        params["config"][telescope]["exposurelist"] = gwemopt.utils.get_exposures(params, params["config"][telescope], params["config"][telescope]["segmentlist"])
        if len(params["config"][telescope]["exposurelist"]) == 0:
            params["config"][telescope]["n_windows"] = 0
            params["config"][telescope]["tot_obs_time"] = 0.0
            continue

        nexp,_ = np.array(params["config"][telescope]["exposurelist"]).shape
        params["config"][telescope]["n_windows"] = nexp
        tot_obs_time = np.sum(np.diff(np.array(params["config"][telescope]["exposurelist"]))) * 86400.
        params["config"][telescope]["tot_obs_time"] = tot_obs_time

    return params

def get_segments(params, config_struct):
    """
    Adding to params something for orbit type (e.g. L1 Halo or file name if we have coords?)

    No other changes needed here. Maybe modifications are not needed to second half
    since this function returns just "segmentlist" if we specify sat_sun_restriction
    to be non-zero. "lat", "long", and "elevs" are used only to check when it's night...
    """

    gpstime = params["gpstime"]
    event_mjd = Time(gpstime, format='gps', scale='utc').mjd

    segmentlist = segments.segmentlist()
    n_windows = len(params["Tobs"]) // 2
    start_segments = event_mjd + params["Tobs"][::2]
    end_segments = event_mjd + params["Tobs"][1::2]
    for start_segment, end_segment in zip(start_segments,end_segments):
        segmentlist.append(segments.segment(start_segment,end_segment))

    # Not sure this next part is really necessary after cutting night time segments calculations
    segmentlist.coalesce()

    return segmentlist

def get_segments_tile(config_struct, radec, segmentlist):
    """
    TO DO::: Need to check that all the ra-dec definitions across lisabeta, gwemopt, and synex
    are coherent... I think there are a fair few slip ups moving between Earth frame, Lisa frame,
    Athena frame, and SSB frame... Maybe we should write some routines to translate between them,
    including calls to lisatools for simplification.

    This function is heavily modified because the functions are just not
    useful as are in "segments.py" for a telescope like Athena.
    """
    # Set defaults for constraints just in case they aren't handled elsewhere.
    if "sat_moon_restriction" not in config_struct:
        config_struct["sat_moon_restriction"]=20.
    elif config_struct["sat_moon_restriction"]==None:
        config_struct["sat_moon_restriction"]=20.

    if "sat_earth_restriction" not in config_struct:
        config_struct["sat_earth_restriction"]=30.
    elif config_struct["sat_earth_restriction"]==None:
        config_struct["sat_earth_restriction"]=30.

    if "sat_sun_restriction" not in config_struct:
        config_struct["sat_sun_restriction"]=45.
    elif config_struct["sat_sun_restriction"]==None:
        config_struct["sat_sun_restriction"]=45.

    # Get constraints for tile proximity to Sun, Earth, and Moon.
    moon_constraint = float(config_struct["sat_moon_restriction"])
    earth_constraint = float(config_struct["sat_earth_restriction"])
    sun_constraint = float(config_struct["sat_sun_restriction"])

    # Unpack orbit stuff times_as_dates
    elapsed_time_from_start=config_struct["orbit_dict"]["elapsed_time_from_start"]
    times_mjd=config_struct["orbit_dict"]["times_mjd"]
    M_f_A_radecs=config_struct["orbit_dict"]["Moon_From_Athena_radecs"]
    E_f_A_radecs=config_struct["orbit_dict"]["Earth_From_Athena_radecs"]
    S_f_A_radecs=config_struct["orbit_dict"]["Sun_From_Athena_radecs"]

    # Angular radii in degrees
    moon_radii=np.arctan(1737400./np.linalg.norm(config_struct["orbit_dict"]["Moon_From_Athena"],axis=0))*180./np.pi
    earth_radii=np.arctan(consts.R_earth.value/np.linalg.norm(config_struct["orbit_dict"]["Earth_From_Athena"],axis=0))*180./np.pi
    sun_radii=np.arctan(consts.R_sun.value/np.linalg.norm(config_struct["orbit_dict"]["Sun_From_Athena"],axis=0))*180./np.pi

    # Angular distances from tile direction and three astrophysical objects
    # Can we pass the 'angular_distance' function an array of directions? I think so if its ndarray...
    Moon_AngDist=angular_distance(M_f_A_radecs[0,:], M_f_A_radecs[1,:], radec.ra.value, radec.dec.value)
    Earth_AngDist=angular_distance(E_f_A_radecs[0,:], E_f_A_radecs[1,:], radec.ra.value, radec.dec.value)
    Sun_AngDist=angular_distance(S_f_A_radecs[0,:], S_f_A_radecs[1,:], radec.ra.value, radec.dec.value)

    # Make the segmentlists for times when we are outside restriction limits
    seglst=[segments.segment(times_mjd[i],times_mjd[i+1]) for i in range(len(times_mjd)-1) if Moon_AngDist[i+1]>=moon_constraint and Moon_AngDist[i+1]>=moon_radii[i+1]]
    moonsegmentlist = segments.segmentlist(seglst)
    moonsegmentlist.coalesce()
    seglst=[segments.segment(times_mjd[i],times_mjd[i+1]) for i in range(len(times_mjd)-1) if Earth_AngDist[i+1]>=earth_constraint and Earth_AngDist[i+1]>=earth_radii[i+1]]
    earthsegmentlist = segments.segmentlist(seglst)
    earthsegmentlist.coalesce()
    seglst=[segments.segment(times_mjd[i],times_mjd[i+1]) for i in range(len(times_mjd)-1) if Sun_AngDist[i+1]>=sun_constraint and Sun_AngDist[i+1]>=sun_radii[i+1]]
    sunsegmentlist = segments.segmentlist(seglst)
    sunsegmentlist.coalesce()

    # Coalesce with detector availabilities
    tilesegmentlistdic = segments.segmentlistdict()
    tilesegmentlistdic["observations"] = segmentlist
    tilesegmentlistdic["moon"] = moonsegmentlist
    tilesegmentlistdic["earth"] = earthsegmentlist
    tilesegmentlistdic["sun"] = sunsegmentlist
    tilesegmentlist = tilesegmentlistdic.intersection(["observations","moon","earth","sun"])
    tilesegmentlist.coalesce()

    return tilesegmentlist

def get_segments_tiles(params, config_struct, tile_struct, verbose=True):

    # Get rough details of sun, earth, moon radecs in ATHENA frame through orbit
    config_struct = get_telescope_orbit(config_struct,verbose=verbose)

    segmentlist = config_struct["segmentlist"]

    if verbose: print("Generating Athena segments for tiles...")

    keys = tile_struct.keys()

    ras = [tile_struct[key]["ra"] for key in keys]
    decs = [tile_struct[key]["dec"] for key in keys]

    # Convert to RA, Dec. -- changing to a class of astropy for angular representations..
    # this could be useful to implement when carrying out orbit calculation...
    radecs = astropy.coordinates.SkyCoord(
            ra=np.array(ras)*u.degree, dec=np.array(decs)*u.degree, frame='icrs')

    if params["doParallel"]:
        tilesegmentlists = Parallel(n_jobs=params["Ncores"])(delayed(get_segments_tile)(config_struct, radec, segmentlist) for radec in radecs)
        for ii,key in enumerate(keys):
            tile_struct[key]["segmentlist"] = tilesegmentlists[ii]
    else:

        for ii,key in enumerate(keys):
            #if np.mod(ii,100) == 0:
            #    print("Generating segments for tile %d/%d"%(ii+1,len(radecs)))
            radec = radecs[ii]

            if params["doMinimalTiling"]:
                if ii == 0:
                    keys_computed = [key]
                    radecs_computed = np.atleast_2d([radec.ra.value, radec.dec.value])
                    tilesegmentlist = get_segments_tile(config_struct, radec, segmentlist)
                    tile_struct[key]["segmentlist"] = tilesegmentlist
                else:
                    seps = angular_distance(radec.ra.value, radec.dec.value,
                                            radecs_computed[:,0],
                                            radecs_computed[:,1])
                    sepmin = np.min(seps)
                    sepamin = np.argmin(seps)
                    if sepmin <= 5.0:
                        key_computed = keys_computed[sepamin]
                        tile_struct[key]["segmentlist"] = copy.deepcopy(tile_struct[key_computed]["segmentlist"])
                    else:
                        keys_computed.append(key)
                        radecs_computed = np.vstack((radecs_computed,[radec.ra.value, radec.dec.value]))
                        tilesegmentlist = get_segments_tile(config_struct, radec, segmentlist)
                        tile_struct[key]["segmentlist"] = tilesegmentlist

            else:
                tilesegmentlist = get_segments_tile(config_struct, radec, segmentlist)
                tile_struct[key]["segmentlist"] = tilesegmentlist

    return tile_struct

def get_telescope_orbit(config_struct,SAVETOFILE=False,verbose=True):
    """
    Overhead function to either load a custom orbit from
    config_struct["orbitFile"]=full_file_path,
    or calculate an orbit approximation from additional parameters
    stored in config struct.

    NB : 'SAVETOFILE' is an override flag that forces overwrite of orbitFile.
    """
    if os.path.isfile(config_struct["orbitFile"]) and not SAVETOFILE:
        with open(config_struct["orbitFile"], 'rb') as f:
            config_struct["orbit_dict"] = pickle.load(f)
    else:
        config_struct=calc_telescope_orbit(config_struct,SAVETOFILE,verbose=verbose)

    return config_struct

def calc_telescope_orbit(config_struct,SAVETOFILE,verbose=True):
    """
    NB: heliospheric long and latitude assumed here. config_struct has a number of
        added parameters for orbit calculation... Consider adding a switch for
        orbit appromitations so we can put in different approximation models.
        custom orbits can be loaded by assigning
        "config_struct["orbitFile"]=full_file_path"
        and calling "get_telescope_orbit(config_struct)"".

    TO DO: Convert SYNEX detector locations from SSB oordinates to geocentric
    coords. think there is a module in ephem to do this. config_struct["OrbitParams"]
    should contain these directly, therefore their calculation should be within the
    SYNEX class generation as a seperate function to create gwemopt-specific utilities.

    Angles measured with reference to plane whose normal is the along Sun-Earth axis,
    and whose x-y axes are defined by the intersection with Sun-Earth ecliptic
    and the normal to the x-axis colinear with terrestial north. This makes the
    Earth->Sun direction the z-axis.
    Note that we assume here the orbit is elliptic (on it's long side) but centred
    on the L2 point. No idea if this is ok but when you look at SOHO's orbit as an
    example this seems to be the case.
    """

    # Starting and end times to calculate orbit over
    # Always start from where phi0 is taken -- start of science measurements
    t = Time(config_struct["gps_science_start"], format='gps', scale='utc')
    date_start = ephem.Date(t.iso)
    if "segmentlist" in config_struct:
        # Use the end of observation segments as end time
        segmentlist = config_struct["segmentlist"]
        date_end = ephem.Date(Time(segmentlist[-1][1], format='mjd', scale='utc').iso)
    else:
        # Use gps start time of science measurements
        date_end = ephem.Date(Time(t.mjd+config_struct["mission_duration"]*364.25, format='mjd', scale='utc').iso)

    # Get number of intervals as function of latency time
    # so we know how things move between tile changes...
    interval_in_days=float(date_end)-float(date_start)
    n_times = int((interval_in_days*24.*60.*60)//config_struct["exposuretime"])

    # Check we aren't calculating some stupid number of points...
    if n_times>5000: n_times=5000

    # Set the times - make sure we get the start and end times too. LISA/Athena sources might only have a few days to work with
    times = np.linspace(0.,interval_in_days,n_times+1)

    # Athena-only motion if not frozen
    if not config_struct["frozenAthena"]:
        L2_phis =np.pi*2.*times/config_struct["period"] + config_struct["phi_0"]*np.pi/180.
        L2_orbit_dist = config_struct["semi_maj"]*(1.-config_struct["eccentricity"]**2)/(1.+config_struct["eccentricity"]*np.cos(L2_phis))
        xs_L2 = L2_orbit_dist*np.cos(L2_phis)
        xs_L2 += config_struct["semi_maj"]-np.amax(xs_L2) # shift so origin is at centre of ellipse and not focal point
        ys_L2 = L2_orbit_dist*np.sin(L2_phis)
        zs_L2 = np.zeros(np.shape(xs_L2))

        # Rotate to periapsis argument
        from scipy.spatial.transform import Rotation as R
        rotation_radians = np.radians(config_struct["ArgPeriapsis"])
        rotation_axis = np.array([0., 0., 1.])
        rotation_vector = rotation_radians * rotation_axis
        rotation = R.from_rotvec(rotation_vector)
        OrbPlane_periapsis = rotation.apply(np.stack((xs_L2,ys_L2,zs_L2), axis = 1))

        # Solution in reference plane (rotated for inclination)
        rotation_axis = np.array([np.cos(config_struct["AscendingNode"]*np.pi/180.), np.sin(config_struct["AscendingNode"]*np.pi/180.), 0.])
        rotation_vector = (config_struct["inc"]*np.pi/180.)*rotation_axis # rotation_radians * rotation_axis
        rotation = R.from_rotvec(rotation_vector)
        RefPlane_orbit_dist = rotation.apply(OrbPlane_periapsis)
    else:
        RefPlane_orbit_dist=np.array([[0.,0.,0.]]*len(times))

    # Solution rotated again to have x along Earth-Sun axis, y tangent to Earth orbit, z directed up from Sun-Earth ecliptic plane.
    RefPlane_orbit_EclFrame = np.transpose(np.array([(-vec[2],-vec[0],vec[1]) for vec in RefPlane_orbit_dist]))

    # Decouple from Earth orbit
    # NOTE : ephem gives helio longitude and latitude relative to Earth for Sun AND Moon classes. Earth dist is in AU.
    # NOTE : hlon in [0,2*pi], hlat in [-pi/2,pi/2]
    # NOTE : RA in [-pi,pi], dec in [-pi/2,pi/2]
    # NOTE : phi in [-pi,pi], theta in [0,pi] AND output of hp.vec2dir is [theta,phi], not [phi,theta].
    Earth_Orbit = np.transpose(np.array([(ephem.Sun(ephem.Date(date_start+dt)).hlon,ephem.Sun(ephem.Date(date_start+dt)).hlat,ephem.Sun(ephem.Date(date_start+dt)).earth_distance) for dt in times])) # Can we optimize this to not init class 'Sun' twice for iteration?
    UnitVec_sun_from_earth = hp.dir2vec(np.pi+Earth_Orbit[0,:]*180./np.pi,-Earth_Orbit[1,:]*180./np.pi,lonlat=True)
    Earth_From_Athena = RefPlane_orbit_EclFrame+UnitVec_sun_from_earth*config_struct["L2_from_Earth"]
    Sun_From_Athena = Earth_From_Athena+UnitVec_sun_from_earth*consts.au.value

    # Get position of moon through orbit
    # NOTE : ephem gives helio longitude and latitude relative to Earth for Sun and Moon classes. Earth dist is in AU.
    Moon_Orbit = np.transpose(np.array([(ephem.Moon(ephem.Date(date_start+dt)).hlon,ephem.Moon(ephem.Date(date_start+dt)).hlat,ephem.Moon(ephem.Date(date_start+dt)).earth_distance) for dt in times])) # Can we optimize this to not init class 'Sun' twice for iteration?
    UnitVec_moon_from_earth = hp.dir2vec(Moon_Orbit[:-1,:]*180./np.pi,lonlat=True)
    Moon_From_Athena = Earth_From_Athena+UnitVec_moon_from_earth*Moon_Orbit[2,:]*consts.au.value

    # Translate everything to theta-phis
    Moon_From_Athena_thetaphi = hp.vec2dir(Moon_From_Athena)
    Earth_From_Athena_thetaphi = hp.vec2dir(Earth_From_Athena)
    Sun_From_Athena_thetaphi = hp.vec2dir(Sun_From_Athena)

    # Change from theta-phis to ra-decs - check we had radians or degrees...
    Moon_From_Athena_radecs=np.empty(np.shape(Moon_From_Athena_thetaphi))
    Earth_From_Athena_radecs=np.empty(np.shape(Earth_From_Athena_thetaphi))
    Sun_From_Athena_radecs=np.empty(np.shape(Sun_From_Athena_thetaphi))
    Moon_From_Athena_radecs[0,:] = np.rad2deg(Moon_From_Athena_thetaphi[1,:]) # phi
    Moon_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Moon_From_Athena_thetaphi[0,:]) # theta
    Earth_From_Athena_radecs[0,:] = np.rad2deg(Earth_From_Athena_thetaphi[1,:])
    Earth_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Earth_From_Athena_thetaphi[0,:])
    Sun_From_Athena_radecs[0,:] = np.rad2deg(Sun_From_Athena_thetaphi[1,:])
    Sun_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Sun_From_Athena_thetaphi[0,:])

    # Gather for return - do we want to put these inside config_struct ?
    times_mjd=[t.mjd+dt for dt in times]
    astrophysical_bodies_from_athena_radecs={
                        "elapsed_time_from_start":times,
                        "times_mjd":times_mjd,
                        "Moon_From_Athena_radecs":Moon_From_Athena_radecs,
                        "Earth_From_Athena_radecs":Earth_From_Athena_radecs,
                        "Sun_From_Athena_radecs":Sun_From_Athena_radecs,
                        "Moon_From_Athena":Moon_From_Athena,
                        "Earth_From_Athena":Earth_From_Athena,
                        "Sun_From_Athena":Sun_From_Athena
                        }

    # Put into config struct to keep things tidy
    config_struct["orbit_dict"]=astrophysical_bodies_from_athena_radecs

    # Write to file
    if SAVETOFILE:
        with open(config_struct["orbitFile"], 'wb') as f:
            pickle.dump(astrophysical_bodies_from_athena_radecs, f)
        if verbose: print("Saved orbit to :",config_struct["orbitFile"])

    return config_struct

def angular_distance(ra1, dec1, ra2, dec2):

    delt_lon = (ra1 - ra2)*np.pi/180.
    delt_lat = (dec1 - dec2)*np.pi/180.
    dist = 2.0*np.arcsin( np.sqrt( np.sin(delt_lat/2.0)**2 + \
         np.cos(dec1*np.pi/180.)*np.cos(dec2*np.pi/180.)*np.sin(delt_lon/2.0)**2 ) )

    return dist/np.pi*180.
