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
from astropy.time import Time, TimeDelta
import astropy.units as u
from joblib import Parallel, delayed
import ephem

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

        # This block of modifications should not be necessary after modifications to "get_segments"
        # seg=segments.segment(params["config"][telescope]["segmentlist"][0][0],params["config"][telescope]["segmentlist"][-1][1])
        # params["config"][telescope]["segmentlist"].clear()
        # params["config"][telescope]["segmentlist"].append(seg)

        params["config"][telescope]["exposurelist"] = gwemopt.utils.get_exposures(params, params["config"][telescope], params["config"][telescope]["segmentlist"])
        if len(params["config"][telescope]["exposurelist"]) == 0:
            params["config"][telescope]["n_windows"] = 0
            params["config"][telescope]["tot_obs_time"] = 0.0
            continue

        nexp, junk = np.array(params["config"][telescope]["exposurelist"]).shape
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

    ###
    #
    # Maybe the next modifications are not needed since it returns just the above
    # segmentlist if we specify sat_sun_restriction to be non-zero.
    # the lat and long and elevs are used only to check when it's night...
    #
    ###

    # Extrace las, longs, elevations for orbit
    lats = [str(l) for l in config_struct["latitude"]] # List of strings instead, update in loop
    longs = [str(l) for l in config_struct["longitude"]] # List of strings instead, update in loop
    elevs = [el for el in config_struct["elevation"]] # List of elevations if we have full orbit info, otherwise a list with just one element, update in loop

    observer = ephem.Observer()
    OrbitalCount = 0
    observer.lat = lats[OrbitalCount] # str(config_struct["latitude"])
    observer.lon = longs[OrbitalCount] # str(config_struct["longitude"])
    # observer.horizon = str(-12.0) # not needed for us
    observer.elevation = elevs[OrbitalCount] # config_struct["elevation"]

    date_start = ephem.Date(Time(segmentlist[0][0], format='mjd', scale='utc').iso)
    date_end = ephem.Date(Time(segmentlist[-1][1], format='mjd', scale='utc').iso)
    observer.date = ephem.Date(Time(segmentlist[0][0], format='mjd', scale='utc').iso)

    sun = ephem.Sun()
    nightsegmentlist = segments.segmentlist()
    while date_start < date_end:
        date_rise = observer.next_rising(sun, start = date_start)
        date_set = observer.next_setting(sun, start = date_start)
        if date_set > date_rise:
            date_set = observer.previous_setting(sun, start = date_start)

        astropy_rise = Time(date_rise.datetime(), scale='utc').mjd
        astropy_set  = Time(date_set.datetime(), scale='utc').mjd

        segment = segments.segment(astropy_set,astropy_rise)
        nightsegmentlist = nightsegmentlist + segments.segmentlist([segment])
        nightsegmentlist.coalesce()

        date_start = date_rise
        observer.date = date_rise

        # Update position of telescope
        OrbitalCount+=1
        observer.lat = lats[OrbitalCount]
        observer.lon = longs[OrbitalCount]
        observer.elevation = elevs[OrbitalCount]

    segmentlistdic = segments.segmentlistdict()
    segmentlistdic["observations"] = segmentlist
    segmentlistdic["night"] = nightsegmentlist

    #load the sun retriction for a satelite
    try:
        sat_sun_restriction = config_struct["sat_sun_restriction"]
    except:
        sat_sun_restriction = 0.0

    #in the case of satellite use don't intersect with night segment and take all observation time available
    if sat_sun_restriction:

        segmentlist.coalesce()

        return segmentlist

    segmentlist = segmentlistdic.intersection(["observations","night"])
    segmentlist.coalesce()

    return segmentlist

def get_segments_tile(config_struct, observatory, radec, segmentlist):
    """
    NB: "observatory" input is not used...
    NB: we assume heliospheric long and latitude...

    This function is heavily modified because the functions are just not
    useful as are in "segments.py" for a telescope like Athena.

    TO DO: Convert SYNEX detector locations from SSB oordinates to geocentric
    coords. think there is a module in ephem to do this. config_struct["OrbitParams"]
    should contain these directly, therefore their calculation should be within the
    SYNEX class generation as a seperate function to create gwemopt-specific utilities.
    """
    # Do we have constraints for Sun, Earth, and / or Moon proximity?
    if "moon_constraint" in config_struct:
        moon_constraint = float(config_struct["moon_constraint"])
    else:
        moon_constraint = 20.0

    if "earth_constraint" in config_struct:
        earth_constraint = float(config_struct["earth_constraint"])
    else:
        earth_constraint = 30.0 # No idea about this number

    if "sun_constraint" in config_struct:
        sun_constraint = float(config_struct["sun_constraint"])
    elif "sat_sun_restriction" in config_struct:
        sun_constraint = float(config_struct["sat_sun_restriction"])
    else:
        # Bottom of page 46 in https://arxiv.org/pdf/1903.04083.pdf states AXIS (NOT Athena...) needs 45 deg at least from Sun for resolution and thermal control...
        sun_constraint = 45.0

    #######
    #
    # Generate orbit... New function to move to Detector class? Include possability to load and save from/to file.
    #
    #######

    # Angles measured with reference to plane whose normal is the along Sun-Earth axis,
    # and whose x-y axes are defined by the intersection with Sun-Earth ecliptic
    # and the normal to the x-axis colinear with terrestial north. This makes the
    # Earth->Sun direction the z-axis.
    # Note that we assume here the orbit is elliptic (on it's long side) but centred on the L2 point. No idea if this is ok but when
    # you look at SOHO's orbit as an example this seems to be the case.
    inc = 70. # deg (raised from reference plane from axis perpendicular 'north' from Eun)Earth ecliptic plane.
    MeanRadius = 750000000. # meters (from earth-orbit normal axis)
    semi_maj = MeanRadius # equivalent to semimaj axis ONLY IF we say we are really orbiting the centre and not the focal point
    L2_from_Earth = 1500000000. # meters (L2 from Earth)
    eccentricity = 0.8
    ArgPeriapsis = 0. # angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
    AscendingNode = 0.
    phi_0 = 0. # initial phase of Athena when measurments start
    period = 180. # In days, for one complete halo orbit about L2

    # Starting and end times to calculate orbit over
    date_start = ephem.Date(Time(segmentlist[0][0], format='mjd', scale='utc').iso)
    date_end = ephem.Date(Time(segmentlist[-1][1], format='mjd', scale='utc').iso)

    # Get number of 6hr in whole observation period
    n_6hr = (float(date_end)-float(date_start))//0.25

    # Set the times and angles to solve orbit for
    times = np.array([0.125+ii*0.25 for ii in range(n_6hr)])
    L2_phis =np.pi*2.*times/period + phi_0
    angles = L2_phis-ArgPeriapsis

    # Rough solution in orbital plane
    L2_orbit_dist = semi_maj*(1.-eccentricity**2)/(1.+eccentricity*np.cos(angles))
    xs_L2 = L2_orbit_dist*np.cos(L2_phis)
    ys_L2 = L2_orbit_dist*np.sin(L2_phis)
    zs_L2 = np.zeros(np.shape(xs_L2))

    # Solution in reference plane
    import scipy.spatial.transform.Rotation.from_rotvec
    rotation_radians = np.radians(inc)
    rotation_axis = np.array([np.cos(AscendingNode), np.sin(AscendingNode), 0.])
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    RefPlane_orbit_dist = rotation.apply(zip(xs_L2,ys_L2,zs_L2))

    # Solution wrt Earth in Earth-Sun ecliptic plane -- x is now tangent to Earth orbit, y along Earth-Sun axis, z directed up out of Sun-Earth ecliptic plane.
    Earth_dist = [(vec[0],L2_from_Earth-vec[2],vec[1]) for vec in RefPlane_orbit_dist]

    # Solution rotated again to have x radially out from sun (through earth), y perpendicular and in eclipticplane, z directed up from Sun-Earth ecliptic plane.
    # Also adding Earth distance from Sun so we are co-rotating the Sun but in heliocentric coords
    Earth_From_Athena = np.transpose(np.array([(vec[1],-vec[0],vec[3]) for vec in Earth_dist]))

    # Decouple from Earth orbit -- ephem uses degrees and dunno what healpy.dir2vec uses.
    # NOTE : ephem does not give helio longitude and latitude relative to Sun for the Sun or Moon classes...
    # they are both relative to Earth...
    Earth_Orbit = np.transpose(np.array([(ephem.Sun(ephem.Date(d0+dt)).hlon,ephem.Sun(ephem.Date(date_start+dt)).hlat,ephem.Sun(ephem.Date(date_start+dt)).earth_distance) for dt in times])) # Can we optimize this to not init class 'Sun' twice for iteration?
    Sun_From_Athena = Earth_From_Athena+hp.dir2vec(Earth_Orbit[:-1,:],lonlat=True)*Earth_Orbit[2,:]

    # Get position of moon through orbit
    # NOTE : ephem does not give helio longitude and latitude relative to Sun for the Sun or Moon classes...
    # they are both relative to Earth...
    Moon_Orbit = np.transpose(np.array([(ephem.Moon(ephem.Date(d0+dt)).hlon,ephem.Moon(ephem.Date(date_start+dt)).hlat,ephem.Moon(ephem.Date(date_start+dt)).earth_distance) for dt in times])) # Can we optimize this to not init class 'Sun' twice for iteration?
    Moon_From_Athena = Earth_From_Athena+hp.dir2vec(Earth_Orbit[:-1,:],lonlat=True)*Earth_Orbit[2,:]

    # Translate everything to theta-phis
    Moon_From_Athena_radecs = hp.vec2dir(Moon_From_Athena)
    Earth_From_Athena_radecs = hp.vec2dir(Earth_From_Athena)
    Sun_From_Athena_radecs = hp.vec2dir(Sun_From_Athena)

    # Change from theta-phis to ra-decs - check we had radians or degrees...
    Moon_From_Athena_radecs[0,:] = np.rad2deg(Moon_From_Athena_radecs[1,:]) # phi
    Moon_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Moon_From_Athena_radecs[0,:]) # theta
    Earth_From_Athena_radecs[0,:] = np.rad2deg(Moon_From_Athena_radecs[1,:])
    Earth_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Moon_From_Athena_radecs[0,:])
    Sun_From_Athena_radecs[0,:] = np.rad2deg(Moon_From_Athena_radecs[1,:])
    Sun_From_Athena_radecs[1,:] = np.rad2deg(0.5*np.pi - Moon_From_Athena_radecs[0,:])

    # Angular distances from tile direction and three objects
    # Can we pass the 'angular_distance' function an array of directions? I think so if its ndarray...
    Moon_AngDist=angular_distance(Moon_From_Athena_radecs[0,:], Moon_From_Athena_radecs[1,:], radec[0], radec[1])
    Earth_AngDist=angular_distance(Earth_From_Athena_radecs[0,:], Earth_From_Athena_radecs[1,:], radec[0], radec[1])
    Sun_AngDist=angular_distance(Sun_From_Athena_radecs[0,:], Sun_From_Athena_radecs[1,:], radec[0], radec[1])

    # See when we are within the restriction limits
    moon_times=times[np.where(Moon_AngDist>=moon_constraint)]
    earth_times=times[np.where(Earth_AngDist>=earth_constraint)]
    sun_times=times[np.where(Sun_AngDist>=sun_constraint)]

    # Make the segmentlists
    seglst=[segments.segment(moon_times[2*ii],moon_times[2*ii+1]) for ii in range(len(moon_times)//2]
    moonsegmentlist = segments.segmentlist(seglst)
    moonsegmentlist.coalesce()
    seglst=[segments.segment(earth_times[2*ii],earth_times[2*ii+1]) for ii in range(len(earth_times)//2]
    earthsegmentlist = segments.segmentlist(seglst)
    earthsegmentlist.coalesce()
    seglst=[segments.segment(sun_times[2*ii],sun_times[2*ii+1]) for ii in range(len(sun_times)//2]
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

def get_segments_tiles(params, config_struct, tile_struct):

    observatory = astropy.coordinates.EarthLocation(
        lat=config_struct["latitude"]*u.deg, lon=config_struct["longitude"]*u.deg, height=config_struct["elevation"]*u.m)

    segmentlist = config_struct["segmentlist"]

    print("Generating segments for tiles...")

    ras = []
    decs = []
    keys = tile_struct.keys()
    for key in keys:
        ras.append(tile_struct[key]["ra"])
        decs.append(tile_struct[key]["dec"])

    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
            ra=np.array(ras)*u.degree, dec=np.array(decs)*u.degree, frame='icrs')

    if params["doParallel"]:
        tilesegmentlists = Parallel(n_jobs=params["Ncores"])(delayed(get_segments_tile)(config_struct, observatory, radec, segmentlist) for radec in radecs)
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
                    tilesegmentlist = get_segments_tile(config_struct, observatory, radec, segmentlist)
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
                        tilesegmentlist = get_segments_tile(config_struct, observatory, radec, segmentlist)
                        tile_struct[key]["segmentlist"] = tilesegmentlist

            else:
                tilesegmentlist = get_segments_tile(config_struct, observatory, radec, segmentlist)
                tile_struct[key]["segmentlist"] = tilesegmentlist

    return tile_struct

def angular_distance(ra1, dec1, ra2, dec2):

    delt_lon = (ra1 - ra2)*np.pi/180.
    delt_lat = (dec1 - dec2)*np.pi/180.
    dist = 2.0*np.arcsin( np.sqrt( np.sin(delt_lat/2.0)**2 + \
         np.cos(dec1*np.pi/180.)*np.cos(dec2*np.pi/180.)*np.sin(delt_lon/2.0)**2 ) )

    return dist/np.pi*180.
