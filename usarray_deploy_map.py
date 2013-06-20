#!/usr/bin/env python

"""
@name       usarray_deploy_map.py
@author     Rob Newman <rlnewman@ucsd.edu>
@created    2005-11-15
@modified   2012-01-09
@notes      1. Was originally a C-shell script, converted to bash, converted 
               to Python as a callable module or stand alone script
            2. Add array of station colors
            3. Add STDIN commands to dynamically plot lat/lon
            4. Dynamic plotting of grid lines
            5. Remove wget request and just use pure Datascope
            6. Add Alaska inset map
            7. Add capability to handle inframet deployment
"""

import sys
import os
from shutil import move
import tempfile
from optparse import OptionParser
from subprocess import call, check_call
# Load datascope functions
sys.path.append(os.environ['ANTELOPE'] + '/data/python/antelope')
import datascope as antdb
from stock import pfupdate, pfget, pfget_arr, epoch2str, epoch, str2epoch
from time import time, gmtime, strftime
import datetime

# return the first year that data is available for a given project
def get_start_year():
  return 2004

def process_command_line(argv):
    """Return a 6-tuple: (verbosity, year, month, maptype, deploytype, size).
    'argv' is a list of arguments, or 'None' for ''sys.argv[1:]''.
    """

    if argv is None:
        argv = sys.argv[1:]

    # Limit the choices to two
    maptypes = ['cumulative', 'rolling']
    deploytypes = ['seismic', 'inframet']

    # Initialize the parser object
    usage = "Usage: %prog [options] YYYY MM"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", action="store_true", dest="verbose", help="verbose output", default=False)
    parser.add_option("-x", action="store_true", dest="debug", help="debug script", default=False)
    parser.add_option("-d", "--deploytype", action="store", type="string", dest="deploytype", help="type of deployment to plot", default="seismic")
    parser.add_option("-t", "--type", action="store", type="string", dest="maptype", help="type of map to plot", default=False)
    parser.add_option("-s", "--size", action="store", type="string", dest="size", help="generate different sizes", default=False)
    options, args = parser.parse_args(argv)

    if options.verbose:
        verbose = True
    else:
        verbose = False

    if options.debug:
        debug = True
    else:
        debug = False

    if options.deploytype:
        deploytype = options.deploytype
        if not deploytype in deploytypes:
            print "Your deployment type ('%s') must be either '%s' or '%s'. Goodbye" % (deploytype, deploytypes[0], deploytypes[1])
            exit()
    else:
        print "You have not defined a deploytype to plot, exiting"
        exit()

    if options.maptype:
        maptype = options.maptype
        if not maptype in maptypes:
            print "Your map type ('%s') must be either '%s' or '%s'. Goodbye" % (maptype, maptypes[0], maptypes[1])
            exit()
        maptype = [maptype]
    else:
        print "You have not defined a maptype to plot, so both types will be created"
        maptype = maptypes

    if options.size:
        size = options.size
    else:
        size = False

    if len(args) < 2:
        today = datetime.date.today()
        m = today.month
        y = today.year
        if m == 1:
          m = 12
          y -= 1
        else:
          m -= 1
        month = m
        year = y
        print "You have not specified a year and/or month."
        print "Using default value of last whole month, which is %d %02d" % (year, month)
    else:
        year = int(args[0])
        month = int(args[1])
        if year < get_start_year():
            print "Your year integer ('%s') must be four characters. Goodbye." % year
            exit()
        if 12 < month < 1:
            print "Bad month number ('%d') specified. Goodbye." % month
            exit()

    return verbose, debug, year, month, maptype, deploytype, size

def generate_times(year, month):
    """Generate start and end time unix timestamps for dbsubsets """

    month=int(month)
    year=int(year)
    next_month = (month + 1) % 12
    next_year = year + 1 if next_month==1 else year
    start_time = str2epoch('%d-%02d-01 00:00' % (year, month))
    end_time = str2epoch('%d-%02d-01 00:00' % (next_year, next_month))
    return start_time, end_time

def generate_inframet_locations(db, mtype, deploytype, year, month, imap=False, verbose=False, debug=False):
    """Generate inframet locations for specific
    periods in time and write out to xy files 
    suitable for GMT
    """
    # Build the Datascope query str. 
    # For some reason this list comprehensions 
    # has to be at the top of a function?
    # Cannot reproduce in independent tests?

    qstr = '|'.join([ '|'.join(v) for k,v in imap.iteritems()])
    start_time, end_time = generate_times(year, month)

    if verbose or debug:
        print "  - generate_inframet_locations(): Infrasound: Searching sitechan table for chans that match: %s" % qstr

    infraptr = antdb.dbopen(db, 'r')

    process_list = [
        'dbopen sitechan',
        'dbjoin deployment',
        'dbjoin site',
        'dbsubset chan=~/(%s)/' % qstr,
        'dbsubset ondate <= %s' % end_time # Remove future deployed stations
    ]

    if mtype == 'rolling':
        process_list.append('dbsubset endtime >= %s' % start_time) # No decommissioned stations for rolling plot
    elif mtype != 'cumulative':
        print "generate_inframet_locations(): Inframet Error: Map type ('%s') is not recognized" % mtype
        exit()

    process_list.append('dbsort sta ondate chan time')

    try:
        infraptr = antdb.dbprocess(infraptr, process_list)
    except Exception,e:
        print "  - generate_inframet_locations(): Dbprocessing failed with exception: %s" % e
    else:
        all_stations = {}

        infra_tmp_all = tempfile.mkstemp(suffix='.xy',
                                         prefix='deployment_list_inframet_ALL_')

        infra_tmp_ncpa = tempfile.mkstemp(suffix='.xy',
                                          prefix='deployment_list_inframet_NCPA_')

        infra_tmp_setra = tempfile.mkstemp(suffix='.xy',
                                           prefix='deployment_list_inframet_SETRA_')

        infra_tmp_mems = tempfile.mkstemp(suffix='.xy',
                                          prefix='deployment_list_inframet_MEMS_')

        file_list = {'complete':infra_tmp_all[1], 'ncpa':infra_tmp_ncpa[1],
                     'setra':infra_tmp_setra[1], 'mems':infra_tmp_mems[1]}

        counter = {'complete':0, 'ncpa':0, 'setra':0, 'mems':0}

        if mtype == 'cumulative':
            infra_tmp_decom = tempfile.mkstemp(
              suffix='.xy',
              prefix='deployment_list_inframet_DECOM_'
            )
            # Add the DECOM by hand as it is a manufactured
            # file, not a snet per se. Call it _DECOM to force
            # it to plot first
            file_list['1_DECOM'] = infra_tmp_decom[1]
            counter['decom'] = 0
        try:
            infraptr_grp = antdb.dbgroup(infraptr, 'sta')
        except Exception,e:
            print "  - generate_inframet_locations(): Dbgroup failed with exception: %s" % e
        else:
            # {{{ Get values into a easily digestible dict
            for i in range(antdb.dbquery(infraptr_grp, antdb.dbRECORD_COUNT)):
                infraptr_grp[3] = i
                sta, [db, view, end_rec, start_rec] = \
                    antdb.dbgetv(infraptr_grp, 'sta', 'bundle')
                all_stations[sta] = {'sensors': {'MEMS':False, 'NCPA':False,
                                                 'SETRA':False},
                                     'location': {'lat':0, 'lon':0}}
                for j in range(start_rec, end_rec):
                    infraptr[3] = j
                    # Cannot use time or endtime as that applies to the station, not to the inframet sensor
                    ondate, offdate, chan, lat, lon = \
                        antdb.dbgetv(infraptr, 'ondate', 'offdate', 'chan', 'lat', 'lon')
                    all_stations[sta]['location']['lat'] = lat
                    all_stations[sta]['location']['lon'] = lon

                    ondate = epoch(ondate)

                    if offdate > 0:
                        offdate = epoch(offdate)
                    else:
                        offdate = 'NULL'

                    if chan == 'LDM_EP':
                        if ondate <= end_time and (offdate == 'NULL' or offdate >= start_time):
                            all_stations[sta]['sensors']['MEMS'] = True
                    elif chan == 'BDF_EP' or chan == 'LDF_EP':
                        if ondate <= end_time and (offdate == 'NULL' or offdate >= start_time):
                            all_stations[sta]['sensors']['NCPA'] = True
                    elif chan == 'BDO_EP' or chan == 'LDO_EP':
                        if ondate <= end_time and (offdate == 'NULL' or offdate > start_time):
                            all_stations[sta]['sensors']['SETRA'] = True
                    else:
                        print "   - ***Channel %s not recognized***" % chan
            # }}}
            if debug:
                print all_stations
            # {{{ Process dict
            for sta in sorted(all_stations.iterkeys()):
                if verbose or debug:
                     print "   - Working on station %s" % sta
                lat = all_stations[sta]['location']['lat']
                lon = all_stations[sta]['location']['lon']
                sensors = all_stations[sta]['sensors']
                if mtype == 'rolling':
                    if sensors['MEMS'] and sensors['NCPA'] and sensors['SETRA']:
                        os.write(infra_tmp_all[0], "%s    %s    # %s \n" % (lat, lon, sta))
                        counter['complete'] += 1
                    elif sensors['MEMS'] and sensors['NCPA']:
                        os.write(infra_tmp_ncpa[0], "%s    %s    # %s \n" % (lat, lon, sta))
                        counter['ncpa'] += 1
                    elif sensors['MEMS'] and sensors['SETRA']:
                        os.write(infra_tmp_setra[0], "%s    %s    # %s \n" % (lat, lon, sta))
                        counter['setra'] += 1
                    elif sensors['MEMS']:
                        os.write(infra_tmp_mems[0], "%s    %s    # %s \n" % (lat, lon, sta))
                        counter['mems'] += 1
                elif mtype == 'cumulative':
                    if not sensors['MEMS'] and not sensors['NCPA'] and not sensors['SETRA']:
                        os.write(infra_tmp_decom[0], "%s    %s    # DECOM %s \n" % (lat, lon, sta))
                        counter['decom'] += 1
                    else:
                        if sensors['MEMS'] and sensors['NCPA'] and sensors['SETRA']:
                            os.write(infra_tmp_all[0], "%s    %s    # %s \n" % (lat, lon, sta))
                            counter['complete'] += 1
                        elif sensors['MEMS'] and sensors['NCPA']:
                            os.write(infra_tmp_ncpa[0], "%s    %s    # %s \n" % (lat, lon, sta))
                            counter['ncpa'] += 1
                        elif sensors['MEMS'] and sensors['SETRA']:
                            os.write(infra_tmp_setra[0], "%s    %s    # %s \n" % (lat, lon, sta))
                            counter['setra'] += 1
                        elif sensors['MEMS']:
                            os.write(infra_tmp_mems[0], "%s    %s    # %s \n" % (lat, lon, sta))
                            counter['mems'] += 1
            # }}}
            os.close(infra_tmp_all[0])
            os.close(infra_tmp_mems[0])
            if mtype == 'cumulative':
                os.close(infra_tmp_decom[0])
            antdb.dbfree(infraptr_grp)
        antdb.dbclose(infraptr)
    return file_list, counter
 
def generate_sta_locations(db, mtype, deploytype, year, month, verbose=False, debug=False):
    """Generate station locations for specific
    periods in time and write out to xy files 
    suitable for GMT
    """
    start_time, end_time = generate_times(year, month)

    # {{{ Get the networks
    snetptr = antdb.dbopen(db, 'r')
    snetptr = antdb.dbprocess(snetptr,
                            ['dbopen site',
                             'dbjoin snetsta',
                             'dbjoin deployment'])
    snetptr = antdb.dbsort(snetptr,'snet', unique=True)
    usnets = []
    try:
        for i in range(antdb.dbquery(snetptr, antdb.dbRECORD_COUNT )):
            snetptr['record'] = i
            mysnet = antdb.dbgetv(snetptr,'snet')[0]
            usnets.append(mysnet)
        antdb.dbclose(snetptr)
    except Exception, e:
        print "generate_sta_locations(): Exception occurred: %s" % e
    # }}}

    # {{{ Define dbops
    process_list = [
        'dbopen site', 
        'dbjoin snetsta', 
        'dbjoin deployment', 
        'dbsubset time <= %s' % end_time
    ]
    dbptr = antdb.dbopen(db, 'r')
    if mtype == 'rolling':
        process_list.append('dbsubset endtime >= %s' % start_time)
    elif mtype != 'cumulative':
        print "generate_sta_locations(): Map type ('%s') is not recognized" % mtype
        exit()
    process_list.append('dbsort snet sta')
    dbptr = antdb.dbprocess(dbptr,process_list)
    # }}}

    file_list = {}
    counter = {}

    dfile = tempfile.mkstemp(suffix='.xy', prefix='deployment_list_DECOM_')
    decom_ptr = dfile[0]
    decom_name = dfile[1]

    if mtype == 'cumulative':
        this_decom_counter = 0

    # {{{ Loop over unqiue snets
    for s in usnets:
        stmp = tempfile.mkstemp(suffix='.xy',
                                prefix='deployment_list_%s_' % s)
        file_ptr = stmp[0]
        file_name = stmp[1]
        if verbose:
            print "generate_sta_locations(): Working on network: %s" % s

        try:
            dbptr_snet = antdb.dbsubset(dbptr, 'snet=~/%s/' % s )
        except Exception, e:
            print "Error occurred: %s" % e
        else:
            this_counter = 0
            for i in range( antdb.dbquery( dbptr_snet,
                                          antdb.dbRECORD_COUNT ) ):
                dbptr_snet[3] = i
                if mtype == 'rolling':
                    sta, lat, lon, snet = antdb.dbgetv(dbptr_snet,'sta', 'lat',
                                                          'lon', 'snet')
                    os.write(file_ptr, "%s    %s    # %s %s\n" % (lat, lon,
                                                                  snet, sta))
                    this_counter = this_counter + 1
                elif mtype == 'cumulative':
                    sta, lat, lon, snet, sta_time, sta_endtime = antdb.dbgetv(
                      dbptr_snet, 'sta', 'lat', 'lon', 'snet', 'time', 'endtime')
                    if sta_endtime >= start_time:
                        os.write(file_ptr, "%s    %s    # %s %s\n" %
                                 (lat, lon, snet, sta))
                        this_counter = this_counter + 1
                    else:
                        os.write(decom_ptr, "%s    %s    # DECOM %s\n" %
                                 (lat, lon, sta))
                        this_decom_counter = this_decom_counter + 1
            counter[s] = this_counter
            os.close(file_ptr)
            file_list[s] = file_name
    if mtype == 'cumulative':
        counter['decom'] = this_decom_counter
    # }}}

    antdb.dbclose(dbptr)

    # Add the DECOM by hand as it is a manufactured 
    # file, not a snet per se. Call it _DECOM to force
    # it plot first
    file_list['1_DECOM'] = decom_name
    os.close(decom_ptr)

    return file_list, counter

def main(argv=None):
    """Main processing script
    for all maps
    """
    print "Start of script at time %s" % strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

    verbose, debug, year, month, maptype, deploytype, size = process_command_line(argv)
        
    if debug:
        print "\n*** DEBUGGING ON ***"
        print "*** No grd or grad files - just single color for speed ***\n"

    common_pf = 'common.pf'
    stations_pf = 'stations.pf'

    print " - Creating **%s** maps" % deploytype
    if verbose:
        print " - Parse configuration parameter file (%s)" % common_pf
        print " - Parse stations parameter file (%s)" % stations_pf

    wet_rgb = '202/255/255'

    pfupdate(common_pf)
    pfupdate(stations_pf)

    dbmaster = pfget(common_pf, 'USARRAY_DBMASTER')
    networks = pfget_arr(stations_pf, 'network')
    infrasound = pfget_arr(stations_pf, 'infrasound')
    colors = pfget_arr(stations_pf, 'colors')
    # Force the tmp dir environmental variable
    tmp = pfget(common_pf, 'TMP')
    os.environ['TMPDIR'] = os.environ['TEMP'] = os.environ['TMP'] = tmp
    gmtbindir = pfget(common_pf, 'GMT_BIN')
    usa_coords = pfget_arr(common_pf, 'USACOORDS')
    ak_coords = pfget_arr(common_pf, 'AKCOORDS')
    web_output_dir = pfget(common_pf, 'CACHE_MONTHLY_DEPLOYMENT')
    web_output_dir_infra = pfget(common_pf, 'CACHE_MONTHLY_DEPLOYMENT_INFRA')
    infrasound_mapping = pfget(common_pf, 'INFRASOUND_MAPPING')
    output_dir = '/var/tmp' # FOR TESTING
    sys.path.append(gmtbindir)

    # Make sure execution occurs in the right directory
    cwd = os.getcwd()
    path_parts = cwd.split('/')
    if path_parts[-1] == 'deployment_history' and path_parts[-2] == 'bin':
        if verbose or debug:
            print ' - Already in the correct current working directory %s' % cwd
    else:
        cwd = os.getcwd() + '/bin/deployment_history'
        if verbose or debug:
            print ' - Changed current working directory to %s' % cwd
        os.chdir(cwd)

    for m in maptype:
        if size == 'wario':
            ps = tempfile.mkstemp(suffix='.ps', prefix='deployment_history_map_%s_%d_%02d_%s_WARIO_' % (deploytype, year, month, m))
            png = 'PNG not created for tiled display wario. Create by hand in Photoshop'
        else:
            ps = tempfile.mkstemp(suffix='.ps', prefix='deployment_history_map_%s_%d_%02d_%s_' % (deploytype, year, month, m))
            if deploytype == 'inframet':
                finalfile = 'deploymap_%s_%d_%02d.%s.png' % (deploytype, year, month, m)
            else:
                finalfile = 'deploymap_%d_%02d.%s.png' % (year, month, m)
            png = '%s/%s' % (output_dir, finalfile)

        if verbose or debug or size:
            print ' - Working on maptype: %s' % m
            print ' - Temp postscript file: %s' % ps[1]
            print ' - Output target: %s' % png

        # Make sure we set some GMT parameters for just this script
        # {{{ GMTSET
        # Plot media
        if size == 'wario':
            paper_orientation = 'landscape'
            paper_media = 'b0'
        else:
            paper_orientation = 'portrait'
            paper_media = 'a1'
        try:
            retcode = check_call("gmtset"+" PAGE_COLOR 255/255/255 PAGE_ORIENTATION %s PAPER_MEDIA %s" % (paper_orientation, paper_media), shell=True)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

        # Basemap Anotation Parameters
        try:
            retcode = check_call("gmtset"+" ANNOT_OFFSET_PRIMARY 0.2c ANNOT_OFFSET_SECONDARY 0.2c LABEL_OFFSET 0.2c", shell=True)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

        # Basemap Layout Parameters
        try:
            retcode = check_call("gmtset"+" FRAME_WIDTH 0.2c MAP_SCALE_HEIGHT 0.2c TICK_LENGTH 0.2c X_AXIS_LENGTH 25c Y_AXIS_LENGTH 15c X_ORIGIN 2.5c Y_ORIGIN 2.5c UNIX_TIME_POS -0.2c/-0.2c", shell=True)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

        # Miscellaneous
        try:
            retcode = check_call("gmtset"+" LINE_STEP 0.025c MEASURE_UNIT inch", shell=True)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e
        # }}} GMTSET

        # Determine region of interest and center of plot
        # The lat and lon padding ensures we get full topo and bathy.
        minlon = int(usa_coords['MINLON'])
        maxlon = int(usa_coords['MAXLON'])
        minlat = int(usa_coords['MINLAT'])
        maxlat = int(usa_coords['MAXLAT'])
        region = '%s/%s/%s/%s' % (minlon, minlat, maxlon, maxlat) + 'r'
        centerlat = (maxlat - minlat)/2 + minlat
        centerlon = (maxlon - minlon)/2 + minlon

        ak_minlon = int(ak_coords['MINLON'])
        ak_maxlon = int(ak_coords['MAXLON'])
        ak_minlat = int(ak_coords['MINLAT'])
        ak_maxlat = int(ak_coords['MAXLAT'])
        ak_region = '%s/%s/%s/%s' % (ak_minlon, ak_minlat, ak_maxlon, ak_maxlat) + 'r'
        ak_centerlat = (ak_maxlat - ak_minlat)/2 + ak_minlat
        ak_centerlon = (ak_maxlon - ak_minlon)/2 + ak_minlon

        if size == 'wario':
            center = '%s/%s/%s' % (centerlon, centerlat, '44') + 'i'
            ak_center = '%s/%s/%s' % (ak_centerlon, ak_centerlat, '10') + 'i'
        else:
            center = '%s/%s/%s' % (centerlon, centerlat, usa_coords['WIDTH']) + 'i'
            ak_center = '%s/%s/%s' % (ak_centerlon, ak_centerlat, ak_coords['WIDTH']) + 'i'

        if verbose or debug:
            print ' - GMT USA region string: %s' % region
            print ' - GMT USA center location string: %s' % center
            print ' - GMT AK region string: %s' % ak_region
            print ' - GMT AK center location string: %s' % ak_center

        if deploytype == 'seismic':
            station_loc_files, counter = generate_sta_locations(dbmaster, m, deploytype, year, month, verbose, debug)
            rgbs = {'1_DECOM':'77/77/77'} # Init with the missing color and force to be first plotted
        elif deploytype == 'inframet':
            station_loc_files, counter = generate_inframet_locations(dbmaster, m, deploytype, year, month, infrasound_mapping, verbose, debug)
            rgbs = {'1_DECOM':'255/255/255'} # Init with the missing color and force to be first plotted

        snets_text = {}
        for key in sorted(station_loc_files.iterkeys()):
            if deploytype == 'seismic':
                if key in networks:
                      color = networks[key]['color']
                      rgbs[key] = colors[color]['rgb'].replace(',', '/')
                      snets_text[key] = networks[key]['abbrev'].replace(' ', '\ ')
            elif deploytype == 'inframet':
                if debug:
                    print "\tWorking on inframet key: %s" % key
                if key in infrasound:
                      color = infrasound[key]['color']
                      rgbs[key] = colors[color]['rgb'].replace(',', '/')
                      snets_text[key] = infrasound[key]['name'].replace(' ', '\ ')
        # Extra key for the decommissioned stations group
        if m == 'cumulative':
            if deploytype == 'seismic':
               color = networks['DECOM']['color']
               rgbs['decom'] = colors[color]['rgb'].replace(',', '/')
            elif deploytype == 'inframet':
               color = infrasound['decom']['color']
               rgbs['decom'] = colors[color]['rgb'].replace(',', '/')

        if verbose or debug:
            print ' - Working on contiguous United States'

        # Create the contiguous United States topography basemap

        # {{{ Contiguous United States

        if debug == True:
            try:
                retcode = check_call("pscoast -R%s -JE%s -Df -A5000 -S%s -G40/200/40 -V -X2 -Y2 -K >> %s" % (region, center, wet_rgb, ps[1]), shell=True)
            except OSError, e:
                print >>sys.stderr, "pscoast for contiguous United States execution failed:", e
        else:
            try:
                retcode = check_call("grdimage data/usa.grd -R%s -JE%s -Cdata/land_ocean.cpt -Idata/usa.grad -V -E100 -X2 -Y2 -K >> %s" % (region, center, ps[1]), shell=True)
            except OSError, e:
                print >>sys.stderr, "grdimage for usa.grd execution failed:", e

        # {{{ START: PLOT TOPO, BATHYMETRY CORRECTLY
        # Plot land areas below sea level correctly

        # {{{ START: Salton Sea co-ords -R-116.8/-115/32/34
        # Define a clip region
        try:
            retcode = check_call("psclip data/saltonsea.xy -R%s -JE%s -V -K -O >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Salton Sea psclip execution failed:", e

        # Make the Salton Sea be 'land-only' and put into the clipping region
        try:
            retcode = check_call("grdimage data/saltonsea.grd -V -R%s -JE%s -Cdata/land_only.cpt -Idata/saltonsea.grad -O -K >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Salton Sea grdimage execution failed:", e

        # Color the Salton Sea Blue
        try:
            retcode = check_call("pscoast -V -R%s -JE%s -C%s -Df -O -K >> %s" % (region, center, wet_rgb, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Salton Sea pscoast execution failed:", e

        # Close psclip
        try:
            retcode = check_call("psclip -C -K -O >> %s" % ps[1], shell=True)
        except OSError, e:
            print >>sys.stderr, "Salton Sea psclip execution failed:", e

        # }}} END: Salton Sea

        # {{{ START: Death Valley co-ords -R
        # Define a clip region
        try:
            retcode = check_call("psclip data/deathvalley.xy -R%s -JE%s -V -K -O >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Death Valley psclip execution failed:", e

        # Make Death Valley be 'land-only' and put into the clipping region
        try:
            retcode = check_call("grdimage data/deathvalley.grd -V -R%s -JE%s -Cdata/land_only.cpt -Idata/deathvalley.grad -O -K >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Death Valley grdimage execution failed:", e

        # Color the Salton Sea Blue
        try:
            retcode = check_call("pscoast -V -R%s -JE%s -C%s -Df -O -K >> %s" % (region, center, wet_rgb, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Death Valley pscoast execution failed:", e

        # Close psclip
        try:
            retcode = check_call("psclip -C -K -O >> %s" % ps[1], shell=True)
        except OSError, e:
            print >>sys.stderr, "Death Valley psclip execution failed:", e

        # }}} END: Death Valley

        # }}} END: PLOT TOPO, BATHYMETRY CORRECTLY

        # {{{ Plot wet areas and coastline
        # Plot wet areas (not coast)
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -W0.5p,%s -S%s -A0/2/4 -Df -O -K >> %s" % (region, center, wet_rgb, wet_rgb, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot coastline in black
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -W0.5p,0/0/0 -Df -O -K >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot major rivers
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -Ir/0.5p,0/0/255 -Df -O -K >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot national (N1) and state (N2) boundaries
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -N1/5/0/0/0 -N2/1/0/0/0 -Df -O -K >> %s" % (region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # }}} Plot wet areas and coastline

        # {{{ Overlay the grid
        try:
            retcode = check_call("psbasemap -X0 -Y0 -R%s -JE%s -V -Bg%swesn -Lf-75/30/36/500k+l -O -K >> %s" % (region, center, usa_coords['GRIDLINES'], ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "psbasemap execution failed:", e
        # Add stations from local text files
        for key in sorted(station_loc_files.iterkeys()):
            if size == 'wario':
                symsize = '0.3'
            else:
                symsize = '0.15'
            if key == 'IU' or key == 'US': 
                # Plots diamond symbols for US backbone stations
                try:
                    retcode = check_call("psxy %s -R -JE -V -Sd%s -G%s -W -L -O -K -: >> %s" % (station_loc_files[key], symsize, rgbs[key], ps[1]), shell=True)
                except OSError, e:
                    print >>sys.stderr, "psxy execution failed:", e
            else:
                try:
                    retcode = check_call("psxy %s -R -JE -V -St%s -G%s -W -L -O -K -: >> %s" % (station_loc_files[key], symsize, rgbs[key], ps[1]), shell=True)
                except OSError, e:
                    print >>sys.stderr, "psxy execution failed:", e
        # }}} Overlay the grid

        # }}} Contiguous United States

        if verbose or debug:
            print ' - Working on Alaska inset'

        # {{{ Alaska

        if debug == True:
            try:
                retcode = check_call("pscoast -R%s -JE%s -Df -A5000 -S%s -G40/200/40 -V -X0.1i -Y0.1i -O -K >> %s" % (ak_region, ak_center, wet_rgb, ps[1]), shell=True)
            except OSError, e:
                print >>sys.stderr, "pscoast for Alaska execution failed:", e
        else:
            try:
                retcode = check_call("grdimage data/alaska.grd -R%s -JE%s -Cdata/land_ocean.cpt -Idata/alaska.grad -V -E100 -X0.1i -Y0.1i -O -K >> %s" % (ak_region, ak_center, ps[1]), shell=True)
            except OSError, e:
                print >>sys.stderr, "grdimage for alaska.grd execution failed:", e

        # {{{ Plot wet areas and coastline
        # Plot wet areas (not coast)
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -W0.5p,%s -S%s -A0/2/4 -Df -O -K >> %s" % (ak_region, ak_center, wet_rgb, wet_rgb, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot coastline in black
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -W0.5p,0/0/0 -Df -O -K >> %s" % (ak_region, ak_center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot major rivers
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -Ir/0.5p,0/0/255 -Df -O -K >> %s" % (ak_region, ak_center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # Plot national (N1) and state (N2) boundaries
        try:
            retcode = check_call("pscoast"+" -V -R%s -JE%s -N1/5/0/0/0 -N2/1/0/0/0 -Df -O -K >> %s" % (ak_region, ak_center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "pscoast execution failed:", e

        # }}} Plot wet areas and coastline

        # {{{ Overlay the grid
        try:
            retcode = check_call("psbasemap -X0 -Y0 -R%s -JE%s -V -Bg%swesn -Lf-145/57/60/500k+l -O -K >> %s" % (ak_region, ak_center, ak_coords['GRIDLINES'], ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "psbasemap execution failed:", e
        # Add stations from local text files
        for key in sorted(station_loc_files.iterkeys()):
            if size == 'wario':
                symsize = '0.3'
            else:
                symsize = '0.15'
            if key == 'IU' or key == 'US': 
                # Plots diamond symbols for US backbone stations
                try:
                    retcode = check_call("psxy %s -R -JE -V -Sd%s -G%s -W -L -O -K -: >> %s" % (station_loc_files[key], symsize, rgbs[key], ps[1]), shell=True)
                except OSError, e:
                    print >>sys.stderr, "psxy execution failed:", e
            else:
                try:
                    retcode = check_call("psxy %s -R -JE -V -St%s -G%s -W -L -O -K -: >> %s" % (station_loc_files[key], symsize, rgbs[key], ps[1]), shell=True)
                except OSError, e:
                    print >>sys.stderr, "psxy execution failed:", e
        # }}} Overlay the grid

        # }}} Alaska

        # Clean up station files
        for key in sorted(station_loc_files.iterkeys()):
            os.unlink(station_loc_files[key])

        if verbose or debug:
            print ' - Working on year and month timestamp'
        # Create the text files of year & month and legend
        time_file = tempfile.mkstemp(suffix='.txt', prefix='year_month_')
        # time_file = "%syear_month.txt" % tmp
        tf = open(time_file[1], 'w')
        tf.write("-75.5    17    20    0    1    BR    %d\ %02d" % (year, month))
        tf.close()

        if verbose or debug:
            print ' - Working on copyright file'
        copyright_file = tempfile.mkstemp(suffix='.txt', prefix='copyright_')
        # copyright_file = "%scopyright.txt" % tmp
        cf = open(copyright_file[1], 'w')
        cf.write("-67    48.7    11    0    1    BR    (c)\ 2004\ -\ %s\ Array\ Network\ Facility,\ http://anf.ucsd.edu" % year)
        cf.close()

        if verbose or debug:
            print ' - Working on snet files'
        snets_file = tempfile.mkstemp(suffix='.txt', prefix='snets_')
        sf = open(snets_file[1], 'w')
        if size == 'wario':
            legend_symsize = '0.3'
        else:
            legend_symsize = '0.15'
        snet_file_txt = 'G 0.06i\n'
        snet_file_txt += 'H 14 Helvetica-Bold Network Legend\n'
        snet_file_txt += 'D 0.06i 1p\n' # A horizontal line
        snet_file_txt += 'N 1\n'
        snet_file_txt += 'V 0 1p\n'
        for key in sorted(snets_text.iterkeys()):
            if key == 'IU' or key == 'US': 
                snet_symbol = 'd'
            else:
                snet_symbol = 't'
            snet_file_txt += 'S 0.1i %s %s %s 0.25p 0.3i %s\ [%s]\n' % (snet_symbol, legend_symsize, rgbs[key], snets_text[key], counter[key])
        if m == 'cumulative':
            snet_file_txt += 'D 0.06i 1p\n' # A horizontal line
            snet_file_txt += 'S 0.1i %s %s %s 0.25p 0.3i Decommissioned\ [%s]\n' % (snet_symbol, legend_symsize, rgbs['decom'], counter['decom'])
        sf.write(snet_file_txt)
        sf.close()

        # Overlay the copyright notice
        if verbose or debug:
            print ' - Overlay the copyright notice'
        try:
            retcode = check_call("pstext %s -R%s -JE%s -V -D0.25/0.25 -S2p/255/255/255 -P -O -K >> %s" % (copyright_file[1], region, center, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Copyright msg plotting error: pstext execution failed:", e
        else:
            os.unlink(copyright_file[1])

        # Overlay the date legend stamp
        if verbose or debug:
            print ' - Overlay the date legend stamp'
        try:
            retcode = check_call("pstext " + time_file[1] + " -R" + region + " -JE" + center + " -V -D0.25/0.25 -W255/255/255o1p/0/0/0 -C50% -P -O -K >> " + ps[1], shell=True)
        except OSError, e:
            print >>sys.stderr, "Time msg plotting error: pstext execution failed:", e
        else:
            os.unlink(time_file[1])

        # Overlay the snet legend
        if verbose or debug:
            print ' - Overlay the snet legend'

        if deploytype == 'seismic':
            legend_width = '2.6'
            legend_height = '2.98'
        elif deploytype == 'inframet':
            legend_width = '4.2'
            legend_height = '1.7'

        try:
            retcode = check_call("pslegend %s -R%s -JE%s -V -D-90/26/%si/%si/TC -F -G255/255/255 -P -O >> %s" % (snets_file[1], region, center, legend_width, legend_height, ps[1]), shell=True)
        except OSError, e:
            print >>sys.stderr, "Network (snet) legend plotting error: pslegend execution failed:", e
        else:
            os.unlink(snets_file[1])

        # Run Imagemagick convert cmd on Postscript output
        if size == 'wario':
            print " - Your file for wario is ready for photoshop and is called %s" % ps[1]
        else:
            if verbose or debug:
                print " - Running Imagemagick's convert command on postscript file %s" % ps[1]
            try:
                retcode = check_call("convert -trim -depth 16 +repage %s %s" % (ps[1], png), shell=True)
            except OSError, e:
                print >>sys.stderr, "Execution failed:", e
            else:
                os.unlink(ps[1])

            if deploytype == 'inframet':
                web_output_dir = web_output_dir_infra

            if verbose or debug:
                print " - Going to move %s to %s/%s" % (png, web_output_dir, finalfile)
            try:
                move(png, '%s/%s' % (web_output_dir, finalfile))
            except OSError, e:
                print >>sys.stderr, "shutil.move failed:", e
            else:
                print " - Your file is ready and is called %s/%s" % (web_output_dir, finalfile)

    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)
