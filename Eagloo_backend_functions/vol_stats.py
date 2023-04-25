#! /usr/bin/env python

#Created by David Shean, modified by G. Pfundstein


#Little utility for volume change analysis from dz rasters
#Input should already be clipped appropriately (e.g. masked except over glaciers)

import sys
import os
import datetime
import argparse

import numpy as np

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import timelib

#def getparser():
#    parser = argparse.ArgumentParser(description="Compute volume/mass change stats from DEM difference")
#    parser.add_argument('fn', type=str, help='Elevation difference filename (dz.tif)')
#    parser.add_argument('-rho', type=float, default=0.917, help='Density for mass change calculation (default: %(default)s)')
#    return parser

def volume_mass_stats(parser,                              # description
                      fn,                                  # dDEM file
                      rho,                                 # density of the ice
                      path_out,                            # path where the files will be saved
                      title,                               # Title 
                      calculate_uncertainty,               # Boolean, would you to like calculate the uncertainty 
                      uncertainty_of_the_volume_change=0,  # mean standard error of the 
                      error_rho_density=0,                 # error of the density
                      error_glacier_outline=0              # error of the glacier outline
                      ):
    
    #Parser : description
    #fn :  , type : str
    #rho : density for mass change calculation, type : float, default=0.917
    #This is mean density for N Cascades snow
    #rho = 0.5
    #Density of pure ice 
   
    #If number is in kg/m³ rather than g/cc
    if rho > 10.:
        rho /= 1000.

    #Clip negative values to 0
    filt = False 

    src_ds = iolib.fn_getds(fn)
    res = geolib.get_res(src_ds, square=True)[0]
    bma = iolib.ds_getma(src_ds)

    #Attempt to extract t1 and t2 from input filename
    ts = timelib.fn_getdatetime_list(fn)
    #Hardcode timestamps
    #ts = [datetime.datetime(2013,9,10), datetime.datetime(2014,5,14)]

    dt_yr = None
    if len(ts) == 2:
        dt = ts[1] - ts[0]
        year = datetime.timedelta(days=365.25)
        dt_yr = dt.total_seconds()/year.total_seconds()

    #Can add filter here to remove outliers, perc_fltr(0.01, 99.9)
    if filt:
        mask = np.ma.getmaskarray(bma)
        bma[bma < 0] = 0
        bma = np.ma.array(bma, mask=mask)

    #Print out stats
    print('\n')
    stats = malib.print_stats(bma)
    print('\n')

    count = stats[0]
    area = res**2*count
    mean = stats[3]
    med = stats[5]

    s_m3 = np.ma.sum(bma)*res**2 
    s_km3 = s_m3/1E9 
    s_mwe = mean*rho
    s_gt = s_km3*rho
    #s_mm = s_gt/374
    #https://climatesanity.wordpress.com/conversion-factors-for-ice-and-water-mass-and-volume/
    s_mm = -s_gt/360
    
    #uncertainty 
    if calculate_uncertainty:
        mass_balance_error=abs(s_m3)*np.sqrt((uncertainty_of_the_volume_change/s_m3)**2+
                            (error_rho_density/rho)**2+
                            (error_glacier_outline/area)**2
                            )
        mass_balance_error_gt= mass_balance_error/1E9
        sea_level_rise_error= mass_balance_error_gt/360
        s_mwe_error=abs(s_mwe)*np.sqrt((uncertainty_of_the_volume_change/s_m3)**2+
                            (error_rho_density/rho)**2+
                            (error_glacier_outline/area)**2
                            )
    #
    if calculate_uncertainty:
        if dt_yr is not None:
            print(" - Period : %s to %s: %0.2f yr  " % (ts[0], ts[1], dt_yr))
            print(" - Volume difference : %0.0f m³ (%0.0f m³/yr)  " % (s_m3, s_m3/dt_yr))
            print(" - Volume difference : %0.3f km³ (%0.3f km³/yr)  " % (s_km3, s_km3/dt_yr))
            print(" - Density : %0.3f ± g/cc " % rho)
            print(" - Mass difference : %0.5f GT (%0.5f GT/yr) " % (s_gt, s_gt/dt_yr))
            print(" - Contribution to sea level rise : %0.6f mm SLR (%0.6f mm/yr) " % (s_mm, s_mm/dt_yr))
            print(" - Meter water equivalent : %0.5f m.w.e. (%0.5f m.w.e./yr) " % (s_mwe, s_mwe/dt_yr))
            open(f"{path_out}/volume_mass_{title}_results.doc ", 'w').close()
            with open(f"{path_out}/volume_mass_{title}_results.doc", "a") as o:
                o.write(" - Period : ** %s to %s: %0.2f yr**\n\n" % (ts[0], ts[1], dt_yr))
                o.write(" - Area : **%0.2f ± %0.2f km² **\n\n" % (area/1E6, error_glacier_outline/1E6 ))
                o.write(" - Volume difference : **%0.0f ± %0.0f m³ (%0.0f ± %0.0f m³/yr)**\n\n" % (s_m3, uncertainty_of_the_volume_change,s_m3/dt_yr,uncertainty_of_the_volume_change/dt_yr))
                o.write(" - Volume difference : **%0.4f ± %0.4f km³ (%0.4f ± %0.4f km³/yr)**\n\n" % (s_km3,uncertainty_of_the_volume_change/1E9, s_km3/dt_yr, uncertainty_of_the_volume_change/(dt_yr*1E9)))
                o.write(" - Density :** %0.3f ± %0.3f g/cc \n\n" % (rho,error_rho_density))
                o.write(" - Mass difference : **%0.4f ± %0.4f GT (%0.4f ± %0.4f GT/yr)**\n\n" % (s_gt, mass_balance_error_gt, s_gt/dt_yr,mass_balance_error_gt/dt_yr))
                o.write(" - Contribution to sea level rise : **%0.6f ± %0.6f mm SLR (%0.6f ± %0.6f mm/yr)**\n\n" % (s_mm, sea_level_rise_error, s_mm/dt_yr,sea_level_rise_error/dt_yr))
                o.write(" - Meter water equivalent :** %0.4f ± %0.4f m.w.e. (%0.4f ± %0.4f m.w.e./yr)**\n\n" % (s_mwe, s_mwe_error,s_mwe/dt_yr,s_mwe_error/dt_yr))
        else:
            print(" - Area: %0.2f km²" % (area/1E6))
            print(" - Volume difference : %0.0f m³" % s_m3)
            print(" - Volume difference : %0.3f km³" % s_km3) 
            print(" - Density: %0.3f g/cc" % rho)
            print(" - Mass difference : %0.3f GT" % s_gt)
            print(" - Contribution to sea level rise : %0.6f mm SLR" % s_mm)
            print(" - Meter water equivalent : %0.3f m.w.e." % s_mwe)
            open(f"{path_out}/volume_mass_{title}_results.doc", 'w').close()
            with open(f"{path_out}/volume_mass_{title}_results.doc", "a") as o:
                o.write(" - Area : **%0.2f ± %0.2f km² **\n\n" % (area/1E6, error_glacier_outline/1E6 ))
                o.write(" - Volume difference : **%0.0f ± %0.3f m³ **\n\n" % (s_m3, uncertainty_of_the_volume_change) )
                o.write(" - Volume difference : **%0.4f ± %0.4f km³ **\n\n" % (s_km3, uncertainty_of_the_volume_change/1E9)) 
                o.write(" - Density:** %0.3f ± %0.3f g/cc** \n\n" % (rho,error_rho_density))
                o.write(" - Mass difference : **%0.4f ± %0.4f GT **\n\n" % (s_gt,mass_balance_error_gt))
                o.write(" - Contribution to sea level rise : **%0.6f ± %0.6f mm SLR** \n\n" % (s_mm, sea_level_rise_error))
                o.write(" - Meter water equivalent : **%0.4f ± %0.4f m.w.e.** \n\n" % (s_mwe, s_mwe_error))
                
    else: # No uncertainties calculate
        if dt_yr is not None:
            print(" - Period : %s to %s: %0.2f yr  " % (ts[0], ts[1], dt_yr))
            print(" - Volume difference : %0.0f m³ (%0.0f m³/yr)  " % (s_m3, s_m3/dt_yr))
            print(" - Volume difference : %0.3f km³ (%0.3f km³/yr)  " % (s_km3, s_km3/dt_yr))
            print(" - Density : %0.3f g/cc " % rho)
            print(" - Mass difference : %0.3f GT (%0.3f GT/yr) " % (s_gt, s_gt/dt_yr))
            print(" - Contribution to sea level rise : %0.6f mm SLR (%0.6f mm/yr) " % (s_mm, s_mm/dt_yr))
            print(" - Meter water equivalent : %0.3f m.w.e. (%0.3f m.w.e./yr) " % (s_mwe, s_mwe/dt_yr))
            open(f"{path_out}/volume_mass_{title}_results.doc ", 'w').close()
            with open(f"{path_out}/volume_mass_{title}_results.doc", "a") as o:
                o.write(" - Period : **%s to %s: %0.2f yr**\n\n" % (ts[0], ts[1], dt_yr))
                o.write(" - Area : **%0.2f km²** \n\n" % (area/1E6))
                o.write(" - Volume difference : **%0.0f m³ (%0.0f m³/yr)**\n\n" % (s_m3, s_m3/dt_yr))
                o.write(" - Volume difference :** %0.3f km³ (%0.3f km³/yr)**\n\n" % (s_km3, s_km3/dt_yr))
                o.write(" - Density :** %0.3f g/cc**" % (rho))
                o.write(" - Mass difference : **%0.3f GT (%0.3f GT/yr)**\n\n" % (s_gt, s_gt/dt_yr))
                o.write(" - Contribution to sea level rise : **%0.6f mm SLR (%0.6f mm/yr)**\n\n" % (s_mm, s_mm/dt_yr))
                o.write(" - Meter water equivalent : **%0.3f m.w.e. (%0.3f m.w.e./yr)**\n\n" % (s_mwe, s_mwe/dt_yr))
        else:
            print(" - Area: %0.2f km²" % (area/1E6))
            print(" - Volume difference : %0.0f m³" % s_m3)
            print(" - Volume difference : %0.3f km³" % s_km3) 
            print(" - Density: %0.3f g/cc" % rho)
            print(" - Mass difference : %0.3f GT" % s_gt)
            print(" - Contribution to sea level rise : %0.6f mm SLR" % s_mm)
            print(" - Meter water equivalent : %0.3f m.w.e." % s_mwe)
            open(f"{path_out}/volume_mass_{title}_results.doc", 'w').close()
            with open(f"{path_out}/volume_mass_{title}_results.doc", "a") as o:
                o.write(" - Area :**%0.2f km² **\n\n" % (area/1E6))
                o.write(" - Volume difference : **%0.0f m³ **\n\n" % s_m3)
                o.write(" - Volume difference :** %0.3f km³ **\n\n" % s_km3) 
                o.write(" - Density: **%0.3f g/cc **\n\n" % (rho))
                o.write(" - Mass difference :** %0.3f GT** \n\n" % (s_gt))
                o.write(" - Contribution to sea level rise : **%0.6f mm SLR** \n\n" % (s_mm))
                o.write(" - Meter water equivalent : **%0.3f m.w.e.** \n\n" % (s_mwe))
        
    print('\n')

