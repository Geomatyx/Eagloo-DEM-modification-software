#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 17:48:29 2022

@author: Guillaume Pfundstein

V2
"""

# =============================================================================
# Process : 
# 1/ Coregistration with ICP + Nuth and Kaab method 
# 2/ Calculates subtraction between two DEMs
# 3/ Computes hypsometric interpolation
# 4/ Calculates the volume of the ddem and mass (linear conversion)

# Files that are computed in the coregistration step are used in subtraction, interpolation, and volumes calculation steps.
# It does not exist any correlation between interpolation and subtraction steps.
# =============================================================================

import os
import geoutils as gu
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import xdem
import rasterio
import rasterio.mask
import geopandas
import fiona
import Eagloo_backend_functions.coreg_xdem as coreg_xdem
import shutil
import time


from rasterio.transform import from_origin
from datetime import datetime
from osgeo import gdal
from Eagloo_backend_functions.vol_stats import volume_mass_stats as vol_stats 
from fpdf import FPDF
from datetime import datetime
from PIL import Image   #to know th size of image (width /height)
import PySimpleGUI as sg
from time import sleep
import multiprocessing

def time_as_int():
    return int(round(time.time() * 100))

def time_elapse(start_time):
    current_time = time_as_int() - start_time 
    

def Treatment_dem(
                  # =============================================================================
                  #Beggining entries
        
                  dem_1,                # - Type: Geotiff - Oldest DEM used for subtraction, file to co-register
                  dem_2 ,               # - Type: Geotiff - (optionnal for coregistration, None if empty) Latest DEM used for subtraction (dem_2-dem_1), dem_2 is coregistered only if it is different of dem_coreg_ref
                  vector_glacier ,      # - Type: Shapefile -  Glacier mask of the year where the glacier is the biggest
                  path_out ,            # - Type:path -  path where files will be saved
                  coordinate_system ,   # - Type: coordinate_system - ex:EPSG:32633 -  coordinate system used during the gdal.Warp process
                  vector_frame ,        # - Type: Shapefile - Shapefile of the global area, this file is utilized for preprocessing data 
                  resolution ,          # - Type: Integer - Define the resolution, this resolution will be applied to each DEM, x and y are identical
                  coregistration ,      # - Type: Boolean - If "True" coregistration is activated
                  dem_coreg_ref,        # - Type: Geotiff - It is necessary for the coregistration, this DEM is used as a reference for the coregistration,  
                  method_coregistration,# - Type: String - One of these choices must be selected  "Nuth and Kaab","ICP", "ICP + Nuth and Kaab"
                  subtraction,          # - Type: Boolean - Calculate layer dem2-dem1 and print a graph of this one
                  use_dem_coreg_sub,    # - Type: Boolean - if you want select the DEMs calculated during the co-registration step
                  interpolation,        # - Type: Boolean - if True interpolation is activated
                  use_dem_coreg_inter,  # - Type: Boolean - if true, the DEM calculated in the Co-registration step will be used in the interpolation
                  method_interpolation, # - Type: String - One of these choices must be selected  "Nuth and Kaab","ICP", "ICP + Nuth and Kaab"
                  dem_inter_ref,        # - Type: Geotiff - DEM used as a reference for the hypsometric interpolation and coregistration
                  ddem_calculated,      # - Type: Geotiff - dDEM previously calculated by the user
                  rho_density,                  # - Type: Number - ice density
                  calculation_mass,     # - Type: Boolean - If true calculation mass step is activated
                  ref_standardization_nonstationarity_error, # - Type: Geotiff - reference DEM used to calculate standardization non stationnary error
                  ddem_standardization_nonstationarity_error, # - Type: Geotiff - ddem for calculating the error
                  unstandardized_nonstationarity_error,  # - Type: Boolean - If true unstandardized nonstationnary error is activated
                  standardization_nonstationarity_error, # - Type: Boolean - If true standardized nonstationnary error is activated
                  standardized_average_shapefile_error,  # - Type: Boolean - If true standardized average error is calculated
                  glacier_outlines_error,  # - Type: Shapefile - Shapefile of the glacier where you want to know the error
                  ddem_volume,          # - Type: Geotiff - dDEM tiff previously calculated by the user 
                  use_ddem_sub,         # - Type: Boolean - if true, calculation mass and volume will be calculated with the dDEM calulated in the subtraction step
                  use_ddem_inter,       # - Type: Boolean - if true, calculation mass and volume will be calculated with the dDEM calulated in the interpolation step
                  error_rho_density,                 # error of the density
                  error_glacier_outline,
                  name_pdf,             # - Type: String - Name of the pdf
                  #End entries
                  # =============================================================================
                  ): 
    # %% Loader
    # layout the form
    
    number_main_process=8 #intialization number of process - generate pdf and initialisation
    for i in (coregistration,subtraction,interpolation,calculation_mass,unstandardized_nonstationarity_error,standardization_nonstationarity_error,standardized_average_shapefile_error):
        if i:
            number_main_process+=4
            
    layout = [[sg.Text('Processing...',justification="center")],
              [sg.Text('', size=(50, 1),key='-text_main_process-',justification="center")],
              [sg.Text('', size=(10, 1),key='-percentage-',justification="center")],
              [sg.ProgressBar(max_value=number_main_process, orientation='h', size=(20, 4), key='-main_process-')],
              [sg.Cancel()]]
    
    # create the form`
    # must finalize since not running an event loop
    window = sg.Window('Eagloo : processing', layout, finalize=True, element_justification='center',)

    # Get the element to make updating easier
    progress_bar = window['-main_process-']
    percentage=window['-percentage-']
    text_progress_bar_main=window['-text_main_process-']
    
    
    progress_bar_value=0
    
    #Beginning update loader
    def loader_in_progress(progress_bar_value):
        progress_bar_value+=1
        progress_bar.update_bar(progress_bar_value) 
        percentage_number=(progress_bar_value/number_main_process)*100
        if percentage_number>100:
            percentage_number=100
        print(percentage_number)
        percentage.update(" {:.0f} % ".format(percentage_number))
        return progress_bar_value
        
        event, values = window.Read(timeout=0)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            window.close()
            return "The function has been stopped"
        
    #end update loader
    
    
    # %% Preprocessing
    
    text_progress_bar_main.update('Preprocessing of your DEMs')
    print("The process of your DEMs is launched")
    #Intialization function
    figure_plot=[None,None, None, None,None] #list return at the end of the algorithm: figure computed during matplotlib process
    
    
    if ddem_volume:
        ddem_volume_out=os.path.basename(ddem_volume).split('.')[0] # remove path and extension
    
    #Creating temporary folder
    directory = "temporary"
    os.makedirs(os.path.join(path_out, directory) , exist_ok="False") 
    path_temporary=os.path.join(path_out, directory)
    
    progress_bar_value=loader_in_progress(progress_bar_value)
    if vector_glacier:
        vector_glacier_out=os.path.basename(vector_glacier).split('.')[0]
        # Beginning homogenization coordinate glacier shapefile
        shapefile_glacier = geopandas.read_file(vector_glacier)
        shapefile_glacier= shapefile_glacier.to_crs(coordinate_system)
        vector_glacier_preprocessed="{}/preprocessed_{}.shp".format(path_temporary,vector_glacier_out)
        shapefile_glacier.to_file(vector_glacier_preprocessed)
        # End homogenization coordinate glacier shapefile
    
    progress_bar_value=loader_in_progress(progress_bar_value)
    
    # Beginning homogenization coordinate frame shapefile
    vector_frame_out=os.path.basename(vector_frame).split('.')[0]  # remove path and extension
    shapefile_frame = geopandas.read_file(vector_frame)
    shapefile_frame= shapefile_frame.to_crs(coordinate_system)
    vector_frame_preprocessed="{}/preprocessed_{}.shp".format(path_temporary,vector_frame_out)
    shapefile_frame.to_file(vector_frame_preprocessed)
    # End homogenization coordinate glacier shapefile
    
    progress_bar_value=loader_in_progress(progress_bar_value)
    
    outlines_glacier = gu.Vector("{}".format(vector_glacier_preprocessed))
    print(f"The glacier shapefile and the frame shapefile were reprojected to - {coordinate_system}")
    
    if dem_1: 
        dem1_out=os.path.basename(dem_1).split('.')[0] 
        dem_1_name="{}/preprocessing_{}.tif".format(path_temporary,dem1_out)
        gdal.Warp(dem_1_name ,dem_1 , cutlineDSName=vector_frame_preprocessed, dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
        dem1 = xdem.DEM("{}".format(dem_1_name))                #opening dem_1_name
        print(f"The earliest DEM was homogenized, resized and the new resolution is {resolution} m")
        
    if dem_2:
        dem2_out=os.path.basename(dem_2).split('.')[0]
        dem_2_name="{}/preprocessing_{}.tif".format(path_temporary,dem2_out)
        gdal.Warp(dem_2_name ,dem_2 , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
        dem2 = xdem.DEM("{}".format(dem_2_name))
        print(f"The latest DEM was homogenized, resized and the new resolution is {resolution} m")

    progress_bar_value=loader_in_progress(progress_bar_value)

 # %% Coregistration
    # =============================================================================
    #Beginning Coregistration ICP + Nuth And Kaab
    # =============================================================================
    if coregistration : 
        # %%% Pretreatment coregistration
        
        
        text_progress_bar_main.update('Co-registration of the earliest DEM')
        try: 
            if dem_2:
                if dem_coreg_ref:
                    if dem_2!=dem_coreg_ref:
                        dem_coreg_ref_out=os.path.basename(dem_coreg_ref).split('.')[0]
                        dem_ref_preprocessing="{}/preprocessing_{}.tif".format(path_temporary,dem_coreg_ref_out)
                        gdal.Warp(dem_ref_preprocessing ,dem_coreg_ref , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)  # !!!warning mask is not apply currently
                        print(f"The reference DEM for coregistration was homogenized, resized and the new resolution is {resolution} m")
                        
                        demref_coreg = xdem.DEM("{}".format(dem_ref_preprocessing))
                    else:
                        demref_coreg=dem2
                        dem_coreg_ref_out=dem2_out
                        dem_ref_preprocessing=dem_2_name
                else:
                    demref_coreg=dem2
                    dem_coreg_ref_out=dem2_out
                    dem_ref_preprocessing=dem_2_name
            else:
                if dem_coreg_ref:
                    dem_coreg_ref_out=os.path.basename(dem_coreg_ref).split('.')[0]
                    dem_ref_preprocessing="{}/preprocessing_{}.tif".format(path_temporary,dem_coreg_ref_out)
                    gdal.Warp(dem_ref_preprocessing ,dem_coreg_ref , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)  # !!!warning mask is not apply currently
                    print(f"The reference DEM for coregistration was homogenized, resized and the new resolution is {resolution} m")
                    
                    demref_coreg = xdem.DEM("{}".format(dem_ref_preprocessing))
        except Exception:
            return 'No DEM reference for the co-registration'
        
                
        inlier_mask = ~outlines_glacier.create_mask(demref_coreg)
        
        #Mask co-registration file with the glacier outlines shapefile
        with fiona.open(vector_glacier_preprocessed, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]
    
        with rasterio.open(dem_ref_preprocessing) as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, invert=True)
            out_meta = src.meta
    
        out_meta.update({"driver": "GTiff",
                          "height": out_image.shape[1],
                          "width": out_image.shape[2],
                          "transform": out_transform})
        coregistration_ref_masked="{}/reference_coregistration_masked_{}.tif".format(path_temporary,dem_coreg_ref_out)
        with rasterio.open(coregistration_ref_masked, "w", **out_meta) as dest:
            dest.write(out_image)
        coreg_ref_masked=xdem.DEM(coregistration_ref_masked)
        print("The reference DEM for coregistration was masked with the glacier outlines shapefiles")
        
        _ = dem1.reproject(demref_coreg) # reproject dem1 on dem2
        
        diff_before_dem1=coreg_ref_masked.data - dem1.data #calculate the difference between the reference DEM and the dem1
        transform_dem_1 = coreg_ref_masked.transform
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        # %%% Select method for coregistration
        
        # select which method will be used for the co-registration
        method=None
        if method_coregistration=="ICP":
            method= coreg_xdem.ICP()
             
        
        elif method_coregistration=="Nuth and Kääb":
            method= coreg_xdem.NuthKaab()
            
        
        elif method_coregistration=="ICP + Nuth and Kääb":
            method= coreg_xdem.ICP() + coreg_xdem.NuthKaab()

            
        elif method_coregistration=="Nuth and Kääb + Deramp":
            method=coreg_xdem.NuthKaab() + coreg_xdem.Deramp(degree=1)
        
        elif method_coregistration=="Deramp + Nuth and Kääb":
            method= coreg_xdem.Deramp(degree=1)+coreg_xdem.NuthKaab()
            
        
        elif method_coregistration=="Vertical offset + ICP + Nuth and Kääb":
             method=coreg_xdem.BiasCorr()+ coreg_xdem.ICP() + coreg_xdem.NuthKaab()
        
        elif method_coregistration=="Deramp(0) + ICP + Nuth and Kääb":
             method=coreg_xdem.Deramp(degree=0)+ coreg_xdem.ICP() + coreg_xdem.NuthKaab()
             
             
        elif method_coregistration=="ICP + Nuth and Kääb + Deramp":
            method= coreg_xdem.ICP() + coreg_xdem.NuthKaab() + coreg_xdem.Deramp(degree=1)
            
        elif method_coregistration=="Deramp + ICP + Nuth and Kääb":
            method =  coreg_xdem.Deramp(degree=1) + coreg_xdem.ICP()  + coreg_xdem.NuthKaab()   
        
        elif method_coregistration=="ICP + Deramp + Nuth and Kääb":
            method= coreg_xdem.ICP() + coreg_xdem.Deramp(degree=1) + coreg_xdem.NuthKaab()
            
             
        assert method!= None , " method_coregistration variable should be set out by choosing following entries : ICP,Nuth and Kaab,ICP + Nuth and Kaab "
        
        print(f"The method for co-registration is :{method_coregistration}")
        
        method.fit(reference_dem=coreg_ref_masked.data, dem_to_be_aligned=dem1.data,inlier_mask=inlier_mask,transform=transform_dem_1,verbose=True)
        dem1_coreg = method.apply(dem=dem1, transform=transform_dem_1) #Launch the co-registration
        offsets_1=method.to_matrix()
        
        teta_x_rotation_1= np.arctan2(offsets_1[2,1],offsets_1[2,2])*(180/np.pi)    #calculate rotation from transformation matrix
        teta_y_rotation_1= np.arctan2(-offsets_1[2,0],np.sqrt(np.square(offsets_1[2,1])+np.square(offsets_1[2,2])))*(180/np.pi)
        teta_z_rotation_1= np.arctan2(offsets_1[1,0],offsets_1[0,0])*(180/np.pi)
        
        if method_coregistration=="ICP":
            centroid_1=method._meta["centroid"]
            residual_1=method._meta["residual"]
        
        
        #Save the file co-registered
        dem1_coreg_path="{}/coregistration_{}_{}_on_{}.tif".format(path_out,method_coregistration,dem1_out,dem_coreg_ref_out)
        dem1_coreg.save(dem1_coreg_path)
    
        diff_after_dem1 = coreg_ref_masked.data - dem1_coreg.data #calculate the difference after the reference DEM and the dem1
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        # %%%  display the graph of the difference before and after coregistration
        
        #Calculate the maximum and minimum of the two graphs
        bounds_min = np.linspace(np.min(diff_after_dem1.squeeze()), 0, 129)
        bounds_max = np.linspace(0, np.max(diff_after_dem1.squeeze()), 129)[1:]
            # the zero is only needed once
            # in total there will be 257 bounds, so 256 bins
        bounds = np.concatenate((bounds_min, bounds_max), axis=None)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        
        plt_extent = [
           demref_coreg.bounds.left,
           demref_coreg.bounds.right,
           demref_coreg.bounds.bottom,
           demref_coreg.bounds.top,
       ] 
        
        figure_plot[0]=plt.figure(num=0,figsize=(9, 6))
       
        plt.subplot2grid((1, 15), (0, 0), colspan=7)
        plt.title("Elevation difference \n before coregistration ")
        plt.imshow(diff_before_dem1.squeeze(), cmap="RdYlBu", norm=norm,extent=plt_extent)
        plt.subplot2grid((1, 15), (0, 7), colspan=7)
        plt.title("Elevation difference \n after coregistration ")
        img = plt.imshow(diff_after_dem1.squeeze(), cmap="RdYlBu", norm=norm,extent=plt_extent )
        plt.subplot2grid((1, 15), (0, 14), colspan=1)
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        
        plt.suptitle('Co-registration, method : {} \n Co-registered DEM :{} \n Reference DEM :{}'.format(method_coregistration,dem1_out,dem_coreg_ref_out))
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.axis("off")
        coreg_1_image="{}/coregistration_{}_{}_ref_{}.png".format(path_out,method_coregistration,dem1_out,dem_coreg_ref_out)
        plt.savefig(coreg_1_image)
        plt.show()
        
        #calculate NMAD before and after coregistration
        NMAD_before_DEM1=xdem.spatialstats.nmad(diff_before_dem1)
        NMAD_after_DEM1=xdem.spatialstats.nmad(diff_after_dem1)
        print(f"Error before coregistration: {NMAD_before_DEM1:.3f} m")
        print(f"Error after coregistration: {NMAD_after_DEM1:.2f} m")

        print ( " \n {} has been co-registered - {} is the reference  \n" .format(dem1_out,dem_coreg_ref_out) )
        
        #Beginning update loader
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        event, values = window.Read(timeout=0)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            window.close()
            return "The function has been stopped"
            
            
        #end update loader
        
        # %%%  Coregistration of the seconf DEM
        if dem_2 == dem_coreg_ref: 
            dem2_coreg=dem2
        
        elif dem_2 and dem_coreg_ref and dem_2 != dem_coreg_ref:   
            # %%%% Prtocessing
            text_progress_bar_main.update('Co-registration of the latest DEM')
            
            _ = dem2.reproject(demref_coreg)
            print("dem_2 is different of dem_coreg_ref, thus dem_2 is co-registering and dem_coreg_ref is the reference")
            diff_before_dem2 = coreg_ref_masked.data - dem2.data
            transform_dem_2 = demref_coreg.transform
            
            method.fit(reference_dem=coreg_ref_masked.data, dem_to_be_aligned=dem2.data,inlier_mask=inlier_mask,transform=transform_dem_2,verbose=True,)
    
            dem2_coreg = method.apply(dem=dem2, transform=transform_dem_2)
            offsets_2=method.to_matrix()
            teta_x_rotation_2= np.arctan2(offsets_2[2,1],offsets_2[2,2]) *(180/np.pi)   #calculate rotation from transformation matrix
            teta_y_rotation_2= np.arctan2(-offsets_2[2,0],np.sqrt(np.square(offsets_2[2,1])+np.square(offsets_2[2,2])))*(180/np.pi)
            teta_z_rotation_2= np.arctan2(offsets_2[1,0],offsets_2[0,0])*(180/np.pi)
            
            if method_coregistration=="ICP" :
                centroid_2=method._meta["centroid"]
                residual_2=method._meta["residual"]
            
            dem2_coreg_path="{}/coregistration_{}_{}_on_{}.tif".format(path_out,method_coregistration,dem2_out,dem_coreg_ref_out)
            dem2_coreg.save(dem2_coreg_path)
    
            diff_after_dem2 = coreg_ref_masked.data - dem2_coreg.data
            # %%%%  display the graph of the difference before and after
           
            bounds_min = np.linspace(np.min(diff_after_dem2.squeeze()), 0, 129)
            bounds_max = np.linspace(0, np.max(diff_after_dem2.squeeze()), 129)[1:]
                # the zero is only needed once
                # in total there will be 257 bounds, so 256 bins
            bounds = np.concatenate((bounds_min, bounds_max), axis=None)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            
            figure_plot[1]=plt.figure(num=1,figsize=(9, 6))
            
            plt.subplot2grid((4, 15), (1, 0), colspan=7)
            plt.title("Elevation difference \nbefore coregistration ")
            plt.imshow(diff_before_dem2.squeeze(), cmap="RdYlBu",norm=norm,extent=plt_extent)
            plt.subplot2grid((4, 15), (1, 7), colspan=7)
            plt.title("Elevation difference \nafter coregistration ")
            img = plt.imshow(diff_after_dem2.squeeze(), cmap="RdYlBu", norm=norm,extent=plt_extent)
            plt.subplot2grid((4, 15), (1, 14), colspan=1)
            cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
            cbar.set_label("Elevation change (m)")
            plt.suptitle('Co-registration, method : {} :\n Co-registered  DEM:{} \n Reference DEM :{}'.format(method_coregistration,dem2_out,dem_coreg_ref_out))
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.axis("off")
            coreg_2_image="{}/coregistration_{}_{}_ref_{}.png".format(path_out,method_coregistration,dem2_out,dem_coreg_ref_out)
            plt.savefig(coreg_2_image)
            plt.show()
           
            print ( "\n {} has been coregistered - {} is the reference \n" .format(dem2_out,dem_coreg_ref_out))
            
            #calculate NMAD before and after coregistration
            NMAD_before_DEM2=xdem.spatialstats.nmad(diff_before_dem2)
            NMAD_after_DEM2=xdem.spatialstats.nmad(diff_after_dem2)
            print(f"Error before coregistration: {NMAD_before_DEM2:.3f} m")
            print(f"Error after coregistration: {NMAD_after_DEM2:.2f} m")
            
        progress_bar_value=loader_in_progress(progress_bar_value)
    #End Coregistration ICP + Nuth And Kaab
    
    # %% Subtraction
    # =============================================================================
    # Beginning subtraction between the 2 DEMs
    # =============================================================================
    if subtraction:
        # %%% Processing
        #if you want to use the DEMs coregistered previously
        text_progress_bar_main.update('Subtraction : Earliest DEM - Latest DEM')
        if use_dem_coreg_sub and coregistration:  
            _ = dem1_coreg.reproject(dem2_coreg)
            ddem_subtraction= dem2_coreg-dem1_coreg
        else:
            ddem_subtraction= dem2-dem1
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        name_ddem="{}/subtraction_{}_{}.tif".format(path_temporary,dem1_out,dem2_out)
        ddem_subtraction.save(name_ddem)
        
        #clip ddem with glacier outlines shapefile
        with fiona.open(vector_glacier_preprocessed, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]
    
        with rasterio.open(name_ddem) as src:
            ddem_out_image, out_transform = rasterio.mask.mask(src, shapes, filled=True)
            ddem_out_meta = src.meta
            
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        ddem_out_meta.update({"driver": "GTiff",
                         "height": ddem_out_image.shape[1],
                         "width": ddem_out_image.shape[2],
                         "transform": out_transform})
        
        
        ddem_subtraction_name="DDEM_subtraction_clipped_{}_{}.tif".format(dem1_out,dem2_out)
        ddem_clipped_name="{}/{}".format(path_out,ddem_subtraction_name)
        with rasterio.open(ddem_clipped_name, "w" ,**ddem_out_meta) as dest:
            dest.write(ddem_out_image)
        ddem_clipped_out=xdem.DEM(ddem_clipped_name)
        #end clipping
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        # %%% Display the graph
        bounds_min = np.linspace(np.min(ddem_clipped_out.data.squeeze()), 0, 129)
        bounds_max = np.linspace(0, np.max(ddem_clipped_out.data.squeeze()), 129)[1:]
            # the zero is only needed once
            # in total there will be 257 bounds, so 256 bins
        bounds = np.concatenate((bounds_min, bounds_max), axis=None)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        
        plt_extent = [
            ddem_clipped_out.bounds.left,
            ddem_clipped_out.bounds.right,
            ddem_clipped_out.bounds.bottom,
            ddem_clipped_out.bounds.top,
        ]  
        
        figure_plot[2]=plt.figure(num=2,figsize=(8, 5))
        plt.suptitle("Subtraction : {} - {} " .format(dem2_out,dem1_out))
        img = plt.imshow(ddem_clipped_out.data.squeeze(), cmap="RdYlBu", norm=norm, extent=plt_extent)
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        subtraction_image="{}/DDEM_subtraction_clipped_{}_{}.png".format(path_out,dem1_out,dem2_out)
        plt.savefig(subtraction_image)
        plt.show()   
        
        print ("\n {} - {} has been processed \n" .format(dem2_out,dem1_out))
        
        #Beginning update loader
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        event, values = window.Read(timeout=0)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            window.close()
            return "The function has been stopped"
            
        #end update loader
        
        
    
    # %% Interpolation
    
    if interpolation:
        # %%% Preprocessing
        text_progress_bar_main.update(f'Interpolation : {method_coregistration}')
        
        if dem_inter_ref:
            dem_inter_ref_out=os.path.basename(dem_inter_ref).split('.')[0]
            #beginning preprocessing
            dem_ref_hypso_name_out="{}/preprocessing_{}.tif".format(path_temporary,dem_inter_ref_out)
            gdal.Warp(dem_ref_hypso_name_out ,dem_inter_ref , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
            demref_hypso = xdem.DEM("{}".format(dem_ref_hypso_name_out), datetime=datetime(2013,8,1))
        else:
            dem_ref_hypso_name_out=dem2_out
            demref_hypso=dem2
            
        src = gdal.Open("{}".format(dem_ref_hypso_name_out))                #open file with gdal - type:gdal.Dataset - this layer us used to get information on rasters (resolution, position of the first pixel)
        gt=src.GetGeoTransform()
        glacier_index_map = outlines_glacier.rasterize(demref_hypso)
        
        inlier_mask = outlines_glacier.create_mask(demref_hypso)
        #end preprocessing
        outlines = {
            datetime(1959, 8, 1): outlines_glacier,
        }
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        if ddem_calculated:
            ddemc_out=os.path.basename(ddem_calculated).split('.')[0] 
            #Preprocessing of the ddem_calculated resolution, resize frame and coordinate system
            ddemc_name="{}/preprocessing_{}.tif".format(path_temporary,ddemc_out)
            gdal.Warp(ddemc_name ,ddem_calculated , cutlineDSName=vector_frame_preprocessed, dstSRS=coordinate_system,  xRes=resolution, yRes=resolution, cropToCutline=True)
            ddemc = xdem.DEM("{}".format(ddemc_name))                #opening dem_1_name
            
            ddem_unclipped= ddemc
            
            name_ddem="{}/subtraction_{}.tif".format(path_temporary,ddemc_out)
            ddem_unclipped.save(name_ddem)
        
        else:
            _ = dem1.reproject(dem2)
            #if you want to use the DEMs coregistered previously
            if use_dem_coreg_inter and coregistration:   
                ddem_unclipped= dem2_coreg-dem1_coreg
            #else your own dem will be used
            else:
                ddem_unclipped= dem2-dem1
            name_ddem="{}/subtraction_{}_{}.tif".format(path_temporary,dem1_out,dem2_out)
            ddem_unclipped.save(name_ddem)
       
        progress_bar_value=loader_in_progress(progress_bar_value)  
      #------------------------------------------------------------------------
      #Beginning clipping
        with fiona.open(vector_glacier_preprocessed, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]
            
        with rasterio.open(name_ddem) as src:
            ddem_out_image, out_transform = rasterio.mask.mask(src, shapes, filled=True)
            ddem_out_meta = src.meta
    
        ddem_out_meta.update({"driver": "GTiff",
                         "height": ddem_out_image.shape[1],
                         "width": ddem_out_image.shape[2],
                         "transform": out_transform})
        
        
        
        if ddem_calculated:
            ddem_interpolation_name="DDEM_subtraction_clipped_{}.tif".format(ddemc_out)
            ddem_clipped_name="{}/DDEM_subtraction_clipped_{}.tif".format(path_temporary,ddemc_out)
            
        else:
            ddem_interpolation_name="DDEM_subtraction_clipped_{}_{}.tif".format(dem1_out,dem2_out)
            ddem_clipped_name="{}/DDEM_subtraction_clipped_{}_{}.tif".format(path_temporary,dem1_out,dem2_out)
        
        with rasterio.open(ddem_clipped_name, "w",**ddem_out_meta) as dest:
            dest.write(ddem_out_image)
            
        #End clipping
        #------------------------------------------------------------------------
        # %%% Processing
        ddem_clipped_out=xdem.DEM(ddem_clipped_name)
        ddem_clipped=xdem.DEM(ddem_clipped_name)
        
        
        ddem_1 =xdem.dDEM(
            ddem_clipped,
            start_time=np.datetime64("1959-08-01"),
            end_time=np.datetime64("2013-08-01")
        )
        ddem_1.data /= (2013-1959)
        
        #choose interpolation method  
        if method_interpolation=="linear interpolation":
            ddem_filled=xdem.volume.linear_interpolation(ddem_clipped.data, max_search_distance=20, extrapolate=False, force_fill=True)
            nodata_value=0.0
        elif method_interpolation=="hypsometric interpolation":
            ddem_filled=xdem.volume.hypsometric_interpolation(voided_ddem=ddem_clipped.data, ref_dem=demref_hypso.data,mask=inlier_mask)
            nodata_value="nan"
        elif method_interpolation=="local hypsometric interpolation":
            ddem_filled=xdem.volume.local_hypsometric_interpolation(voided_ddem=ddem_clipped.data, ref_dem=demref_hypso.data, mask=inlier_mask, min_coverage=0.5,nodata=1, plot=False)
            nodata_value=1
        
        ddem_filled_correction=ddem_filled.squeeze()/0.01851849713  
        
        # ddem_filled_correction=xdem.volume.local_hypsometric_interpolation(voided_ddem=ddem_1.data ,ref_dem=demref_hypso.data, mask=inlier_mask)
        
        pixelsizex= gt[1]
        pixelsizey=-gt[5]
        coordinatex=gt[0]
        coordinatey=gt[3]
        
        transform_interpolation = from_origin(coordinatex, coordinatey, pixelsizex, pixelsizey)
        
        if ddem_calculated:
            interpolation_name="{}/interpolation_{}_{}.tif".format(path_out,method_interpolation,ddemc_out)
        else:
            interpolation_name="{}/interpolation_{}_{}_{}.tif".format(path_out,method_interpolation,dem1_out,dem2_out)
        new_dataset= rasterio.open(interpolation_name, 'w',height=ddem_filled_correction.data.shape[0], width=ddem_filled_correction.data.shape[1], driver='GTiff', count=1, dtype=ddem_filled_correction.dtype, crs=coordinate_system, transform=transform_interpolation, nodata=nodata_value) 
        new_dataset.write(ddem_filled_correction,1)
        new_dataset.close()
           
        progress_bar_value=loader_in_progress(progress_bar_value)   
        
        print(gt)
        print("nterpolation has been processed.")
        
        # %%% display the graph between and before interpolation
        
        bounds_min = np.linspace(np.min(ddem_filled_correction.squeeze()), 0, 129)
        bounds_max = np.linspace(0, np.max(ddem_filled_correction.squeeze()), 129)[1:]
            # the zero is only needed once
            # in total there will be 257 bounds, so 256 bins
        bounds = np.concatenate((bounds_min, bounds_max), axis=None)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
           
        plt_extent = [
            ddem_clipped_out.bounds.left,
            ddem_clipped_out.bounds.right,
            ddem_clipped_out.bounds.bottom,
            ddem_clipped_out.bounds.top,
        ]         
    
        figure_plot[3]=plt.figure(num=3,figsize=(8, 5))
        if ddem_calculated:
            plt.suptitle('Interpolation, method: {} \n dDEM = {} ;\n Reference DEM: {}'.format(method_interpolation,ddemc_out,dem_inter_ref_out))
            
        else:
            plt.suptitle('Interpolation, method: {} \n dDEM = {}-{} ;\n Reference DEM: {}'.format(method_interpolation,dem2_out,dem1_out,dem_inter_ref_out))
        plt.subplot2grid((1, 15), (0, 0), colspan=7)
        plt.title("Before interpolation ")
        plt.imshow(ddem_clipped_out.data.squeeze(), cmap="RdYlBu", norm=norm, extent=plt_extent)
        plt.subplot2grid((1, 15), (0, 7), colspan=7)
        plt.title("After interpolation")
        img = plt.imshow(ddem_filled_correction.squeeze(), cmap="RdYlBu", norm=norm, extent=plt_extent)
        plt.subplot2grid((1, 15), (0, 14), colspan=1)
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if ddem_calculated:
            interpolation_image="{}/interpolation_{}_ddem_{}_ref_{}.png".format(path_out,method_interpolation,ddemc_out,dem_inter_ref_out)
            
        else:
            interpolation_image="{}/interpolation_{}_ddem_={}-{}_ref_{}.png".format(path_out,method_interpolation,dem1_out,dem2_out,dem_inter_ref_out)
        plt.savefig(interpolation_image)
        plt.show()
        
    
        # graph of the interpolation
        mask = outlines_glacier.create_mask(ddem_1)
        
        ddem_bins = xdem.volume.hypsometric_binning(ddem_1.data[mask], demref_hypso.data[mask])/0.01851849713 
        stds = xdem.volume.hypsometric_binning(ddem_1.data[mask], demref_hypso.data[mask], aggregation_function=np.std)
        
        figure_plot[4]=plt.figure(num=4,figsize=(8, 8))
        plt.grid(zorder=0)
        plt.plot(ddem_bins["value"], ddem_bins.index.mid, linestyle="--", zorder=1)
        plt.barh(
            y=ddem_bins.index.mid,
            width=stds["value"],
            left=ddem_bins["value"] - stds["value"] / 2,
            height=(ddem_bins.index.left - ddem_bins.index.right) * 1,
            zorder=2,
            edgecolor="black",
        )
        
        for bin in ddem_bins.index:
            plt.vlines(ddem_bins.loc[bin, "value"], bin.left, bin.right, color="black", zorder=3)
        plt.xlabel("Elevation change (m / a)")
        plt.twiny()
        plt.barh(
            y=ddem_bins.index.mid,
            width=ddem_bins["count"] / ddem_bins["count"].sum(),
            left=0,
            height=(ddem_bins.index.left - ddem_bins.index.right) * 1,
            zorder=2,
            alpha=0.2,
        )
        plt.xlabel("Normalized area distribution (hypsometry)") 
        plt.ylabel("Elevation (m a.s.l.)")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if ddem_calculated:
            graph_hypso_image="{}/graph_interpolation_{}_ddem_{}_ref_{}.png".format(path_out,method_interpolation,ddemc_out,dem_inter_ref_out)
            
        else:
            graph_hypso_image="{}/graph_interpolation{}_{}_{}_ref_{}.png".format(path_out,method_interpolation,dem1_out,dem2_out,dem_inter_ref_out)
        plt.savefig(graph_hypso_image)
        plt.show()
        
        #Beginning update loader
        progress_bar_value=loader_in_progress(progress_bar_value)
        
        event, values = window.Read(timeout=0)
        if event == sg.WIN_CLOSED or event == 'Cancel':
            window.close()
            return "The function has been stopped"
            
        #end update loader
        
    
#%% Standardization error
    uncertainty_of_the_volume_change=0 #initialization if errors is unwanted
    if standardization_nonstationarity_error or unstandardized_nonstationarity_error or standardized_average_shapefile_error:
        
        #%%% Preprocessing
        
        text_progress_bar_main.update('Errors calculation')
        # We start by estimating the non-stationarities and deriving a terrain-dependent measurement error as a function of both
        # slope and maximum curvature, as shown in the  :ref:`sphx_glr_auto_examples_plot_nonstationarity_error.py` example.
        
        # Load the data
        
        if ref_standardization_nonstationarity_error:
            #beginning preprocessing
            dem_ref_error_name_out=os.path.basename(ref_standardization_nonstationarity_error).split('.')[0]
            dem_ref_error__out="{}/preprocessing_{}.tif".format(path_temporary,dem_ref_error_name_out)
            gdal.Warp(dem_ref_error__out ,ref_standardization_nonstationarity_error , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
            #end preprocessing
            ref_dem_error = xdem.DEM(f"{dem_ref_error__out}")
        elif dem2:
            ref_dem_error=dem2
            dem_ref_error_name_out=dem2_out
        
        
        
        if ddem_standardization_nonstationarity_error:
            #beginning preprocessing
            ddem_error_out=dem_inter_ref_out=os.path.basename(ddem_standardization_nonstationarity_error).split('.')[0]
            ddem_error_name_out="{}/preprocessing_{}.tif".format(path_temporary,ddem_error_out)
            gdal.Warp(ddem_error_name_out ,ddem_standardization_nonstationarity_error , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
            #end preprocessing
            dh_error = xdem.DEM(f"{ddem_error_name_out}")
        
        else: 
            dh_error=dem2-dem1
            
        src = gdal.Open(f"{dem_1_name}")               #open file with gdal - type:gdal.Dataset - this layer us used to get information on rasters (resolution, position of the first pixel)
        gt=src.GetGeoTransform()
        
        mask_glacier = outlines_glacier.create_mask(dh_error)
        
        
        #%%% Calculate spatial error without standardization
        slope, planc, profc = \
                            xdem.terrain.get_terrain_attribute(dem=ref_dem_error.data,
                                               attribute=['slope', 'planform_curvature', 'profile_curvature'],
                                               resolution=ref_dem_error.res)
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        # Remove values on unstable terrain
        dh_arr = dh_error.data[~mask_glacier]
        slope_arr = slope[~mask_glacier]
        planc_arr = planc[~mask_glacier]
        profc_arr = profc[~mask_glacier]
        maxc_arr = np.maximum(np.abs(planc_arr),np.abs(profc_arr))
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        # Remove large outliers
        dh_arr[np.abs(dh_arr) > 4 * xdem.spatialstats.nmad(dh_arr)] = np.nan
        
        # Define bins for 2D binning
        custom_bin_slope = np.unique(np.concatenate([np.nanquantile(slope_arr,np.linspace(0,0.95,20)),
                                                     np.nanquantile(slope_arr,np.linspace(0.96,0.99,5)),
                                                     np.nanquantile(slope_arr,np.linspace(0.991,1,10))]))
        
        custom_bin_curvature = np.unique(np.concatenate([np.nanquantile(maxc_arr,np.linspace(0,0.95,20)),
                                                     np.nanquantile(maxc_arr,np.linspace(0.96,0.99,5)),
                                                     np.nanquantile(maxc_arr,np.linspace(0.991,1,10))]))
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        if unstandardized_nonstationarity_error:
            text_progress_bar_main.update('Unstdandardized non stationnarity errror calculation')
            
            df = xdem.spatialstats.nd_binning(values=dh_arr, list_var=[slope_arr, maxc_arr], list_var_names=['slope', 'maxc'],
                                      statistics=['count', np.nanmedian, xdem.spatialstats.nmad],
                                      list_var_bins=[custom_bin_slope,custom_bin_curvature])
            xdem.spatialstats.plot_2d_binning(df, 'slope', 'maxc', 'nmad', 'Slope (degrees)', 'Maximum absolute curvature (100 m$^{-1}$)', 'NMAD of dh (m)', scale_var_2='log', vmin=2, vmax=10)
            # Estimate an interpolant of the measurement error with slope and maximum curvature
            slope_curv_to_dh_err = xdem.spatialstats.interp_nd_binning(df, list_var_names=['slope', 'maxc'], statistic='nmad', min_count=10)
            maxc = np.maximum(np.abs(profc), np.abs(planc))
            
            unstandardized_NMAD_slope_maxc_out=f"{path_out}/unstandardized_NMAD_slope_maxc.png"
            plt.savefig(unstandardized_NMAD_slope_maxc_out)
            
            # Estimate a measurement error per pixel
            dh_err = slope_curv_to_dh_err((slope, maxc))
            
            plt.figure(figsize=(8, 5))
            
            plt.imshow(dh_err.squeeze(), cmap="Reds", vmin=0, vmax=30)
            cbar = plt.colorbar()
            cbar.set_label('Elevation measurement error (m)')
            unstandardized_elevation_difference_out=f"{path_out}/unstandardized_elevation_difference.png"
            plt.suptitle("Unstandardized spatial distribution \n of the elevation measurement error over the area \n based on the measurement error dependent on slope and maximum curvature \n\n " )
            plt.tight_layout()
            plt.savefig(unstandardized_elevation_difference_out)
            plt.show()
            
            #Beginning update loader
            progress_bar_value+=1
            progress_bar.update_bar(progress_bar_value) 
            
            event, values = window.Read(timeout=0)
            if event == sg.WIN_CLOSED or event == 'Cancel':
                window.close()
                return "The function has been stopped"
                
            #end update loader
            progress_bar_value=loader_in_progress(progress_bar_value)
            #%%% Calculate spatial error with standardization 
        if standardization_nonstationarity_error or standardized_average_shapefile_error:
            text_progress_bar_main.update('Standardized non stationnarity errors calculation')
            # Perform 2D binning to estimate the measurement error with slope and maximum curvature
            df = xdem.spatialstats.nd_binning(values=dh_arr, list_var=[slope_arr, maxc_arr], list_var_names=['slope', 'maxc'],
                                              statistics=['count', np.nanmedian, np.nanstd],
                                              list_var_bins=[custom_bin_slope,custom_bin_curvature])
            
            progress_bar_value=loader_in_progress(progress_bar_value)
            
            # Estimate an interpolant of the measurement error with slope and maximum curvature
            slope_curv_to_dh_err = xdem.spatialstats.interp_nd_binning(df, list_var_names=['slope', 'maxc'], statistic='nanstd', min_count=30)
            maxc = np.maximum(np.abs(profc), np.abs(planc))
            
            # Estimate a measurement error per pixel
            dh_err = slope_curv_to_dh_err((slope, maxc))
            
            progress_bar_value=loader_in_progress(progress_bar_value)
            #------------------------------------------------------------
            # Using the measurement error estimated for each pixel, we standardize the elevation differences by applying
            # a simple division:
            
            z_dh = dh_error.data/dh_err
            
            
            # We remove values on glacierized terrain and large outliers.
            z_dh.data[mask_glacier] = np.nan
            z_dh.data[np.abs(z_dh.data)>4] = np.nan
            
            
            
            #------------------------------------------------------------
            
            
            # We perform a scale-correction for the standardization, to ensure that the standard deviation of the data is exactly 1.
            standard_deviation_before_scale_correction=np.nanstd(z_dh.data)
            print('Standard deviation before scale-correction: {:.1f}'.format(standard_deviation_before_scale_correction))
            scale_fac_std = np.nanstd(z_dh.data)
            z_dh = z_dh/scale_fac_std
            standard_deviation_after_scale_correction=np.nanstd(z_dh.data)
            print('Standard deviation after scale-correction: {:.1f}'.format(standard_deviation_after_scale_correction))
            
            plt.figure(figsize=(8, 5))
            plt_extent = [
                ref_dem_error.bounds.left,
                ref_dem_error.bounds.right,
                ref_dem_error.bounds.bottom,
                ref_dem_error.bounds.top,
            ]
            ax = plt.gca()
            outlines_glacier.ds.plot(ax=ax, fc='none', ec='tab:gray')
            ax.plot([], [], color='tab:gray', label='Glacier outlines')
            plt.imshow(z_dh.squeeze(), cmap="RdYlBu", vmin=-3, vmax=3, extent=plt_extent)
            cbar = plt.colorbar()
            cbar.set_label('Standardized elevation differences (m)')
            plt.legend(loc='lower right')
            standardized_elevation_difference_clipped_out=f"{path_out}/standardized_elevation_difference_clipped_out.png"
            plt.tight_layout()
            plt.savefig(standardized_elevation_difference_clipped_out)
            plt.show()
            
            progress_bar_value=loader_in_progress(progress_bar_value)
            # Now, we can perform an analysis of spatial correlation as shown in the :ref:`sphx_glr_auto_examples_plot_vgm_error.py`
            # example, by estimating a variogram and fitting a sum of two models.
            df_vgm = xdem.spatialstats.sample_empirical_variogram(
                values=z_dh.data.squeeze(), gsd=dh_error.res[0], subsample=50, runs=30, n_variograms=10, random_state=42)
            
            fun, params = xdem.spatialstats.fit_sum_model_variogram(['Sph','Sph'], empirical_variogram=df_vgm)
            xdem.spatialstats.plot_vgm(df_vgm, xscale_range_split=[100, 1000, 10000], list_fit_fun=[fun],
                                       list_fit_fun_label=['Standardized double-range variogram'])
            variance_sample_count_standardization=f"{path_out}/variance_sample_count_standardization.png"
            plt.savefig(variance_sample_count_standardization)
            
            
            # With standardized input, the variogram should converge towards one. With the input data close to a stationarity
            # variance, the variogram will be more robust as it won't be affected by changes in variance due to terrain- or
            # instrument-dependent variability of measurement error. The variogram should only capture changes in variance due to
            # spatial correlation.
            
            #------------------------------------------------------------
            #calculate the error map with standardization
            
            #Beginning update loader
            progress_bar_value=loader_in_progress(progress_bar_value)
            
            event, values = window.Read(timeout=0)
            if event == sg.WIN_CLOSED or event == 'Cancel':
                window.close()
                return "The function has been stopped"
                
            #end update loader
            
            #%%% Calculate error with standardization on a glacier
           
            dh_z_err = z_dh.data[~mask_glacier]
            
            # Remove large outliers
            dh_z_err[np.abs(dh_z_err) > 4 * xdem.spatialstats.nmad(dh_z_err)] = np.nan
            
            
            
            df_slope_maxc = xdem.spatialstats.nd_binning(values=dh_z_err, list_var=[slope_arr, maxc_arr], list_var_names=['slope', 'maxc'],
                                      statistics=['count', np.nanmedian, xdem.spatialstats.nmad],
                                      list_var_bins=[custom_bin_slope,custom_bin_curvature])
            xdem.spatialstats.plot_2d_binning(df_slope_maxc, 'slope', 'maxc', 'nmad', 'Slope (degrees)', 'Maximum absolute curvature (100 m$^{-1}$)', 'NMAD of dh (m)', scale_var_2='log', vmin=0, vmax=1.5)
            curvature_slope_nmad_out=f"{path_out}/curvature_slope_nmad_out.png"
            plt.savefig(curvature_slope_nmad_out)
            
            
            
            # calculate the layer to know the exact error everywhere on a map - destandardize
            slope_curv_to_dh_err = xdem.spatialstats.interp_nd_binning(df_slope_maxc, list_var_names=['slope', 'maxc'], statistic='nmad', min_count=30)
            maxc = np.maximum(np.abs(profc), np.abs(planc))
            dh_err_2 = slope_curv_to_dh_err((slope, maxc))
            dh_err_final=np.abs(dh_err_2*dh_err*scale_fac_std)
                
            if standardization_nonstationarity_error:
                progress_bar_value=loader_in_progress(progress_bar_value)
                text_progress_bar_main.update('Standardized non stationnarity errors calculation')
                progress_bar_value=loader_in_progress(progress_bar_value)
                plt.figure(figsize=(8, 5))
                plt.suptitle("Destandardized spatial distribution \n of the elevation measurement error over the area \n based on the measurement error dependent on slope and maximum curvature \n\n " )
                plt.imshow(dh_err_final.squeeze(), cmap="Reds", vmin=0, vmax=30, extent=plt_extent)
                cbar = plt.colorbar()
                cbar.set_label('Elevation measurement error (m)')
                
                progress_bar_value=loader_in_progress(progress_bar_value)
                
                standardized_elevation_difference_out=f"{path_out}/standardized_elevation_difference_out.png"
                plt.tight_layout()
                plt.savefig(standardized_elevation_difference_out)
                plt.show()
                
                pixelsizex= gt[1]
                pixelsizey=-gt[5]
                coordinatex=gt[0]
                coordinatey=gt[3]
                
                transform_standardized_nonstationnarity_error = from_origin(coordinatex, coordinatey, pixelsizex, pixelsizey)
                
                standardized_spatial_error_name="{}/standardized_spatial_error_{}_{}.tif".format(path_out,dem1_out,dem_ref_error_name_out)
                new_dataset= rasterio.open(standardized_spatial_error_name, 'w',height=dh_err_final.shape[1], width=dh_err_final.shape[2], driver='GTiff', count=1, dtype=dh_error.data.squeeze().dtype, crs=coordinate_system, transform=transform_standardized_nonstationnarity_error, nodata="nan") 
                new_dataset.write(dh_err_final.squeeze(),1)
                new_dataset.close()
                
                
                progress_bar_value=loader_in_progress(progress_bar_value)
            #------------------------------------------------------------
            
            
            # **How to use this standardized spatial analysis to compute final uncertainties?**
            #
            # Let's take the example of two glaciers of similar size: local_glacierbreen and Medalsbreen, which are respectively
            # north and south-facing. The south-facing Medalsbreen glacier is subject to more sun exposure, and thus should be
            # located in higher slopes, with possibly higher curvatures.
            # Beginning homogenization coordinate glacier shapefile
            if standardized_average_shapefile_error:
                
                if glacier_outlines_error:
                    shapefile_glacier_errors = geopandas.read_file(glacier_outlines_error)
                    name_shapefile_glacier_errors_out=os.path.basename(glacier_outlines_error).split('.')[0]
                    shapefile_glacier_errors= shapefile_glacier_errors.to_crs(coordinate_system)
                else: #if glacier_outlines_error is empty, the glacier outline, the vector_glacier is used
                    shapefile_glacier_errors= shapefile_glacier
                    name_shapefile_glacier_errors_out=vector_glacier_out
                
                area_shapefile_glacier_errors=np.sum(shapefile_glacier_errors.area)#calculate the area - use to calculate the number of effective samples
                print(area_shapefile_glacier_errors)
                vector_glacier_errors_preprocessed="{}/preprocessed_glacier_errors_preprocessed_{}.shp".format(path_temporary,name_shapefile_glacier_errors_out)
                shapefile_glacier_errors.to_file(vector_glacier_errors_preprocessed)
                # End homogenization coordinate glacier shapefile
                
                local_glacier_shp = gu.Vector(vector_glacier_errors_preprocessed)
                local_glacier_mask = local_glacier_shp.create_mask(dh_error)
                
                progress_bar_value=loader_in_progress(progress_bar_value)
                
                plt.figure(figsize=(8, 5))
                ax = plt.gca()
                plt_extent = [
                    ref_dem_error.bounds.left,
                    ref_dem_error.bounds.right,
                    ref_dem_error.bounds.bottom,
                    ref_dem_error.bounds.top,
                ]
                plt.imshow(slope.squeeze(), cmap="Blues", vmin=0, vmax=40, extent=plt_extent)
                cbar = plt.colorbar(ax=ax)
                cbar.set_label('Slope (degrees)')
                local_glacier_shp.ds.plot(ax=ax, fc='none', ec='tab:olive', lw=2)
                plt.plot([], [], color='tab:gray', label=f'Glacier outlines,{name_shapefile_glacier_errors_out}')
                plt.legend(loc='lower left')
                slope_out=f"{path_out}/slope_out_{name_shapefile_glacier_errors_out}.png"
                plt.tight_layout()
                plt.savefig(slope_out)
                plt.show()
                
                progress_bar_value=loader_in_progress(progress_bar_value)
                
                average_slope_local_glacier=np.nanmean(slope[local_glacier_mask])
                print('Average slope of local glacier: {:.1f}'.format(average_slope_local_glacier))
                
               
                # We calculate the number of effective samples for each glacier based on the variogram
                local_glacier_neff = xdem.spatialstats.neff_circ(area_shapefile_glacier_errors,[(params[0], 'Sph', params[1]),
                                                                   (params[2], 'Sph', params[3])])
                
               
                print('Number of effective samples of the local glacier: {:.1f}'.format(local_glacier_neff))
                
                
                # Due to the long-range spatial correlations affecting the elevation differences, both glacier have a similar, low
                # number of effective samples. This transcribes into a large standardized integrated error.
                
                local_glacier_z_err = 1/np.sqrt(local_glacier_neff)
                
                progress_bar_value=loader_in_progress(progress_bar_value)
                
                print('Standardized integrated error of local glacier: {:.3f}'.format(local_glacier_z_err))
                
                
                # Finally, we destandardize the spatially integrated errors based on the measurement error dependent on slope and
                # maximum curvature. This yields the uncertainty into the mean elevation change for each glacier.
                
                # Destandardize the uncertainty
                local_glacier_dh_err_integrated = np.nanmean(dh_err_final[local_glacier_mask])*local_glacier_z_err*scale_fac_std
                
                # Derive mean elevation change
                local_glacier_dh = np.nanmean(dh_error.data[local_glacier_mask])
                
                #mean standardized error on the glacier
                local_glacier_dh_err_mean=np.nanmean(dh_err_final[local_glacier_mask])
                
                # Plot the result
                plt.figure(figsize=(8, 5))
                plt.suptitle(f"Destandardized spatial integrated errors \n based on the measurement error dependent on slope and maximum curvature \n glacier:{name_shapefile_glacier_errors_out}\n\n" )
                ax = plt.gca()
                plt.imshow(dh_error.data.squeeze(), cmap="RdYlBu", vmin=-50, vmax=50, extent=plt_extent)
                cbar = plt.colorbar(ax=ax)
                cbar.set_label('Elevation differences (m)')
                local_glacier_shp.ds.plot(ax=ax, fc='none', ec='tab:olive', lw=2)
                
                plt.plot([],[], color='tab:olive', label=f' Glacier selected : {name_shapefile_glacier_errors_out}')
                
                ax.text(local_glacier_shp.ds.centroid.x.values[0], local_glacier_shp.ds.centroid.y.values[0]-1500,
                        '{:.2f} \n$\\pm$ {:.2f}'.format(local_glacier_dh, local_glacier_dh_err_integrated), color='tab:grey', fontweight='bold',
                        va='top', ha='center', fontsize=12)
                
                plt.legend(loc='lower left')
                elevation_difference_out=f"{path_out}/elevation_difference_out_{name_shapefile_glacier_errors_out}.png"
                plt.tight_layout()
                plt.savefig(elevation_difference_out)
                plt.show()
                
                uncertainty_of_the_volume_change = local_glacier_dh_err_integrated*area_shapefile_glacier_errors
                
                # Because of slightly higher slopes and curvatures, the final uncertainty for Medalsbreen is larger by about 10%.
                # The differences between the mean terrain slope and curvatures of stable terrain and those of glaciers is quite limited
                # on Svalbard. In high moutain terrain, such as the Alps or Himalayas, the difference between stable terrain and glaciers,
                # and among glaciers, would be much larger.
                
                #Beginning update loader
                progress_bar_value=loader_in_progress(progress_bar_value)
                
                event, values = window.Read(timeout=0)
                if event == sg.WIN_CLOSED or event == 'Cancel':
                    window.close()
                    return "The function has been stopped"
                  
                #end update loader
                print(uncertainty_of_the_volume_change)
    # %% Calculation mass
    # =============================================================================
    # Beginning volume and mass statistics
    # =============================================================================
    if calculation_mass: 
        text_progress_bar_main.update('Mass and volume calculation')
        
        if ddem_volume:
            vol_stats('mass',ddem_volume,rho_density,path_temporary,'own',standardized_average_shapefile_error,uncertainty_of_the_volume_change,error_rho_density, error_glacier_outline)
            
        if interpolation and use_ddem_inter: 
            print("volume and mass stats calculated with the dDEM interpolated")
            vol_stats('mass',interpolation_name,rho_density,path_temporary,'interpolation', standardized_average_shapefile_error,uncertainty_of_the_volume_change,error_rho_density, error_glacier_outline)
            
        if subtraction and use_ddem_sub:
            print("volume and mass stats calculated with the dDEM subtacted")
            vol_stats('mass',ddem_clipped_name,rho_density,path_temporary,'subtraction',standardized_average_shapefile_error, uncertainty_of_the_volume_change,error_rho_density, error_glacier_outline)  
        
        progress_bar_value=loader_in_progress(progress_bar_value)
        progress_bar_value=loader_in_progress(progress_bar_value)
        progress_bar_value=loader_in_progress(progress_bar_value)
        
         
        
    # end volume and mass statistics
      
   #%% Generate PDF
    now=datetime.now()
    date=f'{now.strftime("%d")} {now.strftime("%B")} {now.strftime("%Y")}'  
    
    text_progress_bar_main.update('Generating PDF report')
    progress_bar_value=loader_in_progress(progress_bar_value)
    class PDF(FPDF):
        #%%% Functions
        def __init__(self):
            super().__init__()
            self.WIDTH = 170
            self.HEIGHT = 180
            self.WIDTH_PAGE=210       #A4 width
            self.position_graph_pdf="left"  #used to display graphs at right or at left
            self.position_y=0
            self.number_title=1  #number of the title
            self.number_subtitle=1 # number of the subtitle
            self.set_margins(left=15,top=15,right=15)
            
        def header(self):
            # Custom logo and positioning
            # Create an `assets` folder and put any wide and short image inside
            # Name the image `logo.png`
            self.image("Eagloo_frontend_functions/eagloo_name_dark.png",10, 10, 35)
            self.set_font('helvetica', 'B', 11)
            self.cell(self.WIDTH -45)
            self.cell(60, 2, f"Report : {name_pdf}", 0, align='R')
            self.ln(20)
            
        def footer(self):
            # Page numbers in the footer
            self.set_y(-15)
            self.set_font('Times', 'I', 8)
            self.set_text_color(128)
            self.cell(0, 10, 'Page ' + str(self.page_no()), 0, align= 'C')
        
        # Style of the PDF title
        def title_pdf(self,text):
            self.ln(5)
            self.set_font('helvetica', 'B', 15)
            self.cell(0,5,f"{text}",align='C')
            self.ln(5)
        
        # Style of the chapter titles
        def chapter_title(self, label):
            # helvetica12
            self.set_font('helvetica', 'B', 14)
            # Background color
            self.set_fill_color(243, 242, 242)
            self.ln(3)
            # Title
            self.multi_cell(0, 12, "{} ".format(label),0,'L',1)
            # Line break
            self.ln(3)
            
        def image_body(self, image1):
            # Determine how many plots there are per page and set positions
            # and margins accordingly
            self.ln(5)
            self.image(image1, (self.WIDTH_PAGE-self.WIDTH)/2, None,self.WIDTH)
            self.ln(5)
            
        def print_raster(self, raster,namegraph,position):
            #entries
            #raster1 :Gtiff
            #postion: "right" or "left"
            #output: raster print on the page
            graph_pdf=xdem.DEM(f"{raster}")
            plt_extent = [
                graph_pdf.bounds.left,
                graph_pdf.bounds.right,
                graph_pdf.bounds.bottom,
                graph_pdf.bounds.top,
            ]
            plt.suptitle(namegraph)
            img = plt.imshow(graph_pdf.data.squeeze(), cmap='turbo', vmin=np.min(graph_pdf.data.squeeze()), vmax=np.max(graph_pdf.data.squeeze()), extent=plt_extent)
            cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
            cbar.set_label("Elevation change (m)")
            path_graph="{}/{}.png".format(path_temporary,namegraph)
            plt.savefig(path_graph)
            plt.show()   
            
            if position=="left":
                self.position_graph_pdf="right"
                self.image(path_graph, 5, None,100)
                self.ratio_image= Image.open(path_graph).size[1]/Image.open(path_graph).size[0]
                self.position_y=self.get_y()-(self.ratio_image*100)
            else:
                self.image(path_graph,105,self.position_y,100)
                self.position_graph_pdf="left"
                self.ln(10)
                
                
                  
        #Style of the intermediary title    
        def text_title(self, title):
            self.ln(5)
            self.set_font('helvetica', '', 14)
            self.multi_cell(0,5,f"{title}" ,'L','L',markdown=True)
            self.ln(10)
            
        #Display a .doc file
        def text_body(self, path):
            # Read text file
            with open(path, 'rb') as fh:
                txt = fh.read().decode('utf-8')
            # Times 12
            self.set_font('Times', '', 12)
            # Output justified text
            self.multi_cell(0, 5, txt,markdown=True)
            # Line break
            self.ln(5)
            # Mention in italics
        
        #style of the paragraph text
        def text_description(self,text):
            self.ln(5)
            self.set_font('Times', '', 12)
            self.multi_cell(0,5,f" {text}",0,'L',markdown=True)
            self.ln(5) 
            
        
        #%%% Write page
        def print_page(self):
            # Generates the report
            self.add_page()
            
            
            
            self.title_pdf("Processing report ")
            self.title_pdf("Eagloo")
            if name_pdf:
                self.title_pdf(f"Name : {name_pdf}")
            self.title_pdf(f"{date}")
            self.ln(10)
            self.title_pdf("____________________________________________")
            self.ln(10)
            
            #determine the number of treatment(s) for the pdf title
            n_t=0
            if coregistration:
                n_t+=1
            if subtraction:
                n_t+=1
            if interpolation:
                n_t+=1
            if calculation_mass:
                n_t+=1
            
            # write the pdf title
            if n_t==1:
                self.title_pdf(" Method of processing :")
            
            if n_t>1:
                self.title_pdf(" Methods of processing :")
                
            if coregistration:
                self.title_pdf(f" Co-registration - method : {method_coregistration}")
                
            if subtraction:
                self.title_pdf("Subtraction")
                
            if interpolation:
                self.title_pdf(f"Interpolation - method : {method_interpolation}")
            
            if standardization_nonstationarity_error or unstandardized_nonstationarity_error or standardized_average_shapefile_error : 
                self.title_pdf("Non-stationarity measurment errors")
                
            if calculation_mass:
                self.title_pdf("Volume and mass calculation")
                
            self.ln(10)
            
            self.image("Eagloo_frontend_functions/eagloo_logo_drawing.png", (self.WIDTH_PAGE-40)/2, None,40)
            self.ln(10)
            self.image("Eagloo_frontend_functions/eagloo_name_dark.png", (self.WIDTH_PAGE-40)/2, None,40)
            self.add_page()
            
            #%%% user entries
            self.chapter_title(f" {self.number_title} - Entry names")
            self.number_title+=1
            if dem_1:
                self.text_description(f" - Earliest DEM : {dem1_out} \n")
                
            if dem_2:
                self.text_description(f" - Latest DEM : {dem2_out} \n")
                
            self.text_description(f" - Contour line shapefile of the biggest glacier(Use for masking) : {vector_glacier_out} \n")
            self.text_description(f" - Path where the files were saved : {path_out} ")
            self.text_description(f" - Coordinate system : {coordinate_system} ")
            self.text_description(f" - Shapefile frame of the study area : {vector_frame_out} ")
            self.text_description(f" - Resolution :  {resolution} m ") 
            
            self.ln(5)
            
            #%%% DEM entries (graph)
            self.chapter_title(f" {self.number_title} -  Digital Elevation Model entries")
            self.number_title+=1 #increment the number of the chapter
            if dem1:
                self.print_raster(dem_1,f" Earliest DEM: {dem1_out}",self.position_graph_pdf)
                
            if dem2:
                self.print_raster(dem_2,f" Latest DEM: {dem2_out}",self.position_graph_pdf)
                
            if coregistration:
                self.print_raster(dem_coreg_ref,f" Reference DEM for the co-registration: {dem_coreg_ref_out}",self.position_graph_pdf)
                
            if interpolation:
                self.print_raster(dem_inter_ref,f" Reference DEM for the interpolation: {dem_inter_ref_out}",self.position_graph_pdf)
                
            if ddem_volume:
                self.print_raster(ddem_volume,f" DDEM to calculate the volume and the mass: {ddem_volume_out}",self.position_graph_pdf)

            self.ln(5)
            self.ln(5) 
            
            #%%% Coregistration
            if coregistration:
                self.chapter_title(f" {self.number_title} - Co-registration of the DEM: {dem1_out}")
                self.number_title+=1 #increment the number of the chapter
                self.text_description(f" - Reference DEM for co-registration : {dem_coreg_ref_out} ")
                self.text_description(f" - Method : {method_coregistration} ")
                self.ln(5) 
                self.text_description(" Translation (sum) :")
                self.text_description(f" - East offset (dx) :  **{offsets_1[0,3] :.2f} m **")
                self.text_description(f" - North offset(dy):  **{offsets_1[1,3] :.2f} m**")
                self.text_description(f" - Vertical offset (dz) :  **{offsets_1[2,3] :.2f} m** ")
                self.ln(5) 
                if method_coregistration=="ICP" or method_coregistration=="ICP + Nuth and Kaab":
                    self.text_description(" ICP rotation :")
                    self.text_description(f" - Rotation around the east axis (x):  **{teta_x_rotation_1 :.2f} ° **") 
                    self.text_description(f" - Rotation around the north axis (y) :  **{teta_y_rotation_1 :.2f} °**")
                    self.text_description(f" - Rotation around the vertical axis (z) :  **{teta_z_rotation_1 :.2f} °**")
                    if method_coregistration=="ICP":
                        self.text_description(f" - Centroid - East coordinate :  **{centroid_1[0]:.2f} m; North coordinate : {centroid_1[1]:.2f} m**")
                        self.text_description(f" - Residual :  **{residual_1 :.2f} m")
                        
                    
                self.text_description(f" - Error before the co-registration(NMAD) :  **{NMAD_before_DEM1 :.2f} m**")
                self.text_description(f" - Error after the co-registration(NMAD) :  **{NMAD_after_DEM1 :.2f} m**")
                self.image_body(coreg_1_image)
                self.add_page()
                
                if dem_2 and dem2!=dem2_coreg:
                    self.chapter_title(f" {self.number_title} - Co-registration of the DEM: {dem2_out}")
                    self.number_title+=1
                    self.text_description(f" - Reference DEM for co-registration : {dem_coreg_ref_out} ")
                    self.text_description(f" - Method : {method_coregistration} ")
                    self.ln(5) 
                    self.text_description(" Translation (sum) :")
                    self.text_description(f" - East offset (dx) :  **{offsets_2[0,3]} m** ")
                    self.text_description(f" - North offset(dy):  **{offsets_2[1,3]} m **")
                    self.text_description(f" - Vertical offset (dz) :  **{offsets_1[2,3]} m **")
                    self.ln(5) 
                    if method_coregistration=="ICP" or method_coregistration=="ICP + Nuth and Kaab":
                        self.text_description(" ICP rotation :")
                        self.text_description(f" - Rotation around the east axis (x):  **{teta_x_rotation_2 :.2f} °** ") 
                        self.text_description(f" - Rotation around the north axis (y):  **{teta_y_rotation_2 :.2f} ° **") 
                        self.text_description(f" - Rotation around the vertical axis (z) :  **{teta_z_rotation_2 :.2f} ° **")
                        if method_coregistration=="ICP":
                            self.text_description(f" - Centroid - East coordinate :  **{centroid_2[0] :.2f} m; North coordinate : {centroid_2[1] :.2f} m**")
                            self.text_description(f" - Residual :  **{residual_2 :.2f}**")
                    self.ln(5)
                    self.text_description(f" - Error before the co-registration(NMAD) :  **{NMAD_before_DEM2 :.2f} m**")
                    self.text_description(f" - Error after the co-registration(NMAD) :  **{NMAD_after_DEM2 :.2f} m**")
                    self.image_body(coreg_2_image)
                    self.add_page()
                    
                    
            #%%% Subtraction
            if subtraction:
                self.chapter_title(f" {self.number_title} - Subtraction : {dem1_out} - {dem2_out}")
                self.number_title+=1
                self.image_body(subtraction_image)
                self.add_page()
            
            #%%% Interpolation
            if interpolation:
                self.chapter_title(f" {self.number_title} - Interpolation")
                self.number_title+=1
                if dem_1 and dem_2: 
                    self.text_description(f" dDEM = {dem2_out} - {dem1_out}")
                    
                if ddem_calculated:
                    self.text_description(f" dDEM = {ddemc_out}")
                    
                self.text_description(f" - Reference DEM for interpolation:{dem_inter_ref_out}")
                self.text_description(f" - Method : {method_interpolation}")   
                
                if use_dem_coreg_inter and coregistration: 
                    self.text_description(" The DEMs co-registered in the co-registration step were used.")
                    
                self.image_body(interpolation_image)
                self.image_body(graph_hypso_image)
                self.add_page()
                
                
            if standardization_nonstationarity_error or unstandardized_nonstationarity_error or standardized_average_shapefile_error : 
                self.chapter_title(f" {self.number_title} - Non-stationarity measurment errors")
                if unstandardized_nonstationarity_error:
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Unstandardized non-stationarity measurment errors - based on the slope and the maximum absolute curvature")
                    self.image_body(unstandardized_NMAD_slope_maxc_out)
                    self.image_body(unstandardized_elevation_difference_out)
                    self.number_subtitle+=1
                
                if standardization_nonstationarity_error:
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Standardized non-stationarity measurment errors - based on the slope and the maximum absolute curvature")
                    self.text_description(f" Standard deviation before scale-correction :  **{standard_deviation_before_scale_correction :.2f}**") 
                    self.text_description(f" Standard deviation after scale-correction :  **{standard_deviation_after_scale_correction :.2f}**")
                    self.image_body(curvature_slope_nmad_out)
                    self.image_body(variance_sample_count_standardization)
                    self.image_body(standardized_elevation_difference_clipped_out)
                    self.image_body(standardized_elevation_difference_out)
                    self.number_subtitle+=1
                    
                if standardized_average_shapefile_error:
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Standardized integrated error and standardized mean error - based on the slope and the maximum absolute curvature")
                    self.text_description(f" Standard deviation before scale-correction :  **{standard_deviation_before_scale_correction :.2f}**") 
                    self.text_description(f" Standard deviation after scale-correction : ** {standard_deviation_after_scale_correction :.2f}**")
                    self.text_description(f" Surface area (glacier:{name_shapefile_glacier_errors_out}) :  **{area_shapefile_glacier_errors :.2f} m²**")
                    self.text_description(f" Average slope (glacier : {name_shapefile_glacier_errors_out}):  **{average_slope_local_glacier :.2f} °**")
                    self.text_description(f" Number of effective samples  (glacier :{name_shapefile_glacier_errors_out}) :  **{local_glacier_neff :.2f}**")
                    self.text_description(f" Standardized integrated error(glacier :{name_shapefile_glacier_errors_out}) :  **{local_glacier_z_err :.4f}**")
                    self.text_description(" Average elevation difference and integrated error (destandardized)(glacier {}) :  **{:.3f} ± {:.3f} m**".format(name_shapefile_glacier_errors_out,local_glacier_dh, local_glacier_dh_err_integrated))
                    self.text_description(f" Uncertainty of the volume change (glacier :{name_shapefile_glacier_errors_out}) :  **{uncertainty_of_the_volume_change :.2f} m³**")
                    self.text_description(f" Standardized mean error over the glacier : **{local_glacier_dh_err_mean :.2f} m**")
                    self.image_body(slope_out) 
                    self.image_body(elevation_difference_out) 
                    self.number_subtitle+=1                            
                
                self.number_title+=1
                self.add_page()
                
            #%%% Calculation mass
            if calculation_mass:
                self.chapter_title(f" {self.number_title} - Volume and mass calculation ")
                self.number_subtitle=1
                if use_ddem_sub and subtraction: 
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Volume and mass calculation of the dDEM calculated in the subtraction step")
                    self.text_description(f" dDEM : {ddem_subtraction_name}")
                    self.ln(10)
                    self.text_body(f"{path_temporary}/volume_mass_subtraction_results.doc")
                    self.number_subtitle+=1
                
                if use_ddem_inter and interpolation:
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Volume and mass calculation of the dDEM calculated in the interpolation step")
                    self.text_description(f" dDEM : {ddem_interpolation_name}")
                    self.ln(10)
                    self.text_body(f"{path_temporary}/volume_mass_interpolation_results.doc")
                    self.number_subtitle+=1
                
                if ddem_volume: 
                    self.text_title(f" {self.number_title}.{self.number_subtitle} - Volume and mass calculation of the dDEM calculated selected by the user")
                    self.text_body(f"{path_temporary}/volume_mass_own_results.doc")
                    self.number_subtitle+=1
                self.number_title+=1
                self.add_page()
                
    #%%% Launch the pdf page
    #compute the pdf page  and write
    progress_bar_value=loader_in_progress(progress_bar_value)
    pdf = PDF()
    
    pdf.print_page()
    progress_bar_value=loader_in_progress(progress_bar_value)
    pdf.output(f'{path_out}/report_{name_pdf}.pdf') 
    
    #Beginning update loader
    progress_bar_value=loader_in_progress(progress_bar_value)
    
    event, values = window.Read(timeout=0)
    if event == sg.WIN_CLOSED or event == 'Cancel':
        window.close()
        return "The function has been stopped"
        
    #end update loader
    
    #delete temporary folder and files
    shutil.rmtree(path_temporary)
    print("Temporary folder has been deleted")
    
    print("All steps are computed.")
    
    window.close() #close the loader window
    
    
    return figure_plot

    
    
    
    # ------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    #test code for developper
    
test=False


if test:
    Treatment_dem(    
        dem_1="Folgefonna/folgefonna_input/1962_DEM_v2.tif",             
        dem_2="Folgefonna/folgefonna_input/2017_DEM_NVE_resized_3.tif" ,            
        vector_glacier="Folgefonna/folgefonna_input/glacier_outlines/1959_glacier_area_outline/Folgefonna1959_utm32.shp",
        path_out="Folgefonna/files_out",
        coordinate_system="EPSG:32633" ,
        vector_frame= "Folgefonna/folgefonna_input/frame_folgefonna/frame_folgefonna.shp",
        resolution=10,
        coregistration=False ,
        dem_coreg_ref="Folgefonna/folgefonna_input/2017_DEM_NVE_resized_3.tif",
        method_coregistration="ICP + Nuth and Kääb",
        subtraction=False, 
        use_dem_coreg_sub=False,
        interpolation=True,
        use_dem_coreg_inter=False,
        method_interpolation="linear interpolation",
        dem_inter_ref="Folgefonna/coregistered/coregistration_2013/2013_coregistered_on_2017.tif", 
        ddem_calculated="Folgefonna/folgefonna_input/1962_DEM_v2.tif",
        unstandardized_nonstationarity_error=False,
        standardization_nonstationarity_error=False, # - Type: Boolean - If true standardization nonstationnary error is activated
        standardized_average_shapefile_error=False,
        ref_standardization_nonstationarity_error="Folgefonna/folgefonna_input/2017_DEM_NVE_resized_3.tiff", # - Type: Geotiff - reference DEM used to calculate standardization non stationnary error
        ddem_standardization_nonstationarity_error=None, # - Type: Geotiff - ddem for calculating the error
        glacier_outlines_error="Folgefonna/folgefonna_input/glacier_outlines/2018/2018_glacier_outlines_polygon.shp",
        calculation_mass=False, 
        ddem_volume=None,
        use_ddem_sub=True,
        use_ddem_inter=True,
        rho_density=float("0.917"),
        error_rho_density=float("0.05"),                 # error of the density
        error_glacier_outline=float("150"),
        name_pdf="interpolation_1962_v2"
        )

    
    
