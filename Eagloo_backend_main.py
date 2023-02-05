# =============================================================================
# Process : 
# 1/ Coregistration with ICP + Nuth and Kaab method 
# 2/ Calculates subtraction between two DEMs
# 3/ Computes hypsometric interpolation
# 4/ Calculates the volume of the ddem and mass (linear conversion)

# Files that are computed in the coregistration step are used in subtraction, interpolation, and volumes calculation steps.
# It doesn't exist any correlation between interpolation and subtraction steps.
# =============================================================================

import os
import geoutils as gu
import matplotlib.pyplot as plt
import numpy as np
import xdem
import rasterio
import rasterio.mask
import geopandas
import fiona


from rasterio.transform import from_origin
from datetime import datetime
from osgeo import gdal
from Eagloo_backend_functions.vol_stats import volume_mass_stats as vol_stats
from fpdf import FPDF
import shutil


def Treatment_dem(dem_1,                #- Type: Geotiff - Older DEM used for subtraction, file to co-register
                  dem_2 ,               #- Type: Geotiff - (optionnal for coregistration, None if empty) Newer DEM used for subtraction (dem_2-dem_1), dem_2 is coregistered only if it is different of dem_coreg_ref
                  vector_glacier ,
                  path_out ,
                  coordinate_system ,
                  vector_frame ,
                  resolution ,
                  coregistration ,
                  dem_coreg_ref,
                  method_coregistration,
                  subtraction, 
                  use_dem_coreg_sub,
                  interpolation,
                  use_dem_coreg_inter,
                  method_interpolation,
                  dem_hypso_ref, 
                  ddem_calculated,
                  calculation_mass, 
                  ddem_volume,
                  use_ddem_sub,
                  use_ddem_inter,
                  
                  ):
    print("The process is launched")
    # # =============================================================================
    # Necessary entries(To fill out):
            # dem_1 - Type: Geotiff - Older DEM used for subtraction, file to co-register
            # dem_2 - Type: Geotiff - (optionnal for coregistration, None if empty) Newer DEM used for subtraction (dem_2-dem_1), dem_2 is coregistered only if it is different of dem_coreg_ref
            # vector_glacier - Type: Shapefile -  Glacier mask of the year where the glacier is the biggest
            
            # path_out - Type:path -  path where files will be saved
            
            # coordinate_system - ex:EPSG:32633 -  coordinate system used during the gdal.Warp process
            # vector_frame: - Type: Shapefile - Shapefile of the global area, this file is utilized for preprocessing data
            # resolution: Type: Integer - Define the resolution, this resolution will be applied to each DEM, x and y are identical
        
    
    # Coregistration:
        # The coregistration step is launched only if coregistration button is ticked
            
            # Entries :
                
                # coregistration - Type: Boolean - If "True" coregistration is activated
                # dem_coreg_ref - Type: Geotiff - It is necessary for the coregistration, this DEM is used as a reference for the coregistration, 
                # method_coregistration - Type: String - One of these choices must be selected  "Nuth and Kaab","ICP", "ICP + Nuth and Kaab"
    # Sutraction between two DEMs
    
        # subtraction=True       # True -> export the layer dem2-dem1 and print a graph of this one
        # use_dem_coreg_sub - Type: Boolean - if you want select the DEMs calculated during the co-registration step

    
    # #Interpolation
    # interpolation=True     # true -> export the hypsometric interpolation of a subtraction ddem(dem2-dem1)
    # dem_hypso_ref="1959_2013/2013_NVE_WGS.tif"   # DEM used as a reference for the hypsometric interpolation and coregistration
    
    # #Volume and mass calculation
    # calculation_mass=True
    # ddem_volume="1959_2013/hypsometric_interpolation_1959_DEM_contour_2017_DEM_NVE_resized_3.tif"
    
    #End commands
    # =============================================================================
    
    #Beginning Homogeneization of resolution, frame and georeferencing
    
    #Intialization function
    # figure_coreg_dem1,figure_coreg_dem2,figure_sub,figure_inter_map,figure_inter_graph=None
    figure_plot=[None,None, None, None,None] #list return at the end of the algorithm: figure computed during matplotlib process
    # beginning remove path and extension, for writing all files  

    
    
    vector_glacier_out=os.path.basename(vector_glacier).split('.')[0]
    vector_frame_out=os.path.basename(vector_frame).split('.')[0]
    # end remove path and extension, for writing all files
    
    
    #Beginning creating temporary folder
    
    directory = "temporary"
    os.makedirs(os.path.join(path_out, directory) , exist_ok="False") 
    path_temporary=os.path.join(path_out, directory)
    #End creating temporary folder
    
    # Beginning homogenization coordinate glacier shapefile
    shapefile_glacier = geopandas.read_file(vector_glacier)
    # change CRS to coordinate_system
    shapefile_glacier= shapefile_glacier.to_crs(coordinate_system)
    # write shp file
    vector_glacier_preprocessed="{}/preprocessed_{}.shp".format(path_temporary,vector_glacier_out)
    shapefile_glacier.to_file(vector_glacier_preprocessed)
    # End homogenization coordinate glacier shapefile
    
    # Beginning homogenization coordinate glacier shapefile
    shapefile_frame = geopandas.read_file(vector_frame)
    # change CRS to coordinate_system
    shapefile_frame= shapefile_frame.to_crs(coordinate_system)
    # write shp file
    vector_frame_preprocessed="{}/preprocessed_{}.shp".format(path_temporary,vector_frame_out)
    shapefile_frame.to_file(vector_frame_preprocessed)
    # End homogenization coordinate glacier shapefile
    outlines_glacier = gu.Vector("{}".format(vector_glacier_preprocessed))
    
    # =============================================================================
    if dem_1: 
        dem1_out=os.path.basename(dem_1).split('.')[0] 
        dem_1_name="{}/preprocessing_{}.tif".format(path_temporary,dem1_out)
        gdal.Warp(dem_1_name ,dem_1 , cutlineDSName=vector_frame_preprocessed, dstSRS=coordinate_system,  xRes=resolution, yRes=resolution, cropToCutline=True)
        dem1 = xdem.DEM("{}".format(dem_1_name))                #opening dem_1_name
        
    if dem_2:
        dem2_out=os.path.basename(dem_2).split('.')[0]
        dem_2_name="{}/preprocessing_{}.tif".format(path_temporary,dem2_out)
        gdal.Warp(dem_2_name ,dem_2 , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
        dem2 = xdem.DEM("{}".format(dem_2_name))
    
    
    
    
    
    print("Preprocessing is done")
    #End Homogeneization of resolution, frame and georeferencing
    
    # =============================================================================
    
    
    #Beginning Coregistration ICP + Nuth And Kaab
    
    if coregistration : 
        dem_coreg_ref_out=os.path.basename(dem_coreg_ref).split('.')[0]
        #beginning homogenization of reference DEM
        dem_ref_preprocessing="{}/preprocessing_{}.tif".format(path_temporary,dem_coreg_ref_out)
        gdal.Warp(dem_ref_preprocessing ,dem_coreg_ref , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)  # !!!warning mask is not apply currently
        #End homogenization of reference DEM
        
        demref_coreg = xdem.DEM("{}".format(dem_ref_preprocessing))
        inlier_mask = ~outlines_glacier.create_mask(demref_coreg)
        #beginning masking co-registration file
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
        #end masking co-registration file
        
        #beginning choose method
        method=None
        if method_coregistration=="ICP":
            method= xdem.coreg.ICP()
        if method_coregistration=="Nuth and Kaab":
            method= xdem.coreg.NuthKaab()
        if method_coregistration=="ICP + Nuth and Kaab":
            method=xdem.coreg.ICP() + xdem.coreg.NuthKaab()
        assert method!= None , " method_coregistration variable should be set out by choosing following entries : ICP,Nuth and Kaab,ICP + Nuth and Kaab "
        #end choose method
        _ = dem1.reproject(demref_coreg)
        diff_before_dem1=coreg_ref_masked.data - dem1.data
        transform_dem_1 = coreg_ref_masked.transform
        
        method.fit(reference_dem=coreg_ref_masked.data, dem_to_be_aligned=dem1.data,inlier_mask=inlier_mask,transform=transform_dem_1,verbose=True,)
    
        dem1_coreg = method.apply(dem=dem1, transform=transform_dem_1)
        # dem_residuals_1=ICPNuthKaab.residuals()
        dem1_coreg_path="{}/coregistration_{}_{}_on_{}.tif".format(path_out,method_coregistration,dem1_out,dem_coreg_ref_out)
        dem1_coreg.save(dem1_coreg_path)
    
        diff_after_dem1 = coreg_ref_masked.data - dem1_coreg.data
    
        # display the graph of the difference before and after
        figure_plot[0]=plt.figure(num=0,figsize=(8, 5))
        plt.suptitle('Co-registration, method : {} \n Co-registered DEM :{} \n Reference DEM :{}'.format(method_coregistration,dem1_out,dem_coreg_ref_out))
        plt.subplot2grid((1, 15), (0, 0), colspan=7)
        plt.title("Before coregistration. ")
        plt.imshow(diff_before_dem1.squeeze(), cmap="coolwarm_r", vmin=np.min(diff_before_dem1.squeeze()), vmax=np.max(diff_before_dem1.squeeze()))
        plt.axis("off")
        plt.subplot2grid((1, 15), (0, 7), colspan=7)
        plt.title("After coregistration. ")
        img = plt.imshow(diff_after_dem1.squeeze(), cmap="coolwarm_r", vmin=np.min(diff_after_dem1.squeeze()), vmax=np.max(diff_after_dem1.squeeze()))
        plt.axis("off")
        plt.subplot2grid((1, 15), (0, 14), colspan=1)
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        plt.axis("off")
        plt.tight_layout()
        coreg_1_image="{}/coregistration_{}_{}_ref_{}.png".format(path_out,method_coregistration,dem1_out,dem_coreg_ref_out)
        plt.savefig(coreg_1_image)
        plt.show()
        
        #calculate NMAD before and after coregistration
        NMAD_before_DEM1=xdem.spatialstats.nmad(diff_before_dem1)
        NMAD_after_DEM1=xdem.spatialstats.nmad(diff_after_dem1)
        print(f"Error before coregistration: {NMAD_before_DEM1:.3f} m")
        print(f"Error after coregistration: {NMAD_after_DEM1:.2f} m")
        
       
        print ( " \n {} has been co-registered - {} is the reference  \n" .format(dem1_out,dem_coreg_ref_out) )
        
        if dem_2 !=None and dem_2 != dem_coreg_ref: 
             
            
            _ = dem2.reproject(demref_coreg)
            print("dem_2 is different of dem_coreg_ref, thus dem_2 is co-registering and dem_coreg_ref is the reference")
            diff_before_dem2 = coreg_ref_masked.data - dem2.data
            transform_dem_2 = demref_coreg.transform
            
            method.fit(reference_dem=coreg_ref_masked.data, dem_to_be_aligned=dem2.data,inlier_mask=inlier_mask,transform=transform_dem_2)
    
            dem2_coreg = method.apply(dem=dem2, transform=transform_dem_2)
            dem2_coreg_path="{}/coregistration_{}_{}_on_{}.tif".format(path_out,method_coregistration,dem2_out,dem_coreg_ref_out)
            dem2_coreg.save(dem2_coreg_path)
    
            diff_after_dem2 = coreg_ref_masked.data - dem2_coreg.data
            
            # display the graph of the difference before and after
            figure_plot[1]=plt.figure(num=1,figsize=(8, 5))
            plt.suptitle('Co-registration, method : {} :\n Co-registered  DEM:{} \n Reference DEM :{}'.format(method_coregistration,dem2_out,dem_coreg_ref_out))
            plt.subplot2grid((1, 15), (0, 0), colspan=7)
            plt.title("Before coregistration. ")
            plt.imshow(diff_before_dem2.squeeze(), cmap="coolwarm_r", vmin=np.min(diff_before_dem2.squeeze()), vmax=np.max(diff_before_dem2.squeeze()))
            plt.axis("off")
            plt.subplot2grid((1, 15), (0, 7), colspan=7)
            plt.title("After coregistration. ")
            img = plt.imshow(diff_after_dem2.squeeze(), cmap="coolwarm_r", vmin=np.min(diff_after_dem2.squeeze()), vmax=np.max(diff_after_dem2.squeeze()))
            plt.axis("off")
            plt.subplot2grid((1, 15), (0, 14), colspan=1)
            cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
            cbar.set_label("Elevation change (m)")
            plt.axis("off")
            plt.tight_layout()
            coreg_2_image="{}/coregistration_{}_{}_ref_{}.png".format(path_out,method_coregistration,dem2_out,dem_coreg_ref_out)
            plt.savefig(coreg_2_image)
            plt.show()
           
            print ( "\n {} has been coregistered - {} is the reference \n" .format(dem2_out,dem_coreg_ref_out))
            
            #calculate NMAD before and after coregistration
            NMAD_before_DEM2=xdem.spatialstats.nmad(diff_before_dem2)
            NMAD_after_DEM2=xdem.spatialstats.nmad(diff_after_dem2)
            print(f"Error before coregistration: {NMAD_before_DEM2:.3f} m")
            print(f"Error after coregistration: {NMAD_after_DEM2:.2f} m")
            
        
    #End Coregistration ICP + Nuth And Kaab
    
    # =============================================================================
    
    # Beginning subtraction between the 2 DEMs
    if subtraction:
        #if you want to use the DEMs coregistered previously
        if use_dem_coreg_sub and coregistration:  
            _ = dem1_coreg.reproject(dem2_coreg)
            ddem_no_interpolated= dem2_coreg-dem1_coreg
        else:
            ddem_no_interpolated= dem2-dem1
        
        name_ddem="{}/subtraction_{}_{}.tif".format(path_temporary,dem1_out,dem2_out)
        ddem_no_interpolated.save(name_ddem)
        
        #clip ddem with glacier outlines shapefile
        with fiona.open(vector_glacier_preprocessed, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]
    
        with rasterio.open(name_ddem) as src:
            ddem_out_image, out_transform = rasterio.mask.mask(src, shapes, filled=True)
            ddem_out_meta = src.meta
    
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
     
      
        figure_plot[2]=plt.figure(num=2,figsize=(8, 5))
        plt.suptitle("Subtraction : {} - {} " .format(dem2_out,dem1_out))
        img = plt.imshow(ddem_clipped_out.data.squeeze(), cmap="coolwarm_r", vmin=np.min(ddem_clipped_out.data.squeeze()), vmax=np.max(ddem_clipped_out.data.squeeze()))
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        plt.axis("off")
        subtraction_image="{}/DDEM_subtraction_clipped_{}_{}.png".format(path_out,dem1_out,dem2_out)
        plt.savefig(subtraction_image)
        plt.show()   
        
        print ("\n {} - {} has been processed \n" .format(dem2_out,dem1_out))
    # End subtraction difference between the 2 DEMs
    
    # =============================================================================
    
    # Beginning Hypsometric interpolation
    
    if interpolation:
        dem_hypso_ref_out=os.path.basename(dem_hypso_ref).split('.')[0]
        #beginning preprocessing
        dem_ref_hypso_name_out="{}/preprocessing_{}.tif".format(path_temporary,dem_hypso_ref_out)
        gdal.Warp(dem_ref_hypso_name_out ,dem_hypso_ref , cutlineDSName=vector_frame_preprocessed,dstSRS=coordinate_system, xRes=resolution, yRes=resolution, cropToCutline=True)
        demref_hypso = xdem.DEM("{}".format(dem_ref_hypso_name_out), datetime=datetime(2013,8,1))
       
        src = gdal.Open("{}".format(dem_ref_hypso_name_out))                #open file with gdal - type:gdal.Dataset - this layer us used to get information on rasters (resolution, position of the first pixel)
        gt=src.GetGeoTransform()
        glacier_index_map = outlines_glacier.rasterize(demref_hypso)
        
        inlier_mask = outlines_glacier.create_mask(demref_hypso)
        #end preprocessing
        outlines = {
            datetime(1959, 8, 1): outlines_glacier,
        }
        
        
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
        
        #end clipping
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
            nodata_value=None
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
            interpolation_name="{}/hypsometric_interpolation_{}.tif".format(path_out,ddemc_out)
        else:
            interpolation_name="{}/hypsometric_interpolation_{}_{}.tif".format(path_out,dem1_out,dem2_out)
        
        new_dataset= rasterio.open(interpolation_name, 'w',height=ddem_filled_correction.data.shape[0], width=ddem_filled_correction.data.shape[1], driver='GTiff', count=1, dtype=ddem_filled_correction.dtype, crs=coordinate_system, transform=transform_interpolation, nodata=nodata_value) 
        new_dataset.write(ddem_filled_correction,1)
        new_dataset.close()
           
        
        print(gt)
        
        print("Hypsometric interpolation is processed.")
    
        figure_plot[3]=plt.figure(num=3,figsize=(8, 5))
        if ddem_calculated:
            plt.suptitle('Interpolation, method: {} \n dDEM = {} ;\n Reference DEM: {}'.format(method_interpolation,ddemc_out,dem_hypso_ref_out))
            
        else:
            plt.suptitle('Interpolation, method: {} \n dDEM = {}-{} ;\n Reference DEM: {}'.format(method_interpolation,dem2_out,dem1_out,dem_hypso_ref_out))
        plt.subplot2grid((1, 15), (0, 0), colspan=7)
        plt.title("Before interpolation ")
        plt.imshow(ddem_clipped_out.data.squeeze(), cmap="coolwarm_r", vmin=np.min(ddem_clipped_out.data.squeeze()), vmax=np.max(ddem_clipped_out.data.squeeze()))
        plt.axis("off")
        plt.subplot2grid((1, 15), (0, 7), colspan=7)
        plt.title("After interpolation")
        img = plt.imshow(ddem_filled_correction.squeeze(), cmap="coolwarm_r", vmin=np.min(ddem_filled_correction.squeeze()), vmax=np.max(ddem_filled_correction.squeeze()))
        plt.axis("off")
        plt.subplot2grid((1, 15), (0, 14), colspan=1)
        cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
        cbar.set_label("Elevation change (m)")
        plt.axis("off")
        plt.tight_layout()
        if ddem_calculated:
            interpolation_image="{}/interpolation_{}_ddem_{}_ref_{}.png".format(path_out,method_interpolation,ddemc_out,dem_hypso_ref_out)
            
        else:
            interpolation_image="{}/interpolation_{}_ddem_={}-{}_ref_{}.png".format(path_out,method_interpolation,dem1_out,dem2_out,dem_hypso_ref_out)
        plt.savefig(interpolation_image)
        plt.show()
        
    
        # graph of the hypsometric interpolation
    
    
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
        plt.tight_layout()
        if ddem_calculated:
            graph_hypso_image="{}/graph_interpolation_{}_ddem_{}_ref_{}.png".format(path_out,method_interpolation,ddemc_out,dem_hypso_ref_out)
            
        else:
            graph_hypso_image="{}/graph_interpolation{}_{}_{}_ref_{}.png".format(path_out,method_interpolation,dem1_out,dem2_out,dem_hypso_ref_out)
        plt.savefig(graph_hypso_image)
        plt.show()
        
    # End Hypsometric interpolation
    
    # =============================================================================
    
    # Beginning volume and mass statistics
    if calculation_mass: 
        if ddem_volume:
            vol_stats('mass',ddem_volume,0.917,path_temporary,'own')
            
        if interpolation and use_ddem_inter: 
            print("volume and mass stats calculated with hypsometric interpolation")
            vol_stats('mass',interpolation_name,0.917,path_temporary,'interpolation')
            
        if subtraction and use_ddem_sub:
            print("volume and mass stats calculated without hypsometric interpolation")
            vol_stats('mass',ddem_clipped_name,0.917,path_temporary,'subtraction')
            
    # end volume and mass statistics
      
    # =============================================================================
       
    # Generate PDF 
    title="Report : "
    if coregistration:
        title = title + f"co-registration, {method_coregistration} "
    if subtraction:
        title=title + " - Subtraction"
    if interpolation:
        title= title +" - Interpolation"
    if calculation_mass:
        title= title +" - Volume and mass calculation"
    
    
       
    class PDF(FPDF):
        def __init__(self):
            super().__init__()
            self.WIDTH = 180
            self.HEIGHT = 180
            self.position_graph_pdf="left"  #use to display graphs at right or at left
            self.y_iterator=60              #use to display the graphs on y axis     
            
        def header(self):
            # Custom logo and positioning
            # Create an `assets` folder and put any wide and short image inside
            # Name the image `logo.png`
            self.image("Eagloo_frontend_functions/eagloo_logo.png", 8, 8, 20)
            self.set_font('Arial', 'B', 11)
            self.cell(self.WIDTH -45)
            self.cell(60, 2, title, 0, 0, 'R')
            self.ln(20)
            
        def footer(self):
            # Page numbers in the footer
            self.set_y(-15)
            self.set_font('Times', 'I', 8)
            self.set_text_color(128)
            self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')
            
        def chapter_title(self, label):
            # Arial 12
            self.set_font('Times', 'B', 15)
            # Background color
            self.set_fill_color(200, 220, 255)
            self.ln(3)
            # Title
            self.cell(0, 6, "{} ".format(label),  0, 1, 'L', 1)
            # Line break
            self.ln(3)
    
        def image_body(self, image1):
            # Determine how many plots there are per page and set positions
            # and margins accordingly
            self.ln(5)
            self.image(image1, None, None,self.WIDTH)
            self.ln(5)
            
        def print_raster(self, raster,namegraph,position,positiony):
            #entries
            #raster1 :Gtiff
            #postion: "right" or "left"
            #output: raster print on the page
            graph_pdf=xdem.DEM(f"{raster}")
            plt.suptitle(namegraph)
            img = plt.imshow(graph_pdf.data.squeeze(), cmap="coolwarm_r", vmin=np.min(graph_pdf.data.squeeze()), vmax=np.max(graph_pdf.data.squeeze()))
            cbar = plt.colorbar(img, fraction=0.4, ax=plt.gca())
            cbar.set_label("Elevation change (m)")
            plt.axis("off")
            path_graph="{}/{}.png".format(path_out,namegraph)
            plt.savefig(path_graph)
            plt.show()   
            if position=="right":
                self.image(path_graph,100, positiony,100)
            else:
                self.image(path_graph, 0, positiony,100)
                
        def position_graph(self):
            if self.position_graph_pdf=="left":
                self.position_graph_pdf="right"
                
            else :
                self.position_graph_pdf="left"
                self.y_iterator= self.y_iterator+50
                self.ln(3)
            
            
        def text_title(self, title):
            self.ln(5)
            self.set_font('Times', 'B', 14)
            self.cell(0,5,f"{title}")
            self.ln(10)
            
        
        def text_body(self, path):
            # Read text file
            with open(path, 'rb') as fh:
                txt = fh.read().decode('latin-1')
            # Times 12
            self.set_font('Times', '', 12)
            # Output justified text
            self.multi_cell(0, 5, txt)
            # Line break
            self.ln(5)
            # Mention in italics
            
        def text_description(self,text):
            self.ln(5)
            self.set_font('Times', '', 12)
            self.cell(0,5,f"{text}")
            self.ln(5) 
               
        def print_page(self):
            # Generates the report
            self.add_page()
            
            #entries
            self.chapter_title("Entries")
            if dem_1:
                self.text_description(f"- Earliest DEM : {dem1_out} \n")
                
            if dem_2:
                self.text_description(f"- Latest DEM : {dem2_out} \n")
                
            self.text_description(f" - Contour line shapefile of the biggest glacier(Use for masking) : {vector_glacier_out} \n")
            self.text_description(f" - Path where the files were saved : {path_out} \n")
            self.text_description(f" - Coordinate system : {coordinate_system} \n")
            self.text_description(f" - Shapefile frame of the study area : {vector_frame_out} \n")
            self.text_description(f" - Resolution : {resolution} \n\n") 
            self.y_position=self.get_y()
            if dem1:
                self.print_raster(dem_1,f"Earliest DEM: {dem1_out}",self.position_graph_pdf,self.y_position+self.y_iterator)
                self.position_graph()
                
            if dem2:
                self.print_raster(dem_2,f"Latest DEM: {dem2_out}",self.position_graph_pdf,self.y_position+self.y_iterator)
                self.position_graph()
                
            if dem_coreg_ref:
                self.print_raster(dem_coreg_ref,f"Reference DEM for the co-registration: {dem_coreg_ref_out}",self.position_graph_pdf,self.y_position+self.y_iterator)
                self.position_graph()
                
            if dem_hypso_ref:
                self.print_raster(dem_hypso_ref,f"Reference DEM for the interpolation: {dem_coreg_ref_out}",self.position_graph_pdf,self.y_position+self.y_iterator)
                self.position_graph()
                
            if ddem_volume:
                self.print_raster(ddem_volume,f"DDEM to calculate the volume and the mass: {os.path.basename(ddem_volume).split('.')[0]}",self.position_graph_pdf,self.y_position+self.y_iterator)
                self.position_graph()
            
            if coregistration:
                self.chapter_title(f" - Co-registration of the DEM: {dem1_out}")
                self.text_description(f" - Reference DEM for co-registration : {dem_coreg_ref_out} ")
                self.text_description(f" - Method : {method_coregistration} ")
                self.text_description(f"- Error before the co-registration(NMAD) : {NMAD_before_DEM1} ")
                self.text_description(f" - Error after the co-registration(NMAD) : {NMAD_after_DEM1} ")
                self.image_body(coreg_1_image)
                
                if dem_2:
                    self.chapter_title(f" - Co-registration of the DEM: {dem2_out}")
                    self.text_description(f"- Reference DEM for co-registration : {dem_coreg_ref_out} ")
                    self.text_description(f" - Method : {method_coregistration} ")
                    self.text_description(f" - Error before the co-registration(NMAD) : {NMAD_before_DEM2} ")
                    self.text_description(f" - Error after the co-registration(NMAD) : {NMAD_after_DEM2} ")
                    self.image_body(coreg_2_image)
            self.ln(10)
            if subtraction:
                self.chapter_title(f"Subtraction : {dem1_out} - {dem2_out}")
                self.image_body(subtraction_image)
            
            if interpolation:
                self.chapter_title("Interpolation")
                if dem_1 and dem_2: 
                    self.text_description(f"dDEM = {dem2_out} - {dem1_out}")
                    
                if ddem_calculated:
                    self.text_description(f"dDEM = {ddemc_out}")
                    
                self.text_description(f" - Reference DEM for interpolation:{dem_hypso_ref_out}")
                self.text_description(f" - Method : {method_interpolation}")   
                
                if use_dem_coreg_inter and coregistration: 
                    self.text_description(" The DEMs co-registered in the co-registration step were used.")
                    
                self.image_body(interpolation_image)
                self.image_body(graph_hypso_image)
                
            if calculation_mass:
                self.chapter_title("Volume and mass calculation ")
                if use_ddem_sub and subtraction: 
                    self.text_title("Volume and mass calculation of the dDEM calculated in the subtraction step")
                    self.text_description(f"dDEM : {ddem_subtraction_name}")
                    self.text_body(f"{path_temporary}/volume_mass_subtraction_results.txt")
                
                if use_ddem_inter and interpolation:
                    self.text_title("Volume and mass calculation of the dDEM calculated in the interpolation step")
                    self.text_description(f"dDEM : {ddem_interpolation_name}")
                    self.text_body(f"{path_temporary}/volume_mass_interpolation_results.txt")
                
                if ddem_volume: 
                    self.text_title("Volume and mass calculation of the dDEM calculated selected by you")
                    self.text_body(f"{path_temporary}/volume_mass_own_results.txt")
                    
    #compute the pdf page  and write
    pdf = PDF()
    pdf.print_page()
    pdf.output(f'{path_out}/Treatments_report.pdf', 'F')
    
    #delete temporary folder and files
    shutil.rmtree(path_temporary)
    print("Temporary folder has been deleted")
    
    print("All steps are computed.")
    
    
    return figure_plot
   
    
    
    
    # ------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    #test code for developper
    
test=False

if test:
    Treatment_dem(    
        dem_1="1959_2013/1959_DEM_contour.tif" ,              #- Type: Geotiff - Older DEM used for subtraction, file to co-register
        dem_2="1959_2013/2017_DEM_NVE_resized_3.tif" ,            #- Type: Geotiff - (optionnal for coregistration, None if empty) Newer DEM used for subtraction (dem_2-dem_1), dem_2 is coregistered only if it is different of dem_coreg_ref
        vector_glacier="1959_2013/Folgefonna1959_utm32.shp",
        path_out="1959_2013/files_out",
        coordinate_system="EPSG:32633" ,
        vector_frame= "1959_2013/frame_reference_homogeneisation_xdem.shp",
        resolution=50 ,
        coregistration=True ,
        dem_coreg_ref="1959_2013/2013_NVE_WGS_co_delete.tif",
        method_coregistration="Nuth and Kaab",
        subtraction=True, 
        use_dem_coreg_sub=True,
        interpolation=True,
        use_dem_coreg_inter=True,
        method_interpolation="linear interpolation",
        dem_hypso_ref="1959_2013/2013_NVE_WGS.tif", 
        ddem_calculated=None,
        calculation_mass=None, 
        ddem_volume=None,
        use_ddem_sub=True,
        use_ddem_inter=True,
        )

    
    