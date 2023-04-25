#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 17:48:29 2022

@author: Guillaume Pfundstein

V2
"""

path_function="Eagloo_frontend_functions"

import PySimpleGUI as sg
import Eagloo_backend_main as tdem
import os.path
import PIL.Image
import io
import base64
import traceback
import sys
import fitz        #display_pdf
import psutil #CPU
import matplotlib  #print
import inspect
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import concurrent.futures
import queue
import threading
import time
from threading import Thread


"""
    calculate :
    co-registration
    subtraction
    interpolation
    volume and mass

    Copyright 2022 
"""

#style
size_name_menu=45
font=("Helvetica", 11, "bold")
sg.set_options(font=font,
               icon=f'{path_function}/eagloo_logo_drawing.png',
               
               )

sg.LOOK_AND_FEEL_TABLE['Native'] = {
    'BACKGROUND': '#5D6D7E',
    'TEXT': '#F2F3F4',
    'INPUT': '#ABB2B9',
    'FONT': 'Arial 20',
    'TEXT_INPUT': '#17202A',
    'SCROLL': '#a5a4a4',
    'BUTTON': ("#ECECEA", "#C19331"),
    'PROGRESS': sg.DEFAULT_PROGRESS_BAR_COLOR,
    'BORDER': 1,
    'SLIDER_DEPTH': 0,
    'PROGRESS_DEPTH': 0,
    'ACCENT1': '#ff5414',
    'ACCENT2': '#33a8ff',
    'ACCENT3': '#dbf0ff',
    


    }
    
sg.theme('Native')

 # %% Functions
# -------------------------------------------------------------------------------
# Beginning definition function

#matplotlib
def convert_to_bytes(file_or_bytes, resize=None):
    '''
    Will convert into bytes and optionally resize an image that is a file or a base64 bytes object.
    Turns into  PNG format in the process so that can be displayed by tkinter
    :param file_or_bytes: either a string filename or a bytes base64 image object
    :type file_or_bytes:  (Union[str, bytes])
    :param resize:  optional new size
    :type resize: (Tuple[int, int] or None)
    :return: (bytes) a byte-string object
    :rtype: (bytes)
    '''
    
    if isinstance(file_or_bytes, str):
        img = PIL.Image.open(file_or_bytes)
        
    else:
        try:
            img = PIL.Image.open(io.BytesIO(base64.b64decode(file_or_bytes)))
        except Exception:
            dataBytesIO = io.BytesIO(file_or_bytes)
            img = PIL.Image.open(dataBytesIO)

    cur_width, cur_height = img.size
    if resize:
        new_width, new_height = resize
        scale = min(new_height/cur_height, new_width/cur_width)
        img = img.resize((int(cur_width*scale), int(cur_height*scale)), PIL.Image.ANTIALIAS)
    bio = io.BytesIO()
    img.save(bio, format="PNG")
    del img
    return bio.getvalue()

# Convert figure matplot lib to canvas
def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

# Delete figure popup matplotlib
def delete_figure_agg(figure_agg):
    figure_agg.get_tk_widget().forget()
    plt.close('all')
    
def message_display(text):
    layout = [[sg.Image(convert_to_bytes(f'{path_function}/eagloo_logo_drawing.png',resize=(70,70)))],
        [sg.Text(f'{text}')],
        [sg.Button("Close", key="-Close-")]]
    window = sg.Window('Eagloo: Message',layout ,finalize=True, element_justification='center')
    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED or event == "-Close-" :
            break
        
    window.close()



if sys.version_info >= (3, 0):
    _thread_target_key = '_target'
    _thread_args_key = '_args'
    _thread_kwargs_key = '_kwargs'
else:
    _thread_target_key = '_Thread__target'
    _thread_args_key = '_Thread__args'
    _thread_kwargs_key = '_Thread__kwargs'

class ThreadWithReturn(Thread):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._return = None

    def run(self):
        target = getattr(self, _thread_target_key)
        if not target is None:
            self._return = target(*getattr(self, _thread_args_key), **getattr(self, _thread_kwargs_key))

    def join(self, *args, **kwargs):
        super().join(*args, **kwargs)
        return self._return

 # %% Layout
# --------------------------------- Define Layout ---------------------------------

#end definition fonction

menu_def = [['File', ['Open', 'Save','Save as', 'Exit']], ['Help', 'About'],  ] 


 # %%% Input layout

tab_main_layout=[
    [sg.Text('', size=(90, 1))],
    [sg.Text('Earliest DEM (.tif) (Optionnal for Co-registration, Interpolation and calculation)', size=(size_name_menu, 2)), sg.Input(key='-DEM1-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Latest DEM (.tif) (Optionnal for Interpolation and Calculation)', size=(size_name_menu, 2)), sg.Input(key='-DEM2-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Glacier outlines(.shp)', size=(size_name_menu, 1)), sg.Input(key='-vector glacier-'), sg.FileBrowse(file_types=(("Shape files", "*.shp"),("ALL Files","*.*"),))],
    [sg.Text('Coordinate system (ex: EPSG:32633)', size=(size_name_menu, 1)), sg.Input(key='-coordinate_system-')],
    [sg.Text('Vector frame shapefile (.shp)', size=(size_name_menu, 1)), sg.Input(key='-vector_frame-'),sg.FileBrowse(file_types=(("Shape files", "*.shp"),("ALL Files","*.*"),))],
    [sg.Text('Desired resolution:(ex:20)', size=(size_name_menu, 1)), sg.Input(key='-resolution-')],
    ]

tab_coregistration_layout=[
    [sg.Text('Co-register DEM 1 and DEM 2 (if filled out) with a reference DEM', size=(90, 3))],
    [sg.Checkbox('Select to co-register', key='-coregistration-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('Reference DEM for coregistration (if empty : latest DEM is used)', size=(size_name_menu, 1)), sg.Input(key='-reference_coreg-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Method for coregistration', size=(size_name_menu, 1)),sg.InputCombo(('Nuth and Kääb', 'ICP', 'ICP + Nuth and Kääb', 'Nuth and Kääb + Deramp',"Deramp + Nuth and Kääb",'Vertical offset + ICP + Nuth and Kääb',"Deramp(0) + ICP + Nuth and Kääb", "ICP + Nuth and Kääb + Deramp", "Deramp + ICP + Nuth and Kääb", "ICP + Deramp + Nuth and Kääb"), size=(35, 1), key='-method_coregistration-')]
    ]

tab_subtraction_layout=[
    [sg.Text('Calculate the following subtraction : dDEM= DEM 2 - DEM 1 , if coregistration is unticked, the entries DEM1 and DEM2 will be used', size=(90, 3))], 
    [sg.Checkbox('Select to subtract', default=False, key='-subtraction-')],
    [sg.Text(size=(90, 1))],
    [sg.Checkbox('Select to use the DEMs calculated in the co-registration steps', default=True, key='-use_dem_coreg_sub-')],
    ]

tab_interpolation_layout=[
    [sg.Text('Interpolate the following subtraction : dDEM= DEM 2 - DEM 1 , if coregistration is ticked, the co-registered DEM will be used', size=(90, 3))], 
    [sg.Checkbox('Select to interpolate', default=False, key='-interpolation-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('Reference DEM for interpolation (if empty : latest DEM is used)', size=(size_name_menu, 2)), sg.Input( justification='right',key='-reference_interpolation-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Method for interpolation', size=(size_name_menu, 1)),sg.InputCombo(('linear interpolation', 'hypsometric interpolation', 'local hypsometric interpolation'), size=(30, 1), key='-method_interpolation-')],
    [sg.Checkbox('Select to use the DEMs calculated in the co-registration step', default=True, key='-use_dem_coreg_inter-')],
    [sg.Text('dDEM already calculated (Optionnal)', size=(size_name_menu, 1)), sg.Input( justification='right',key='-ddem_calculated-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    ]

tab_errors_layout=[
    [sg.Text('Calculate non-stationarity of elevation measurment errors ', size=(100, 2))],
    [sg.Checkbox('Select to calculate unstandardized non-stationarity errors', default=False, key='-unstandardized_non_stationarity_errors-')],
    [sg.Checkbox('Select to calculate standardized non-stationarity errors', default=False, key='-standardized_non_stationarity_errors-')],
    [sg.Text('Reference DEM for this step (if empty : latest DEM is used)', size=(size_name_menu, 1)), sg.Input(key='-reference_standardized_non_stationarity_errors-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*")))],
    [sg.Text('dDEM for this step (Optionnal)', size=(size_name_menu, 1)), sg.Input(key='-ddem_standardized_non_stationarity_errors-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*")))],
    [sg.Checkbox('Select to calculate integrated and mean standardized error', default=False, key='-standardized_average_shapefile_errors-')],
    [sg.Text('Shapefile of the glacier where you want to know the integrated and mean error (if empty : glacier outline is used)', size=(size_name_menu, 2)), sg.Input(key='-outlines_standardized_non_stationarity_errors-'), sg.FileBrowse(file_types=(("Shape files", "*.shp"),("ALL Files","*.*")))],   
    ]

tab_calculation_layout=[
    [sg.Text('Calculate the volume and mass of one or several dDEM(s) previously calculated ', size=(100, 3))],
    [sg.Checkbox('Select to calculate volume and mass', default=False, key='-calculation_volume-')],
    [sg.Text('Ice density', size=(10, 1)), sg.Input(default_text=0.917,key='-rho_density-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('If the calculation volume entry is ticked (above), you can calculate volume and mass with your own dDEM layer, please fill out the entry below: ', size=(90, 2))],
    [sg.Text('Select your own DEM reference for subtraction (optionnal)', size=(size_name_menu, 2)), sg.Input(key='-ddem_volume-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Checkbox('Select to use the dDEM calculated in the subtraction step', default=True, key='-use_ddem_sub-')],
    [sg.Checkbox('Select to use the dDEM calculated in the interpolation step', default=True, key='-use_ddem_inter-')],
    [sg.Text('If several entries above are filled up, then all will be calculated', size=(90, 1))],
    [sg.Checkbox('Calculate error', default=False, key='-calculate_error-'),sg.Text('   Ice density error (ρ:kg/m³) :', size=(28, 1),key="-text_ice_density_error-"), sg.Input(key='-rho_error_density-',size=(8,1)),sg.Text('    Glacier outline error (m²) :', size=(25, 1),key="-text_glacier_outline_error-"), sg.Input(key='-glacier_outline_error-',size=(5,1))],
    ]

tab_out_entries_layout=[
    [sg.Text('Location where the files will be saved ', size=(size_name_menu, 1)), sg.Input(key='-path_out-'), sg.FolderBrowse("Select Path")],
    [sg.Text('Name of the PDF report', size=(size_name_menu, 1)), sg.Input(key='-name_pdf-')],
    ]


# %%% script layout

script_layout=[
    [sg.Text('Script output', size=(20, 1)),sg.VSeparator(), sg.Text('', size=(10, 1), key='-CPU-')],
    [sg.Output(size=(88, 20), font='Courier 10')],
    [sg.Button("Display output folder",visible=False, key="-display_output_folder-"),sg.Button("Display PDF file",visible=False, key="-display_pdf-")],
    ]

tab_function=[
    [sg.TabGroup([[sg.Tab('Required inputs', tab_main_layout), sg.Tab('Co-registration', tab_coregistration_layout), sg.Tab('Subtraction', tab_subtraction_layout), sg.Tab('Interpolation', tab_interpolation_layout),sg.Tab('Errors', tab_errors_layout), sg.Tab('Mass and volume calculation', tab_calculation_layout),sg.Tab('Output', tab_out_entries_layout)]],key="-tab_entries-")],
    [sg.Submit("Execute"),sg.Button("Next",enable_events=True) ] 
    ] 


# %%% Tab displayer


left_col_image = [
            [sg.Text('Folder'), sg.In(size=(25,1), enable_events=True ,key='-FOLDER-'), sg.FolderBrowse()],
            [sg.Text('Select one image to display it',key='-text_image-')],
            [sg.Listbox(values=[], enable_events=True, size=(40,5),key='-FILE LIST-')],
            [sg.Text('Resize to'), sg.In(key='-W-', size=(5,1)), sg.In(key='-H-', size=(5,1))]
            ]


images_col_image = [[sg.Text('Images')],
              [sg.Text(size=(40,1), key='-TOUT-')],
              [sg.Image(key='-IMAGE-')]]

left_col_graph = [
            [sg.Text('Select one graph to display it',key='-text_graph-')],
            [sg.Listbox(values=[], enable_events=True, size=(40,5),key='-GRAPH_LIST-')],
           ]

images_col_graph = [
                [sg.Text('Graph')],
              [sg.Text(size=(40,1), key='-TOUT-')],
              [sg.Canvas( key='-CANVAS-')]
              ]    

tab_image=[[sg.Column(left_col_image, element_justification='left'),sg.VSeperator(),sg.Column(images_col_image, element_justification='center')]]
tab_graph=[[sg.Column(left_col_graph, element_justification='left'),sg.VSeperator(),sg.Column(images_col_graph, element_justification='center')]]
image_eagloo=[[sg.Image(convert_to_bytes(f'{path_function}/eagloo_logo.png',resize=(300,300)))]]


# %%% Main tab
# beginning main layout
layout = [
    [sg.Menu(menu_def, tearoff=False, key='MENU_BAR')],    
    [sg.Column(tab_function, element_justification='c'), sg.VSeperator(),sg.Column(script_layout, element_justification='center')],
    [sg.HSeparator()], 
    [sg.TabGroup([[sg.Tab('Eagloo', image_eagloo,  key='-image_eagloo-',element_justification='c'),sg.Tab('Image', tab_image,key='-image_tab-',visible=False), sg.Tab('Graph', tab_graph,key='-graph_tab-',visible=False)]],key="-tab_output-",size=(1500,500))],
    ]
   


 # %% Backend of the frontend
window = sg.Window('Eagloo DEM Modification Software',layout , resizable=True, finalize=True, element_justification='center', location=(0,0))
# ---===--- Loop taking in user input and using it to call scripts --- #


tab_keys_out = ('-image_tab-','-graph_tab-')
list_tab=['Required inputs', 'Co-registration','Subtraction','Interpolation','Mass and volume calculation','Errors','Output']

figure_agg = None
calculate_volume_change_error=False
try:
      
    print(f"Running python : {sys.version}--version :{sys.version_info}")
    # image_folder=False
    
    while True:
        
        # if image_folder==False:
        #     image_printed=convert_to_bytes(f'{path_function}/eagloo_logo.png', resize=(150,150))
        #     window['-TOUT-'].update('Select a file from the list to display here:')
        #     window['-IMAGE-'].update(data=image_printed)    
        
        cpu_percent = psutil.cpu_percent(interval=1)
        window['-CPU-'].update(f'CPU : {cpu_percent:02.0f}%')
        
        event, values = window.read()
        
            
        if event == 'EXIT'  or event == sg.WIN_CLOSED:
            break # exit button clicked
    
        elif event in (sg.WIN_CLOSED, 'Exit'):
            break
        #Beginning menu
        elif event == 'Open':
            filename = sg.popup_get_file('Will not see this message', no_window=True)
            sg.popup('You selected:', filename)

        elif event == 'Save as':
            filename = sg.popup_get_file('Name of the file',save_as=True)
            sg.popup('You saved:', filename)
        elif event == 'About':
             sg.Popup('Eagloo DEM','Version 1.0')
        
        
       
        
        if event=="Next":
            index=list_tab.index(window.Element('-tab_entries-').get())
            # window.Element('-tab_entries-').selectab()
            if index<6:
                window['-tab_entries-'].Widget.select(index+1)
            else:
                window['-tab_entries-'].Widget.select(0)

            

 # %%% Execute Treatment_dem
        if event=="Execute":
            if values["-path_out-"]=="":
                window['-tab_entries-'].Widget.select(6)
                message_display(' You have to select a location to save your files.')
                
                
            else:
                
                if values["-standardized_average_shapefile_errors-"] or values["-calculate_error-"]:
                    calculate_volume_change_error=True
                
                try: 
                    error_rho_density=float(values["-rho_error_density-"])
                except: 
                    error_rho_density=0
                try:
                    error_glacier_outline=float(values["-glacier_outline_error-"])
                except: 
                    error_glacier_outline=0
                try:
                    
                    plot_function=tdem.Treatment_dem(
                    dem_1=f'{values["-DEM1-"]}', 
                    dem_2=f'{values["-DEM2-"]}',
                    vector_glacier=f'{values["-vector glacier-"]}',
                    path_out=f'{values["-path_out-"]}', 
                    coordinate_system=f'{values["-coordinate_system-"]}',
                    vector_frame=f'{values["-vector_frame-"]}',
                    resolution=float(values["-resolution-"]),
                    coregistration=values["-coregistration-"],
                    dem_coreg_ref=f'{values["-reference_coreg-"]}',
                    method_coregistration=f'{values["-method_coregistration-"]}',
                    subtraction=values["-subtraction-"],
                    use_dem_coreg_sub=values["-use_dem_coreg_sub-"],
                    interpolation=values["-interpolation-"],
                    use_dem_coreg_inter=values["-use_dem_coreg_inter-"],
                    method_interpolation=f'{values["-method_interpolation-"]}',
                    dem_inter_ref=f'{values["-reference_interpolation-"]}',
                    ddem_calculated=f'{values["-ddem_calculated-"]}',
                    rho_density=float(values["-rho_density-"]),
                    calculation_mass=values["-calculation_volume-"],
                    unstandardized_nonstationarity_error=values["-unstandardized_non_stationarity_errors-"],
                    standardization_nonstationarity_error=values["-standardized_non_stationarity_errors-"],
                    standardized_average_shapefile_error=calculate_volume_change_error,
                    ref_standardization_nonstationarity_error=f'{values["-reference_standardized_non_stationarity_errors-"]}',
                    ddem_standardization_nonstationarity_error=f'{values["-ddem_standardized_non_stationarity_errors-"]}',
                    glacier_outlines_error=f'{values["-outlines_standardized_non_stationarity_errors-"]}',
                    ddem_volume=f'{values["-ddem_volume-"]}',
                    use_ddem_sub=values["-use_ddem_sub-"],
                    use_ddem_inter=values["-use_ddem_inter-"],  
                    error_rho_density=error_rho_density,
                    error_glacier_outline=error_glacier_outline,
                    name_pdf=values["-name_pdf-"]
                    ) 
                    
                      
                    print(f'{plot_function}') 
                     
                    
                    window["-display_pdf-"].update(visible=True) #display the button to display PDF file
                    window["-display_output_folder-"].update(visible=True)
                      
                    # Print graph in image window and popup
                    #conditions to display functions
                    Name_plot=[]
                    try:
                        for i in range(0,len(plot_function)):
                            
                            if i==0 and plot_function[i]:
                                Name_plot.append('Co-registration DEM 1')
                                index_display=i
                            elif i==1 and plot_function[i]:
                                Name_plot.append('Co-registration DEM 2')
                                index_display=i
                            elif i==2 and plot_function[i]:
                                Name_plot.append('Subtraction DEM 2 - DEM 1')
                                index_display=i
                            elif i==3 and plot_function[i]:
                                Name_plot.append('Interpolation')
                                index_display=i
                            elif i==4 and plot_function[i]:
                                Name_plot.append('Graph interpolation')
                                index_display=i
                    
                        window["-GRAPH_LIST-"].update(Name_plot) 
                        try:
                            #display the first graph generated
                            window[tab_keys_out[1]].update(visible=True)
                            fig=plot_function[index_display]
                            plt.show()
                            draw_figure(window['-CANVAS-'].TKCanvas, fig)
                            window[tab_keys_out[1]].select()
                            #display image tab
                            window['-image_tab-'].update(visible=True) 
                            #Hide eagloo image
                            window['-image_eagloo-'].update(visible=False) 
                            
                            #update image to display
                            folder=values["-path_out-"]
                            try:
                                file_list = os.listdir(folder)         # get list of files in folder
                            except:
                                file_list = []
                            fnames = [f for f in file_list if os.path.isfile(
                                os.path.join(folder, f)) and f.lower().endswith((".png", ".jpg", "jpeg", ".tiff", ".bmp"))]
                            window['-FILE LIST-'].update(fnames)
                            window['-FOLDER-'].update(folder)
                        except:
                            pass
                        sg.SystemTray.notify('Eagloo', 'Your DEMs have been processed')
                    except Exception as E:
                        message_display(f' Error {E}')
                        pass 
                    
                except Exception as E:
                    message_display(f' Error {E}')
                    pass 
        # %%% display pdf and output buttons    
        
        #display pdf button
        if event=="-display_pdf-":
            path=f'{values["-path_out-"]}'
            path=path.replace("/","\\")
            path=path.lstrip("\\")
            name_pdf=f'report_{values["-name_pdf-"]}.pdf'
            subprocess.Popen(["explorer.exe",'\\{}\"{}"'.format(path,name_pdf)])
            
        #display output folder button
        if event=="-display_output_folder-":
            path=f'{values["-path_out-"]}'
            path=path.replace('/','\\')
            path=path.lstrip("\\")
            print(path)
            subprocess.Popen(["explorer.exe", '\\' + r''+ path])
            
        
  # %%% Update graph list  
    
        #Last step to print figure image and popup figure
        if event=="-GRAPH_LIST-":
            if figure_agg:
                # ** IMPORTANT ** Clean up previous drawing before drawing again
                delete_figure_agg(figure_agg)
            
            try:
                # get first listbox item chosen (returned as a list)
                choice = values['-GRAPH_LIST-'][0]
                if choice=='Co-registration DEM 1':
                    fig=plot_function[0]
                    plt.show()
                    figure_agg = draw_figure(
                        window['-CANVAS-'].TKCanvas, fig)  # draw the figure
                elif choice=='Co-registration DEM 2':
                    fig=plot_function[1]
                    plt.show()
                    figure_agg = draw_figure(
                        window['-CANVAS-'].TKCanvas, fig)  # draw the figure
                elif choice=='Subtraction DEM 2 - DEM 1':
                     fig=plot_function[2]
                     plt.plot()
                     plt.show()
                     figure_agg = draw_figure(
                          window['-CANVAS-'].TKCanvas, fig)  # draw the figure+
                elif choice=='Interpolation':
                    fig=plot_function[3]
                    plt.show()
                    figure_agg = draw_figure(
                        window['-CANVAS-'].TKCanvas, fig)  # draw the figure
                elif choice=="Graph interpolation":
                    fig=plot_function[4]
                    plt.show()
                    figure_agg = draw_figure(
                        window['-CANVAS-'].TKCanvas, fig)  # draw the figure
                
            except Exception as E:
                message_display(f' Error {E}')  
                pass 
         
        # %%% Update image list       
        #window to print image
        
        if event == '-FOLDER-':                         # Folder name was filled in, make a list of files in the folder
            folder = values['-FOLDER-']
            try:
                file_list = os.listdir(folder)         # get list of files in folder
            except:
                file_list = []
            fnames = [f for f in file_list if os.path.isfile(
                os.path.join(folder, f)) and f.lower().endswith((".png", ".jpg", "jpeg", ".tiff", ".bmp"))]
            window['-FILE LIST-'].update(fnames)
        elif event == '-FILE LIST-':    # A file was chosen from the listbox
            try:
                filename = os.path.join(values['-FOLDER-'], values['-FILE LIST-'][0])
                window['-TOUT-'].update(filename)
                if values['-W-'] and values['-H-']:
                    new_size = int(values['-W-']), int(values['-H-'])
                else:
                    new_size = None
                image_printed=convert_to_bytes(filename, resize=new_size)
                image_folder=True
                window['-IMAGE-'].update(data=image_printed)
                
            except Exception as E:
                print(f'Error: {E}')
                
                pass        # something weird happened making the full filename
        else:
             image_folder=False
        
        # Print graph with matplotlib
           
except Exception as E:
    message_display(f' Error {E}') 
    pass 

window.close()
