#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 17:48:29 2022

@author: guillaume
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
import psutil #CPU

import matplotlib       #print
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
import inspect
import PySimpleGUI as sg


        



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
sg.set_options(font=font)

sg.LOOK_AND_FEEL_TABLE['Native'] = {
    'BACKGROUND': '#5D6D7E',
    'TEXT': '#F2F3F4',
    'INPUT': '#ABB2B9',
    'FONT': 'Arial 20',
    'TEXT_INPUT': '#17202A',
    'SCROLL': '#a5a4a4',
    'BUTTON': ("#efe", "#E59866"),
    'PROGRESS': sg.DEFAULT_PROGRESS_BAR_COLOR,
    'BORDER': 1,
    'SLIDER_DEPTH': 0,
    'PROGRESS_DEPTH': 0,
    'ACCENT1': '#ff5414',
    'ACCENT2': '#33a8ff',
    'ACCENT3': '#dbf0ff'}

sg.theme('Native')



#definition function

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
        except Exception as e:
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


# --------------------------------- Define Layout ---------------------------------



#end definition fonction

menu_def = [['File', ['Open', 'Save','Save as', 'Exit']], ['Help', 'About'],  ] 
 

tab_main_layout=[
    [sg.Text('Please fill out the following entries:', size=(90, 2))],
    [sg.Text('Earliest DEM (.tif) (Optionnal for Co-registration, Interpolation and calculation)', size=(size_name_menu, 2)), sg.Input(key='-DEM1-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Latest DEM (.tif) (Optionnal for Interpolation and Calculation)', size=(size_name_menu, 2)), sg.Input(key='-DEM2-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Glacier outlines(.shp)', size=(size_name_menu, 1)), sg.Input(key='-vector glacier-'), sg.FileBrowse(file_types=(("Shape files", "*.shp"),("ALL Files","*.*"),))],
    [sg.Text('Location where files will be saved ', size=(size_name_menu, 1)), sg.Input(key='-path_out-'), sg.FolderBrowse("Select Path")],
    [sg.Text('Coordinate system (ex: EPSG:32633)', size=(size_name_menu, 1)), sg.Input(key='-coordinate_system-')],
    [sg.Text('Vector frame shapefile (.shp)', size=(size_name_menu, 1)), sg.Input(key='-vector_frame-'),sg.FileBrowse(file_types=(("Shape files", "*.shp"),("ALL Files","*.*"),))],
    [sg.Text('Desired resolution:(ex:20)', size=(size_name_menu, 1)), sg.Input(key='-resolution-')],
    
    ]
tab_coregistration_layout=[
    [sg.Text('Co-register DEM 1 and DEM 2 (if filled out) with a reference DEM', size=(90, 3))],
    [sg.Checkbox('Do you wish to co-register?', key='-coregistration-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('If the co-registration entry is ticked (above), please fill out the entries below:')],
    [sg.Text('Reference DEM for coregistration', size=(size_name_menu, 1)), sg.Input(key='-reference_coreg-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Method for coregistration', size=(size_name_menu, 1)),sg.InputCombo(('Nuth and Kaab', 'ICP', 'ICP + Nuth and Kaab'), size=(20, 1), key='-method_coregistration-')]
    ]

tab_subtraction_layout=[
    [sg.Text('Calculate the following subtraction : dDEM= DEM 2 - DEM 1 , if coregistration is unticked, the entries DEM1 and DEM2 will be used', size=(90, 3))], 
    [sg.Checkbox('Do you wish to subtract?', default=False, key='-subtraction-')],
    [sg.Text(size=(90, 1))],
    [sg.Checkbox('Do you wish to use the DEMs calculated in the co-registration steps?', default=True, key='-use_dem_coreg_sub-')],
    
    ]
tab_interpolation_layout=[
    [sg.Text('Interpolate the following subtraction : dDEM= DEM 2 - DEM 1 , if coregistration is ticked, the co-registered DEM will be used', size=(90, 3))], 
    [sg.Checkbox('Do you wish to interpolate?', default=False, key='-interpolation-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('If the interpolation entry is ticked (above), please fill out the entries below:', size=(90, 2))],
    [sg.Text('Reference DEM for interpolation', size=(size_name_menu, 2)), sg.Input( justification='right',key='-reference_interpolation-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Text('Method for interpolation', size=(size_name_menu, 1)),sg.InputCombo(('linear interpolation', 'hypsometric interpolation', 'local hypsometric interpolation'), size=(30, 1), key='-method_interpolation-')],
    [sg.Checkbox('Do you wish to use the DEMs calculated in the co-registration step?', default=True, key='-use_dem_coreg_inter-')],
    [sg.Text('dDEM already calculated (Optionnal)', size=(size_name_menu, 1)), sg.Input( justification='right',key='-ddem_calculated-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    
    ]
tab_calculation_layout=[
    [sg.Text('Calculate the volume and mass of one or several dDEM(s) previously calculated ', size=(90, 3))],
    [sg.Checkbox('Do you wish to calculate volume and mass?', default=False, key='-calculation_volume-')],
    [sg.Text(size=(90, 1))],
    [sg.Text('If the calculation volume entry is ticked (above), you can calculate volume and mass with your own dDEM layer, please fill out the entry below: ', size=(90, 2))],
    [sg.Text('Select DEM reference for subtraction (optionnal, necessary if you want to calucalte volume and mass)', size=(size_name_menu, 2)), sg.Input(key='-ddem_volume-'), sg.FileBrowse(file_types=(("TIF files", "*.tif"),("ALL Files","*.*"),))],
    [sg.Checkbox('Do you wish to use the dDEM calculated in the subtraction step?', default=True, key='-use_ddem_sub-')],
    [sg.Checkbox('Do you wish to use the dDEM calculated in the interpolation step?', default=True, key='-use_ddem_inter-')],
    [sg.Text('If several entries above are filled up, then all will be calculated', size=(90, 1))],
    ]

tab_function=[
    [sg.TabGroup([[sg.Tab('Required inputs', tab_main_layout, tooltip='tip'), sg.Tab('Co-registration', tab_coregistration_layout), sg.Tab('Subtraction', tab_subtraction_layout), sg.Tab('Interpolation', tab_interpolation_layout), sg.Tab('Mass and volume calculation', tab_calculation_layout)]])],
    [sg.Submit("Execute"), ]
    ]

script_layout=[
    [sg.Text('Script output', size=(20, 1)),sg.VSeparator(), sg.Text('', size=(10, 1), key='-CPU-')],
    [sg.Output(size=(88, 20), font='Courier 10')],
    [sg.Cancel(), sg.Button('EXIT')],
    ]

# First the window layout...2 columns

left_col_image = [[sg.Text('Select one graph to display it')],
            [sg.Listbox(values=[], enable_events=True, size=(40,5),key='-GRAPH_LIST-')],
            [sg.Checkbox('Display interactive graph', default=False,enable_events=True, key='-interactive_graph-')],
            [sg.HSeparator()],
            [sg.Text('Folder'), sg.In(size=(25,1), enable_events=True ,key='-FOLDER-'), sg.FolderBrowse()],
            [sg.Listbox(values=[], enable_events=True, size=(40,5),key='-FILE LIST-')],
            [sg.Text('Resize to'), sg.In(key='-W-', size=(5,1)), sg.In(key='-H-', size=(5,1))]]

# For now will only show the name of the file that was chosen
images_col_image = [[sg.Text('Images')],
              [sg.Text(size=(40,1), key='-TOUT-')],
              [sg.Image(key='-IMAGE-')]]
images_col_graph = [[sg.Text('Graph')],
              [sg.Text(size=(40,1), key='-TOUT-')],
              [sg.Canvas( key='-CANVAS-')]]
    


tab_image_graph=[
    [sg.TabGroup([[sg.Tab('Image', images_col_image, tooltip='tip', element_justification="center"), sg.Tab('Graph', images_col_graph, element_justification="center")]],size=(1000,500) )]
    ]
# =============================================================================
# beginning main layout
layout = [
    
    [sg.Menu(menu_def, tearoff=False, key='MENU_BAR')],    
    [sg.Column(tab_function, element_justification='c'), sg.VSeperator(),sg.Column(script_layout, element_justification='c')],
    [sg.HSeparator()],
   
    [sg.Column(left_col_image, element_justification='left'), sg.VSeperator(),sg.Column(tab_image_graph, element_justification='center')],
    
    
    ]
    
    
#end main layout
# =============================================================================



window = sg.Window('Eagloo: Treatments for DEM',layout ,icon=f'{path_function}/eagloo_logo.png', resizable=True, finalize=True, element_justification='center', location=(0,0))
# ---===--- Loop taking in user input and using it to call scripts --- #
window.Maximize()
figure_agg = None
try:
    


       
    print(f"Running python : {sys.version}--version :{sys.version_info}")
    image_folder=False
    
    
    while True:
        
        
        if image_folder==False:
            image_printed=convert_to_bytes(f'{path_function}/eagloo_logo.png', resize=(150,150))
            window['-TOUT-'].update('Select a file from the list to display here:')
            window['-IMAGE-'].update(data=image_printed)    
        
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
        elif event==('-interactive_graph-'):
            if values['-interactive_graph-']==True:
                matplotlib.use('Tkagg')
        
            
        
        #End Menu
        
        
        
        
        
        # launch treatment dem function
        
        
        if event=="Execute":
            
            try:
                
                plot_function=tdem.Treatment_dem(
                dem_1=f'{values["-DEM1-"]}', 
                dem_2=f'{values["-DEM2-"]}',
                vector_glacier=f'{values["-vector glacier-"]}',
                path_out=f'{values["-path_out-"]}', 
                coordinate_system=f'{values["-coordinate_system-"]}',
                vector_frame=f'{values["-vector_frame-"]}',
                resolution=values["-resolution-"],
                coregistration=values["-coregistration-"],
                dem_coreg_ref=f'{values["-reference_coreg-"]}',
                method_coregistration=f'{values["-method_coregistration-"]}',
                subtraction=values["-subtraction-"],
                use_dem_coreg_sub=values["-use_dem_coreg_sub-"],
                interpolation=values["-interpolation-"],
                use_dem_coreg_inter=values["-use_dem_coreg_inter-"],
                method_interpolation=f'{values["-method_interpolation-"]}',
                dem_hypso_ref=f'{values["-reference_interpolation-"]}',
                ddem_calculated=f'{values["-ddem_calculated-"]}',
                calculation_mass=values["-calculation_volume-"],
                ddem_volume=f'{values["-ddem_volume-"]}',
                use_ddem_sub=values["-use_ddem_sub-"],
                use_ddem_inter=values["-use_ddem_inter-"],
                ) 
                
                
                print(f'{plot_function}') 
                
                # Print graph in image window and popup
                #conditions to display functions
                Name_plot=[]
                try:
                    for i in range(0,len(plot_function)):
                        
                        if i==0 and plot_function[i]:
                            Name_plot.append('Co-registration DEM 1')
                        elif i==1 and plot_function[i]:
                            Name_plot.append('Co-registration DEM 2')
                        elif i==2 and plot_function[i]:
                            Name_plot.append('Subtraction DEM 2 - DEM 1')
                        elif i==3 and plot_function[i]:
                            Name_plot.append('Interpolation')
                        elif i==4 and plot_function[i]:
                            Name_plot.append('Graph interpolation')
                            
                   
                    window["-GRAPH_LIST-"].update(Name_plot) 
                except Exception as E:
                    sg.Popup(f' Error : {E}')  
                    pass 
                
            except Exception as E:
                sg.Popup(f' Error : {E}')  
                pass 
            
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
                sg.Popup(f' Error : {E}')  
                pass 
         
            
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
    sg.Popup(f' Error : {E}')  
    pass 
   
    


window.close()