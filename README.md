<head>
  <link rel="stylesheet" type="text/css" href="/styles_readme.css">
</head>

<h1 align="center">DEM modification software</h1>

<p align="center">
<img src="/Eagloo_frontend_functions/eagloo_logo.png" width="200px" height="200px" style="text-align:center;">
</p>

## 1 | Introduction
This software is designed for simplified some Digital Elevation Models (DEMs) modelisation. It is made to co-register, subtract, interpolate, calculate mass and volume, and automatically generate a synthetic PDF. Eagloo is a tool for researchers, scientists, and anyone working with DEMs.

## 2 | Purposes of this software  
- Create a user interface that allows the use of specific functions for processing and analyzing digital elevation models (DEMs) without any knowledge of code.
- Improve the efficiency of DEM processing by providing an easier method for selecting input data using a pipeline approach.
- Automatically generate reports.

## 3 | External Open Source ressources used
Various external libraries and frameworks were utilized, with the most frequently used resources presented below. A comprehensive list of all libraries used can be found at the bottom of the page.

### Frontend
The user interface was developed using the [PySimpleGUI](https://www.pysimplegui.org/en/latest/)(under GNU GPL license) library.

### Processing data
To modify certain input and output data, the [GDAL](https://github.com/OSGeo/gdal)(under MIT license) framework was utilized.

### Co-registration, Interpolation, Errors
The functions used for co-registration, interpolation, and calculating errors of the digital elevation models (DEMs) are from the [xDEM](https://github.com/GlacioHack/xdem)(under MIT license).

### PDF 
To automatically generate a PDF report, the [FPDF](http://www.fpdf.org/)(no usage restriction) library was utilized. 

## 4 | Features
Eagloo offers a variety of features to help you analyse your DEMs, including:

### - Co-registration: 
Eagloo offers multiple co-registration methods, including Nuth and Kääb, Iterative Closest Point (ICP), and ICP + Nuth and Kääb.

### - Subtraction: 
Easily subtract one DEM from another to analyze changes in terrain over time.

### - Interpolation: 
Eagloo's interpolation capabilities include:
  1. Linear interpolation
  2. Hypsometric interpolation
  3. Local hypsometric interpolation

### - Mass and volume calculation: 
Calculate the mass and volume of specific areas within your DEMs.

### - Synthetic PDF generation: 
Eagloo can automatically generate a synthetic PDF to help you quickly and easily visualize your analysis results.

## 5 | Installation
Currently, Eagloo is only available on Linux.
  1. Clone this repository: git clone https://github.com/Geomatyx/Eagloo_software.git
  2. Install dependencies, right-click on "installing_dependencies.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")
  3. Open the software, right-click on "Eagloo_software.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")

## 6 | Contributing
To contribute to this project, please follow these steps:

  1. Fork this repository
  2. Create a new branch: git checkout -b feature-branch
  3. Make your changes and commit them: git commit -am 'Add some feature'
  4. Push to the branch: git push origin feature-branch
  5. Create a pull request 

## 7 | License
This repository is open-source and is licensed under the MIT license. See the [LICENSE](/LICENSE) file for details.

## 8 | Recognition of Open Source use

In the Eagloo Repository these external open source packages were used at least one time. 
If you use Open Source software in your project, be sure and supply information about the packages you used.

    
    attrs
    Brotli
    ConfigParser
    cryptography
    Cython
    dl
    docutils
    Fiona
    fpdf
    GDAL
    geopandas
    geoutils
    HTMLParser
    importlib_metadata
    ipython
    ipywidgets
    Jinja2
    keyring
    lockfile
    lxml
    matplotlib
    numpy
    opencv_python
    ordereddict
    pandas
    Pillow
    protobuf
    psutil
    pygeotools
    pyOpenSSL
    PySimpleGUI
    pytransform3d
    railroad
    rasterio
    richdem
    scikit_image
    scipy
    Sphinx
    toml
    tornado
    tqdm
    xdem








