<head>
  <link rel="stylesheet" type="text/css" href="/styles_readme.css">
</head>

<h1 align="center">DEM modification software</h1>

<p align="center">
<img src="/Eagloo_frontend_functions/eagloo_logo.png" height="300px" style="text-align:center;">
</p>
  
## 1 | Introduction
Eagloo is a package designed for simplified the modification and analysis of Digital Elevation Models (DEMs). It is made to co-register, subtract, interpolate, compute errors, calculate mass and volume, and automatically generate a synthetic PDF.
It is a tool for researchers, scientists, and anyone working with DEMs.
  
## 2 | Features
This software offers a variety of features to help you analyse your DEMs, including:
  
### - Co-registration: 
Eagloo offers multiple co-registration methods, including Nuth and Kääb, Iterative Closest Point (ICP), and ICP + Nuth and Kääb.

### - Subtraction: 
Easily subtract one DEM from another to analyze changes in terrain over time.

### - Interpolation: 
Eagloo's interpolation capabilities include:
  1. Linear interpolation
  2. Hypsometric interpolation
  3. Local hypsometric interpolation

### - Errors:
Errors can be calculated through those different methods: 
  1. Unstandardized non-stationnarity errors 
  2. Standardized non-stationnarity errors
  3. Integrated and mean standardized error

### - Mass and volume calculation: 
Calculate the mass and volume of specific areas within your DEMs.

### - Synthetic PDF generation: 
Eagloo can automatically generate a synthetic PDF to help you quickly and easily visualize your analysis results.
You can click [here](/Eagloo_frontend_functions/github_assets/example_report_2013_2017_nuth_kaab.pdf) to see an example of a pdf.

## 3 | Presentation of the Graphical User Interface
<p align="center">
<img src="/Eagloo_frontend_functions/github_assets/eagloo_frontend.jpg" height="500px" style="text-align:center;">
</p>

## 4 | Installation
Eagloo is only available on Linux.
  1. Clone this repository: git clone https://github.com/Geomatyx/Eagloo_software.git
  2. Install dependencies, right-click on "installing_dependencies.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")
  3. Open the software, right-click on "Eagloo_software.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")

## 5 | Work carried out 
The work is consisting to:
- Create a user interface that allows the use of specific functions for processing and analyzing digital elevation models (DEMs) without any knowledge of code
- Improve the efficiency of DEM processing by using a pipeline approach
- Create of a function for automatically generate reports

## 6 | External Open Source ressources used
Various external libraries and packages were utilized, with the most frequently used resources presented below. A comprehensive list of all libraries used can be found [here](#9 | Recognition of Open Source use).

### Frontend
The user interface was developed using the [PySimpleGUI] (https://www.pysimplegui.org/en/latest/) package (under GNU GPL license).

### Processing data
To modify certain input and output data, the [GDAL](https://github.com/OSGeo/gdal) framework was utilized (under MIT license).

### Co-registration, Interpolation, Errors
The functions used for co-registration, interpolation, and calculating errors of the digital elevation models (DEMs) are from the [xDEM](https://github.com/GlacioHack/xdem) package (under MIT license).

### PDF 
To automatically generate a PDF report, the [FPDF](http://www.fpdf.org/) library was utilized(no usage restriction). 

## 7 | Recognition of Open Source use

In the Eagloo Repository these external open source packages were used at least one time. 

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
    
## 8 | License
This repository is open-source and is licensed under the MIT license. See the [LICENSE](/LICENSE) file for details.

## 9 | Contributing
To contribute to this project, please follow these steps:

  1. Fork this repository
  2. Create a new branch: git checkout -b feature-branch
  3. Make your changes and commit them: git commit -am 'Add some feature'
  4. Push to the branch: git push origin feature-branch
  5. Create a pull request 






