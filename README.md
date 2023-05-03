<head>
  <link rel="stylesheet" type="text/css" href="/styles_readme.css">
</head>

<h1 align="center">DEM modification software</h1>

<p align="center">
<img src="/Eagloo_frontend_functions/eagloo_logo.png" height="300px" style="text-align:center;">
</p>
<br><br><br>

## 1 | Introduction
Eagloo is a package is designed to simplify the modification and analysis of Digital Elevation Models (DEMs). It is made to co-register, subtract, interpolate, compute errors, calculate mass and volume, and automatically generate a synthetic PDF.
It is a tool for researchers, scientists, and anyone working with DEMs.
<br><br><br>

## 2 | Features
This software offers a variety of features to help you analyse your DEMs, including:
  <br><br> 
  
### 2.1 | Co-registration: 
Eagloo offers multiple co-registration methods, including Nuth and Kääb, Iterative Closest Point (ICP), and ICP + Nuth and Kääb.
<br><br>

### 2.2 | Subtraction: 
Easily subtract one DEM from another to analyze changes in terrain over time.
<br><br>

### 2.3 | Interpolation: 
Eagloo's interpolation capabilities include:
  1. Linear interpolation
  2. Hypsometric interpolation
  3. Local hypsometric interpolation
<br><br>

### 2.4 | Errors:
Errors can be calculated through the following methods: 
  1. Unstandardized non-stationnarity errors 
  2. Standardized non-stationnarity errors
  3. Integrated and mean standardized error
<br><br>

### 2.5 | Mass and volume calculation: 
Calculate the mass and volume of specific areas within your DEMs.
<br><br>

### 2.6 | Synthetic PDF generation: 
Eagloo can automatically generate a synthetic PDF to help you efficiently visualize results of your analysis.
You can click [here](/Eagloo_frontend_functions/github_assets/example_report_2013_2017_nuth_kaab.pdf) to see an example of a pdf report.
<br><br><br>

## 3 | Presentation of the Graphical User Interface
<p align="center">
<img src="/Eagloo_frontend_functions/github_assets/eagloo_frontend.jpg" height="500px" style="text-align:center;">
</p>
<br><br><br>

## 4 | Installation
Eagloo is only available on Linux.
  1. Clone this repository: git clone https://github.com/Geomatyx/Eagloo_software.git
  2. Install dependencies, right-click on "installing_dependencies.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")
  3. Open the software, right-click on "Eagloo_software.sh" and click on "run as program" (if "run as program" is not displayed: Preferences -> Permissions -> "Allow executing file as program")
<br><br><br>

## 5 | Work carried out 
The work which was carried out there is:
- Create a user interface that allows the use of specific functions for processing and analyzing digital elevation models (DEMs) without any knowledge of code
- Improve the efficiency of DEM processing by using a pipeline approach
- Create of a function for automatically generate reports
<br><br><br>

## 6 | External Open Source ressources used
Various external libraries and packages were utilized, with the most frequently used resources presented below. A comprehensive list of all libraries used can be found [here](#7--recognition-of-open-source-use).
<br><br>

### 6.1 | Frontend
The user interface was developed using the [PySimpleGUI](https://www.pysimplegui.org/en/latest/) package (under GNU GPL license).
<br><br>

### 6.2 | Processing data
To modify certain input and output data, the [GDAL](https://github.com/OSGeo/gdal) framework was utilized (under MIT license). <br>
Source: <br>
GDAL/OGR contributors (2022). GDAL/OGR Geospatial Data Abstraction
  software Library. Open Source Geospatial Foundation. URL https://gdal.org <br>
  [![Zenodo](https://zenodo.org/badge/doi/10.5281/zenodo.5884351.svg)](https://zenodo.org/record/5884351)
<br><br>

### 6.3 | Co-registration, Interpolation, Errors
The functions used for co-registration, interpolation, and calculating errors of the digital elevation models (DEMs) are from the [xDEM](https://github.com/GlacioHack/xdem) package (under MIT license).<br>
Source: <br>
xDEM: [![Zenodo](https://zenodo.org/badge/doi/10.5281/zenodo.4809697.svg)](https://zenodo.org/record/4809698)

Related study:

- **Coregistration**:
  - Horizontal shift from aspect/slope relationship of *[Nuth and Kääb (2011)](https://doi.org/10.5194/tc-5-271-2011)*,
  - Iterative closest point (ICP) of *[Besl and McKay (1992)](http://dx.doi.org/10.1109/34.121791)*,
- **Bias correction**:
  - Along-track multi-sinusoidal noise by basin-hopping of *[Girod et al. (2017)](https://doi.org/10.3390/rs9070704)*,
- **Uncertainty analysis**:
  - Heteroscedasticity and multi-range correlations from stable terrain of *[Hugonnet et al. (2022)](https://doi.org/10.1109/JSTARS.2022.3188922)*,
- **Terrain attributes**:
  - Slope, aspect and hillshade of either *[Horn (1981)](http://dx.doi.org/10.1109/PROC.1981.11918)* or *[Zevenbergen and Thorne (1987)](http://dx.doi.org/10.1002/esp.3290120107)*,
  - Profile, plan and maximum curvature of *[Zevenbergen and Thorne (1987)](http://dx.doi.org/10.1002/esp.3290120107)*,
  - Topographic position index of *[Weiss (2001)](http://www.jennessent.com/downloads/TPI-poster-TNC_18x22.pdf)*,
  - Terrain ruggedness index of either *[Riley et al. (1999)](http://download.osgeo.org/qgis/doc/reference-docs/Terrain_Ruggedness_Index.pdf)* or *[Wilson et al. (2007)](http://dx.doi.org/10.1080/01490410701295962)*,
  - Roughness of *[Dartnell (2000)](http://dx.doi.org/10.14358/PERS.70.9.1081)*,
  - Rugosity of *[Jenness (2004)](https://doi.org/10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2)*,
  - Fractal roughness of *[Taud et Parrot (2005)](https://doi.org/10.4000/geomorphologie.622)*.

<br><br>

### 6.4 | PDF 
To automatically generate a PDF report, the [FPDF](http://www.fpdf.org/) library was utilized(no usage restriction). 
<br><br><br>

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
<br><br><br>   

## 8 | License
This repository is open-source and is licensed under the MIT license. See the [LICENSE](/LICENSE) file for details.
<br><br><br>

## 9 | Contributing
To contribute to this project, please follow these steps:

  1. Fork this repository
  2. Create a new branch: git checkout -b feature-branch
  3. Make your changes and commit them: git commit -am 'Add some feature'
  4. Push to the branch: git push origin feature-branch
  5. Create a pull request 






