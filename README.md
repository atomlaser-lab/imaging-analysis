# Imaging Analysis

This repository contains code for analyzing absorption images of cold atomic samples.  Using the standard method of illuminating an atomic sample with near-resonant light and recording the intensity pattern and then comparing that to a reference image with no atoms present, it constructs a map of optical depth (OD) and uses that for determinining properties of the sample.

## Use

The analysis code is broken into four main classes which can be invoked using a wrapper file.  These classes are
  - `RawImageData`: describes the raw image data, and has methods for loading that data from a directory
  - `AtomImageConstants`: describes various constants for defining an absorption image, such as detuning, pixel size, etc.
  - `AtomCloudFit`: mainly provides methods for fitting image data to various distributions.
  - `AbsorptionImage`: the over-arching class that describes the image.  Has three "constant" properties `raw`, `constants`, and `fitdata` which are each instances of `RawImageData`, `AtomImageConstants`, and `AtomCloudFit`, respectively.  Other properties describe the image and the sample, such as peak OD, number of atoms, temperature, etc.

Below is a bare-bones example of how one can read an image set corresponding to a single absorption image and analyze that image:
```
c = AbsorptionImage;
c.raw.load('filenames','last','directory','.','index',1);
c.constants.set('tof',30e-3,'exposuretime',30e-6,'detuning',0);
c.fitdata.set('roirow',[400,600],'roicol',[400,600],'fittype','gauss2d');
c.makeImage();
c.fit();
c.plotAllData([0,3]);
```
The first line creates a blank `AbsorptionImage` object `c`.  The second line loads the last set of images in the nominated directory (here the current directory `.`); last in this context is ordered according to time.  The index specifies how many sets to go back from the last image: a value of 1 is the last image set taken, a value of 2 is the second-to-last image, and so on.

The third line sets some constants associated with the image; in this case, we set the time-of-flight `tof` to 30 ms, the exposure time to 30 us, and the detuning to 0 MHz.  The fourth line sets the "region-of-interest" (ROI) to be the rectangular region with corners at (400,400), (400,600), (600,400), and (600,600).  The fit type here is set to a 2D Gaussian distribution.

The fifth line creates the absorption image from the raw data based on the values in the property `c.constants`.  The sixth line fits the absorption data in the ROI using the given fit type.  The last line plots the image, the ROI, and the x and y marginal distributions on a figure along with marginal distributions for the fit.  Information about the atomic sample, such as temperature, number of atoms, position, etc, are stored as properties in `c`.

## Description

### RawImageData Class

If we consider an absorption image to be comprised of a set of raw images, this class handles the loading of that set of raw images.  It has only three properties: the `directory` from which we should read image files, the `files` themselves represented as structures with file names, sizes, and other information, and finally the image data as a 3D array where each "slice" of the 3D array is one of the images in the set.

The user will primarily use this class through one of two methods: `load()` or the static method `loadImageSets()`.  Both of these methods will load raw images based on the variable argument list. The difference between the two methods is that the "dynamic" method `load()` applies only to the `RawImageData` object that called it, and can only load one image set at a time, while the static method can load multiple image sets and creates as many `RawImageData` objects as necessary. Below is the difference:
```
raw = RawImageData;                                                     %Create a blank instance of the RawImageData class
raw.load('filenames','last','index',1);                                 %Load the last image set in the default directory
raws = RawImageData.loadImageSets('filenames','last','index',1:5);      %Load the last 5 image sets in the default directory and store each in its one RawImageData object.
```
In the last line, `raws` is a 5x1 vector of `RawImageData` objects corresponding to the last five image sets, while `raw` is a single instance of `RawImageData` corresponding to the last image set.

Valid name/value pairs for either `load()` or `loadImageSets()` are
  - `filenames` or `files`: file names for the image sets or the keyword `'last'`.  If you use `'last'`, you should specify an index detailing how many image sets ago to load (last image, second-to-last image, etc).  If using file names, the file names should be a cell array of file names.  If you are using `loadImageSets()` with multiple image sets, then you should specify `filenames` as a cell array of cell arrays, where each cell in the parent array is a cell array of file names corresponding to the files in each image set.
  - `directory`: the directory from which to load data.  If not given, it defaults to `RawImageData.DEFAULT_DIRECTORY`.
  - `index` or `idx`: the index to use when paired with `'last'`.  This must be a scalar when using `load()` but can be a vector with `loadImageSets()`.
  - `length` or `len`: the number of image files per image set.  The default is 2 (an image with atoms and an image without atoms), but you can load more images if that is necessary.
  - `dims`: the dimensions of the images to load.  Useful when loading cropped images
  - `datatype`: the data type/format that the images are stored in.  Valid values are `'mono8'`, `'mono16'`, or `'raw8'`, corresponding to the same values on the Point Grey Grasshopper cameras.

Once the images are loaded, information about the files used can be accessed using the `files` property.  The image data, which consists of the integer representation of the voltage read from each pixel, is accessible using the `images` property.

### AtomImageConstants Class

Associated with the atomic sample and the absorption image are many parameters/constants, and these are captured by the `AtomImageConstants` class.  One can set these constants either at the construction of a particular instance of the class or afterwards.  For example, the three methods below can all be used for setting property values:
```
c = AtomImageConstants('Rb87','freqs',2*pi*[50,50,75]);             %Set atom type to Rb87 and then set property values at construction
c.set('detuning',0,'tof',30e-3);                                    %Set property values using the set() function
c.polarizationCorrection = 1.5;                                     %Set property value directly
```
When property values are set at construction, the first argument must be the atom type.  Currently, the only allowed value is `'Rb87'`.  When this is used, default values are set for many other properties such as the `mass`, `gamma` (linewidth in MHz), and saturation intensity, to name a few.  Users should read through the property declarations in the class definition for all property values.

### AtomCloudFit Class

We get parameters describing the atomic sample by fitting to the measured spatial distribution, and the fitting procedure is handled by the `AtomCloudFit` class.  Fitting starts by defining a "region-of-interest" or ROI which is a rectangular area defined by two 2-element vectors which give the start and finish of the rows and columns to include.  Additionally, one can set the step size to reduce the number of points used for fitting -- especially useful for 2D fits that can otherwise take a long time.  One can set these parameters using either the syntax
```
f = AtomCloudFit;       %Create a blank object
f.roiRow = [100,500];   %Set the start and end values for the rows in the ROI
f.roiCol = [250,750];   %Set the start and end values for the columns in the ROI
f.roiStep = [1,2];      %Set the step size for the row values (first element) and the column values (second element)
```
or
```
f = AtomCloudFit('roirow',[100,500],'roicol',[250,750],'roistep',[1,2]);
```
or
```
f = AtomCloudFit;
f.set('roirow',[100,500],'roicol',[250,750],'roistep',[1,2]);
```
All three methods are equivalent.  In the above examples, I have set the ROI to be a region bounded by the rows 100 to 500 and the columns 250 to 750, inclusive.  The image data will be sampled at every row and every second column.  When specifying the step size, if you use a single value that is assumed to apply to both rows and columns.

One can set the fit type/fit function using one of the below syntaxes
```
f.fittype = 'gauss2d';
%OR
f.set('fittype','gauss2d');
```
The valid fit types are
  - `'none'`: No fit is applied
  - `'sum'`: No fit is applied, but the image is summed within the ROI and converted into a number of atoms
  - `'gauss1d'`: Integrates the image data separately along the x and y axes and fits each marginal distribution to a 1D Gaussian with an offset and a linear variation in the backgroud.  The functional form is `z = A*exp(-(x-x0).^2/(2*s^2))+z0+linx*(x-x0)`.
  - `'tf1d'`: Fits the 1D marginal distributions to a 1D Thomas-Fermi profile of the form `z = A*(1-((x-x0)/s)^2).^2.*(abs(x-x0) < s) + z0 + linx*(x-x0)`.
  - `'2comp1d'`: Fits a two-component model comprised of a 1D Gaussian and a 1D Thomas-Fermi distribution
  - `'gauss2d'`: Fits a 2D Gaussian distribution to the image data with a uniform offset and linear variations along x and y in the background
  - `'tf2d'`: Fits a 2D Thomas-Fermi distribution to the image data with a uniform offset and linear variations in the background.  The form is `z = A*(1-((x-x0)/sx)^2-((y-y0)/sy)^2).^1.5.*abs((((x-x0)/sx)^2+((y-y0)/sy)^2) < 1) + z0 + linx*(x-x0) + liny*(y-y0)`.
  - `'2comp2d'`: Fits a 2D two-component distribution to the image data

With an ROI and a fit type defined, we now create our fit objects -- the x and y position vectors and the image data, as well as marginal distributions as needed -- and fit the data.  Assuming that the position vectors are labelled `x` and `y` and that the image data is `image` then we use
```
f.makeFitObjects(x,y,image);
f.fit();
```
Optionally, you can specify the fit type as the first argument in `f.fit(fittype)` where `fittype` is one of the allowed values above.  The function `f.makeFitObjects()` creates position vectors `f.x` and `f.y` corresponding to the position vectors restricted to the ROI, and marginal distributions `f.xdata` and `f.ydata`.  The 2D image data restricted to the ROI is stored in `f.image`.  Assuming that the fit converges, the resulting parameters are returned as a `CloudParameter` object in the property `f.params`.  The fit function that was used is returned in `f.fitfunc` and the marginal x and y distributions for the fitted data are in `f.xfit` and `f.yfit`.  The object `CloudParameter` has properties
  - `gaussAmp`: A one or two element vector of the Gaussian fit amplitudes
  - `becAmp`: A one or two element vector of the BEC fit amplitudes
  - `offset`: A one or two element vector of the uniform offsets
  - `lin`: A two element vector of the linear background variation coefficients
  - `pos`: A two element vector of the center position of the cloud in [x,y] format
  - `gaussWidth`: A two element vector of the Gaussian widths
  - `becWidth`: A two element vector of the BEC widths
  - `cloudAngle`: A one element vector of the rotation angle of the cloud for special 2D Gaussian fits
The property `f.params` is used for extracting the temperature and number of atoms from the fit.


### AbsorptionImage Class

This is the main, over-arching class that describes the entire absorption image and atomic sample.  It has properties `x` and `y` that describe the position on the image, a "raw" optical depth (OD) corresponding to the logarithm of the ratio of input to output light `image`, and a corrected version of that OD in `imageCorr`.  The corrections that are applied are to account for saturation of the OD due to, for instance, off-resonant light, saturation intensity corrections, and corrections due to the polarization of the light.  This corrected optical depth can be integrated over the ROI and the result divided by the absorption cross section and multiplied by a detuning factor to get the total number of atoms in the ROI.

The class has three "immutable" properties `raw`, `constants`, and `fitdata` corresponding to instances of `RawImageData`, `AtomImageConstants`, and `AtomCloudFit` classes.  These properties can only be set at the creation of the instance which means that the data type can never be changed, but since the classes are passed by reference we can still change the internal properties.

Finally, the `AbsorptionImage` class has the following properties for describing the atomic sample:
  - `N`: the number of atoms extracted from the fit
  - `Nsum`: the number of atoms extracted by integrating over the ROI
  - `pos`: the center position of the sample extracted from the fit
  - `gaussWidth`: the width of the fitted Gaussians
  - `T`: the temperature of the clouds extracted from the Gaussian widths, the time of flight, and the trap frequencies
  - `peakOD`: the peak OD in the ROI
  - `PSD`: the phase-space density of the sample
  - `becFrac`: the fraction of atoms in the sample that are in the condesate, as extracted from the fits.  Always 0 for the fit type `sum`.
  - `becWidth`: the width of the condensate

We start by creating an instance of the `AbsorptionImage` class.  If we already have instances of `RawImageData`, `AtomImageConstants`, and `AtomCloudFit` that we want to use, we can call
```
c = AbsorptionImage(raw,constants,fitdata);
```
where `raw`, `constants`, and `fitdata` are instances of `RawImageData`, `AtomImageConstants`, and `AtomCloudFit, respectively.  We can also create a blank object
```
c = AbsorptionImage;
```
instead.  If we then copy properties from these other class instances using
```
c.raw.copy(raw);
c.constants.copy(constants);
c.fitdata.copy(fitdata);
```
Note that we can't use `c.raw = raw` since the `raw` property in `AbsorptionImage` has the `SetAccess = immutable` flag set.  Besides, since `RawImageData` is passed by reference, this kind of assignment does not create a copy of `raw`, so if it gets changed someone else it will affect the current instance of `AbsorptionImage`.

Assuming that we have properly set the immutable properties `raw`, `constants`, and `fitdata`, we can create the optical depth maps `c.image` and `c.imageCorr` using the `makeImage()` method
```
c.makeImage();
%OR
c.makeImage(imgIndexes);
```
The first way of calling `makeImage()` assumes that each absorption image is a set of two raw images where the first image is with the atoms and the second image is without the atoms.  The second way of calling the method assumes that the first element in `imgIndexes` points to the raw image that is with the atoms and the second element points to the raw image that is without the atoms.  This can be useful if one is taking multiple images with atoms but only one reference image if, for instance, one is doing multi-state imaging.

With the OD maps in hand, we can now fit our data using
```
c.fit();
%OR
c.fit('method',calcmethod,'fittype',fittype);
```
where the first call type assumes that one wants to use the fit type set in `fitdata`, whereas the second method will overwrite that method.  The parameter `'method'` refers to how the number of atoms is calculated when using 1D distributions.  If set to `'y'`, it calculates the number of atoms from the y marginal distribution, and similarly for `'x'`.  If set to `'xy'` or `'yx'` it calculates from a geometric mean of the two results.

The method `c.fit()` also sums over the region of interest to get the number of atoms `Nsum` (as compared to the fitted number of atoms `N`).  If a 2D fit is used, the uniform offset from the fit is subtracted from the OD data before summation, and this can greatly improve the reliability of calculating interferometer phases in a model-independent way.

Finally, one can plot the absorption image using various plot functions.  Users should look through the plotting functions for more specifics, but one can plot the x or y marginal distributions, the absorption data, or all of the data.

