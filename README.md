# Imaging Analysis

This repository contains code for analyzing absorption images of cold atomic samples.  Using the standard method of illuminating an atomic sample with near-resonant light and recording the intensity pattern and then comparing that to a reference image with no atoms present, it constructs a map of optical depth (OD) and uses that for determinining properties of the sample.

## Use

The analysis code is broken into six main classes which can be invoked using a wrapper file.  These classes are
  - `RawImageData`: describes the raw image data, and has methods for loading that data from a directory
  - `AtomImageConstants`: describes various constants for defining an absorption image, such as detuning, pixel size, etc.
  - `CloudParameters` : defines the properties used in fits and is used in defining bounds and starting points for the fits
  - `AtomCloudFit`: mainly provides methods for fitting image data to various distributions.
  - `AtomCloud` : defines the properties of the particular cloud of atoms under analysis. Containts instances of `AtomImageConstants` and `AtomCloudFit`.
  - `AbsorptionImage`: the over-arching class that describes the image.  Has two immutable properties `raw` and `constants` which are each instances of `RawImageData` and `AtomImageConstants`, respectively.  Also contains an *array* of `AtomCloud` objects, each describing a cloud of atoms in the image.

Below is a bare-bones example of how one can read an image set corresponding to a single absorption image containing a single atom cloud and analyze that image:
```
img = AbsorptionImage;
img.raw.load('filenames','last','directory','.','index',1);
img.constants.set('tof',30e-3,'exposuretime',30e-6,'detuning',0);
img.setClouds(1);
img.clouds.fitdata.set('roirow',[400,600],'roicol',[400,600],'fittype','gauss2d');
img.makeImage();
img.fit();
img.plotAllData([0,3]);
```
The first line creates a blank `AbsorptionImage` object `img`.  The second line loads the last set of images in the nominated directory (here the current directory `.`); last in this context is ordered according to time.  The index specifies how many sets to go back from the last image: a value of 1 is the last image set taken, a value of 2 is the second-to-last image, and so on.

The third line sets some constants associated with the image; in this case, we set the time-of-flight `tof` to 30 ms, the exposure time to 30 us, and the detuning to 0 MHz.  Other properties are retained as default values.  The fourth line defines the number of atom clouds in the image; in this case it is 1.  The fifth line sets the "region-of-interest" (ROI) to be the rectangular region with corners at (400,400), (400,600), (600,400), and (600,600).  The fit type here is set to a 2D Gaussian distribution.

The sixth line creates the absorption image from the raw data based on the values in the property `img.constants`.  The seventh line fits the absorption data in the ROI using the given fit type.  The last line plots the image from an OD of 0 to an OD of 3, the ROI, and the x and y marginal distributions on a figure along with marginal distributions for the fit.  Information about the atomic sample, such as temperature, number of atoms, position, etc, are stored as properties in `c`.

## Description

### RawImageData Class

If we consider an absorption image to be comprised of a set of raw images, this class handles the loading of that set of raw images.  It has only three public properties: the `directory` from which we should read image files, the `files` themselves represented as structures with file names, sizes, and other information, and finally the image data as a 3D array where each "slice" of the 3D array is one of the images in the set.  An additional property which has protected 'SetAccess' is a status property, which can be used to flag errors in loading images.

The user will primarily use this class through one of two methods: `load()` or the static method `loadImageSets()`.  Both of these methods will load raw images based on the variable argument list. The difference between the two methods is that the "dynamic" method `load()` applies only to the `RawImageData` object that called it, and can only load one image set at a time, while the static method can load multiple image sets and creates as many `RawImageData` objects as necessary. Below is the difference:
```
raw = RawImageData;                                                     %Create a blank instance of the RawImageData class
raw.load('filenames','last','index',1);                                 %Load the last image set in the default directory
raws = RawImageData.loadImageSets('filenames','last','index',1:5);      %Load the last 5 image sets in the default directory and store each in its one RawImageData object.
```
In the last line, `raws` is a 5x1 vector of `RawImageData` objects corresponding to the last five image sets, while `raw` is a single instance of `RawImageData` corresponding to the last image set.

Valid name/value pairs for either `load()` or `loadImageSets()` are
  - `filenames` or `files`: file names for the image sets or the keyword `'last'`.  If you use `'last'`, you should specify an index detailing how many image sets ago to load (last image, second-to-last image, etc).  If using file names, the file names should be a cell array of file names or a cell array of file structures output from `dir()`.  If you are using `loadImageSets()` with multiple image sets, then you should specify `filenames` as a cell array of cell arrays, where each cell in the parent array is a cell array of file names corresponding to the files in each image set.
  - `directory`: the directory from which to load data.  If not given, it defaults to `RawImageData.DEFAULT_DIRECTORY`.
  - `index` or `idx`: the index to use when paired with `'last'`.  This must be a scalar when using `load()` but can be a vector with `loadImageSets()`.
  - `length` or `len`: the number of image files per image set.  The default is 2 (an image with atoms and an image without atoms), but you can load more images if that is necessary.  Currently, the `AbsorptionImage.makeImage()` method can handle three raw images, where the first image is the image with atoms, the second image is the image without atoms, and the last image is a dark image to correct for camera background.
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
Optionally, you can specify the fit type as the first argument in `f.fit(fittype)` where `fittype` is one of the allowed values above.  The function `f.makeFitObjects()` creates position vectors `f.x` and `f.y` corresponding to the position vectors restricted to the ROI, and marginal distributions `f.xdata` and `f.ydata`.  The 2D image data restricted to the ROI is stored in `f.image`.  Assuming that the fit converges, the resulting parameters are returned as a `CloudParameters` object in the property `f.params`.  The object `CloudParameters` has properties
  - `gaussAmp`: A one or two element vector of the Gaussian fit amplitudes
  - `becAmp`: A one or two element vector of the BEC fit amplitudes
  - `offset`: A one or two element vector of the uniform offsets
  - `lin`: A two element vector of the linear background variation coefficients
  - `pos`: A two element vector of the center position of the cloud in [x,y] format
  - `gaussWidth`: A two element vector of the Gaussian widths
  - `becWidth`: A two element vector of the BEC widths
  - `cloudAngle`: A one element vector of the rotation angle of the cloud for special 2D Gaussian fits
The property `f.params` is used for extracting the temperature and number of atoms from the fit.  

The marginal x and y distributions for the fitted data are in `f.xfit` and `f.yfit`.  Fit residuals are returned in the `f.residuals` property which is a 2D array for 2D fits and a structure with fields `x` and `y` for 1D fits.  Fit residuals are calculated as `f.residuals = f.image - (fit function)`, so in principle the fit function can be reconstructed from the residuals and input data.  Similarly, the background of the fit is returned in the property `f.bg` in the same manner for 2D and 1D fits, and this background field can be used to calculate the number of atoms in a cloud by directly summing over the OD map.

When preparing for the fit, the `AtomCloudFit` class attempts to guess upper and lower bounds for the fit as well as a starting point.  In most situations these guesses work well, but in some situations the user may want to set their own constraints and starting values.  The `AtomCloudFit` class has the properties `lb`, `ub`, and `guess`, each of which are an instance of `CloudParameters`, which are the lower bounds, upper bounds, and starting value, respectively.  If the user sets the value of a property in one of these objects, then that value is used in place of the automatically generated bound/guess.  For instance, suppose the user wants to more tightly constrain the width of a BEC fit.  They would set
```
f.lb.set('becwidth',650e-6*[1,1]);    %Set lower bound for x and y widths to be 650 um
f.ub.set('becwidth',850e-6*[1,1]);    %Set the upper bound for the x and y widths to be 850 um
f.guess.set('becwidth',750e-6*[1,1]); %Set the starting values for the x and y widths to be 750 um
```
These bounds and guesses are then used in place of the automatically generated ones.  Any combination of values can be set.  The `AtomCloudFit` class keeps track of which constraints to override by assuming that if the corresponding property in the user-defined constraint is set to the empty vecot `[]` (the default value) then it should use the automatically generated constraint.

### AtomCloud Class

This class describes a particular atom cloud in an absorption image, as for atom interferometry there can be multiple atom clouds in each image.  It contains properties `fitdata` and `constants` which are instances of `AtomCloudFit` and `AtomImageConstants`, respectively.  Note that `constants` here is the *same* object as in the parent `AbsorptionImage` class, as these constants are obviously shared with each atom cloud in an image.  Additionally, the `AtomCloud` class has the following properties for describing an atom cloud:
  - `N`: the number of atoms extracted from the fit
  - `Nsum`: the number of atoms extracted by integrating over the ROI
  - `pos`: the center positions of the sample extracted from the fit
  - `gaussWidth`: the widths of the fitted Gaussians
  - `T`: the temperature of the clouds extracted from the Gaussian widths, the time of flight, and the trap frequencies
  - `peakOD`: the peak OD in the ROI
  - `PSD`: the phase-space density of the sample
  - `becFrac`: the fraction of atoms in the sample that are in the condesate, as extracted from the fits.  Always 0 for the fit type `sum`.
  - `becWidth`: the widths of the condensate

### AbsorptionImage Class

This is the main, over-arching class that describes the entire absorption image and the atomic samples.  It has properties `x` and `y` that describe the position on the image, a "raw" optical depth (OD) corresponding to the logarithm of the ratio of input to output light `ODraw`, and a corrected version of that OD in `ODcorr`.  The corrections that are applied are to account for saturation of the OD due to, for instance, off-resonant light, saturation intensity corrections, and corrections due to the polarization of the light.  This corrected optical depth can be integrated over the ROI and the result divided by the absorption cross section and multiplied by a detuning factor to get the total number of atoms in the ROI.

The class has two "immutable" properties `raw` and `constants` corresponding to instances of `RawImageData` and `AtomImageConstants` classes.  These properties can only be set at the creation of the instance which means that the data type can never be changed, but since the classes are passed by reference we can still change the internal properties.  Finally, the class contains an *array* of instances of `AtomCloud` corresponding to the atom clouds within the image.

We start by creating an instance of the `AbsorptionImage` class.  If we already have instances of `RawImageData` and `AtomImageConstants` that we want to use, we can call
```
img = AbsorptionImage(raw,constants);
```
where `raw` and `constants` are instances of `RawImageData` and `AtomImageConstants`, respectively.  We can also create a blank object
```
img = AbsorptionImage;
```
instead.  If we then copy properties from these other class instances using
```
img.raw.copy(raw);
img.constants.copy(constants);
```
Note that we can't use `c.raw = raw` since the `raw` property in `AbsorptionImage` has the `SetAccess = immutable` flag set.  Besides, since `RawImageData` is passed by reference, this kind of assignment does not create a copy of `raw`, so if it gets changed someone else it will affect the current instance of `AbsorptionImage`.

Assuming that we have properly set the immutable properties `raw` and `constants`, we can create the optical depth maps `img.ODraw` and `img.ODcorr` using the `makeImage()` method
```
img.makeImage();
%OR
img.makeImage(imgIndexes);
```
The first way of calling `makeImage()` assumes that each absorption image is a set of two raw images where the first image is with the atoms and the second image is without the atoms.  The second way of calling the method assumes that the first element in `imgIndexes` points to the raw image that is with the atoms and the second element points to the raw image that is without the atoms.  This can be useful if one is taking multiple images with atoms but only one reference image if, for instance, one is doing multi-state imaging.

Next, we need to indicate how many atom clouds are in our image, and then define the fitting parameters for these clouds.  Supposing that we have only one cloud, we would use
```
img.setClouds(1);
img.clouds.fitdata.set('roiRow',[400,700],'roiCol',[500,800],'fittype','gauss2d');
```
where we have set the fitting region to be the rows `[400,700]` and the columns `[500,800]`, and the fit type is a 2D Gaussian.

With the OD maps in hand, we can now fit our data using
```
c.fit();
```
which automatically creates the fit objects for every cloud and then performs the nominated fits.  The fitting procedure also sums the image over the region of interest after correcting for the fitted background of the absorption image.  This result, stored in `AtomCloud.Nsum`, can be compared to the fitted value `AtomCloud.N`, but be warned that this comparison is only reliable for 2D fits.

Finally, one can plot the absorption image using various plot functions.  Users should look through the plotting functions for more specifics, but one can plot the x or y marginal distributions, the absorption data, or all of the data.  For most purposes, the method `AbsorptionImage.plotAllData(ODrange)` will suffice, where `ODrange` is a vector of colour axis limits to apply.  If more than one cloud has been defined, then the marginal distributions for each cloud will be plotted on the same axes.



