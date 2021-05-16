# Imaging Analysis

This repository contains code for analyzing absorption images of cold atomic samples.  Using the standard method of illuminating an atomic sample with near-resonant light and recording the intensity pattern and then comparing that to a reference image with no atoms present, it constructs a map of optical depth (OD) and uses that for determinining properties of the sample.

## Use

The analysis code is broken into four main classes which can be invoked using a wrapper file.  These classes are
  - `RawImageData`: describes the raw image data, and has methods for loading that data from a directory
  - `AtomImageConstants`: describes various constants for defining an absorption image, such as detuning, pixel size, etc.
  - `AtomCloudFit`: mainly provides methods for fitting image data to various distributions.
  - `AbsorptionImage`: the over-arching class that describes the image.  Has three "constant" properties `raw`, `constants`, and `fitdata` which are each instances of `RawImagData`, `AtomImageConstants`, and `AtomCloudFit`, respectively.  Other properties describe the image and the sample, such as peak OD, number of atoms, temperature, etc.

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


