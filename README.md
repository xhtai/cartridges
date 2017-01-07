<!-- README.md is generated from README.Rmd. Please edit that file -->
cartridges
==========

This package contains utilities to read and plot cartridge case images. It pre-processes these images and implements an algorithm to compare the images. It produces a similarity score for each pairwise comparison and computes the probability of obtaining a higher score by chance. This is work in progress.

Images must be in the standard format as released by the National Institute of Standards and Technology (NIST). Such data are available from the [NIST Ballistics Toolmark Research Database](https://tsapps.nist.gov/NRBTD), and in particular we are developing methods for images of breech face impressions using 2D ring light. This methodology has been tested on images from the NBIDE study in the database. An example image from the database is below.

![](README-unnamed-chunk-2-1.png)

Background
----------

When a gun is fired, it leaves marks on the bottom surface of the cartridge, and these marks are thought to be unique to the gun, in other words each gun is thought to produce unique marks. Here the marks that we are interested in are the breechface marks, which are marked out in the image below. These marks are caused by the bottom surface of the cartridge pressing against the breech block of the gun during the firing process.

![](README-BF.png)

Since the marks are thought to be unique, when cartridge cases are collected from crime scenes, they can be compared to cartridge cases that have been collected previously, to see if they come from the same gun as something that has been seen before. Because of the large number of cartridge cases being collected, rather than comparing these physical cartridge cases, we are interested in using automated algorithms to compare these images of the bottom surface of the cartridge cases.

Given a new image and a database of images, this package implements tools to assess how similar the new image is to the images in the database, producing a similarity metric. To quantify the uncertainty in this comparison procedure, and also attach more meaning to the similarity scores, we also propose a method of computing the probability of obtaining a larger score by chance.

Description of method
---------------------

These are the steps for one pairwise comparison. There are 4 pre-processing steps before we compute the two measures that we are interested in.

1.  Automatically select breechface marks
2.  Level image
3.  Remove circular symmetry
4.  Outlier removal and filtering
5.  Maximize correlation by translations and rotations
6.  Compute probability of obtaining a higher score by chance

Steps 2, 4 and 5 are implementations of work that has been done by Vorburger and co-authors (NISTIR, 2007), and Roth, Carriveau, Liu, Jain (IEEE, 2015). We illustrate each of these steps using an example image.

### Step 1

The first step is to automatically select the breechface marks, and this is broken into two steps, finding the primer region and then removing the firing pin impression. A rough schema for finding the primer region is as follows.

![](README-primer.png)

To remove the firing pin impression, we use use a similar set of steps, but we start with an edge detector to first identify the firing pin region. This would work even if the firing pin impression isn't circular.

![](README-FP1.png)

Now we perform a second pass where we apply an edge detector again, to try to remove any remaining marks that we might have missed the first time. This is necessary for some images where the firing pin impression might not be as highly contrasted with the surrounding breechface impression. In this particular example a second pass is not actually required, since the entire firing pin impression has been removed. Instead, several small marks on the breechface were picked out and removed. This is an unintended consequence but one that we live with for now. It turns out that we still achieve good results for most images in later comparisons, but this is an area for further improvement (TODO: in particular, in the NBIDE data set, NBIDE009, NBIDE110, NBIDE128 have poorer performance due to this step).

![](README-FP2.png)

### Step 2

Step 2 is to level the image. The reason for this step is that the base of the cartridge case may not be level, and may instead be tilted slightly on a plane. Images of such a surface may have differences in brightness that are planar in nature. Here we fit a plane that captures these differences, and taking the residuals ensures that the resulting image is free from such differences in brightness.

In this example the original image is slightly darker in the bottom left corner and brighter on the top right, and we can see this in the left panel in the figure below. The residuals on the right panel are free from any such effects. We take the residuals for further processing.

![](README-unnamed-chunk-3-1.png)

### Step 3

Step 3 is to remove any circular symmetry. The reason for this step is that apart from the base of the cartridge not being level, there could also be differences in depth that are circular in nature, for example the surface may slope inwards towards the center. This would cause differences in brightness that are circular in nature, for example the center of the image may be darker than the edges. Like in the previous step, we fit a model that captures this circular symmetry, and then take the residuals, so that the residuals would be free from any circular symmetry.

The model that we are using is a linear combination of circularly symmetric basis functions. This model assumes that pixels located the same distance from the center of the image take the same value, and the first few matrices in the basis are given in the figure below, where each figure in the panel represents one matrix. Each matrix takes the value 1 for pixels that are the same distance from the center, and zero otherwise. Basis are enumerated from center outwards.

![](README-basis.png)

We represent these matrices as circularly symmetric basis functions, which take an ij coordinate as an input and return the value 0 or 1. An image can then be decomposed as follows:

![](README-equation.png)

![](README-bulletPoints.png)

To fit this model, the coefficient of each basis function is the mean of pixel values of pixels in that basis function. The fitted coefficients for our example image are in the figure below.

![](README-basisCoefs.png)

Because of the large number of basis functions, with each only containing only a few pixels, the variance of the coefficients is very high. We fit a local smoother through the coefficients to get a smoothed circularly symmetric model. The fitted model and the residuals are below. These residuals are free from both planar bias from the previous step, and circular symmetry.

![](README-removeCircular.png)

### Step 4

The last pre-processing step is outlier removal and filtering. Outliers are removed and inpainted so that they won't affect the similarity scores being computed, and filtering highlights some features of the image. This methodology was described by Vorburger and co-authors (NISTIR, 2007) and implemented in MATLAB by Roth, Carriveau, Liu, Jain (IEEE, 2015).

After all pre-processing, we get the following image.

![](README-processed.png)

### Step 5

After pre-processing, step 5 involves computing a similarity metric. Again, this step was described by Vorburger and co-authors (NISTIR, 2007) and implemented in MATLAB by Roth, Carriveau, Liu, Jain (IEEE, 2015).

The similarity metric is the correlation between the two images, and this is known in the literature as the maximum cross-correlation function, or *C**C**F*<sub>*m**a**x*</sub>. We first compute the cross-correlation function for each rotation angle:

![](README-CCFequation.png)

where *I*<sub>1</sub> and *I*<sub>2</sub> are the two images, *i* indexes the rows and *j* indexes the columns, and *d**x* and *d**y* represent translations. The *C**C**F* returns a matrix of correlation values, where each entry corresponds to a particular translation, and we store the maximum correlation. Repeating for many rotation angles, we obtain *C**C**F*<sub>*m**a**x*</sub>. Since this is a correlation, it takes values between -1 and 1, and can be interpreted as the maximum correlation between two images after lining them up correctly.

We compare our example image against another image from the NBIDE study, which was obtained using the same gun. We obtain a similarity score of .37, with a rotation angle of −15<sup>∘</sup>, meaning that the second image is rotated 15<sup>∘</sup> counter-clockwise. Plotting the two images with the second correctly rotated, we notice that the breechface marks are now lined up nicely.

![](README-comparison.png)

### Step 6

The last step is to convert each similarity score into a statement of probability. In particular we compute the probability of obtaining a higher score by chance. These probabilities attach meaning to the scores, and also serve as a measure of uncertainty for this procedure. Our proposed method is as follows.

We assume that all *C**C**F*<sub>*m**a**x*</sub> values for non-matches are drawn from the same distribution, and given such a distribution, we compare each newly computed score against this distribution, and compute the right tail proportion. This value is the probability of observing a higher *C**C**F*<sub>*m**a**x*</sub> by chance.

In reality, we do not have access to such a distribution. What we might have is a known database, where we are able to compute all pairwise non-matching scores. These form a sample from the unknown distribution. For example, using the NBIDE data set, we have a total of 108 images from 12 different guns. Doing all pairwise comparisons within the database, we have a total of 10692 non-match scores, which would form a sample from the unknown population of non-match scores. We can then compare .37 against this distribution.

Installation
------------

If you do not have R installed, visit <https://www.r-project.org/>. You will also need to install the `devtools` package from CRAN using `install.packages("devtools")`.

To install this package:

``` r
library(devtools)
devtools::install_github("xhtai/cartridges")
```

`EBImage` is required and is hosted on Bioconductor. To install, use

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

After installing the `cartridges` package, load it using

``` r
library(cartridges)
```

For full package functionality, the following packages are suggested: `fields`, `imager`, `purrr`, `dplyr` and `plyr`. You can install these using `install.packages()`. For `imager`, if you are using a Linux machine you might need `libX11` and `libfftw3` (`sudo apt-get install libfftw3-dev libX11-dev`). For more details see <https://github.com/dahtah/imager>.

#### Difficulties using devtools

If you are using a Linux machine, you might have some difficulties as `devtools` has a number of dependencies (e.g. `libssl-dev`, `libcurl4-gnutls-dev`). If you are unable to install these, you can download the [tarball](https://github.com/xhtai/cartridges/tarball/master) directly and install using `install.packages(file.choose(), repos=NULL)` or `R CMD INSTALL`. Before doing this, you will also need to make sure that you have the following R packages installed: `EBImage`, `methods`, `Matrix` and `magrittr`. Then install `cartridges` using

``` r
install.packages(file.choose(), repos = NULL)
```

A window will pop-up and all you have to do is select the location of the .tar.gz file.

If you are using a Windows machine, you can do the same using a [.zip file](https://github.com/xhtai/cartridges/zipball/master).

Functions available
-------------------

-   `readCartridgeImage`: to read in a raw image and obtain a matrix of pixel values
-   `plotImage`: produces a plot from a matrix of values. Values may be pixel values, residuals, etc.
-   `allPreprocess`: performs steps 1-4 above
-   `calculateCCFmax`: step 5 above
-   `computeProb`: step 6 above

The remaining functions perform intermediate steps.

-   Step 1: `findPrimer`, `findFP`
-   Step 2: `centerBFprimer`, `levelBF`
-   Step 3: `removeCircular`. Some more general functions for the circularly symmetric model are `getBasisFunctions`, `statisticsByBasisFunction`, `fitBasis`, and `getFitted`.
-   Step 4: `cropBorders`, `outlierRejection`, `inpaint_nans`, `gaussianFilter`
-   Step 5: `comparison`

Help files with examples for each of these functions can be accessed using `help(functionName)` or `?functionName`, e.g. `help(plotImage)`.

Example
-------

In illustrating the steps above we processed an example image, taken from the NBIDE study in the NIST database. This is a 2D breechface ring light image, and the original filename in the NIST download is "NBIDE R BF 118.png". The cartridge being imaged is a test fire from a Ruger gun (gun 3 in the study), using PMC ammunition. The raw data can be accessed from this package using `system.file("extdata", "NBIDE R BF 118.png", package="cartridges")`.

We can read in and plot the image as follows:

``` r
exampleImage <- readCartridgeImage(system.file("extdata", "NBIDE R BF 118.png", 
    package = "cartridges"))
plotImage(exampleImage, type = "original")
```

We can perform all pre-processing using the following code. Note that this could take around 10 minutes or more, depending on your set-up. The processed image is available within the package, and can be accessed using `preprocessedExample`.

``` r
processedExample <- allPreprocess(system.file("extdata", "NBIDE R BF 118.png", 
    package = "cartridges"))
```

Now, a second processed image is available as `preprocessedExample2`, and this was produced using the same gun, so the similarity score should be high. To compare these two images, we use

``` r
calculateCCFmax(processedExample, processedExample2)
```

This should take about a minute and produce a score of .37. Finally to determine the probability of obtaining a higher score by chance, we will need to obtain a known database of non-match scores. This was described earlier. Given such a set of scores, we can run

``` r
computeProb(0.37, knownScores)
```

to obtain the required probability.

If we are interested in the results from each of the pre-processing steps, we can run each step manually using the following code. Because some of these steps take a while to run, the following `.rda` files are available in the package to illustrate the results from intermediate steps: `primerExample`, `FPexample`, `removedExample`, and `inpaintedExample`.

``` r
# step 1
primerExample <- findPrimer(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"))
FPexample <- findFP(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"), 
    primer = primerExample)

# step 2
centeredExample <- centerBFprimer(FPexample, primerExample)
leveledExample <- levelBF(centeredExample$centeredBF)

# step 3
removedExample <- removeCircular(leveledExample)

# step 4
croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
outlierNAexample <- outlierRejection(croppedExample)
inpaintedExample <- inpaint_nans(outlierNAexample)

nonBF <- is.na(croppedExample)
processedExample <- gaussianFilter(inpaintedExample, nonBF)
```

Further work
------------

There are possible improvements to both the methodology and the code. Some of these are:

-   The second pass for removing the firing pin impression results in some valid areas being removed. In particular in the NBIDE data set, performance on images 9, 110 and 128 can be improved.
-   Code speed-ups are possible, especially for steps 4 and 5, which were optimized for MATLAB
-   Test on more data.

Credits
-------

This is work with William F. Eddy, with advice from Xiaoyu Alan Zheng, who also maintains the NIST database. Joseph Roth provided MATLAB code which we translated for steps 4 and 5. Max Mitchell suggested code speed-ups.

License
-------

The `cartridges` package is licensed under GPLv3 (<http://www.gnu.org/licenses/gpl.html>).

References
----------
