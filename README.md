Draw-EM Segmentation Software
==========================================

Draw-EM (Developing brain Region Annotation With Expectation-Maximization) is a package of [MIRTK](https://github.com/BioMedIA/MIRTK) developed by Antonios Makropoulos and the [BioMedIA](https://biomedia.doc.ic.ac.uk/) research group. 
It provides a collection of command-line tools as well as pipelines for the segmentation of developing brain MR images.


Installation
------------

Draw-EM is part of MIRTK. 
In order to compile it as part of MIRTK you need to do the following steps:
- download (clone) Draw-EM inside the Packages folder of your MIRTK directory
- enable compile of the package by setting the CMake flag "MODULE_DrawEM" of MIRTK to "ON" (using cmake or ccmake)
- build MIRTK

See the [installation instructions](https://mirtk.github.io/install.html) 
for a step-by-step guide on how to install MIRTK.

The segmentation pipeline further requires the following:
- [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) installed
- The atlases required by Draw-EM need to be downloaded from [here](https://www.doc.ic.ac.uk/~am411/atlases-DrawEM.html) and extracted inside the Draw-EM directory.
- The N4 bias field correction from ITK is included in the ThirdParty/ITK folder. If the executable does not work, you will need to compile ITK and replace the ThirdParty/ITK/N4 binary with the N4BiasFieldCorrectionImageFilter binary


Run
---

The segmentation pipeline can be run with the following script:

pipelines/neonatal-pipeline-v1.1.sh 

The script requires the T2 image and the age at scan of the subject to be segmented (as first and second argument respectively).
Run the script without arguments for a detailed list of options.


License
-------

Draw-EM is distributed under the terms of the Apache License Version 2.
See the accompanying [license file](LICENSE.txt) for details. The license enables usage of
Draw-EM in both commercial and non-commercial applications, without restrictions on the
licensing applied to the combined work.

Draw-EM uses third-party software, namely the "ITK: The Insight Toolkit for Segmentation and Registration".
ITK is distributed under the Apache License Version 2.
Specifically, the N4 bias field correction by Tustison et al. is included (http://www.insight-journal.org/browse/publication/640).
The covered file (N4) and license (LICENSE) can be found in ThirdParty/ITK.


Citation and acknowledgements
-----------------------------

In case you found Draw-EM useful please give appropriate credit to the software.

Publication:

A. Makropoulos et al. Automatic whole brain MRI segmentation of the developing neonatal brain, IEEE TMI, 2014
