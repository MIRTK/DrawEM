Draw-EM Segmentation Software
==========================================

Draw-EM (Developing brain Region Annotation With Expectation-Maximization) is a package of MIRTK developed by Antonios Makropoulos and the [BioMedIA](https://biomedia.doc.ic.ac.uk/) research group. 
It provides a collection of command-line tools as well as pipelines for the segmentation of developing brain MR images.


Installation
------------

Draw-EM is part of MIRTK. 
Installation is enabled by setting the CMake flag "MODULE_DrawEM" of MIRTK to "ON"

The atlases required by Draw-EM need to be downloaded from [here](https://www.doc.ic.ac.uk/~am411/atlases-DrawEM.html) and extracted inside the Draw-EM directory.

See the [installation instructions](https://mirtk.github.io/install.html) 
for a step-by-step guide on how to install the MIRTK.


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

In the event you found Draw-EM useful, please consider giving appropriate credit to the software.

Publication:

A. Makropoulos et al. Automatic whole brain MRI segmentation of the developing neonatal brain, IEEE TMI, 2014
