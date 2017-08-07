Draw-EM Fetal Segmentation Software
==========================================

This is a pilot version for the fetal brain segmentation.
It is based on the dhcp branch.

The new files with respect to the dhcp branch are:
- pipelines/fetal-pipeline.sh
- scripts/fetal-tissue-priors.sh


Installation
------------

See installation of the dhcp branch.
- Additional atlases are required for the fetal segmentation which need to be downloaded from [here](https://www.doc.ic.ac.uk/~am411/atlases-DrawEM.html) and extracted inside the Draw-EM directory.


Run
---

The segmentation pipeline can be run with the following script:

pipelines/fetal-pipeline.sh 

The script requires the T2 image and the age at scan of the subject to be segmented (as first and second argument respectively).
Run the script without arguments for a detailed list of options.


