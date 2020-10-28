# Draw-EM Segmentation Software

![segmentation image](segmentation.png)

Draw-EM (Developing brain Region Annotation With Expectation-Maximization) is a package of [MIRTK](https://github.com/BioMedIA/MIRTK) developed by Antonios Makropoulos and the [BioMedIA](https://biomedia.doc.ic.ac.uk/) research group. 
It provides a collection of command-line tools and pipelines for the segmentation of developing brain MR images.

Draw-EM is used as part of the [dHCP structural pipeline](https://github.com/BioMedIA/dhcp-structural-pipeline) for the structural analysis (segmentation and surface extraction) of the neonatal brain.


## Dependencies
### FSL

The segmentation pipeline uses
[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL). 
See the [installation instructions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) for FSL.


## Installation

Draw-EM is part of MIRTK.

You can install Draw-EM using the provided setup file ([setup.sh](setup.sh)):
```
$ ./setup.sh -DWORKBENCH_install=0 -DSPHERICALMESH_install=0 [-j <num_cores>] 
```

where `num_cores` the number of CPU cores used to compile the software.

## Running the pipeline

The segmentation pipeline can be run as follows:

mirtk neonatal-segmentation <subject_T2> <age_at_scan>

```
Arguments:
  <subject_T2>                  Nifti Image: The T2 image of the subject to be segmented.
  <age_at_scan>                 Integer: Subject age in weeks. This is used to select the appropriate template for the initial registration. 
			        					If the age is <28w or >44w, it will be set to 28w or 44w respectively.
Options:
  -a / -atlas  <atlasname>      Atlas used for the segmentation, options: ALBERT, MCRIB (default: ALBERT)
  -ta / -tissue-atlas  <atlasname>  Atlas used to compute the GM tissue probability, options: neonatal
  -d / -data-dir  <directory>   The directory used to run the script and output the files.
  -c / -cleanup  <0/1>          Whether cleanup of temporary files is required (default: 1)
  -p / -save-posteriors  <0/1>  Whether the structures' posteriors are required (default: 0)
  -t / -threads  <number>       Number of threads (CPU cores) allowed for the registration to run in parallel (default: 1)
  -v / -verbose  <0/1>          Whether the script progress is reported (default: 1)
  -h / -help / --help           Print usage.
```

The produced segmentations will be stored in the `segmentations` folder:
- *\<subject\>_all_labels.nii.gz*: output segmentation with all the labels used by Draw-EM. This file includes some labels that do not exist in the atlases that are helpful for segmentation (e.g. the ALBERTs do not include division of labels into grey matter and white matter, this is done automatically).
- *\<subject\>_labels.nii.gz*: output label segmentation
- *\<subject\>_tissue_labels.nii.gz*: output tissue segmentation
- *\<subject\>_brain_mask.nii.gz*: output brain mask
- *\<subject\>_L_white.nii.gz*, *\<subject\>_R_white.nii.gz*: Left and right mask useful for white matter surface reconstruction (includes white matter and deep grey matter).
- *\<subject\>_L_pial.nii.gz*, *\<subject\>_R_pial.nii.gz*: Left and right mask useful for pial surface reconstruction (includes grey matter, white matter and deep grey matter).

Labels and lookup tables for all the *labels.nii.gz files are provided in the [label_names](label_names) folder according to the atlas used.


## Atlases
The following atlases can be used for segmentation.

#### ALBERTs
The ALBERTs atlases are described in Gousias et al. "Magnetic resonance imaging of the newborn brain: Manual segmentation of labelled atlases in term-born and preterm infants",  NeuroImage, 2012.

They define 50 labels in total, with 16 cortical labels per hemisphere and 16 sub-cortical structures.

The ALBERTs atlases can be used by providing the `-a ALBERT` argument when running the pipeline.

Labels of produced segmentations are provided in the [label_names/ALBERT](label_names/ALBERT) folder.

#### M-CRIB 2.0
The M-CRIB 2.0 atlases are described in Alexander et al. "Desikan-Killiany-Tourville Atlas Compatible Version of M-CRIB Neonatal Parcellated Whole Brain Atlas: The M-CRIB 2.0", Front. Neurosci., 2019.

They define 94 labels in total, with 31 cortical labels per hemisphere (consistent with the Desikan-Killiany-Tourville adult cortical atlas) and 23 sub-cortical structures.

The M-CRIB 2.0 atlases can be used by providing the `-a MCRIB` argument when running the pipeline.

Labels of produced segmentations are provided in the [label_names/MCRIB](label_names/MCRIB) folder.

## License

Draw-EM is distributed under the terms of the Apache License Version 2.
See the accompanying [license file](LICENSE.txt) for details. The license enables usage of
Draw-EM in both commercial and non-commercial applications, without restrictions on the
licensing applied to the combined work.

Draw-EM uses the atlases described in the previous section.

Each atlas is covered by a separate license:
- ALBERTs: [license](label_names/ALBERT/LICENSE.txt)
- M-CRIB 2.0: [license](label_names/MCRIB/LICENSE.txt)

## Releases 
- v1.3: allow segmentation using the M-CRIB 2.0 atlases
- v1.2.1: Corpus Callosum segmentation improvement
- v1.2: dHCP segmentation pipeline, method improvements described in [2]: multi-channel registration, modelling of hyper and hypo-intensities.
- v1.1: initial code release, method described in [1].


## Citation and acknowledgements

In case you found Draw-EM useful please give appropriate credit to the software.

Publications:

1. A. Makropoulos et al. *"Automatic whole brain MRI segmentation of the developing neonatal brain"*, IEEE TMI, 2014
2. A. Makropoulos, E. C. Robinson et al. *"The Developing Human Connectome Project: a Minimal Processing Pipeline for Neonatal Cortical Surface Reconstruction"*, NeuroImage, 2018
