# HUP interictal iEEG atlas

This repository contains code for analyzing short clips of interictal intracranial data, mapping quantitative abnormalities, and predicting epileptogenicity

## Data

Due to the size of the dataset, it is not stored locally in this repository, but can be found at the following link: https://discover.pennsieve.io/datasets/179
Our project also uses the MNI open iEEG atlas dataset, which can be found at the following link: https://mni-open-ieegatlas.research.mcgill.ca/

## Requirements

In addition to matlab, users will need BrainNet viewer to generate renderings: https://www.nitrc.org/projects/bnv/

## Usage

For most general research purposes, using the filtered but otherwise unprocessed data found on Pennsieve Discover may be most useful. 

To illustrate the usage of our methods for mapping quantitative abnormalities, please see the demo script for a simple application of the atlas to a 'test' patient. The end product of this demo is a rendering of predicted abnormalities, similar to those seen in Figure 7 of the paper. 


## Contributing
Please contact the corresponding authors if you would like to contribute data to the atlas. Be sure to see the methods of our paper (Bernabei et al., Brain 2022) and that of Frauscher et al., Brain 2018 for details.

