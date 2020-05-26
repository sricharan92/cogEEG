# EEG Preprocessing pipeline for CogLab

This is a repository created for EEG Preprocessing for Cognition Lab, IISc. The aim is for people using both systems (BioSemi/EGI) to be able to use this common pipeline so that there is homogeneity in the analyses in lab. 

## Getting Started

There are a few toolboxes that are required by this pipeline: 
1. Fieldtrip toolBox (version?)
2. NoiseTools (do we require this?)

The absolute paths for these need to be set to the variables 'ftPath', 'ntPath'. 

Please set the variables like 'subjects', 'orig_fs', 'resample_fs', 'events_req', 'trial_end' etc. in `preprocessing.m`

## Current order of steps

### EGI

* Load the data
* Divide the data