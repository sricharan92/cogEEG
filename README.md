# EEG Preprocessing pipeline for CogLab

This is a repository created for EEG Preprocessing for Cognition Lab, IISc. The aim is for people using both systems (BioSemi/EGI) to be able to use this common pipeline so that there is homogeneity in the analyses in lab. 

## Getting Started

There are a few toolboxes that are required by this pipeline: 
1. Fieldtrip toolBox (version?)
2. NoiseTools (do we require this?)

The absolute paths for these need to be set to the variables 'ftPath', 'ntPath'. 

Please set important variables in `preprocessing.m` eg: 'subjects', 'orig_fs', 'resample_fs', 'events_req', 'trial_end', 'bp_freq' ... 

Importantly, please set 'homeDir' as well. This is the directory which should have the raw data in '/Raw/'

The structure of data inside this should be
```
- Raw
-- subjectXX
---- Task (if exists) OR this could be the last level
------ if Task, this will be the last level 
```
In the last level, the raw data files exist. Could be `.mff` or `.bdf` according to your system of acquisition. 

## Current order of steps

### EGI

* Load the data
* Divide the data