# fMRI-phase-toolbox
Toolbox for processing and analysis of spatiotemporal image phase fluctuations during fMRI scans. Presented at ISMRM 2023.
Main script: `/pipeline/fmri_bg_phase.m`

## Requirements:

1. Christian Kames QSM toolbox (QSM.m) available here:
https://github.com/kamesy/QSM.m

2. 4D fMRI Data in either: 

- NIFTI format, magnitude and phase saved independently
- Par/Rec format, but you then must have MRecon access
 

3. Reasonable Computing Power

## Pipeline:

1. Data import: Read in Magnitude image, Phase image and Physiological log files
2. Phase unwrapping 
3. Background field removal (isolation)
4. Reshape data and perform BSS to extract respiratory trace
5. Comparison and Plotting
