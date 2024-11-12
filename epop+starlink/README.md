Motivation: 
`Bright unintended electromagnetic radiation from second-generation Starlink satellites', Bassa et al. 2024, A&A

e-POP, following Bernhardt et al. 2023, JPP
https://epop.phys.ucalgary.ca/

Flash RRI conjunctions obtained from
https://epop-data.phys.ucalgary.ca/flash_rri_conjunctions.csv
gives datetime of conjunction and NORAD CAT ID of object involved

RRI data obtained from https://epop-data.phys.ucalgary.ca/
These are HDF5 files! 

FFT processing following 
https://epop.phys.ucalgary.ca/wp-content/uploads/2024/10/CASSIOPE-RRI-tutorial-v1.2.ipynb

starlink-id.ipynb: identifies using space-track which RRI conjunctions recorded in supplied .csv were Starlink satellites, outputs .pkl with conjunction datetime, Starlink NORAD CAT ID and name

epop_data.ipynb: obtains RRI lv1 files for all times given in the RRI conjunctions list as .zip files

epop_data_fft.ipynb: extracts these .zip files and generates both dipole spectrograms for each file

epop_fft_plot.ipynb: for list of Starlink conjunctions, generates TCA, minimum distance, and plots spectrogram for time interval of interest