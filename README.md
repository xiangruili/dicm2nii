# INSTALL BIDS converter
Copy paste in your Matlab command line:
````matlab
% download
websave('dicm2nii.zip','https://github.com/tanguyduval/dicm2nii/archive/master.zip')
unzip dicm2nii.zip
% add to Matlab path
addpath('dicm2nii-master')
% set preferences for BIDS
setpref('dicm2nii_gui_para','rstFmt',3)
setpref('dicm2nii_gui_para','save_json',true)
% run
dicm2nii
````

If you don't have Matlab, download compiled versions:
https://github.com/tanguyduval/dicm2nii/releases

# DICOM to NIfTI conversion, DICOM and NIfTI tools, NIfTI visualization (version 2019.01.12)

# dicm2nii
Convert DICOM into NIfTI. It can also convert PAR/XML/REC, HEAD/BRIK, MGZ and BrainVoyager files into NIfTI.

# nii_tool
Create, load, save NIfTI file. Support both version 1 and 2 NIfTI, and variety of data type.

# nii_viewer
Visualize NIfTI. Can also visualize any file convertible to NIfTI by dicm2nii.

# nii_moco
Perform motion correction on a NIfTI.

# nii_stc
Perform slice timing correction on a NIfTI.

# nii_xform
Transform a NIfTI into different resolution, or into a template space.

# dicm_hdr, dicm_img, dicm_dict
Read DICOM header and image, independent of Matlab Image Processing Toolbox. 

# rename_dicm, sort_dicm, anonymize_dicm
DICOM tools performing the tasks as indicated by the name.
