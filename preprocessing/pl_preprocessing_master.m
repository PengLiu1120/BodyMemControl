%% 3T MRI image preprocessing
% This is the master script for 3T MRI image data preprocessing
% Including both structural and functional images
%% ........................................................................functions
% *step 1: pl_struc_preprocessing*
% structural preprocessing of T1-weighted images
% *step 2: pl_fmap_VDM_calculation*
% calculating field maps using resting state image as reference
% *step 3* pl_fsl_skull_strpping
% functional images skull stripping using fsl
% can be performed before or in parallel with step 1 and step 2
% *step 4: pl_modulation_func_preprocessing*
% TNT func preprocessing
% realignment & unwarp
% coregistration
% normalisation
% smoothing
% *step 5: pl_recall_func_preprocessing*
% recall func preprocessing
% realignment & unwarp
% coregistration
% normalisation
% smoothing
% *step 6: pl_localiser_func_preprocessing*
% localiser func preprocessing
% realignment & unwarp
% coregistration
% normalisation
% smoothing
%% ........................................................................credits
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@uni-tuebingen.de
% Last updated 19 May 2025 by P.Liu
%% ........................................................................tidy up
clear all
close all
clc
%% ........................................................................set defaults
% .........................................................................specify data path
dir_base = '/mnt/volume/tacmem_workspace/derivatives/data';

% .........................................................................specify subjects
sub = {'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'};

% .........................................................................specify parameters
anat = '/anat';
fmap = '/fmap';
func = '/func';

% .........................................................................
T1w = '_T1w.nii,1';
skull_stripped = 'skull_stripped_bias_corrected_T1w.nii';
T1 = '_T1w.nii';

% .........................................................................
phase = '_phasediff.nii,1';
magnitude = '_magnitude1.nii,1';
epi = '_task-TNT_bold.nii,1';
vdm = '_phasediff.nii';

% .........................................................................
mod = '_task-TNT_bold_skull_stripped.nii';
recall = '_task-recall_bold_skull_stripped.nii';
localizer = '_task-localizer_bold_skull_stripped.nii';

% .........................................................................switch for functions
% .........................................................................struc preprocessing
switch_struc = true;

% .........................................................................fieldmap calculation
switch_fmap = true;

% .........................................................................func preprocessing
switch_mod_func = true;
switch_recall_func = true;
switch_localiser_func = true;

%% ........................................................................preprocessing pipeline
for i_sub = 1:size(sub,2)

    curr_sub = sub{i_sub};

    % .....................................................................struc preprocessing
    if switch_struc

        pl_struc_preprocessing(dir_base, curr_sub, anat, T1w);

    end

    % .....................................................................fieldmap calculation
    if switch_fmap

        pl_fmap_VDM_calculation(dir_base, curr_sub, anat, fmap, func, phase, magnitude, epi, T1w);

    end

    % .....................................................................modulation preprocessing
    if switch_mod_func

        pl_modulation_func_preprocessing(dir_base, curr_sub, func, anat, fmap, mod, skull_stripped, vdm, T1);

    end

    % .....................................................................recall preprocessing
    if switch_recall_func

        pl_recall_func_preprocessing(dir_base, curr_sub, func, anat, fmap, recall, skull_stripped, vdm, T1);

    end

    % .....................................................................localiser preprocessing
    if switch_localiser_func

        pl_localiser_func_preprocessing(dir_base, curr_sub, func, anat, fmap, localizer, skull_stripped, vdm, T1);

    end

end