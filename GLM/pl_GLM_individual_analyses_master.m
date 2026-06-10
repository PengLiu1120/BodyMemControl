%% 3T MRI General Linear Model (GLM) individual analysis
% This is the master script for 3T MRI GLM individual analysis for visualisation in Freesufer
% It includes the GLM analyses for TNT phase, recall phase and functional localiser
% To run this analysis, .dur and .ons files must be provided
%% ........................................................................credits
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@uni-tuebingen.de
% Last updated 18th Nov 2025 by P.Liu
%% ........................................................................tidy up
clear all
close all
clc
%% ........................................................................set defaults
% .........................................................................specify data path
dir_base = '/mnt/volume/tacmem_workspace/derivatives/data';
dir_res = '/mnt/volume/tacmem_workspace/derivatives/results/individual_analyses';

% .........................................................................specify subjects
sub = {'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'};

% .........................................................................specify parameters
func = '/func';
beh = '/beh';
modfolder = '/TNT';
recallfolder = '/recall';
localiserfolder = '/localiser';
prefix = 'u';
mod = '_task-TNT_bold_skull_stripped.nii';
recall = '_task-recall_bold_skull_stripped.nii';
localiser = '_task-localizer_bold_skull_stripped.nii';

% .........................................................................specify conditions
mod_cond = {'Think', 'NoThink'};
recall_cond = {'Think_Recall', 'NoThink_Recall', 'Baseline_Recall'};
localiser_cond = {'Left_Hand', 'Right_Hand','Left_Foot', 'Right_Foot'};

% .........................................................................specify contrasts
mod_contrast = {'Think>NoThink', 'NoThink>Think'};
mod_con  = [1 -1; -1 1];

recall_contrast = {'Think>Baseline', 'NoThink>Baseline', 'Think>NoThink', 'NoThink>Think'};
recall_con  = [1 0 -1; 0 1 -1; 1 -1 0; -1 1 0];

localiser_contrast = {'Left_Hand', 'Right_Hand','Left_Foot', 'Right_Foot', 'Hand', 'Foot'};
localiser_con = [3 -1 -1 -1; -1 3 -1 -1; -1 -1 3 -1; -1 -1 -1 3; 1 1 -1 -1; -1 -1 1 1];

% .........................................................................switch for functions
% .........................................................................modulation
switch_mod = false;

% .........................................................................recall
switch_recall = true;

% .........................................................................localiser
switch_localiser = false;

%% ........................................................................GLM subject analyses
for i_sub = 1:size(sub,2)

    curr_sub = sub{i_sub};

    % .....................................................................modulation GLM 
    if switch_mod

        pl_modulation_GLM_analysis(dir_base, dir_res, curr_sub, beh, func, prefix, modfolder, mod, mod_cond, mod_contrast, mod_con);

    end

    % .....................................................................recall GLM
    if switch_recall

        pl_recall_GLM_analysis(dir_base, dir_res, curr_sub, beh, func, prefix, recallfolder, recall, recall_cond, recall_contrast, recall_con);

    end

    % .....................................................................localiser GLM
    if switch_localiser

        pl_localiser_GLM_analysis(dir_base, dir_res, curr_sub, beh, func, prefix, localiserfolder, localiser, localiser_cond, localiser_contrast, localiser_con);

    end

end