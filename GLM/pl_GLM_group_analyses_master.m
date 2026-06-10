%% 3T MRI General Linear Model (GLM) group analysis
% This is the master script for 3T MRI GLM group analysis
% It includes the GLM analyses for TNT phase, recall phase and functional localiser
% To run this analysis, GLM subject analysis must be performed beforehand
%% ........................................................................credits
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@uni-tuebingen.de
% Last updated 26th Feb 2026 by P.Liu
%% ........................................................................tidy up
clear all
close all
clc
%% ........................................................................set defaults
% .........................................................................specify data path
dir_base = '/mnt/volume/tacmem_workspace/derivatives/results/group_analyses';
data_base = '/mnt/volume/tacmem_workspace/derivatives/results/subject_analyses';
beh_base = '/mnt/volume/tacmem_workspace/derivatives/results/group_behavioural_results';

% .........................................................................specify subjects
sub = {'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'};

% .........................................................................specify parameters
localiser = '/localiser';
mod = '/TNT';
recall = '/recall';
covari_recall = '/recall_covariate';
body_parts = {'/left_hand', '/right_hand','/left_foot', '/right_foot', '/hand', '/foot'};
body_part_names = {'left hand', 'right hand','left foot', 'right foot', 'hand', 'foot'};
mod_conditions = {'/think-nothink', '/nothink-think'};
mod_names = {'think-nothink', 'nothink-think'};
recall_conditions = {'/think-baseline', '/nothink-baseline', '/think-nothink', '/nothink-think'};
recall_names = {'think-baseline', 'nothink-baseline', 'think-nothink', 'nothink-think'};

% .........................................................................specify contrasts
body_part_contrast = {'con_0001.nii', 'con_0002.nii','con_0003.nii', 'con_0004.nii', 'con_0005.nii', 'con_0006.nii'};
mod_contrast = {'con_0001.nii', 'con_0002.nii'};
recall_contrast = {'con_0001.nii', 'con_0002.nii','con_0003.nii', 'con_0004.nii'};

% .........................................................................read accuracy
think_accuracy = load(fullfile(beh_base, 'think_accuracy.acc'));
baseline_accuracy = load(fullfile(beh_base, 'baseline_accuracy.acc'));
nothink_accuracy  = load(fullfile(beh_base, 'nothink_accuracy.acc'));

think_baseline_acc = think_accuracy - baseline_accuracy;
nothink_baseline_acc = nothink_accuracy - baseline_accuracy;
think_nothink_acc = think_accuracy - nothink_accuracy;
nothink_think_acc = think_accuracy - nothink_accuracy;

acc_names = {think_baseline_acc, nothink_baseline_acc, think_nothink_acc, nothink_think_acc};

% .........................................................................switch for functions
% .........................................................................TNT
switch_mod = false;

% .........................................................................recall
switch_recall = false;

% .........................................................................recall covariate
switch_covari_recall = true;

% .........................................................................localiser
switch_localiser = false;

%% ........................................................................GLM group analyses
% .........................................................................TNT GLM
if switch_mod

    for i_mod = 1:size(mod_conditions,2)

        curr_mod = mod_conditions{i_mod};
        curr_contrast = mod_contrast{i_mod};
        curr_name = mod_names{i_mod};

        pl_modulation_group_analysis(dir_base, data_base, mod, curr_mod, curr_contrast, curr_name, sub)

    end

end

% .........................................................................recall GLM
if switch_recall

    for i_recall = 1:size(recall_conditions,2)

        curr_recall = recall_conditions{i_recall};
        curr_contrast = recall_contrast{i_recall};
        curr_name = recall_names{i_recall};

        pl_recall_group_analysis(dir_base, data_base, recall, curr_recall, curr_contrast, curr_name, sub)

    end

end

% .........................................................................covariate recall GLM
if switch_covari_recall

    for i_covari_recall = 1:size(recall_conditions,2)

        curr_covari_recall = recall_conditions{i_covari_recall};
        curr_contrast = recall_contrast{i_covari_recall};
        curr_name = recall_names{i_covari_recall};
        curr_acc = acc_names{i_covari_recall};

        pl_recall_covariate_group_analysis(dir_base, data_base, recall, covari_recall, curr_covari_recall, curr_contrast, curr_name, sub, curr_acc)

    end

end

% .........................................................................localiser GLM
if switch_localiser

    for i_body = 1:size(body_parts,2)

        curr_body = body_parts{i_body};
        curr_contrast = body_part_contrast{i_body};
        curr_name = body_part_names{i_body};

        pl_localiser_group_analysis(dir_base, data_base, localiser, curr_body, curr_contrast, curr_name, sub)

    end

end