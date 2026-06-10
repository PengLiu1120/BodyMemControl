%% 3T MRI Psychophysiological Interactions (PPI) group analysis
% This is the script to run Psychophysiological Interactions (PPI) analysis
% To run this analysis, VOI must be provided
% To run this analysis, PPI subject analysis must be performed beforehand
%% ........................................................................credits
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@uni-tuebingen.de
% Last updated 16th Mar 2026 by P.Liu
%% ........................................................................tidy up
clear all
close all
clc
%% ........................................................................set defaults
% .........................................................................specify data path
dir_base = '/mnt/volume/tacmem_workspace/derivatives/results/group_PPI_analyses';
data_base = '/mnt/volume/tacmem_workspace/derivatives/results/subject_PPI_analyses';
beh_base = '/mnt/volume/tacmem_workspace/derivatives/results/group_behavioural_results';

% .........................................................................specify subjects
sub = {'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'};
%'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30',,'sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'

% .........................................................................specify parameters
mod = '/TNT';
recall = '/recall';
covari_recall = '/recall_covariate';

% .........................................................................specify conditions
mod_contrast = {'think_nothink', 'nothink_think'};
recall_contrast = {'think_baseline', 'nothink_baseline', 'think_nothink', 'nothink_think'};

% .........................................................................specify contrasts
TNT_VOI = {'TNT_r_insula', 'TNT_l_dlPFC', 'TNT_l_IPL', 'TNT_r_precuneus', 'TNT_l_precuneus'};
NTT_VOI = {'NTT_r_dlPFC', 'NTT_r_vlPFC'};
recall_TB_VOI = {'TB_r_SI', 'TB_r_M1', 'TB_l_IPL', 'TB_r_IPL'};

% .........................................................................switch for functions
% .........................................................................TNT
switch_mod = true;

% .........................................................................recall
switch_recall = false;

%% ........................................................................group PPI GLM analyses
% .........................................................................TNT PPI GLM
if switch_mod

    for i_mod = 1:size(TNT_VOI,2)

        curr_contrast = mod_contrast{1};
        curr_VOI = TNT_VOI{i_mod};
        
        pl_modulation_PPI_group_analysis(dir_base, data_base, mod, curr_contrast, curr_VOI, sub)

    end

    for i_mod = 1:size(NTT_VOI,2)

        curr_contrast = mod_contrast{2};
        curr_VOI = NTT_VOI{i_mod};

        pl_modulation_PPI_group_analysis(dir_base, data_base, mod, curr_contrast, curr_VOI, sub)

    end

end

% .........................................................................recall GLM
if switch_recall

    for i_recall = 1:size(recall_TB_VOI,2)

        curr_contrast = recall_contrast{1};
        curr_VOI = recall_TB_VOI{i_recall};

        pl_recall_PPI_group_analysis(dir_base, data_base, recall, curr_contrast, curr_VOI, sub)

    end

end