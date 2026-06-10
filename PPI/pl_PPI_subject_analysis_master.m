%% 3T MRI Psychophysiological Interactions (PPI) subject analysis
% This is the script to run Psychophysiological Interactions (PPI) analysis
% for each contrast of each subject
% To run this analysis, VOI must be provided
%% ........................................................................credits
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@uni-tuebingen.de
% Last updated 11th Mar 2026 by P.Liu
%% ........................................................................tidy up
clear all
close all
clc
%% ........................................................................set defaults
% .........................................................................specify data path
dir_base = '/mnt/volume/tacmem_workspace/derivatives/data';
dir_res = '/mnt/volume/tacmem_workspace/derivatives/results/subject_analyses';
dir_PPI_res = '/mnt/volume/tacmem_workspace/derivatives/results/subject_PPI_analyses';

% .........................................................................specify subjects
sub = {'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'};
%'sub-03','sub-05','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-20','sub-24','sub-27','sub-28','sub-29','sub-30','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-39','sub-40','sub-42','sub-43','sub-44','sub-47'

% .........................................................................specify parameters
func = '/func';
modfolder = '/TNT';
recallfolder = '/recall';
spmmat = '/SPM.mat';
prefix = 'swu';
mod = '_task-TNT_bold_skull_stripped.nii';
recall = '_task-recall_bold_skull_stripped.nii';

mod_contrast = {'think_nothink', 'nothink_think'};
mod_weight  = {[1 1 1; 2 1 -1], [1 1 -1;2 1 1]};

TNT_VOI = {'TNT_r_insula', 'TNT_l_dlPFC', 'TNT_l_IPL', 'TNT_r_precuneus', 'TNT_l_precuneus'};
TNT_VOI_coord = {[36 22 -4], [-46 28 30], [-38 -46 46], [18, -60, 24], [-14 -62 26]};
TNT_VOI_radi = {6, 6, 6, 6, 6};

NTT_VOI = {'NTT_r_dlPFC', 'NTT_r_vlPFC'};
NTT_VOI_coord = {[34 32 46], [54 26 2]};
NTT_VOI_radi = {6, 6};

recall_contrast = {'think_baseline'};
recall_weight = {[1 1 1; 3 1 -1]};

recall_TB_VOI = {'TB_r_SI', 'TB_r_M1', 'TB_l_IPL', 'TB_r_IPL'};
recall_TB_VOI_coord = {[30, -28, 46], [22, -26, 66], [-50, -38, 26], [56, -34, 48]};
recall_TB_VOI_radi = {6, 6, 6, 6};

% .........................................................................switch for functions
% .........................................................................modulation PPI
switch_mod_PPI = false;

% .........................................................................modulation PPI GLM
switch_mod_PPI_GLM = false;

% .........................................................................recall PPI
switch_recall_PPI = true;

% .........................................................................recall PPI GLM
switch_recall_PPI_GLM = true;

%% ........................................................................subject PPI GLM analyses
for i_sub = 1:size(sub,2)

    curr_sub = sub{i_sub};

    for i_VOI = 1:size(TNT_VOI,2)

        curr_VOI = TNT_VOI{i_VOI};
        curr_VOI_coord = TNT_VOI_coord{i_VOI};
        curr_VOI_radi = TNT_VOI_radi{i_VOI};

        curr_contrast = mod_contrast{1};
        curr_weight = mod_weight{1};
        PPI_name = [curr_VOI '_' curr_contrast];
        regressor_3 = ['Psych_' curr_contrast];
        PPI_result = ['PPI_' curr_VOI '_' curr_contrast];

        % .................................................................TNT PPI
        if switch_mod_PPI

            pl_modulation_PPI_analysis(dir_res, curr_sub, modfolder, spmmat, curr_VOI, curr_VOI_coord, curr_VOI_radi, curr_weight, PPI_name);

        end

        % .................................................................TNT PPI GLM
        if switch_mod_PPI_GLM

            pl_modulation_PPI_GLM_analysis(dir_base, dir_res, dir_PPI_res, curr_sub, PPI_result, func, prefix, modfolder, curr_VOI, curr_contrast, mod, regressor_3)

        end

    end

    for i_VOI = 1:size(NTT_VOI,2)

        curr_VOI = NTT_VOI{i_VOI};
        curr_VOI_coord = NTT_VOI_coord{i_VOI};
        curr_VOI_radi = NTT_VOI_radi{i_VOI};

        curr_contrast = mod_contrast{2};
        curr_weight = mod_weight{2};
        PPI_name = [curr_VOI '_' curr_contrast];
        regressor_3 = ['Psych_' curr_contrast];
        PPI_result = ['PPI_' curr_VOI '_' curr_contrast];

        % .................................................................TNT PPI
        if switch_mod_PPI

            pl_modulation_PPI_analysis(dir_res, curr_sub, modfolder, spmmat, curr_VOI, curr_VOI_coord, curr_VOI_radi, curr_weight, PPI_name);

        end

        % .................................................................TNT PPI GLM
        if switch_mod_PPI_GLM

            pl_modulation_PPI_GLM_analysis(dir_base, dir_res, dir_PPI_res, curr_sub, PPI_result, func, prefix, modfolder, curr_VOI, curr_contrast, mod, regressor_3)

        end

    end

    for i_VOI = 1:size(recall_TB_VOI,2)

        curr_VOI = recall_TB_VOI{i_VOI};
        curr_VOI_coord = recall_TB_VOI_coord{i_VOI};
        curr_VOI_radi = recall_TB_VOI_radi{i_VOI};

        curr_contrast = recall_contrast{1};
        curr_weight = recall_weight{1};
        PPI_name = [curr_VOI '_' curr_contrast];
        regressor_3 = ['Psych_' curr_contrast];
        PPI_result = ['PPI_' curr_VOI '_' curr_contrast];

        % .................................................................recall PPI
        if switch_recall_PPI

            pl_recall_PPI_analysis(dir_res, curr_sub, recallfolder, spmmat, curr_VOI, curr_VOI_coord, curr_VOI_radi, curr_weight, PPI_name);

        end

        % .................................................................recall PPI GLM
        if switch_recall_PPI_GLM

            pl_recall_PPI_GLM_analysis(dir_base, dir_res, dir_PPI_res, curr_sub, PPI_result, func, prefix, recallfolder, curr_VOI, curr_contrast, recall, regressor_3)

        end

    end

end