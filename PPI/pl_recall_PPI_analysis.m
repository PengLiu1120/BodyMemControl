function pl_recall_PPI_analysis(dir_res, sub, folder, spmmat, VOI, VOI_coord, VOI_radi, weight, PPI_name)
% .........................................................................
% inputs
% .........................................................................
% data directory
% subject
% parameters
% .........................................................................
% outputs
% .........................................................................
% TNT PPI
% .........................................................................
% Written by P.Liu
% Last Updated 11th Mar 2026
%% ........................................................................TNT GLM analysis
% .........................................................................subject spm.mat path
dir_sub_spm = fullfile(dir_res, folder, sub, spmmat);

% .........................................................................spm batch initiating
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% .........................................................................spm PPI specification
matlabbatch{1}.spm.util.voi.spmmat = {dir_sub_spm};
matlabbatch{1}.spm.util.voi.adjust = 0;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = VOI;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = VOI_coord;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = VOI_radi;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';

matlabbatch{2}.spm.stats.ppi.spmmat = {dir_sub_spm};
matlabbatch{2}.spm.stats.ppi.type.ppi.voi(1) = cfg_dep('Volume of Interest:  VOI mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','voimat'));
matlabbatch{2}.spm.stats.ppi.type.ppi.u = weight;
matlabbatch{2}.spm.stats.ppi.name = PPI_name;
matlabbatch{2}.spm.stats.ppi.disp = 1;

% .........................................................................run the batch
spm_jobman('run', matlabbatch);