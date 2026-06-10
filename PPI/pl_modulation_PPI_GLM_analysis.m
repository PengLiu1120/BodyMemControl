function pl_modulation_PPI_GLM_analysis(dir_base, dir_res, dir_PPI_res, sub, PPI_result, func, prefix, folder, VOI, contrast, TNT, regressor_3)
% .........................................................................
% inputs
% .........................................................................
% data directory
% subject
% session
% parameters
% .........................................................................
% outputs
% .........................................................................
% modulation PPI GLM
% .........................................................................
% Written by P.Liu
% Last Updated 11th Mar 2026
%% ........................................................................modulation PPI GLM analysis
% .........................................................................subject result path
dir_sub_res = fullfile(dir_PPI_res, folder, sub, contrast, VOI);

% .........................................................................load subject PPI result
PPI_result_path = fullfile(dir_res, folder, sub, [PPI_result '.mat']);
load(PPI_result_path);

%  ........................................................................subject data
func_file = fullfile(dir_base, sub, func, [prefix, sub, TNT]);

% .........................................................................spm batch initiating
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% .........................................................................spm PPI GLM specification
matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_sub_res};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('expand',func_file));

matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-Interaction';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;

matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'BOLD';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = regressor_3;
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;
%%
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PPI-Interaction';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

% .........................................................................run the batch
spm_jobman('run', matlabbatch);