function pl_modulation_GLM_analysis(dir_base, dir_res, sub, beh, func, prefix, folder, mod, cond_names, contrast, c_con)
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
% modulation GLM
% .........................................................................
% Written by P.Liu
% Last Updated 26th May 2025
%% ........................................................................modulation GLM analysis
% .........................................................................subject onset path
dir_sub_beh = fullfile(dir_base, sub, beh);

% .........................................................................subject result path
dir_sub_res = fullfile(dir_res, folder, sub);

%  ........................................................................subject data
func_file = fullfile(dir_base, sub, func, [prefix, sub, mod]);

% .........................................................................spm batch initiating
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% .........................................................................spm GLM specification
matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_sub_res};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('expand',func_file));

for con=1:size(cond_names,2)

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).name = cond_names{con};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).onset = (load(fullfile(dir_sub_beh,[cond_names{con} '.ons'])));
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).duration = (load(fullfile(dir_sub_beh,[cond_names{con} '.dur'])));
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(con).orth = 1;

end

matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
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

for c=1:size(contrast,2)

    matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = contrast{c};
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = c_con(c,:);
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'none';

end

matlabbatch{3}.spm.stats.con.delete = 0;

% .........................................................................run the batch
spm_jobman('run', matlabbatch);