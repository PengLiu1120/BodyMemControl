function pl_localiser_group_analysis(dir_base, data_base, localiser, body_part, contrast, name, sub)
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
% localiser GLM group result
% .........................................................................
% Written by P.Liu
% Last Updated 1st Oct 2025
%% ........................................................................localiser GLM group analysis
% .........................................................................specify condition
dir_body_part = fullfile(dir_base, localiser, body_part);

dir_contrast = {};

for i = 1:numel(sub)
    sub_id = sub{i};
    sub_dir = fullfile(data_base, localiser, sub_id, contrast);
    scan_string = [sub_dir, ',1'];

    dir_contrast{i, 1} = scan_string;
end

% .........................................................................spm batch initiating
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% .........................................................................spm GLM specification
matlabbatch{1}.spm.stats.factorial_design.dir = {dir_body_part};

matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = dir_contrast;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = name;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

% .........................................................................run the batch
spm_jobman('run', matlabbatch);