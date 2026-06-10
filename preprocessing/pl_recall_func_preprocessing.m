function pl_recall_func_preprocessing(datadir, sub, func, anat, fmap, recall, skull_stripped, vdm, T1)
% .........................................................................
% inputs
% .........................................................................
% skull stripped focus map bold data
% .........................................................................
% outputs
% .........................................................................
% realigned coregistered normalised smoothed focus map bold data
% .........................................................................
% Written by P.Liu
% Last Updated 25 Nov 2024 by P.Liu
%% ........................................................................focus map preprocessing
func_file = fullfile(datadir, sub, func, [sub, recall]);
anat_file = fullfile(datadir, sub, anat, skull_stripped);
fmap_file = fullfile(datadir, sub, fmap, ['vdm5_sc', sub, vdm]);
def_field = fullfile(datadir, sub, anat, ['y_', sub, T1]);

% .........................................................................processing parameters
realign_quality = 0.9;
realign_sep = 4;
realign_fwhm = 5;
smoothing_kernel = [4 4 4]; % FWHM in mm
normalization_voxel_size = [2 2 2];
bounding_box = [-78 -112 -70; 78 76 85];

% .........................................................................spm batch initiating
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% .........................................................................realign and unwarp
matlabbatch{1}.spm.spatial.realignunwarp.data.scans = cellstr(spm_select('expand',func_file));
matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = {fmap_file};
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = realign_quality;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = realign_sep;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = realign_fwhm;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

% .........................................................................coregistration
matlabbatch{2}.spm.spatial.coreg.estimate.ref = {anat_file};
matlabbatch{2}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{2}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% .........................................................................normalization
matlabbatch{3}.spm.spatial.normalise.write.subj.def = {def_field};
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = bounding_box;
matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = normalization_voxel_size;
matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 7;
matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';

% .........................................................................smoothing
matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{4}.spm.spatial.smooth.fwhm = smoothing_kernel;
matlabbatch{4}.spm.spatial.smooth.dtype = 0;
matlabbatch{4}.spm.spatial.smooth.im = 0;
matlabbatch{4}.spm.spatial.smooth.prefix = 's';

% .........................................................................run the batch
spm_jobman('run', matlabbatch);

end