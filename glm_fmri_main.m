
% glm_fmri_main

% Main wrapper script for running a glm on fmri data.  This assumes files
% are organized according to the structure set u

% Before using this script, regressor time series should be generated and
% saved in subjDir/regs directory

% There are four main parts:

%      1) Define relevant subjects, scan runs, stims, etc. to model
%      2) Create regressor time series and a glm design matrix (glm_fmri_mat)
%      3) Fitting the glm model to data (glm_fmri_fit)
%      4) exploring and visualizing model fit statistics (glm_fmri_vis, not yet implemented)

% kelly May 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% give subject ID and relevant file names

% *filenames should be relative to the subject's main directory*

subject = 'sa10';

funcFile = 'afni/all_scaled.nii';

maskFile = 'ROIs/habenula_func_def.nii.gz'; % leave as '' to not use a mask

runs = [1:3]; % scan runs to include

outFileSuffix = ['hab'];

irf = 'spm_hrf'; %options are 'tent' , 'spm_hrf', or 'per_trial'

%%%%% Regressors of interest %%%%%

% files should contain M x n regressor values where M is the # of
% volumes acquired for all runs and n is the # of regressors
% modeled for each stim
regFiles = {['regs/juiceslide_reg_',irf,'.1D'],...
    ['regs/neutralslide_reg_',irf,'.1D'],...
    ['regs/shockslide_reg_',irf,'.1D']};

stims = {'juice','neutral','shock'}; % should correspond to regFiles


%%%%% Baseline/Nuisance regressors %%%%%

motionRegsFile = 'afni/vr.1D'; % leave as '' to not include

% number of polynomial baseline regressors to include for each scan run
nPolyRegs = 3; % 1st is constant, 2nd is linear, 3rd is quadratic, etc.

censor_first_trs = 3; % # of first scans to censor; use [] to censor none


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHOULDNT HAVE TO EDIT BELOW HERE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Design matrix

expPaths = getSAPaths(subject);

[X, regLabels, regIndx] = glm_fmri_mat(expPaths, subject, runs, stims, regFiles,...
    motionRegsFile, nPolyRegs, censor_first_trs);

%% Fit the model to data

func = readFileNifti(funcFile);

% if maskFile is specified, mask the 1st vol of the fmri data w/zeros outside the mask
if (maskFile)
    mask=readFileNifti(fullfile(expPaths.subj, maskFile));
else
    mask.data = ones(size(func.data));
end

% [B, seB, tB, pB, Rsq] = glm_fmri_fit(func.data, X, regIndx, mask.data);
% B = glm_fmri_fit(func.data,X, regIndx, mask.data);

% USE:
% stats = glm_fmri_stats(Y,X,B,regIndx)
% to get Rsq, F, and other stats of model fit

ts = roi_mean_ts(func.data,mask.data);
stats = glm_fmri_stats(ts,X,regIndx);





% %% save results
% 
% cd(expPaths.results);
% 
% fprintf('\nsaving results for subject %s...\n\n', subject);
% 
% switch irf
%     
%     case 'tent'
%         
%         % for tent irfs, its useful to save nifti files w/beta estimates for each
%         % stim, e.g.,:
%         
%         func.dim(4) = (c-b)./(n-1); % so outfiles will have correct time gap btwn tents in the header
%         
%         outNii = makeGlmNifti(func, ['neutralB_irf',outFileSuffix], B(:,:,:,nIndx));
%         writeFileNifti(outNii);
%         
%         outNii = makeGlmNifti(func,['juiceB_irf',outFileSuffix], B(:,:,:,jIndx));
%         writeFileNifti(outNii);
%         
%         outNii = makeGlmNifti(func,['shockB_irf',outFileSuffix], B(:,:,:,sIndx));
%         writeFileNifti(outNii);
%         
%     case 'spm_hrf'
%         
%         % for regressors convolved w hrf, may be more useful to save model
%         % fit stats for all stim in one nifti file, e.g.:
%         
%         outNii = makeGlmNifti(func,['glm_matlab',outFileSuffix],...
%             B(:,:,:,nIndx), tB(:,:,:,nIndx), pB(:,:,:,nIndx),...
%             B(:,:,:,jIndx), tB(:,:,:,jIndx), pB(:,:,:,jIndx),...
%             B(:,:,:,sIndx), tB(:,:,:,sIndx), pB(:,:,:,sIndx));
%         
%         writeFileNifti(outNii);
%         
% end % irf
% 
% 
% 
% 
% 

