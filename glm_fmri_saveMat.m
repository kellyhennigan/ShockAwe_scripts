% make a glm design matrix and save it out as a mat file along with
% regressor labels and an index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files, params, etc

subjects = getSASubjects('fmri');
% subjects={'sa26'};
% subjects = getSASubjects('hab');

runs = [1:3]; % scan runs to include

% irf = 'tent'; % spm_hrf, per_trial, or tent
% irfParamsStr = '_0_10_6';

% irf = 'spm_hrf';
irf='per_trial';
irfParamsStr = '_6_16_1_1_10_0_32';
% irfParamsStr = '_8_16_1_1_10_0_32';

outMatName = ['glm_',irf, irfParamsStr];
% outMatName = 'glm_habPPI.mat';

%%%%% Regressors of interest %%%%%

regFiles = {['regs/juiceslide_',irf, irfParamsStr],...
    ['regs/neutralslide_',irf, irfParamsStr],...
    ['regs/shockslide_',irf, irfParamsStr]};

stims = {'juice','neutral','shock'}; % should correspond to regFiles


%%%%% Baseline/Nuisance regressors %%%%%

motionRegsFile = 'afni/vr.1D'; % leave as '' to not include

nPolyRegs = 3; % 1st is constant, 2nd is linear, 3rd is quadratic, etc.

censor_first_trs = 3; % # of first scans to censor; use [] to censor none

TRreg = 1; % 1 to include TRs as a regressor
if TRreg
    outMatName = [outMatName,'wTRs'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subject loop

for s = 1:length(subjects)
%     for s=11
    subject = subjects{s};
    
    if strcmp(subject,'sa26')
        runs = [1,3];
    else
        runs = [1:3];
    end
    
    expPaths = getSAPaths(subject);
    
    %% Design matrix
    
    [X, regLabels, regIndx] = glm_fmri_mat(expPaths, subject, runs, stims, regFiles,...
        motionRegsFile, nPolyRegs, censor_first_trs, TRreg);
    
    %% save it
    
    cd(expPaths.design_mats); 
    save([outMatName,'.mat'],'X','regLabels','regIndx');
    
    fprintf(['\n\nfinished subject', subject]);
    
end % subjects