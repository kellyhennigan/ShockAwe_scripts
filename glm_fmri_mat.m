function [X, regLabels, regIndx] = glm_fmri_mat(expPaths, subject, runs, ...
    stims, regFiles, motionRegsFile, nPolyRegs, censor_first_trs, TRreg)

% creates a design matrix for specified regressors of interest and baseline
% regressors
%

%%%%%%% INPUTS:

% expPaths - paths to subject's directories
% subject - subject id string
% runs - scan runs to include
% stims - cell array of stim labels corresponding to regFiles
% regFiles - file names containing regressor time series
% motionRegsFile - name of file containing motion regressors to include
% nPolyRegs - number of baseline regressors to include.  1 is a constant, 2
%   for a linear trend, 3 is quadratic, etc.
% censor_first_trs - integer specifying the # of volumes to censor at the
%   beginning of each run. Giving [] or 0 means none will be excluded.


%%%%%%% OUTPUTS:

% X - glm design matrix with a column for each regressor and nScans rows
% regLabels - cell array of string labels corresponding to each reg column
% regIndx - vector w/ 1,2, etc. for each stim specified by stims and 0 for
%   all other (baseline) regs

% also plots the design matrix in a couple figures

% kelly May 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


regLabels = {}; % regressor labels
regIndx = [];   % 1,2, etc. for stim of interest; baseline regs will be 0

cd(expPaths.subj);

fprintf('\n\nMaking design matrix for subject %s...\n', subject);


%% get regressors of interest  

for k = 1:length(stims)
    stimRegs{k} = dlmread(regFiles{k});
    n_regs = size(stimRegs{k},2); % # of regs for this stim
    regLabels = [regLabels, repmat({stims{k}}, [1 n_regs])];
    regIndx = [regIndx, repmat(ones.*k,1,n_regs)];
end

% concatenate across stims
stimRegs = horzcat(stimRegs{:});


%% get motion regs

% assumes motion reg file is from afni output: 9 columns in the following
% order:
%                     n = sub-brick index
%                     roll  = rotation about the I-S axis }
%                     pitch = rotation about the R-L axis } degrees CCW
%                     yaw   = rotation about the A-P axis }
%                       dS  = displacement in the Superior direction  }
%                       dL  = displacement in the Left direction      } mm
%                       dP  = displacement in the Posterior direction }
%                    rmsold = RMS difference between input brick and base brick
%                    rmsnew = RMS difference between output brick and base brick

% ****** so only use columns 2-7 as regressors *******

if motionRegsFile
    motionRegs = dlmread(motionRegsFile);
    motionRegs = motionRegs(:,2:7);  % assumes afni format of motion regs file
    n_regs = size(motionRegs,2); % # of motion regs
    regLabels = [regLabels, repmat({'motion'}, [1 n_regs])];
    regIndx = [regIndx, repmat(0,1,n_regs)];
end


%% baseline polynomial regressors & censor first trs

nPCols = nPolyRegs*numel(runs);   % # of regressor columns for nPolyRegs/run

for j = 1:numel(runs)                 % runs
    
    % get scan triggers to check # of vols/run
    trigFile = fullfile(expPaths.slicetimes, [subject,'_', num2str(runs(j)), '_slicetimes.txt']);
    trigs = dlmread(trigFile,'',1,0);
    
    % polynomial drift regressors
    baseRegs{j,1} = zeros(numel(trigs), nPCols);
    colIndx = (j-1)*nPolyRegs+1;
    for b = 0:nPolyRegs-1
        baseRegs{j,1}(:,colIndx) = (1:numel(trigs)).^b;
        colIndx = colIndx+1;
    end
    
    % censor first trs of each run?
    if censor_first_trs
        censorRegs = zeros(numel(trigs),numel(runs));
        censorRegs(1:censor_first_trs,j) = 1;
        baseRegs{j,1} = [baseRegs{j,1},censorRegs];
    end
      
end

regLabels = [regLabels, repmat({'poly_base'}, [1 nPCols])];
regIndx = [regIndx, repmat(0,1,nPCols)];
if censor_first_trs
    regLabels = [regLabels, repmat({'censor_first_trs'}, [1 numel(runs)])];
    regIndx = [regIndx, repmat(0,1,numel(runs))];
end

% concatenate across runs
baseRegs = vertcat(baseRegs{:});


%%  concatenate all regressors

X = [stimRegs motionRegs baseRegs];

%% include TRs as a regressor? 

if TRreg
    
    TRs = [];
    for j = 1:numel(runs)
        trigFile = fullfile(expPaths.slicetimes, [subject,'_', num2str(runs(j)), '_slicetimes.txt']);
        trigs = dlmread(trigFile,'',1,0);
        runTRs = trigs(2:end)-trigs(1:end-1);
        TRs = [TRs;0;runTRs];
    end
    
    X(:,end+1) = TRs;
    regLabels{end+1}='TRreg';
    regIndx(end+1) = 0;
    
end


%% plot design matrix

% figure
% imagesc(X, [0 1])
% title('GLM Design Matrix')
% % set(gca,'XTickLabel',regLabels)
% colormap(gray)
% 

fprintf('\nDone with design matrix.\n\n');

end




