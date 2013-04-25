function all_runs = glm_fmri_makeRegs(subj, stim, runs, catRuns)
%
% creates a separate regressor for each trial
%
% event times are convolved with spm's gamma function hrf
%
% INPUTS:
% subj - string that's the name of the subject's directory (e.g., 'sa01')
% stim - string to find and match in stimfile w/onset times
% runs - integers specifying which scan runs to include (e.g., [1:3])
% catRuns - 0 to not concatenate runs, 1 to do so (default is to
% concatenate)
%
% OUTPUTS:,
% all_runs - regressor time series in units of TRs
% also saves out a txt file for each run & all runs unless catRuns is set to 0
%
% NOTE: User should edit first section below to specify desired irf
% parameters, etc.
%
% Kelly, August 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT THIS PART AS NEEDED

irf = 'per_trial'; % current options are 'spm_hrf', 'tent', or 'per_trial'
% this will be used as a suffix for the out files

% irf parameters
switch irf
    
    case {'spm_hrf','per_trial'}
        
        stim_duration =2;  % duration of stim (in seconds) to model; set to 0 for instantaneous impulse
        sample_rate = 0.1;             % sample rate in seconds (upsampled)
        params = [6 16 1 1 10 0 32];         % set parameters for hrf; defaults are: P = [6 16 1 1 6 0 32];
        hrf = spm_hrf(sample_rate,params);  % spm's hrf
        
        
    case 'tent'
        
        b = 0;  %  beginning of time window to model after stim onset
        c = 10;  % end of time window to model after stim onset
        n = 6;   % number of tent functions to model for each event of interest
        params = [b c n];
        % note: timegap btwn regressors = (c-b)./(n-1); this value should be >= TR
        
end

param_str = sprintf('_%d', params);
outFileSuffix = ['_',irf,param_str];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shouldn't have to edit below

% if runs argument not given, process all 3 runs
if (~exist('runs', 'var'))
    runs = [1:3];
end

% concatenate runs by default
if (~exist('catRuns', 'var'))
    catRuns = 1;
end

% relevant directories
saPaths = getSAPaths(subj);

% save regs files in directory defined by saPaths.regs
if (~exist(saPaths.regs, 'dir'))
    mkdir(saPaths.regs);
end

%%%%%%%%%%%
%% get stim times

all_runs = [];

for j = runs
    
    stimFilePath = fullfile(saPaths.stimtimes, [subj,'_',num2str(j),'_stimtimes.csv']);
    trigsFile = fullfile(saPaths.slicetimes,[subj,'_',num2str(j),'_slicetimes.txt']);
    trigs = dlmread(trigsFile,'',1,0);
    
    % get stim onset times
    onsets = getStimOnsets(stimFilePath, stim);
    
    
    %% convolve stim events with hrf, or interpolate for tent irf
    
    switch irf
        
        case 'spm_hrf'
            
            % convert trigger times and onsets from seconds to sample rate units
            trigs = ceil(trigs./sample_rate);
            onsets = round(onsets ./ sample_rate);
            
            t = zeros(trigs(end),1); % define regressor time series
            t(onsets) = 1; % 1 when stim occurs
            
            if stim_duration~=0
                ups_stim_dur = (stim_duration/sample_rate)-1;
                for q = onsets'
                    t(q:q+ups_stim_dur) = 1;
                end
            end
            
            % convolve upsampled time series with hrf
            reg_ts = conv(t, hrf);
            
            % scaled to peak at 1
            reg_ts = (reg_ts./max(reg_ts));
            
            % convert time series into units of TRs
            reg_tr = matchTime(trigs, reg_ts);
            stim_reg = reg_tr';
            
        case 'per_trial'
            
            % convert trigger times and onsets from seconds to sample rate units
            trigs = ceil(trigs./sample_rate);
            onsets = round(onsets ./ sample_rate);
            
            for i = 1:length(onsets)
                
                t = zeros(trigs(end),1);  % define regressor time series
                t(onsets(i)) = 1; % 1 when stim occurs
                
                if stim_duration~=0
                    ups_stim_dur = (stim_duration/sample_rate)-1;
                    t(onsets(i):onsets(i)+ups_stim_dur) = 1;
                end
                
                % convolve upsampled time series with hrf
                reg_ts = conv(t, hrf);
                
                % scaled to peak at 1
                reg_ts = (reg_ts./max(reg_ts));
                
                % convert time series into units of TRs
                reg_tr = matchTime(trigs, reg_ts);
                stim_reg(:,i) = reg_tr';
                
            end
            
        case 'tent'
            
            stim_reg = afniTent(trigs, onsets, params);
            
    end
    
    %% save in regs directory
    
    cd(saPaths.regs);

    switch irf
        
        case {'spm_hrf','tent'}
            outFile = [stim,'_',num2str(j),outFileSuffix];
%             dlmwrite(outFile, stim_reg);
            all_runs = [all_runs; stim_reg];
            
        case 'per_trial'
            
            if j ==1
                all_runs = stim_reg;
            else
                next_stim_reg = [zeros(size(stim_reg,1),size(all_runs,2)), stim_reg];
                all_runs = [all_runs,zeros(size(all_runs,1),size(stim_reg,2))];
                all_runs = [all_runs; next_stim_reg];
            end
    end
    
    clear trigs onsets stim_reg
    
end % runs

% save concatenated runs?
if catRuns==1
    allOutFName = [stim,outFileSuffix];
    dlmwrite(allOutFName, all_runs);
end

fprintf(['\nSaved regressor time series for ',stim,'\n']);

% figure
% imagesc(all_runs)
% colormap(gray)
% title([subj,' ',stim,' ',irf])







