function cardiacCorrection_cv(subj, blocks)

% Kimberlee D'Ardenne 7/2005
% corrects intensity values of cardiac gated functional data
% subj is subject ID
% fname is EPI 4-D nifti file
%
% edited by Kelly for Shock & Awe exp files format
% start from any directory
% no need to give file paths (just the file name)

% different from cardiacCorrection in that this solves for the value of t1
% that minimizes the coefficient of variation in A, rather than the
% variance
%
% coefficient of variation of time series A = sqrt(variance(A))./mean(A)
%
%%%%%%%%%%%%%%%%%%%%%%

% get all of subject's directories
saPaths = getSAPaths(subj);

if (~exist(saPaths.cardiacCorr, 'dir'))
    mkdir(saPaths.cardiacCorr);
end

figure; 
title(subj)
hold on

% block loop
for b = 1:numel(blocks)
    
    blockNum = blocks(b);
    
    fprintf('\nWorking on subject %s block %d...\n\n', subj, blockNum);
    
    % epi file
    epi_filepath = fullfile(saPaths.raw_func, ['a_epi', num2str(blockNum), '.nii.gz']);
    f = readFileNifti(epi_filepath);
    n_sl = f.dim(3);  % get number of slices
    
    % trigger times
    trig_filepath = fullfile(saPaths.slicetimes, [subj, '_',num2str(blockNum),'_slicetimes.txt']);
    trigs = dlmread(trig_filepath,'',1,0);
    
    % make sure # of functional volumes and trigger times make sense
    trigs = compareNVolsToTrigs(f.dim(4), trigs);
    
    t_n = trigs(2:end)-trigs(1:end-1); % TRs
    t_n = t_n .* 1000;  % msec
    t_avg = mean(t_n);
    
    %% Estimate T1 values by minimizing the coefficient of variation of A and correct signal with t1 estimate
    
    % for reference, Zhang et al.'s (2006) averaged estimated T1 values were:
    % CSF = 2668 +/- 119 msec
    % gray matter = 1325 +/- 10 msec
    % white matter = 910 +/- 10 msec
   
    sig = f.data(:,:,:,2:end); % don't correct 1st volume acquired
    sig_corr = nan(size(sig));  % the corrected signal
    t1_est = ones(f.dim(1:3)); % estimated t1 values
    
    % define t1 search parameters and options for fminsearch
    t1_0 = 1000;
    t1_min = 200;
    t1_max = 4000;
    options = optimset('MaxFunEvals',1000,'MaxIter',1000); 
    
    for sl = 1:n_sl
        fprintf('\nWorking on slice %d of %d...\n\n', sl, n_sl);
        for i = 1:f.dim(1)
            for j = 1:f.dim(2)
                this_sig = squeeze(sig(i,j,sl,:));    % time series for this voxel
                if (mean(this_sig) < 100)   % then this is a voxel outside the brain, or at least not a voxel of interest
                    this_t1 = 10;
                else
                    % minimize coefficient of variation of signal magnitude based on value of T1
                    this_t1 = fminsearchbnd(@(t1) cvfunc(t1, this_sig, t_n), t1_0, t1_min, t1_max, options);   % @(t1) specifies to return the value for t1 that minimizes the function cvfunc
                end
                if (this_t1 < 250)  % semi-arbitrary value for assuming the voxel is outside of the brain
                    this_t1 = 10;   % if voxel is outside the brain, set the t1 value to 10
                end
                t1_est(i,j,sl) = this_t1;
                
                this_A = this_sig ./ (1 - exp (-t_n ./ this_t1) );  % calculate A (max amp of signal w/out t1 weighting) w/ t1 estimate
                sig_corr(i,j,sl,:) = this_A .* (1 - exp(-t_avg ./ this_t1 ) );   % corrected signal = A * (1 - exp( -TRave / t1 ) )
            end
        end
    end % sl
    
    % plot estimated t1 maps
    subplot(numel(blocks),1,2*b-1)
    imagesc(t1_est(:,:,round(n_sl/2)))
    title([subj, ' block ', num2str(blockNum)])
    subplot(numel(blocks),1,2*b)
    hist(reshape(t1_est,1,[]))
    xlabel('estimated t1 values')
    
    % save estimated t1 map
    cd(saPaths.cardiacCorr);
    t1_file = ['T1_', num2str(blockNum), '_cv.mat'];
    save(t1_file, 't1_est');
    
    % save corrected signal as a nifti
    f_out = f;
    f_out.fname = ['corr_', num2str(blockNum),'_cv.nii.gz'];
    f_out.data(:,:,:,2:end) = sig_corr;
    writeFileNifti(f_out);
    
end % blocks loop

function cv_A = cvfunc(t1, this_sig, t_n)
A = this_sig ./ (1 - exp(-t_n / t1 ));
var_A = var(A);
cv_A = sqrt(var_A)./mean(A);
