function M = afniTent(trigs,onsets,irf_params)

% Estimate Impulse response function using a continuous piecewise linear regression

% based on AFNI's tent function:
% TENT(b,c,n): A tent function deconvolution model,
% ranging between times s + b and s + c after each stimulus time s,
% with n basis functions (and n regression parameters per voxel).
% A tent function is just the colloquial term for a linear B-spline.
% That is tent(x) = max(0, |1-x| )

% INPUTS:

% trigs - volume acquisition times
% onsets - vector of stim onset times (relative to trigs times)

% irf_params - 1x3 vector [b,c,n] where:
    % b -  beginning of time window to model relative to stim onset
    % c -  end of time window to model after stim onset
    % n -  number of tent regressors per event
  
% OUTPUTS: 

% M - regressor matrix of numel(trigs) rows and nKnot columns

% note that for the first and last 'knot', only right and left half of the
% tent regressors will be estimated respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = irf_params(1);
c = irf_params(2);
n = irf_params(3);

L = (c-b)./(n-1);
knot_times = b:L:c;

M = zeros(numel(trigs), n); % matrix for event of interest regressors

for m = 1:numel(onsets) % onsets
    
    stim_onset = onsets(m);
    t = trigs - stim_onset; % trig time relative to stim_onset
    
    for p = 1:numel(knot_times)      % knot_times
    
        knot_time = knot_times(p);
        x = (t - knot_time)./ L;    % tent(x)
        if knot_time==b
            indx = find(x > 0 & x < 1);
        elseif knot_time==c
            indx = find(x > -1 & x < 0);
        else
            indx = find(x > -1 & x < 1);
        end
        regVals = 1 - abs(x(indx));
        M(indx,p)= regVals;
        
    end % knot_times
    
end % onsets

end