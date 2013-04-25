function out = matchTime(slicetimes, reg_ts)

% match scan times to points in upsampled event HRF time course
% for GLM on cardiac gated data

% slicetimesFileName is txt file with volume acquisition times
% regs_ts is regressor time series in .1 sec units

% based on Kim D'Ardenne's matchTime script (10/2011)


%%%%%%%%%%%%%%%%

% omitFirstVols = 3; % omit the X dummy scans of each run



out = zeros(size(slicetimes'));

for i = 1:length(slicetimes)
    n = slicetimes(i);
    if (n > numel(reg_ts))
        out(i) = 0;
    else
        out(i) = reg_ts(n);
    end
end

