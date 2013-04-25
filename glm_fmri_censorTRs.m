function [data,X,regLabels,regIndx] = glm_fmri_censorTRs(data,X,regLabels, regIndx)
%
% this function takes in a 4D fMRI dataset or a single time series, a glm model, and a regressor 
% labels as arguments and returns the data and model with the TRs 
% marked in the regIndx as 'censor TRs' removed and the regressor columns 
% of the model w/censor TRs removed.

%%%%%%%%%%%%%%%%%%%%%% INPUTS:

% data - 4D dataset with time series in the 4th d or a single time series
% as a column vector
% X - design matrix of time x regressors
% regLabels - 1 x N cell vector array where the # of elements = the number 
%           of columns in the design matrix X

%%%%%%%%%%%%%%%%%%%%%% OUTPUTS:

% data - data w/volumes to censor removed
% X - design matrix w/rows to censor removed & censor TR columns removed

% kelly Nov 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

censor_reg_idx = strmatch('censor',regLabels);
[censor_vol_idx,~] = find(X(:,censor_reg_idx));

if (size(X,1)==size(data,1))
    data(censor_vol_idx,:)=[];
elseif (size(X,1)==size(data,4))
    data(:,:,:,censor_vol_idx)=[];
else 
    error('rows of the design matrix must match either the 1st or 4th dim of the data');
end
X(censor_vol_idx,:)=[];
X(:,censor_reg_idx)=[];
regLabels(censor_reg_idx)=[];
regIndx(censor_reg_idx)=[];


