function stats = glm_fmri_stats(Y,X,regIndx)

% this function fits a model (X) to data (Y) and returns a bunch of stats 
% on the model fit. It's somewhat redundant with the glm_fmri_fit() function, 
% but this returns more stats and runs slower. 

%%%%%%%%% INPUTS:

% Y - vector of matrix of time series
% X - glm model
% regIndx - vector w/length = # of regressors, 0 for baseline and >=1 for 
%       regressors of interest

%%%%%%%%% OUTPUTS:

% a stats struct with the following fields:

% 'B'         - regressor weights estimated using OLS
% 'seB'       - standard error for each beta in B
% 'tB'        - t-stat for each regressor coefficient in B
% 'pB'        - corresponding p-values
% 'df_tB'     - df for t-stats

% 'Rsq'       - Rsquared of the full model vs baseline
% 'Ffull'     - F-stat of the full model vs baseline
% 'pFull'     - p-value of the Fstat
% 'df_Ffull   - df for F-stat (df0-df,df) 

% 'Yhat'      - full model's prediction of time series Y
% 'residuals' - Y-Yhat, the error of the model's prediction

% and if regIndx contains a value > 1,

% 'Fpartial'  - F-stat for the each set of parameters specified by regIndx
% 'pFpartial' - coresponding p-values for partial F-stats
% 'Rsqpartial'- Rsq for each set of parameters specified by RegIndx


% kelly May 2012; revised on 4/12/13 to take more than 1 time series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,p] = size(X);  % number of time points & model parameters
df = N - p;         % df = # of data points - parameters estimated
p0 = length(regIndx(regIndx==0)); % # of baseline params
df0 = N - p0;       % df for baseline model

% check to make sure # of X rows = # of Y rows
t_dim = find(size(Y)==N); % which dimension is time?
if t_dim==2
    Y = Y'; % put time in the first d
elseif t_dim~=1
    error('model & data dimensions don''t match, its ambiguous which dim time is, or data has more than two dimensions\n\n');
end

% fprintf('\nfitting the model...\n\n');

%% fit model to data
    
meanY = repmat(mean(Y),[N,1]); % make matrix w/vals of mean(Y) from each column
SStot = sum( (Y-meanY).^2 );

B = pinv(X'*X)*X'*Y;            % fit model using OLS

Yhat = X*B;                     % model's prediction for Y

SSfull = sum( (Yhat - meanY).^2); % SS explained by full model

residuals = Y - Yhat;      % residuals

SSe = sum( (residuals).^2 );     % residual SS, aka squared error

MSe = SSe./df;                  % mean sq error, aka error variance

seY = sqrt(MSe);                % standard error (sd) of the estimate; sd of Y predicted from X


%% Baseline model stats

Xbase = X(:,regIndx==0); % baseline regressors only

Bbase = pinv(Xbase'*Xbase)*Xbase'*Y;

Yhat0 = Xbase * Bbase; % baseline model's prediction for Y

SSbase = sum( (Yhat0 - meanY).^2); % SS explained by baseline model

p0 = p - length(regIndx(regIndx==0));
df0 = N - p0;


%% Coefficient of Multiple Determination, Rsquared

% indicates how well the full model fits the data
% Rsquared is the proportion of the variation in the data (above
% baseline) that is explained by the full regression model

Rsq_full = 100 * (1 - (SSe ./ SStot)); % Rsquared for the full model 
Rsq_roi = 100 * (1 - (SSe + SSbase) ./ SStot); % Rsquared for regressors of interest only


%% F-test for significance of the full model compared to baseline model

%     H0: data = noise          (baseline model)
%     H1: data = noise + signal (full model)

SSroi = SStot - (SSe + SSbase);
MSe_regsOI = SSroi ./ (df0 - df);  % variance explained by regressors of interest

Ffull = MSe_regsOI ./ MSe; %  F(df0-df, df) distribution under the null

pFfull = 1 - fcdf(abs(Ffull), df0-df, df);    % p-value for f-stat


%% Partial F-stats to test the significance of individual stimuli w/ multiple parameters

% Fpartial = MS difference btwn full and reduced model / MSe of full model, which is:

% (SSe_reduced-SSe) / (df_reduced-df )
% ________________________________________
%             SSe / df

% Does an individual stimulus modeled with a set of parameters explain a
% significant amount of variance in the data?  Compare full model SSe to
% SSe of a model without those parameters.
%
% calculate and return Rsquared for the set of parameters, too.  This might be helpful
% for viewing tent regressors.
%

% F-partial stats if there are sets of regressors as grouped by regIndx
if (length(unique(regIndx(regIndx~=0))) < length(regIndx(regIndx~=0)))
    
    Fpartial_regs = unique(regIndx(regIndx~=0));
    
    i = 1;
    for g = Fpartial_regs
        
        red_idx = find(regIndx~=g);
        X_red = X(:,red_idx);
        B_red = pinv(X_red'*X_red)*X_red'*Y;

        Yhat_red = X_red * B_red;               % reduced model's Y estimate,
        
        SS_red = sum( (Yhat_red - meanY).^2); % SS explained by baseline model

        SSe_red = sum( (Y - Yhat_red).^2 );              % squared error,
        
        df_red = N - p + length(red_idx);    % and df
        
        MSe_diff = (SSe_red - SSe) ./ (df_red - df); % variance explained by the g set of regressors
        
        Fpartial(i,:) = MSe_diff ./ MSe;                  % distributed as F(df_red-df, df)
        
        pFpartial(i,:) = 1 - fcdf(abs(Fpartial(i,:)), df0-df, df); % p-value for f-stats
        
        Rsqpartial(i,:) = 100 * (1 - (SSe / SSe_red ));
        
        i=i+1;
    end 
end


%% t-stats to test the significance of each parameter
% 
% % indicates the significance of each parameter in explaining the
% % variation in the data; to check more than one parameter, use
% % partial-F statistic
% 
% % t-stat for significance of each parameter
% 
% s2matrix = MSe.*pinv(X'*X);     % variance-covariance matrix for the regressors
% s2B = diag(s2matrix);           % each regressor's variance
% 
% seB = sqrt(s2B);                % standard error (sd) of each estimated beta
% 
% tB = B ./ seB;                  % distributed as t on N-p df
% 
% pB = 1 - tcdf(abs(tB), df);     % p-value for t-statistics tB


% fprintf('\ndone fitting model.\n\n');


%% return stats

stats = struct();

stats.B = B;                    % fit coefficients (given as input)

% stats.seB = seB;                % standard error of coefficients (sd of estimate)
% stats.tB = tB;                  % t-stat for each regressor estimate
% stats.df_tB = N-2;              % df for t-stat
% stats.pB = pB;                  % p value for t-stat

stats.Rsq_full = Rsq_full;      % Rsquared for full model
stats.Rsq_roi = Rsq_roi;      % Rsquared for regressors of interest
stats.Ffull = Ffull;            % F-stat for full model vs baseline
stats.df_Ffull = [df0-df,df];   % df for F-stat
stats.pFfull = pFfull;          % p value for F-stat

stats.Yhat = Yhat;              % full model's prediction
stats.residuals = residuals;    % residuals

if exist('Fpartial_regs','var')
    stats.Fpartial = Fpartial;
    stats.pFpartial = pFpartial;
    stats.Rsqpartial = Rsqpartial;
end




