function fig = glm_fmri_plotFit(ts,X,regIndx)

% this function plots a time series, model fit, and regressors of interest 
% as specified w/by nonzero values in regIndx

%%%%% INPUTS

% ts      - raw fmri time series (from a voxel or an roi)
% X       - glm model; # of rows should match ts; 1 column for each regressor
% regIndx - index specifying regressors of interest (>0) and nuisance
%           regressors (0)

%%%%% OUTPUTS 

% a figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cb=cbrewer('qual', 'Set1', 8);
% colormap(cb)

colors = getSAColors();

B = pinv(X'*X)*X'*ts;

TRs = 1:length(ts);

Yhat = X*B;

j=sum(X(:,regIndx==1),2).*mean(B(regIndx==1)); % juice regressor * beta
n=sum(X(:,regIndx==2),2).*mean(B(regIndx==2)); % neutral ""
s=sum(X(:,regIndx==3),2).*mean(B(regIndx==3)); % shock ""


fig = figure
set(gcf,'Position',[1 300 1360 384])
plot(TRs,ts','k')
hold on
plot(TRs,Yhat','b') % total model fit
plot(TRs,j+mean(ts),'color',colors(1,:)) % juice
plot(TRs,n+mean(ts),'color',colors(2,:)) % neutral
plot(TRs,s+mean(ts),'color',colors(3,:)) % shock
legend('time series','model fit','juice','neutral','shock')
% subplot(2,1,2)
% set(gcf,'DefaultAxesColorOrder',colors)
% h=plot(TRs,X(:,regIndx~=0)')
% legend('juice','neutral','shock')






