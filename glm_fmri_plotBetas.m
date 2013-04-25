function figH = glm_fmri_plotBetas(betas, stims, irf)

% plots the average beta values for different conditions

% INPUTS:
%   betas: cell array of betas to plot; each cell is a group of regressors
% a row for each subject

%   irf: string identifying the regressor type (as of now this needs to be
% either 'tent' or 'spm_hrf'

% p_vals (optional) - idx w/same # of columns as betas indicating which regressors to plot a '*' over

% OUTPUTS:
%   figH - handle of the plotted figure

% kjh Oct 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p_vals','var')
    p_vals = [];
end

colors = getSAColors();

N = size(betas{1},1);   % number of subjects


figH = figure; hold on

switch irf
    
    case 'tent'
        %         t = 0:2:10; % assume this for tent regs
        for c = 1:length(betas)
            errorbar(mean(betas{c},1),std(betas{c},1)./sqrt(N),'color',colors(c,:))
        end
        yL = ylim;
        legend(stims);
        hold off

%         if ~isempty(p_vals)
%             if(find(p_vals<=.001))
%                 idx=find(p_vals<.001);
%                 for i=1:length(idx)
%                     text(idx(i),yL(2)+range(yL)*.1,'***','FontSize',24,'HorizontalAlignment','center');
%                 end
%             end
%             if(find(p_vals>.001 & p_vals<=.01))
%                 idx = find(p_vals>.001 & p_vals<=.01);
%                 for i=1:length(idx)
%                     text(idx(i),yL(2)+range(yL)*.1,'**','FontSize',24,'HorizontalAlignment','center');
%                 end
%             end
%             if(find(p_vals>.01 & p_vals<=.05))
%                 idx = find(p_vals>.01 & p_vals<=.05);
%                 for i=1:length(idx)
%                     text(idx(i),yL(2)+range(yL)*.1,'*','FontSize',24,'HorizontalAlignment','center');
%                 end
%             end
%         end
   
        
    case {'spm_hrf','afni'}
        for c = 1:length(betas)
            bar(c,mean(betas{c}))
            plot([c,c],[mean(betas{c})+std(betas{c})./sqrt(N),mean(betas{c})-std(betas{c})./sqrt(N)])
        end
        colormap(summer)
        set(gca,'XTickLabel',['.';'j';'.';'n';'.';'s';'.'])
        hold off
        
    case 'per_trial'
        for c = 1:length(betas)
            errorbar(1:length(betas{c}),nanmean(betas{c}),nanstd(betas{c})./sqrt(N),'color',colors(c,:))
        end
        legend('juice','neutral','shock');
        hold off
        
end
title(['mean irf']);
ylabel('Percent signal change')






