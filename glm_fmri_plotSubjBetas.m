function figH = glm_fmri_plotSubjBetas(betas, irf)

% INPUTS:
%   betas: cell array of betas to plot; each cell is a group of regressors
% a row for each subject

%   irf: string identifying the regressor type (as of now this needs to be
% either 'tent' or 'spm_hrf'

% OUTPUTS:
%   figH - handle of the plotted figure

% kjh Oct 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = getSAColors();
N = size(betas{1},1);

figH = figure; hold on

for s = 1:N
    
    subplot(3,ceil(N/3),s); hold on
    
    switch irf
        
        case 'tent'
%             t = 0:2:10; % assume this for tent regs
            for c = 1:length(betas)
                plot(betas{c}(s,:),'.-','color',colors(c,:))
            end
            hold off
            
        case {'spm_hrf','afni'}
            for c = 1:length(betas)
                bar(c,betas{c}(s));
            end
            hold off
            
        case 'per_trial'
            for c = 1:length(betas)
                plot(betas{c}(s,:),'.-','color',colors(c,:))
            end
            hold off
    end
    
end % subjects






