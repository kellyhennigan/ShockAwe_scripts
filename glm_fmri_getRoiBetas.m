function B = glm_fmri_getRoiBetas(roiFiles,roiStrs,irf)
% get an roi's betas from a nifti volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files, params, etc

clear B
cd '~/ShockAwe/data/'

if ~iscell(roiFiles)
    roiFiles = {roiFiles};
end

if ~iscell(roiStrs)
    roiStrs = {roiStrs};
end

nR = length(roiFiles); % number of rois

subjects = getSASubjects('fmri');

stims = {'juice','neutral','shock'}; % should correspond to regFiles

outDir = '/home/kelly/ShockAwe/data/betas';

cd(['/home/kelly/ShockAwe/data/results_',irf]);

if strcmp(irf,'per_trial')
oCount = zeros(N,length(stims),nR); % outlier beta count (only used for per trial irf)
 olVals = [];
 olZVals = [];
end

B = [];
 
% stim loop
for c=1:length(stims)
    
    
    if strcmp(irf,'per_trial')
        roiBetas{1} = nan(N,22); roiBetas = repmat(roiBetas,nR,1); % nan vals for betas
    end
    
    
    %% subject loop
    for s = 1:length(subjects)
        
        subject = subjects{s};
        
        fprintf(['\nWorking on subject ',subject,'...\n\n']);
        
        f=dir([subject,'*',stims{c},'*']);
        nii = readFileNifti(f.name);
        allBetas = reshape(nii.data,prod(nii.dim(1:3)),[]);
        if strcmp(irf,'spm_hrf')
            allBetas = allBetas(:,1); % take only betas, not t vals
        end
        
        for j = 1:nR
            
            roi = readFileNifti(['/home/kelly/ShockAwe/data/ROIs_tlrc/',roiFiles{j}]);
            roiIdx = find(roi.data); % get roi coords index
            
            sub_voxBetas = allBetas(roiIdx,:);
            
            if strcmp(irf,'per_trial')
                zscoreBs = zscore(sub_voxBetas);
                oidx = find(abs(zscoreBs)>3); 
                olVals = [olVals;sub_voxBetas(oidx)];
                olZVals = [olZVals;zscoreBs(oidx)];
                sub_voxBetas(oidx)=nan; % exclude outlier voxel values
                oCount(s,c,j)=oCount(s,c,j)+length(oidx);
            end
            
            sub_roiBetas = nanmean(sub_voxBetas);
            
            roiBetas{j}(s,1:length(sub_roiBetas)) = sub_roiBetas;
            
        end % roiFiles
        
    end % subjects
    
    for j = 1:nR
        outFName = fullfile(outDir,[roiStrs{j},'_',stims{c},'_',irf,'_betas']);
        dlmwrite(outFName,roiBetas{j});
    end
%     
    B = [B, roiBetas]; % all roi betas
    
    clear roiBetas
    
end % stims

