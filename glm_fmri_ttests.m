% % glm_fmri_ttests

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

% irf = 'tent';
% irf = 'habPPI';
% irf = 'spm_hrf';
irf = 'per_trial';

inDir = ['~/ShockAwe/data/results_',irf,'_wo_outliers_8j'];  % relative to subject directory
outDir = ['~/ShockAwe/data/ttests_',irf,'_wo_outliers_8j'];  % relative to subject directory
saveFiles = 1; % if 1, save out nifti and afni files

% maskFile = '~/ShockAwe/data/ROIs_tlrc/group_mask.nii.gz';
maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/mb_group_mask_tlrc.nii.gz';

fileStr = 'shock_corrZ_rlPutamen_p.005';  % get all files from a directory w/this fileStr in file name
volIdx = 1;   %  run ttest on volIdx volume of 4d nifti files

% if fileStr2 is provided, do a ttest for betas of fileStr 1 vs fileStr 2
fileStr2 = 'neutral_corrZ_rlPutamen_p.005';


% get outFile Name
if(fileStr2)
    outNameStr = [fileStr,'-',fileStr2];
else
    outNameStr = [fileStr];
end
if ischar(volIdx)
    outNameStr = [outNameStr,'_vol',volIdx];
elseif (volIdx~=1)
    volIStr = sprintf('%d',volIdx);
    outNameStr = [outNameStr,'_vol',volIStr];
end
outNameStr = 'shock-neutral_corrZ_rlPutamen_p.005';
%% get files

mask = readFileNifti(maskFile);
mask.data = double(mask.data);
idx = find(mask.data);

cd(inDir);

f = dir(['sa*',fileStr,'*nii*']);
N = length(f)
% i=1
for i = 1:N
    nii = readFileNifti(f(i).name);
    subData = reshape(nii.data,prod(nii.dim(1:3)),[]); % each volume of voxels is in a row now
    subData = subData(idx,:);  % get only voxels within the mask
%         d(:,:,i) = subData;
    if(strcmp(volIdx,'peak'))
        [peakB, peakBIdx]=max(subData(:,3:5),[],2);  % only look for peak between 4-8 sec after stim onset
        all_peakBIdx(:,i) = peakBIdx+2; % add 2 to index to make up for not using the first 2 regressors at time 0 or 2
        d(i,:) = peakB;
    else
        d(i,:) = mean(subData(:,volIdx),2); % take mean over columns in case there's more than 1 volIdx
    end
end

if(fileStr2)
    f = dir(['*',fileStr2,'*nii*']);
    if length(f)~=N
        error('the number of cond1 and cond2 files need to be the same!');
    end
    for i = 1:N
        nii = readFileNifti(f(i).name);
        subData = reshape(nii.data,prod(nii.dim(1:3)),[]); % each volume of voxels is in a row now
        subData = subData(idx,:);  % get only voxels within the mask
        if(strcmp(volIdx,'peak'))
            thisIdx=sub2ind(size(subData),(1:length(idx))',all_peakBIdx(:,i)); % get this sub's index for peak vol
            d2(i,:) = subData(thisIdx);
        else
            d2(i,:) = mean(subData(:,volIdx),2); % take mean over columns in case there's more than 1 volIdx
        end
    end
end



%% ttests

if(fileStr2)
    meanStat = mean(d-d2);
    [~,p,~,stats] = ttest(d,d2);
else
    meanStat = mean(d);
    [~,p,~,stats] = ttest(d);
end
tstat = stats.tstat;

tstatMap = mask.data;
tstatMap(idx) = tstat;

meanStatMap = mask.data;
meanStatMap(idx) = meanStat;

pMap = mask.data;
pMap(idx) = p;


%% save out nifti file & convert to afni file
if(saveFiles)
    cd(outDir);
    outNii = makeGlmNifti(mask,outNameStr,meanStatMap,tstatMap);
    writeFileNifti(outNii);
    %      system(sprintf('3dbucket -prefix %s -fbuc %s',outNameStr,[outNameStr,'.nii.gz']));
    %     system(sprintf('3drefit -view ''tlrc'' %s+orig.',outNameStr));
    %     system(sprintf('3drefit -sublabel 0 %s -sublabel 1 %s -substatpar 1 fitt %d %s+tlrc.',...
    %         ['mean_',outNameStr],['Tstat_',outNameStr],N-1,outNameStr));
end

