% make a design matrix, fit the model to data and save out the results

% handy if the same model/fit are being used a lot to avoid having to
% re-estimate the model every time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

subjects = getSASubjects('fmri');

irf='per_trial';
% irf = 'spm_hrf';
irf_param_str = '_6_16_1_1_10_0_32';
% irf = 'tent';
% irf_param_str = '_0_10_6';
% irf = 'habPPI';

gSpace = 'tlrc';  % '' for native space

matName = ['design_mats/glm_',irf,irf_param_str,'wTRs.mat'];
% matName = 'design_mats/glm_habPPI.mat';

funcFile = ['afni/all_scaled_s',gSpace,'.nii.gz'];

maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/Striatum.nii.gz';
% maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/midbrain_func_tlrc.nii.gz';
% maskFilePath = 'ROIs/DA_func.nii.gz';

% outDir = ['/home/kelly/ShockAwe/data/results_',irf];  % relative to subject directory
% outDir = 'results';  % within subject directory
outDir = ['/home/kelly/ShockAwe/data/bs_wo_outliers'];  % relative to subject directory



stim = {'juice','neutral','shock'};
% stim = {'hab_ts','habXjuice','habXneutral','habXshock'}; % should correspond to regFiles

if(strcmp(irf,'per_trial'))
    omitOutlierTrials = 1; % 1 to omit outlier trials, 0 otherwise
    omitCount = zeros(length(subjects),length(stim));    % subjects as rows, stim as columns, value denotes # of omitted trials
    
    seedStr = 'rVTA_svc';
    cd ~/ShockAwe/data/betas
    for c=1:length(stim)
%         if strmatch(stim{c},'juice')
%             seedBetas{c} = dlmread([seedStr,'_',stim{c},'_per_trial_8_betas']);
%         else
        seedBetas{c} = dlmread([seedStr,'_',stim{c},'_per_trial_betas']);
%         end
    end
    
    seedStr2 = 'llatSN_svc';
    for c=1:length(stim)
%         if strmatch(stim{c},'juice')
%             seedBetas{c} = dlmread([seedStr,'_',stim{c},'_per_trial_8_betas']);
%         else
            seedBetas2{c} = dlmread([seedStr2,'_',stim{c},'_per_trial_betas']);
%         end
    end

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subject loop

for s =4:length(subjects)
    
    subject = subjects{s};
    if ~(strcmp(subject,'sa32') && strcmp(irf,'per_trial'))
        
        expPaths = getSAPaths(subject);
        cd(expPaths.subj);
        
        % get design matrix & data
        load(matName);
        func = readFileNifti(funcFile);
        mask=readFileNifti(maskFile);
        %         mask=readFileNifti(maskFilePath);
        
        % remove TRs to censor from the data & the associated rows from the design matrix
        [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx);
        
        % fit model to data
        [B,~,tB,~,~] = glm_fmri_fit(data,X,regIndx,mask.data);
        
        cd(outDir);
        
        switch irf
            
            case 'spm_hrf'
                
                for c = 1:length(stim)
                    outName = [subject,'_',stim{c},'_betaT'];
                    outNii = makeGlmNifti(mask,outName,B(:,:,:,[strmatch(stim{c},regLabels)]),tB(:,:,:,[strmatch(stim{c},regLabels)]));
                    writeFileNifti(outNii);
                end
                
            case 'per_trial'
                
                % get habenula per trial betas for correlation
                dim=size(B);
                for c = 1:length(stim)
                    seed_betas=seedBetas{c}(s,~isnan(seedBetas{c}(s,:)))';
                    seed_betas2=seedBetas2{c}(s,~isnan(seedBetas2{c}(s,:)))';
                    
                    vox_betas=B(:,:,:,[strmatch(stim{c},regLabels)]);
                    vox_betas2=reshape(vox_betas,prod(dim(1:3)),[])';
                    
                    
                    %% note: this is for using 2 seed betas; uncomment the code below to use just one seed
                    if omitOutlierTrials % then omit betas > 3 sds from beta series mean
                   
                        omitTIdx=[find(abs(zscore(seed_betas))>3);find(abs(zscore(seed_betas2))>3)];
                        
                        omitCount(s,c)= omitCount(s,c)+length(unique(omitTIdx));
                        seed_betas(unique(omitTIdx))=[];
                        seed_betas2(unique(omitTIdx))=[];
                        vox_betas2(unique(omitTIdx),:)=[];
                        
                        % if a voxel has an outlier, exclude that trial,
                        % but only for that voxel
                    end
                    
                        r=corr(seed_betas,vox_betas2);
                        r2=corr(seed_betas2,vox_betas2);
                
                    
                    
%                     if omitOutlierTrials % then omit betas > 3 sds from beta series mean
%                         %                          [omitTIdx,~]=find(abs(zscore([seed_betas,vox_betas2]))>3);
%                         
%                         % if seed beta is an outlier, exclude that trial
%                         omitTIdx=find(abs(zscore(seed_betas))>3);
%                         
%                         omitCount(s,c)= omitCount(s,c)+length(unique(omitTIdx));
%                         seed_betas(unique(omitTIdx))=[];
%                         vox_betas2(unique(omitTIdx),:)=[];
%                         
%                         % if a voxel has an outlier, exclude that trial,
%                         % but only for that voxel
%                         r=corr(seed_betas,vox_betas2);
%                         
%                         [idxR,idxC]=find(abs(zscore(vox_betas2))>3);
%                         
%                         if ~isempty(idxR)
%                             
%                             for v = 1:length(idxR)
%                                 this_seed_betas = seed_betas;
%                                 this_seed_betas(idxR(v))=[];
%                                 
%                                 this_vox_betas = vox_betas2(:,idxC(v));
%                                 this_vox_betas(idxR(v))=[];
%                                 
%                                 r(idxC) = corr(this_seed_betas,this_vox_betas);
%                             end
%                         end
%                         
%                     else
%                         
%                         r=corr(seed_betas,vox_betas2);
%                         
%                     end
                    r3=reshape(r,dim(1:3));
                    r32=reshape(r2,dim(1:3));
                    
                    outName = [subject,'_',stim{c},'_corr_',seedStr];
                    outNii = makeGlmNifti(mask,outName,r3);
                    writeFileNifti(outNii);
                    
                    z = .5.*log((1+r3)./(1-r3));
                    outNameZ = [subject,'_',stim{c},'_corrZ_',seedStr];
                    outNiiZ = makeGlmNifti(mask,outNameZ,z);
                    writeFileNifti(outNiiZ);
                    
                    outName = [subject,'_',stim{c},'_corr_',seedStr2];
                    outNii = makeGlmNifti(mask,outName,r32);
                    writeFileNifti(outNii);
                    
                    z2 = .5.*log((1+r32)./(1-r32));
                    outNameZ = [subject,'_',stim{c},'_corrZ_',seedStr2];
                    outNiiZ = makeGlmNifti(mask,outNameZ,z2);
                    writeFileNifti(outNiiZ);
                    
                    
                end
                
            case 'tent'
                
                for c = 1:length(stim)
                    outName = [subject,'_',stim{c},'_irf'];     % juice
                    outNii = makeGlmNifti(mask,outName,B(:,:,:,[strmatch(stim{c},regLabels)]));
                    writeFileNifti(outNii);
                end
                
            case 'habPPI'
                
                for c = 1:length(stim)
                    outName = [subject,'_',stim{c},'_betaT'];
                    outNii = makeGlmNifti(mask,outName,B(:,:,:,[strmatch(stim{c},regLabels)]),tB(:,:,:,[strmatch(stim{c},regLabels)]));
                    writeFileNifti(outNii);
                end
        end
        
        clear expPaths func mask
        
        fprintf(['\n\nfinished subject', subject]);
        
    end % if not sa32 & per_trial
    
end % subjects
