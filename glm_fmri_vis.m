% script for checking out voxelwise fmri results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all

p = [.001 .0025 .005 .01 .02 .05]
t_1tail=tinv(1-p,17) % one tailed
t_2tail=tinv(1-(p/2),17) % one tailed

%% specify inputs & parameters

t1File = '/home/kelly/ShockAwe/data/TT_N27.nii.gz';

niiDir = '/home/kelly/ShockAwe/data/ttests_spm_hrf'
niiFile = 'shock-neutral.nii.gz';

volIdx = 2; % which volume of the nifti file?

maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/mb_group_mask_tlrc.nii.gz';

smooth = []; % if smooth, this is a 1 x 3 vector for a smoothing kernel size, eg, [3 3 3]

v_thresh = 3.24; % voxel value threshold; abs()<v_thresh will be set to zero
cl_thresh = 6;       % cluster size threshold

lm = linspace(1,0,64)';
cmap = [zeros(64,1),lm,ones(64,1)];
cmap = [cmap; autumn(64)];  % colormap
% cmap = autumn(64);

c_range = [-7 7]; % range for cmap; vals > v_range(2) will be plotted as v_range(2)

plane=1;

acpcSlices=[];

saveFigs=0;
% outName = 'sn_peakSTN_coronal_y-18.eps'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get files

t1 = readFileNifti(t1File);

nii=readFileNifti(fullfile(niiDir,niiFile));
nii.data = nii.data(:,:,:,volIdx);
nii.data = double(nii.data);

if (maskFile)
    mask = readFileNifti(maskFile);
    mask.data=double(mask.data);
    nii.data=nii.data.*mask.data;
end

%% threshold nii.data & get cluster info

if (smooth)
    nii.data = smooth3(nii.data,'gaussian',smooth); % smooth data
end

nii.data(abs(nii.data)<v_thresh)=0; % set vals < v_thresh to 0

if (cl_thresh)
    [C,nii]=nii_cluster(nii,cl_thresh,26); % cluster threshold
    length(C)
    max([C(:).n])
end


%% which acpcSlices?

ind = find(nii.data);
[i j k]=ind2sub(size(nii.data),ind);
coords = mrAnatXformCoords(nii.qto_xyz,[i j k]);
x=coords(:,1); y = coords(:,2); z = coords(:,3);


%% plot overlay on structural image
%
if (saveFigs)
    cd ~/Dropbox/SA_writeup/figures/
    h=figure(1)
    saveas(h,outName,'epsc')
end

[imgRgb, maskRgb,h,sl] = plotOverlayImage(nii,t1,cmap,c_range,plane,acpcSlices);
