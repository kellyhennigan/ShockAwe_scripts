function [imgRgb, maskRgb,h,acpcSlices,resNiiVol] = plotOverlayImage(nii,t1,cmap,c_range,plane,acpcSlices)
%
% overlay rgb maps onto gray scaled anatomical background
%
% % INPUTS:
%
%     nii - nii file with vol = nii.data(:,:,:,1) used for overlay
%     t1 - anatomical nii file to be grayscaled background image
%     cmap - colormap for overlay; default is autumn
%     c_range - clip range for overlay;
%         abs(vol) < c_range(1) will be set to zero
%         abs(vol) > c_range(2) will be colored the same as c_range(2)
%     plane - 1 for sag, 2 for coronal, 3 for axial images; default is 3
%     acpcSlices - slices to plot in acpc coords; default is the mode
%         activation slice
%
% % OUTPUTS:
%     imgRgb - cell array of M x N x 3 arrays with r,g,b values along
%              the 3rd dimension
%     maskRgb - " " binary arrays with 1s for voxels included in imgRgb
%              overlay and 0s otherwise
%     h - cell array of figure handles
%     acpcSlices - acpcSlices that correspond to images & figures
%
%
% kjh March 2013
% based heavily on Bob's mrAnatOvelayMontage() function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol = nii.data(:,:,:,1); % overlay volume
xform = nii.qto_xyz;     % it's img to acpc xform
t1Vol = double(t1.data); % structural volume for background
t1Xform = t1.qto_xyz;    % it's img to acpc xform

% if an argument isn't given then assign default values
if(~exist('cmap','var') || isempty(cmap)) cmap = autumn; end

if(~exist('c_range','var') || isempty(c_range))
    idx = find(vol);
    c_range(1) = min(vol(idx));
    c_range(2) = max(vol(idx));
end

if(~exist('plane','var') || isempty(plane))
    plane = 3; % default to axial slices
end

if (isempty(acpcSlices) || ~exist('acpcSlices','var'))
    [vi, vj, vk] = ind2sub(size(vol),find(vol));
    acpcCoords = mrAnatXformCoords(xform,[vi vj vk]);
    acpcSlices = round(unique(acpcCoords(:,plane)))
end

if size(acpcSlices,1)>1
    acpcSlices=acpcSlices'; % make acpcSlices a row vector
end

lm=linspace(c_range(1),c_range(2),length(cmap))'; % values to color mapping

%% threshold vol for overlay

% vol(abs(vol)<c_range(1))=0;
% vol(abs(vol)>c_range(2))=c_range(2);


%% resample vol to be same size as t1Vol w/same dimensions

bb = t1Xform*[1,1,1,1;[size(t1Vol),1]]'; % bounding box for resampled vol
bb = bb(1:3,:)';
interp = [0 0 0 0 0 0]; % [1 1 1 0 0 0] for trilinear interpolation
% interp = [0 0 0 7 7 7]; % [1 1 1 0 0 0] for trilinear interpolation

vol = mrAnatResliceSpm(vol, inv(xform), bb, t1.pixdim, interp, false);

resNiiVol = vol; % to return as output
%% reorient t1Vol and vol depending on slice orientation

switch plane
    case 1 % sagittal; eyes looking left
        t1Vol = permute(t1Vol,[3 2 1]); t1Vol = flipdim(t1Vol,1);
        vol = permute(vol,[3 2 1]); vol = flipdim(vol,1);
        imgCoords = mrAnatXformCoords(t1.qto_ijk,[acpcSlices;ones(1,length(acpcSlices));ones(1,length(acpcSlices))]);
        slStr = 'x = ';
        
    case 2 % coronal; top of the head is up
        t1Vol = permute(t1Vol,[3 1 2]); t1Vol = flipdim(t1Vol,1); 
        vol = permute(vol,[3 1 2]); vol = flipdim(vol,1); 
        imgCoords = mrAnatXformCoords(t1.qto_ijk,[ones(1,length(acpcSlices));acpcSlices;ones(1,length(acpcSlices))]);
        slStr = 'y = ';
        
    case 3 % axial; eyes looking up
        t1Vol = permute(t1Vol,[2 1 3]);  t1Vol = flipdim(t1Vol,1);
        vol = permute(vol,[2 1 3]); vol = flipdim(vol,1);
        imgCoords = mrAnatXformCoords(t1.qto_ijk,[ones(1,length(acpcSlices));ones(1,length(acpcSlices));acpcSlices]);
        slStr = 'z = ';
end
imgSlices = imgCoords(:,plane);

%% make rgb overlay images

for i = 1:length(imgSlices)
    
    % anatomical background image (grayscaled)
    t1Sl=t1Vol(:,:,imgSlices(i));
    t1Sl = t1Sl./max(max(t1Sl));  % scale vals from 0 to 1
    t1Rgb = repmat(t1Sl,[1,1,3]); % repping this means rgb values will be gray scaled
    
    % binary mask of overlay rgb image
    sl=vol(:,:,imgSlices(i));
    olMask = zeros(size(sl)); olMask(abs(sl)>0) = 1; % make a mask of rgb image
    olMask = repmat(olMask,[1,1,3]); % mask of rgb image
    idx_ol = find(olMask); % index for all vals included in rgb overlay
    
    % overlay rgb image w/rgb values from the colormap
    olRgb = olMask;
    [~,ci]=pdist2(lm,sl(abs(sl)>0),'euclidean','smallest',1); % ci is the colormap index for that voxel
    olRgb(idx_ol) = cmap(ci,:);   % fills in rgb values for each overlay voxel
    
    % rgb mask & rgb image w/anatomical background & rgb overlay
    maskRgb{i} = olMask;
    imgRgb{i} = t1Rgb;
    imgRgb{i}(idx_ol) = olRgb(idx_ol);
    
    scSize = get(0,'ScreenSize');
    h{i} = figure;
    pos = get(gcf,'Position');
    set(gcf,'Position',[scSize(3)-pos(3), scSize(4)-pos(4), pos(3), pos(4)])
    image(imgRgb{i})
    axis equal; axis off;
    set(gca,'Position',[0,0,1,1]);   
%     text(size(sl,2)-20,size(sl,1)-20,[slStr,num2str(acpcSlices(i))],'color',[1 1 1])
    
end % slices

%% plot a colorbar as a separate figure

figure
pos = get(gcf,'Position');
set(gcf,'Position',[scSize(3)-pos(3), 60, pos(3), 100])
image(1:length(lm))
colormap(cmap)
xticks = get(gca,'XTick');
xticks = interp1([1 length(cmap)],c_range,xticks);
xticks=num2str(xticks',2);
set(gca,'XTickLabel',xticks)
set(gca,'YTickLabel',[])
set(gca,'box','off');
set(gcf,'Color','w','InvertHardCopy','off','PaperPositionMode','auto');
title('Colorbar')

end


