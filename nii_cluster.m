function [C,nii] = nii_cluster(nii, cl_thresh,cl_connectivity)
% [C,nii] = nii_cluster(nii)
%
% Determines location of clusters of non-zero voxels in
% nii.data(:,:,:,1) and returns a summary of cluster information and voxel
% indices in C.
%
% If cl_thresh is given, clusters w/fewer voxels than that number will be
% removed from the nii file and not included in the C struct w/cluster
% info.
%
% cl_connectivity is a parameter for clustering  - use 6 for faces
% touching, 18 for edges touching, or 26 for corners touching
%
% Written by Samuel McClure, 20??
% modified Feb. 2004
%
% modified by Kelly for nifti files, Oct 2012; to include a
% cluster threshold, March 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nii.data = nii.data(:,:,:,1); % cluster data from the first volume only
voxDim = nii.pixdim(1:3); % voxel dimensions

if (~exist('cl_thresh','var') || isempty(cl_thresh)) cl_thresh=0; end

if (~exist('cl_connectivity','var') || isempty(cl_connectivity)) cl_connectivity=26; end

switch cl_connectivity
    
    case 6  % voxel faces must be touching
        cl_cutoff = max(voxDim);
        
    case 18 % voxels can share faces or edges
        cl_cutoff = sqrt((max(voxDim).^2).*2); % voxels with touching corners are clustered together 

    case 26 % voxels can share corners, edges, or faces
        cl_cutoff = sqrt((max(voxDim).^3).*2); % voxels with touching corners are clustered together

end

% define struct w/ cluster info
C = struct();
C.ind = [];
C.vals = [];
C.n = [];
C.meanCoord = [];
C.maxT = [];
C.maxTCoord = [];

% cluster threshold 
mask=zeros(size(nii.data));
maskIdx = find(abs(nii.data)>0);
mask(maskIdx)=1;
mask = bwareaopen(mask,cl_thresh,cl_connectivity);
nii.data = nii.data.*mask;

% get index for remaining voxels; don't include NaNs
ind = find(abs(nii.data)>0);

if (length(ind)>1)

[ni nj nk]=ind2sub(size(nii.data),ind);
pts = mrAnatXformCoords(nii.qto_xyz,[ni nj nk]); % acpc Coords

% let matlab find the point clusters
pd = pdist(pts);
% dendrogram(linkage(pd))
cl = cluster(linkage(pd),'cutoff',cl_cutoff,'criterion','distance');

num_cl = max(cl);

for i = 1 : num_cl
    
    cl_ind = find(cl == i);
        
    these_pts=pts(cl_ind,:);
    C(i).ind = ind(cl_ind);
    C(i).vals = nii.data(ind(cl_ind));
    C(i).n = size(cl_ind, 1);
    C(i).meanCoord = round(mean(these_pts));
    
    [~,maxTind] = max(abs(C(i).vals));
    C(i).maxT = C(i).vals(maxTind);
    C(i).maxTCoord = these_pts(maxTind,:);
    
end

% if isempty(fields(C))
%     disp(['no clusters passed threshold']);
% else
%     disp(['returning struct C with info on ',num2str(length(C)),' clusters...']);
% end
%  

% else
%     disp('there are less than 2 nonzero vals in nii.data, returning cluster info as an empty struct')

end 

end % function
