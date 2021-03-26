function centroid = atlas_centroids(atlas_filepath)
% this function reads in an atlas nifti file, get all unique labels, and
% output the centroid for each label

% read in atlas
V=niftiinfo(atlas_filepath); % get header
atlas = niftiread(V); % get 3D matrix

% get distance matrix for all electrodes and atlas ROIs
labels = unique(atlas(find(atlas)));
for i = 1:length(labels)
    idx = find(atlas==labels(i));
    [x,y,z] = ind2sub(size(atlas),idx);
    centroid(i,:) = [mean(x) mean(y) mean(z)];
end