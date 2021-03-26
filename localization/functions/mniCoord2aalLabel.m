function [electrode_regions,distance_matrix,atlas_info] = mniCoord2Label(coords_filepath,atlas_name)
% labels = mniCoord2Label(coords_csv)
% takes in a csv files where the first 3 columns are (x,y,z) MNI
% coordinates, outputs the AAL atlas names for each coordinate. Uses AAL
% 116 as the default.
%
% Test:
%   [regions,D,atlas] = mniCoord2Label('test_input.csv')
%
% Input:
%   csv_filepath (str): filepath to coordinate files organized (Nx3)
%
% Output:
%   electrode_regions (cell): each cell row corresponds to the matching csv
%       file row
%   distance_matrix (float): NxN matrix of distance between coordinates. 
%       Coordinates from atlas listed first, coordinates from csv second.
%   atlas_info (cell): 3 columns [region abr., region name, intensity in 
%       atlas
%
% Thomas Campbell Arnold
% tcarnold@seas.upenn.edu
% 3/13/2019
%  
% Updates:
%   3/22/2019 - added nearest neighbor output for zero labels
% 

%coords_filepath = 'electrode_coordinates_mni.csv'
%atlas_filepath = 'AAL116.nii'

% read in default AAL116 if not specified
if ~exist('atlas_filepath')
    atlas_name = 'AAL116';
end

% read in AAL116 image
atlas_filepath = [atlas_name,'.nii'];
[mni_coords, mni_labels, NN_flag] = nifti_values(coords_filepath,atlas_filepath);

% get AAL atlas info
fileID = fopen([atlas_name,'.txt']);
atlas_info = textscan(fileID,'%s %s %s');

% Get label names (ex: HIP 4101 & 4102)
for i = 1:length(mni_labels)
[C,ia,ib(i)] = intersect(mni_labels(i),str2num(cell2mat(atlas_info{1,3})));
end
[mni_coords, mni_labels, NN_flag] = nifti_values(coords_filepath,atlas_filepath); % get atlas labels
electrode_regions = [atlas_info{1,2}(ib) num2cell(mni_labels') num2cell(NN_flag)]; % output AAL116 roi names and atlas labels

% get distance matrix for all electrodes and atlas ROIs
centroid = atlas_centroids(atlas_filepath);
distance_matrix = dist([centroid; mni_coords]');

return;