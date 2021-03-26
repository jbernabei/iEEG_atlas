function [mni_coords, mni_labels, NN_flag] = nifti_values(coords_filepath,atlas_filepath)
% [mni_labels] = nifti_values(cords_filepath,atlas_filepath)
% takes in a csv file where the first 3 columns are (x,y,z) MNI
% coordinates and an atlas file in mni space and outputs the atlas labels
% for the csv coordiantes.
%
% Test:
%   [mni_labels] = nifti_values('test_input.csv','AAL116.nii')
%
% Input:
%   csv_filepath (str): filepath to coordinate file organized (Nx3)
%   atlas_filepath (str): filepath to atlas file in MNI space (.gz unzip)
%
% Output:
%   mni_labels (double): atlas values at each given coordinate
%   NN_flag (double): flags coordinates where nearest neighbor was
%       indicated due to atlas value 0
%
% Thomas Campbell Arnold
% tcarnold@seas.upenn.edu
% 3/22/2019

% read in AAL image
V=niftiinfo(atlas_filepath); % get header
atlas = niftiread(V); % get 3D matrix
T=V.Transform.T; % get transformation matrix
T=T'; % transpose transformation matrix

% get mni coordinates from csv and transfer to matlab coordinate space
%csv_coords = dlmread(coords_filepath);
csv_coords = coords_filepath; %csv_coords(:,1:3);
mni_coords = mni2cor(csv_coords,T);

NN_flag = zeros(size(mni_coords,1), 2); % variable indicating no label initial found in atlas

% get AAL label based on coordinates
for i = 1:size(mni_coords,1)
%     i
%     if mni_coords(i,1) <= 0
%         mni_coords(i,1) = 1;
%     end
    [mni_coords(i,1), mni_coords(i,2), mni_coords(i,3)];
    mni_labels(i) = atlas(mni_coords(i,1), mni_coords(i,2), mni_coords(i,3));
    radius = 0; % initial radius size
    while mni_labels(i) == 0 % get mode of cubes until answer is achieved
        NN_flag(i,1) = 1; % store value indicating label not localizaed
        radius = radius + 1; % increase radius
        [my_sphere,coords] = gen_sphere(radius); % get coordinates to sample
        x = coords.x + mni_coords(i,1); 
        y = coords.y + mni_coords(i,2);
        z = coords.z + mni_coords(i,3);
        size(atlas);
        try
            ind = sub2ind(size(atlas),x,y,z); % convert to indices
            sphere_vals = atlas(ind); % sample sphere
            [j,k,vals] = find(sphere_vals); % find nonzero and non-NaN values
            if vals % if there are nonzero values, update label
                mni_labels(i) = mode(vals(:)); % get mode
                NN_flag(i,2) = radius; % store distance to NN
            end
        catch
            mni_labels(i) = 0; % get mode
            NN_flag(i,2) = radius; % store distance to NN
        end
        
        if radius>=6
            break
        end
        
    end
end

end

function [my_sphere,coords] = gen_sphere(radius)
    dim = 2*radius + 1; % size of box containing sphere
    mid = radius + 1; % midpoint of the box ([mid mid mid] is centroid)
    my_sphere = zeros(dim,dim,dim); % build box
    for i = 1:numel(my_sphere) % loop through each location in box
        [X,Y,Z] = ind2sub([dim,dim,dim],i); % get coordinates
        D = dist([mid mid mid; X Y Z]'); % distance from centroid
        if D <= radius % if less than radius, put in box
            my_sphere(i) = 1;
        end
    end
    
    % get coordinates and offset by center voxel
    [coords.x,coords.y,coords.z] = ind2sub(size(my_sphere), find(my_sphere));
    coords.x  = coords.x - mid;
    coords.y  = coords.y - mid;
    coords.z  = coords.z - mid;
end
