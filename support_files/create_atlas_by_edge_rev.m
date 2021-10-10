function [mean_conn, std_conn, num_samples, sem_conn, raw_atlas] = create_atlas_by_edge_rev(all_conn, all_roi, all_resect, region_list, threshold)
% [mean_conn, std_conn] = create_atlas_by_edge(all_conn, all_roi, all_resect, region_list)
% takes in an array of conectivity structs, an array of 3D mni coordinate
% arrays, an array of resected electrode vectors, and a vector containing
% all mni labels, and outputs two matrices representing the average and
% standard deviation of connections between all regions.
%
% Input:
%   all_conn (cell): cell array containing patient connectivity structs in
%   order
%   all_roi (cell): cell array containing regions of interest corresponding
%   to each electrode for each patient in order
%   all_resect (cell): cell array containing patient resected electrode
%   arrays in order
%   region_list (double): array containing all region labels
%   band (int): frequency band to be used
%   threshold (int): minimum sample size required for edges to be
%   incorporated into the atlas. Default value = 1
%
% Output:
%   mean_conn (double): (i,j) matrix of mean connectivity strengths between 
%   region_list(i) and region_list(j)
%   std_conn (double): (i,j) matrix of standard deviations of connectivity 
%   strengths between region_list(i) and region_list(j)
%   num_samples (double): a matrix where entry (i,j) is the number of
%   samples used to calculate the edge weight between region_list(i) and region_list(j)
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 7/30/2020

% helper function
get_data = @(x) x(triu(true(size(x)),1));

% sets default threshold value if none is given
if ~exist('threshold','var'), threshold = 1; end

% get number of patients
num_patients = length(all_conn);

% get number of regions
num_regions = length(region_list);

% initialize output arrays
mean_conn = NaN(num_regions);
std_conn = NaN(num_regions);
sem_conn = NaN(num_regions);
num_samples = NaN(num_regions);
raw_atlas = NaN(num_regions,num_regions,3576); % change this number

% load in relevant matrices and resected electrodes for each patient
band_matrices = cell(1,num_patients);
resect_booleans = cell(1,num_patients);

for p = 1:num_patients
    % get resected electrodes for the patient
    res_elec_inds = all_resect{p};
    % orient vector vertically
    res_elec_size = size(res_elec_inds);
    if res_elec_size(1) == 1
        res_elec_inds = res_elec_inds.';
    end

    % calculate logical with resected indices
    try
        resect_booleans{p} = cat(2,accumarray(res_elec_inds,1).',...
        zeros(1,length(all_roi{p})-max(res_elec_inds)));
    catch any_error
        resect_booleans{p} = zeros(1,length(all_roi{p}));
    end

    % extract band data
    band_matrix{p} = all_conn{p};
end

% double loop through regions
for i = 1:num_regions % first region

    for j = 1:num_regions % second region
        
        % set up array to hold all edge values
        edge_values_for_region_pair = cell(num_patients,1);
        
        % counter for number of data points for this region pair
        sample_counter = 0;
        
        for p = 1:num_patients
            
            % get electrodes contained within first region
            first_reg_elec = (all_roi{p} == region_list(i) & ~resect_booleans{p});
        
            % extract rows corresponding to electrodes in the first region
            current_matrix = band_matrix{p};
            first_reg_strengths = current_matrix(first_reg_elec,:);
            
            % get electrodes contained within second region
            second_reg_elec = (all_roi{p} == region_list(j) & ~resect_booleans{p});

            % extract connection strengths between the two regions
            patient_strengths = first_reg_strengths(:,second_reg_elec);
            
            if i == j % loop case
                % get only the upper triangular data since it's symmetric
                
                % I edited this to quantify self connections (within ROI) - John
                same_roi_data = current_matrix(first_reg_elec, second_reg_elec);
                patient_strengths = same_roi_data(:); %get_data(patient_strengths);
                
            else % assuming regions are spatially disjoint
                patient_strengths = patient_strengths(:);
            end
            
            % add number of edges to the counter
            sample_counter = sample_counter + length(patient_strengths);
            
            % store electrode edges for this region pair in this patient
            edge_values_for_region_pair{p} = patient_strengths;
            
        end
        
        % record mean, std, and sem if the number of edges exceeds the threshold
        if sample_counter >= threshold
            % take the average and standard deviation of all edge values for
            % the region pair and place both in the corresponding output arrays
            
        
            mean_conn(i,j) = mean(cell2mat(edge_values_for_region_pair));
            std_conn(i,j) = std(cell2mat(edge_values_for_region_pair));
            sem_conn(i,j) = std_conn(i,j)/sqrt(sample_counter);
            num_edge_pairs = length(cell2mat(edge_values_for_region_pair));
            num_current_edges = sum(isnan(raw_atlas(i,j,:)));
            raw_atlas(i,j,((num_current_edges+1):(num_current_edges+num_edge_pairs))) = cell2mat(edge_values_for_region_pair);
            
        end
        
        % record the number of edges for this region pair
        num_samples(i,j) = sample_counter;
        
        if i==j
           mean_conn(i,j) = NaN;
           std_conn(i,j) = NaN;
           sem_conn(i,j) = NaN;
           raw_atlas(i,j,:) = NaN;
        end
        
    end
    
end

raw_atlas(raw_atlas==0) = NaN;

% symmetrize all output matrices
symmetrize = @(x) triu(x) + tril(x.',-1);
mean_conn = symmetrize(mean_conn);
std_conn = symmetrize(std_conn);
sem_conn = symmetrize(sem_conn);
num_samples = symmetrize(num_samples);

end
