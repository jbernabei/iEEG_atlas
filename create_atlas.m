function [median_conn, std_conn, num_samples, sem_conn] = create_atlas(all_conn, all_roi, all_resect, region_list, band, threshold)
% [mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
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
% 7/5/2020

% sets default threshold value if none is given
if ~exist('threshold','var'), threshold = 1; end

% get number of patients
num_patients = length(all_conn);

% initialize output array to store connection strengths
median_conn = NaN(90,90,num_patients);

for p = 1:num_patients
    
    % get electrode regions for the patient
    patient_electrode_regions = all_roi{p};
    % get resected electrodes for the patient
    res_elec_inds = all_resect{p};
    % orient vector vertically
    res_elec_size = size(res_elec_inds);
    if res_elec_size(1) == 1
        res_elec_inds = res_elec_inds.';
    end

    % calculate logical with resected indices
    try
        resect_boolean = cat(2,accumarray(res_elec_inds,1).',...
        zeros(1,length(patient_electrode_regions)-max(res_elec_inds)));
    catch any_error
        resect_boolean = zeros(1,length(patient_electrode_regions));
    end

    % get connection strengths between each pair of electrodes
    adj_matrix = all_conn{p};
    % extract band data
    band_matrix = adj_matrix(band).data;

    % double loop through regions
    for i = 1:90 % first region
        
        % get electrodes contained within first region
        if isempty(resect_boolean)
            p
        end
        first_reg_elec = (patient_electrode_regions == region_list(i) & ~resect_boolean);
        
        % extract rows corresponding to electrodes in the first region
        first_reg_strengths = band_matrix(first_reg_elec,:);
        
%         % calculate connections within the region
%         patient_strengths = first_reg_strengths(:,first_reg_elec);
%         same_elec_mask = eye(size(patient_strengths));
%         patient_strengths(same_elec_mask==1) = NaN;
%         patient_strength = nanmedian(patient_strengths);
%         %nanmedian(patient_strengths(triu(true(size(patient_strengths)),1)));
%         %nanmean(patient_strengths(triu(true(size(patient_strengths)),1)));
%         
%         % add to output array if the region contains any electrodes
%         if ~isnan(patient_strength)
%             median_conn(i,i,p) = patient_strength;
%         else
%             % skip calculations for this region pair if the first region
%             % does not contain any electrodes
%             continue
%         end

        for j = 1:90 % second region (used to be (i+1):90)
            
            % get electrodes contained within second region
            second_reg_elec = (patient_electrode_regions == region_list(j) & ~resect_boolean);

            % extract connection strengths between the two regions
            patient_strengths = first_reg_strengths(:,second_reg_elec);
            
            % average connection strengths between the two regions
            patient_strength = nanmedian(patient_strengths(:));
            
            % add to output array
            median_conn(i,j,p) = patient_strength;
            
            if i==j
               median_conn(i,j,p) = NaN;
            end
        end
    end
end

% take the standard deviation element-wise to get the std matrix
std_conn = std(median_conn,0,3,'omitnan');

% get the number of samples available for each edge
num_samples = sum(~isnan(median_conn),3);

% get the sem
sem_conn = std_conn./sqrt(num_samples);
sem_conn(isinf(sem_conn)) = NaN;

% divide out the number of patients element-wise to get the median matrix
median_conn = median(median_conn,3,'omitnan');

% remove edge values and standard deviations for edges with sample sizes
% less than the threshold
median_conn(num_samples < threshold) = NaN;
std_conn(num_samples < threshold) = NaN;
sem_conn(num_samples < threshold) = NaN;

% symmetrize all output matrices
% symmetrize = @(x) triu(x) + tril(x.',-1);
% mean_conn = symmetrize(mean_conn);
% std_conn = symmetrize(std_conn);
% sem_conn = symmetrize(sem_conn);
% num_samples = symmetrize(num_samples);

end
