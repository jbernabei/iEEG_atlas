function [zscores, pt_zscores, pt_ind,roi_ind] = test_univariate_atlas(mean_vec, std_vec, feat_vals, subject_vec, region_list, channel_roi,type)
    % mean_vec: N_roi x 1 vector of means of each ROI
    % std_vec: N_roi x 1 vector of std of each ROI
    % feat_vals: N_channel x 1 feature values of each channel we want to test
    % region_list: N_roi x 1 vector of ROI integer indices
    % channel ROI: N_channel x 1 vctor of ROI integer indices

    
    zscores = NaN(size(feat_vals));
    pt_zscores = [];
    pt_ind = [];
    roi_ind = [];
    
    for r = 1:length(region_list)
        % select region int and find which channels belong to it
        this_region = region_list(r);
        all_region_channels = find(channel_roi==this_region);
        
        % find feature vals in this region
        region_feat_vals = feat_vals(all_region_channels);
        
        % get region mean and std val
        region_mean_val = mean_vec(r);
        region_std_val = std_vec(r);
        
        zscores(all_region_channels,1) = (region_feat_vals - region_mean_val)./region_std_val;
        
    end
    if strcmp(type,'patient')
        for i = 1:max(subject_vec)
            clear region_abnormality
            this_pt_zscores = zscores(find(subject_vec==i));
            this_pt_roi = channel_roi(find(subject_vec==i));
            this_pt_unique_roi = unique(this_pt_roi);
            for r = 1:length(this_pt_unique_roi)
                this_unique_roi = this_pt_unique_roi(r);
                region_abnormality(r,1) = max(this_pt_zscores(this_pt_roi==this_unique_roi));
                pt_ind = [pt_ind;i];
                roi_ind = [roi_ind;r];
            end
            pt_zscores = [pt_zscores;region_abnormality];
        end
    end
    
    
end