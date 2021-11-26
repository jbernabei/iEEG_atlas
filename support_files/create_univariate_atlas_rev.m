function [mean_vec,median_vec,std_vec,feat_vals] = create_univariate_atlas_rev(raw_data,region_list,channel_roi,method,freq_band)
    
    % raw_data is the underlying data, which should be 12000 x N_channels 
    % in this project, given that we are using 60 seconds of data
    % downsampled to 200 Hz
    
    % region_list is a vector of all possible atlas ROI numerical labels
    
    % channel_roi is an N_channel x 1 vector of each channel's atlas ROI


    % method is a string with the following possibilities:
    % 'entropy' for shannon wavelet entropy
    % 'bandpower' for spectral density

    
    mean_vec = NaN(length(region_list),1);
    median_vec = NaN(length(region_list),1);
    std_vec = NaN(length(region_list),1);

     if strcmp(method,'entropy')
%         
%         % filter in order to extract desired frequency band
%         if strcmp(freq_band,'delta')
%             % bandpass filter 0.5-4 Hz
%             [b1,a1] = butter(4,[0.5 4]/(200/2));
%             filt_data = filter(b1, a1, raw_data);
%             
%         elseif strcmp(freq_band,'theta')
%             % bandpass filter 4-8 Hz
%             [b2,a2] = butter(4,[4 8]/(200/2));
%             filt_data = filter(b2, a2, raw_data);
%             
%         elseif strcmp(freq_band,'alpha')
%             % bandpass filter 8-13 Hz
%             [b3,a3] = butter(4,[8 13]/(200/2));
%             filt_data = filter(b3, a3, raw_data);
%             
%         elseif strcmp(freq_band,'beta')
%             % bandpass filter 13-30 Hz
%             [b4,a4] = butter(4,[13 30]/(200/2));
%             filt_data = filter(b4, a4, raw_data);
%             
%         elseif strcmp(freq_band,'gamma')
%             % bandpass filter 30-80 Hz
%             [b5,a5] = butter(4,[30 80]/(200/2));
%             filt_data = filter(b5, a5, raw_data);
%             
%         end
        
        % loop through time and channels to create entropy matrix
        for i = 1:60
            for j = 1:size(raw_data,2)
                start_inds = (i-1)*200+1;
                end_inds = 200*i;
                all_wentropy(i,j) = wentropy(raw_data((start_inds:end_inds),j),'shannon');
            end
        end
        
        % take median of each channel's entropy across time
        feat_vals = median(real(log(-1*all_wentropy)))';
        
    elseif strcmp(method,'bandpower')
        
        % we compute a normalized broadband spectral decomposition first
        [pxx,f] = pwelch(raw_data,400,200,[0.5,2:80],200);
        pxx_norm = real(normalize(pxx,'norm',1));
        
        % extract spectral density in desired frequency band
        if strcmp(freq_band,'delta')
            % power in 1-4 Hz band;
            feat_vals = median(pxx_norm(1:4,:))';
        elseif strcmp(freq_band,'theta')
            feat_vals = median(pxx_norm(4:8,:))';
        elseif strcmp(freq_band,'alpha')
            feat_vals = median(pxx_norm(8:13,:))';
        elseif strcmp(freq_band,'beta')
            feat_vals = median(pxx_norm(13:30,:))';
        elseif strcmp(freq_band,'gamma')
            feat_vals = median(pxx_norm(30:80,:))';
        end
         
    end
    
    % map to atlas
    for r = 1:length(region_list)
        % select region int and find which channels belong to it
        this_region = region_list(r);
        all_region_channels = find(channel_roi==this_region);
        
        % find feature vals in this region
        region_feat_vals = feat_vals(all_region_channels);
        
        % assign mean, median, std
        mean_vec(r) = nanmean(region_feat_vals);
        median_vec(r) = nanmedian(region_feat_vals);
        std_vec(r) = nanstd(region_feat_vals);
    end
end