function [pt_adj] = compute_FC(iEEG_signals, sample_rate, window_len, type, band)

% iEEG signals must be a N-by-P matrix where N is number of samples and P
% is number of channels

% window_len is in seconds

% % generate common average reference and re-reference
% common_avg = mean(iEEG_signals);
% car_iEEG = iEEG_signals - common_avg;
% 
% % 50 Hz IIR notch
% wo1 = 50/(sample_rate/2);  
% bw1 = wo1/35;
% [b1,a1] = iirnotch(wo1,bw1);
% 
% % 60 Hz IIR notch
% wo2 = 60/(sample_rate/2);  
% bw2 = wo2/35;
% [b2,a2] = iirnotch(wo2,bw2);
% 
% % 1 - 70 Hz
% % fourth order filter
% [b3,a3] = butter(4,[1 70]/(sample_rate/2));
% 
% % apply filters
% signal1 = filter(b1, a1, car_iEEG);
% signal2 = filter(b2, a2, signal1);
% signal3 = filter(b3, a3, signal2);

% Pearson correlation
if strcmp(type,'corr')
num_windows = floor(size(iEEG_signals,1)./(sample_rate.*window_len));
pt_adj = zeros(size(iEEG_signals,2),size(iEEG_signals,2),num_windows);
window_size = window_len.*sample_rate;

    % filter in order to extract desired frequency band
    if strcmp(band,'delta')
        % bandpass filter 0.5-4 Hz
        [b1,a1] = butter(4,[0.5 4]/(200/2));
        iEEG_signals = filter(b1, a1, iEEG_signals);

    elseif strcmp(band,'theta')
        % bandpass filter 4-8 Hz
        [b2,a2] = butter(4,[4 8]/(200/2));
        iEEG_signals = filter(b2, a2, iEEG_signals);

    elseif strcmp(band,'alpha')
        % bandpass filter 8-13 Hz
        [b3,a3] = butter(4,[8 13]/(200/2));
        iEEG_signals = filter(b3, a3, iEEG_signals);

    elseif strcmp(band,'beta')
        % bandpass filter 13-30 Hz
        [b4,a4] = butter(4,[13 30]/(200/2));
        iEEG_signals = filter(b4, a4, iEEG_signals);

    elseif strcmp(band,'gamma')
        % bandpass filter 30-80 Hz
        [b5,a5] = butter(4,[30 80]/(200/2));
        iEEG_signals = filter(b5, a5, iEEG_signals);

    end



for i = 1:num_windows
    start_ind = (i-1).*window_size+1;
    stop_ind = i.*window_size;
    pt_adj(:,:,i) = corr(iEEG_signals([start_ind:stop_ind],:));
    
end

pt_adj = median(pt_adj,3);

elseif strcmp(type,'coh')
    
    num_chs = size(iEEG_signals,2);
    pt_adj = zeros(num_chs,num_chs);
    for i = 1:num_chs
        for j = 1:num_chs
            [cxy,f] = mscohere(iEEG_signals(:,i),iEEG_signals(:,j),200,0,[1:80],200);
            if strcmp(band,'delta')
                pt_adj(i,j) = mean(cxy(1:4));
            elseif strcmp(band,'theta')
                pt_adj(i,j) = mean(cxy(5:8));
            elseif strcmp(band,'alpha')
                pt_adj(i,j) = mean(cxy(9:13));
            elseif strcmp(band,'alpha/theta')
                pt_adj(i,j) = mean(cxy(5:13));
            elseif strcmp(band,'beta')
                pt_adj(i,j) = mean(cxy(14:30));
            elseif strcmp(band,'gamma')
                pt_adj(i,j) = mean(cxy(31:end));
            end
        end
    end
    
    
end

end