%% Set up workspace
clear all

% set up path - change iEEG_atlas_path to where you download the repository
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas_dev';

% other paths may stay the same
metadata = readtable('data/atlas_metadata_final.xlsx');
custom_atlas = readtable('data/custom_atlas.xlsx');

% anatomical atlas information
region_list = zeros(1,90); % for the 90 AAL regions we will be using
region_names = cell(1,90);
fi = fopen("localization/AAL116_WM.txt");
for j = 1:90
    label = split(fgetl(fi));
    region_list(j) = str2double(label{3});
    region_names{j} = label{2};
end

% extract metadata
age_onset = metadata{:,11};
age_surgery = metadata{:,12};
HUP_outcome = floor(metadata{:,3});

% extract outcome
for s = 1:length(age_surgery)
    engel_scores = metadata{s,3:5};
    engel_scores = rmmissing(engel_scores);
    early_outcome(s) = floor(engel_scores(1)); % early outcome is 6 months
    late_outcome(s) = floor(engel_scores(end)); % late outcome is latest time point up to 24 months
end

therapy_field = metadata{:,6};
implant_field = metadata{:,7};
target_field = metadata{:,8};
laterality_field = metadata{:,9};
lesion_field = metadata{:,10};
gender_field = metadata{:,13};

lesion_field = strcmp(lesion_field,'Lesional'); % 1 for lesional
therapy_field = strcmp(therapy_field,'Ablation'); % 1 ffor ablation

load('data/MNI_atlas_new.mat')
load('data/HUP_atlas_May25.mat')

all_pts = [Patient; (patient_no+106)]; % combine patient numbers for MNI & HUP

HUP_outcome_all = zeros(length(soz_ch),1);
for i = 1:60
    HUP_outcome_all(find(patient_no==i)) = late_outcome(i); % we use outcome as late outcome
end

clear NodesLeft
clear NodesLeftInflated
clear NodesRightInflated
clear NodesRight
clear NodesRegionLeft
clear NodesRegionRight
clear Data_N2
clear Data_N3
clear Data_R

% some colors for plotting, etc
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color3 = [255, 248, 209]./255; 
color4 = [103 55 155]/255;
color6 = [78 172 91]/255;
color7 = [0.9290, 0.6940, 0.1250];
color8 = [0.8500, 0.3250, 0.0980];

% assign wake clips from MNI and HUP into a final data matrix
all_wake_data = [Data_W,wake_clip];
    
%% perform localization based on AAL116 atlas (which has an internal WM mask also)
all_coords = [ChannelPosition;mni_coords];

try [coords, all_roi, NN_ind] = nifti_values(all_coords,'localization/AAL116_WM.nii');
catch ME
    for i = 1:size(all_coords,1)
        if mod(i,100)
            fprintf('%d\n',i)
        end
        try [mni_coords1, mni_roi, NN_ind] = nifti_values(all_coords(i,:),'localization/AAL116_WM.nii');
        catch ME
            mni_roi = NaN; 
        end
        all_roi(i,1) = mni_roi;
    end
end

clear NN_ind
clear coords

% assign AAL regions into new, aggregate regions of custom atlas. (116 -> 40 regions)
region_matrix = custom_atlas{:,2:4};
for c = 1:length(all_roi)
    old_ind = find(region_list==all_roi(c));
    try [row,col,V] = find(region_matrix==old_ind);
        try new_roi(c,1) = row;
        catch anyerror
        new_roi(c,1) = 0;
        end
    catch anyerror
        new_roi(c,1) = 0;
    end

end
%% define normal MNI/HUP and abnormal channels (SOZ/high-spike/RZ)
% set threshold on spike counts
spike_thresh = 24; % this is empirical, 1 spike/hour

% find indices for which spikes are greater than threshold
spike_ind = [spike_24h>spike_thresh];

% define all abnormal channels
abnormal_ch = find([resected_ch+spike_ind+soz_ch]>0)+1772;

% define all seizure onset indices
soz_ch_inds = find(soz_ch)+1772;

% define normal HUP channels
normal_HUP_ch = find([resected_ch+spike_ind+soz_ch]==0)+1772;

% define normal MNI channels
normal_MNI_ch = [1:1772]';

% define all normal channels
all_normal_ch = [normal_MNI_ch;normal_HUP_ch];

% define exclusive irritative zone
EIZ_ch = find(spike_ind - spike_ind.*soz_ch)+1772;

% split EIZ into good vs poor
EIZ_good_ch = find((spike_ind - spike_ind.*soz_ch).*(HUP_outcome_all==1))+1772;
EIZ_poor_ch = find((spike_ind - spike_ind.*soz_ch).*(HUP_outcome_all>1))+1772;

% split SOZ into good vs poor
soz_good_ch = find(soz_ch.*(HUP_outcome_all==1))+1772;
soz_poor_ch = find(soz_ch.*(HUP_outcome_all>1))+1772;

% define non-EIZ non-SOZ poor outcome channels
poor_out_ch = find([resected_ch+spike_ind+soz_ch+(HUP_outcome_all==1)]==0)+1772;

%% count HUP and MNI in original atlas and then new atlas
mni_roi_old = all_roi(1:1772);
hup_roi_old = all_roi(1773:end);

for i = 1:45
    roi_count(i,1) = sum(mni_roi_old==region_list(2*i-1))+sum(mni_roi_old==region_list(2*i));
    roi_count(i,2) = sum(all_roi(normal_HUP_ch)==region_list(2*i-1))+sum(all_roi(normal_HUP_ch)==region_list(2*i));
    this_roi_name = split(string(region_names{2*i}),'_R'); % extract roi name without laterality
    split_roi_name_old{i,1} = this_roi_name(1);
end

for i = 1:20
    roi_count_new_normal(i,1) = length(intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]));
    roi_count_new_normal(i,2) = length(intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]));
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R'); % extract roi name without laterality
    split_roi_name{i,1} = this_roi_name(1);
end

new_counts = [roi_count_new_normal(:,1),(roi_count_new_normal(:,2))];

% Supplemental figure
figure(1);clf;
subplot(1,2,1)
barh(roi_count)
set(gca,'TickLabelInterpreter', 'none');
xlabel('Number of nodes per region')
title('AAL regions')
yticks([1:45])
yticklabels(split_roi_name_old)
subplot(1,2,2)
barh(new_counts)
set(gca,'TickLabelInterpreter', 'none');
title('Composite regions')
yticks([1:21])
yticklabels(split_roi_name)
xlabel('Number of nodes per region')

%% do renderings of original and new atlas
random_colors = rand([40,3]);

color_matrix = ones(117,3);

for i = 1:40
    this_region = custom_atlas{i,2:4};
    for r = 1:length(rmmissing(this_region))
        color_matrix(this_region(r),:) = random_colors(i,:);
    end
end

% writes a file that can be used in BrainNet viewer to visualize 
dlmwrite('atlas_render_color.txt',color_matrix(1:90,:),'delimiter',' ','precision',5)

figure(1);clf;
imagesc(reshape(random_colors,40,1,3))
set(gca,'TickLabelInterpreter', 'none');
yticks([1:40])
yticklabels(custom_atlas{:,1})

%% calculate power spectrum
all_wake_data = [Data_W,wake_clip];

[pxx,f] = pwelch(all_wake_data,400,200,[0.5:0.5:80],200); % original was 200/100
pxx_norm = pxx./sum(pxx);

%% Analyze power spectral density between HUP/MNI, Normal/irritative zone/SOZ

freq_inds = [1,  8; % 0.5 - 4 Hz
             9,  16; % 4 - 8 Hz
             17, 26; % 8 - 13 Hz
             27, 60; % 13 - 30 Hz 
             61, 160]; % 30 - 80 Hz


figure(1);clf % supplemental
figure(2);clf % supplemental
for i = 1:20
    % extract indices that correspond to a given ROI for HUP and MNI data
    % (combining data from L and R hemispheres)
    roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_EIZ = intersect([EIZ_good_ch;EIZ_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_soz = intersect([soz_good_ch;soz_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    
    % median power spectral density across channels
    median_mni = median(pxx_norm(:,roi_mni),2);
    median_hup = median(pxx_norm(:,roi_hup),2);
    median_normal = median(pxx_norm(:,[roi_mni;roi_hup]),2);
    median_EIZ = median(pxx_norm(:,roi_EIZ),2);
    median_soz = median(pxx_norm(:,roi_soz),2);
    
    for j = 1:5
        HUP_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_hup))';
        MNI_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_mni))';
        normal_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],[roi_mni;roi_hup]))';
        EIZ_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_EIZ))';
        soz_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_soz))';
        
        [p,h] = ranksum(HUP_feat,MNI_feat); % HUP vs MNI
        pval_matrix_norm(i,j) = p;
        
        [p,h] = ranksum(normal_feat,EIZ_feat); % normal vs EIZ
        pval_matrix_EIZ(i,j) = p;
        
        [p,h] = ranksum(normal_feat,soz_feat); % normal vs SOZ
        pval_matrix_soz(i,j) = p;
        
        [p,h] = ranksum(EIZ_feat,soz_feat); % EIZ vs SOZ
        pval_matrix_EIZ_vs_soz(i,j) = p;
        
    end

    % find upper and lower bounds (here 25th and 75th percentiles)
    prct_25_75_mni = prctile(pxx_norm(:,roi_mni)',[25,75]);
    prct_25_75_hup = prctile(pxx_norm(:,roi_hup)',[25,75]);
    prct_25_75_normal = prctile(pxx_norm(:,[roi_mni;roi_hup])',[25,75]);
    prct_25_75_EIZ = prctile(pxx_norm(:,roi_EIZ)',[25,75]);
    prct_25_75_soz = prctile(pxx_norm(:,roi_soz)',[25,75]);

    x_axis = [0.5:0.5:80];
    
    figure(1)
    % 3 x 7 plot for 21 regions (aggregating across hemispheres)
    subplot(4,5,i)
    hold on
    plot([4,4],[0,0.4],'k-')
    plot([8,8],[0,0.4],'k-')
    plot([13,13],[0,0.4],'k-')
    plot([30,30],[0,0.4],'k-')
    plot(x_axis,median_mni)
    patch([x_axis fliplr(x_axis)], [prct_25_75_mni(1,:) fliplr(prct_25_75_mni(2,:))], color1, 'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x_axis,median_hup)
    patch([x_axis fliplr(x_axis)], [prct_25_75_hup(1,:) fliplr(prct_25_75_hup(2,:))], color6, 'FaceAlpha',0.5, 'EdgeColor','none')
    
    % add greek letters for frequency bands
    txt_d = '\delta';
    txt_t = '\theta';
    txt_a = '\alpha';
    txt_b = '\beta';
    txt_g = '\gamma';
    txt_sig = '*';
    text(2,0.15,txt_d,'FontSize', 12);
    text(5.5,0.15,txt_t,'FontSize', 12);
    text(9.5,0.15,txt_a,'FontSize', 12);
    text(18,0.15,txt_b,'FontSize', 12);
    text(40,0.15,txt_g,'FontSize', 12);
    if pval_matrix_norm(i,1)<(0.05/100)
        text(2,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_norm(i,2)<(0.05/100)
        text(5.5,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_norm(i,3)<(0.05/100)
        text(9.5,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_norm(i,4)<(0.05/100)
        text(18,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_norm(i,5)<(0.05/100)
        text(40,0.12,txt_sig,'FontSize', 20);
    end
    
    % set axes, labels, subgraph titles
    set(gca, 'XScale','log', 'XTick',[0.5,4,8,13,30,80], 'XTickLabel',[0.5,4,8,13,30,80])
    ylabel('Spectral Density')
    xlabel('Frequency (Hz)')
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R'); % extract roi name without laterality
    split_roi_name{i,1} = this_roi_name(1);
    title(sprintf('%s',this_roi_name(1)), 'Interpreter', 'none') % use that as title
    xlim([0.5,80]) % xlimit 0.5-80 Hz
    ylim([0 0.2]) % ylimit 0-0.2 (arbitrary, just to visualize density)
    hold off
    
    figure(2)
    subplot(4,5,i)
    hold on
    plot([4,4],[0,0.4],'k-')
    plot([8,8],[0,0.4],'k-')
    plot([13,13],[0,0.4],'k-')
    plot([30,30],[0,0.4],'k-')
    plot(x_axis,median_normal)
    patch([x_axis fliplr(x_axis)], [prct_25_75_normal(1,:) fliplr(prct_25_75_normal(2,:))], color1, 'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x_axis,median_EIZ,'Color',color7)
    patch([x_axis fliplr(x_axis)], [prct_25_75_EIZ(1,:) fliplr(prct_25_75_EIZ(2,:))], color7, 'FaceAlpha',0.5, 'EdgeColor','none')
    try plot(x_axis,median_soz,'Color',color2)
    patch([x_axis fliplr(x_axis)], [prct_25_75_soz(1,:) fliplr(prct_25_75_soz(2,:))], color2, 'FaceAlpha',0.5, 'EdgeColor','none')
    catch anyerror
    end
    
    % add greek letters for frequency bands
    text(2,0.15,txt_d,'FontSize', 12);
    text(5.5,0.15,txt_t,'FontSize', 12);
    text(9.5,0.15,txt_a,'FontSize', 12);
    text(18,0.15,txt_b,'FontSize', 12);
    text(40,0.15,txt_g,'FontSize', 12);
    if pval_matrix_EIZ(i,1)<(0.05/100)
        text(2,0.12,txt_sig,'FontSize', 20); 
    end
    if pval_matrix_EIZ(i,2)<(0.05/100)
        text(5.5,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_EIZ(i,3)<(0.05/100)
        text(9.5,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_EIZ(i,4)<(0.05/100)
        text(18,0.12,txt_sig,'FontSize', 20);
    end
    if pval_matrix_EIZ(i,5)<(0.05/100)
        text(40,0.12,txt_sig,'FontSize', 20);
    end
    
    if pval_matrix_soz(i,1)<(0.05/100)
        text(2,0.17,txt_sig,'FontSize', 20);
    end
    if pval_matrix_soz(i,2)<(0.05/100)
        text(5.5,0.17,txt_sig,'FontSize', 20);
    end
    if pval_matrix_soz(i,3)<(0.05/100)
        text(9.5,0.17,txt_sig,'FontSize', 20);
    end
    if pval_matrix_soz(i,4)<(0.05/100)
        text(18,0.17,txt_sig,'FontSize', 20);
    end
    if pval_matrix_soz(i,5)<(0.05/100)
        text(40,0.17,txt_sig,'FontSize', 20);
    end
    
    % set axes, labels, subgraph titles
    set(gca, 'XScale','log', 'XTick',[0.5,4,8,13,30,80], 'XTickLabel',[0.5,4,8,13,30,80])
    ylabel('Spectral Density')
    xlabel('Frequency (Hz)')
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R'); % extract roi name without laterality
    split_roi_name{i,1} = this_roi_name(1);
    title(sprintf('%s',this_roi_name(1)), 'Interpreter', 'none') % use that as title
    xlim([0.5,80]) % xlimit 0.5-80 Hz
    ylim([0 0.2]) % ylimit 0-0.2 (arbitrary, just to visualize density)
    hold off
    
end

%% calculate univariate z-scores of neural activity
% power spectrum and entropy abnormality

all_f = {'delta','theta','alpha','beta','gamma'}; % frequency bands
all_m = {'bandpower','entropy'}; % metrics

c = 0; % c for column in z score matrix
for m = 1:2 % m for different metrics
    for f = 1:5 % f for frequency bands
        c = c+1;
        
        % first create the univariate atlas, taking only normal channels
        % the univariate atlas is mean/median + stdev of a given feature
        % across ROI
        [mean_vec,median_vec,std_vec,feat_vals] = create_univariate_atlas(all_wake_data(:,find(all_normal_ch)),[1:40],new_roi(all_normal_ch),sprintf('%s',all_m{m}),sprintf('%s',all_f{f}));
        
        % we call this function again on all of the data, but are only 
        % using it to extract all of the feature values on each channel
        [~,~,~,feat_vals] = create_univariate_atlas(all_wake_data,[1:40],new_roi,sprintf('%s',all_m{m}),sprintf('%s',all_f{f}));
        
        % now we test the univariate atlas to generate z scores
        [zscores, pt_zscores, pt_ind, roi_ind] = test_univariate_atlas(median_vec, std_vec, feat_vals, all_pts, [1:40], new_roi,'patient');
        
        % store the z scores
        univariate_zscores(:,c) = zscores;
        all_pt_zscores(:,c) = pt_zscores;
        
    end
end
%save('univariate_zscores.mat','univariate_zscores')

%% wavelet entropy validation across all regions
for i = 1:60
    for j = 1:size(all_wake_data,2)
    start_inds = (i-1)*200+1;
    end_inds = 200*i;
    all_wentropy(i,j) = (wentropy(all_wake_data((start_inds:end_inds),j),'shannon'));
    end
end

all_mean_wentropy = log(-1*median(all_wentropy));


figure(1);clf;
for i = 1:20
    % extract indices that correspond to a given ROI for HUP and MNI data
    % (combining data from L and R hemispheres)
    roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_EIZ = intersect([EIZ_good_ch;EIZ_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_soz = intersect([soz_good_ch;soz_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    
    entropy_normal = all_mean_wentropy([roi_mni;roi_hup]);
    entropy_EIZ = all_mean_wentropy(roi_EIZ);
    entropy_soz = all_mean_wentropy(roi_soz);
    
    pval_EIZ(i,1) = ranksum(entropy_normal,entropy_EIZ);
    pval_soz(i,1) = ranksum(entropy_normal,entropy_soz);
    pval_EIZ_vs_soz(i,1) = ranksum(entropy_EIZ,entropy_soz);

    % median power spectral density across channels
    %entropy_mni = median(pxx_norm(:,roi_mni),2);
    %median_hup = median(pxx_norm(:,roi_hup),2);
    plot_cell{1,1} = entropy_normal;
    plot_cell{1,2} = entropy_EIZ;
    plot_cell{1,3} = entropy_soz;
    
    subplot(4,5,i)
    hold on
    violin(plot_cell,'xlabel',{'','',''},'facecolor',[color1;color7;color2],'mc',[],'medc','k');%,'edgecolor','k');
    legend('off')
    txt_sig1 = '+';
    txt_sig2 = '*';
    txt_sig3 = '#';
    ylim([7 23])
    
    if pval_EIZ(i,1) < (0.05./20)
        plot([1,2],[20,20],'k-','LineWidth',2)
        text(1.45,20.25,txt_sig2,'FontSize', 20);
    end
    if pval_soz(i,1) < (0.05./20)
        plot([1,3],[21.75,21.75],'k-','LineWidth',2)
        text(2,22,txt_sig2,'FontSize', 20);
    end
    if pval_EIZ_vs_soz(i,1) < (0.05./20)
        plot([2,3],[20,20],'k-','LineWidth',2)
        text(2.5,20.25,txt_sig2,'FontSize', 20);
    end
    
    ylabel('Log -Entropy')
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R'); % extract roi name without laterality
    split_roi_name{i,1} = this_roi_name(1);
    title(sprintf('%s',this_roi_name(1)), 'Interpreter', 'none') % use that as title
    
    hold off
   
end



%% create atlas data - > edge-wise and patient-wise should be options
% for HUP data, aggregate electrodes to ignore 
% these should be poor-outcome/RZ/SOZ/>1 spike per hour

conn_band = {'delta','theta','alpha','beta','gamma'};
conn_type = {'corr','coh'};

% loop through types of connectivity and frequency bands
for c = 1:2
    for f = 1:5
        
        % loop through subjects
        for s = 1:max(all_pts)
            s 
            
            % extract subject data from all the data (columns are channels)
            pt_data = all_wake_data(:,[all_pts==s]);

            % call functional connectivity code (200 is Hz, 1 is window size in seconds)
            [pt_adj] = compute_FC(pt_data, 200, 1, conn_type{c},conn_band{f});

            % assign median adjacency matrix into cell structure
            all_pt_adj{s} = pt_adj;

        end

        % save adjacency matrices
        save(sprintf('adj_matrices/all_adj_%s_%s_May25.mat',conn_type{c},conn_band{f}),'all_pt_adj')

    end
end

%% assign atlas locations per patient into cell structure
for s = 1:max(all_pts)
    pt_loc{s} = new_roi(all_pts==s)';
end

% indices of all abnormal channels
all_abnormal_ch = zeros(5203,1);
all_abnormal_ch(abnormal_ch) = 1;

for s = 1:166
    
    % check abnormal channels for HUP patients
    if s>106
        
        % abnormal channel indices
        all_abn{s} = find(all_abnormal_ch(all_pts==s));
        if isempty(all_abn{s})
            all_abn{s} = [];
        end
    else
        % if MNI patient, no abnormal channels by definition
        all_abn{s} = [];
    end
end

%% functional cconnectivity atlas construction & z-score code

conn_band = {'delta','theta','alpha','beta','gamma'};
conn_type = {'corr','coh'};

% which feature
feat = 0;

% loop through connectivity and frequency
for c = 1:2
    for f = 1:5
        
        feat = feat+1;
        a = 0
        
        load(sprintf('adj_matrices/all_adj_%s_%s_May25.mat',conn_type{c},conn_band{f}));
        
        if c==2
            for s = 1:166
            all_pt_adj{s} = log(all_pt_adj{s});
            end
        end

        [conn_edge, std_edge, samples_edge, sem_edge, raw_atlas_edge] = create_atlas_by_edge(all_pt_adj, pt_loc, all_abn, [1:40]', 2);
        
        % call for MNI patients
        [conn_edge_mni, std_edge_mni, samples_edge_mni, ~, raw_atlas_edge_mni] = create_atlas_by_edge(all_pt_adj(1:106), pt_loc(1:106), all_abn(1:106), [1:40]', 1);
        
        % call for HUP patients
        [conn_edge_hup, std_edge_hup, samples_edge_hup, ~, raw_atlas_edge_hup] = create_atlas_by_edge(all_pt_adj(107:166), pt_loc(107:166), all_abn(107:166), [1:40]', 1);

        this_feat_conn(feat).mni = conn_edge_mni;
        this_feat_conn(feat).hup = conn_edge_hup;
        
        % loop through all patients
        for s = 1:166 % change this back
            a = a+1;

            % do nodal space
            try [native_adj_scores, corr_val] = test_native_adj(all_pt_adj{s}, pt_loc{s}, conn_edge, std_edge, [1:40]');

            bivariate_native(feat).subj(a).data = native_adj_scores;

            catch anyerror
                bivariate_native(feat).subj(a).data = NaN(length(pt_loc{s}));
            end

        end
    end
    
end

%% correlate connectivity between MNI and HUP
for feat = 1:10
    
    conn_edge_mni = this_feat_conn(feat).mni;
    conn_edge_hup = this_feat_conn(feat).hup;
    
conn_edge_all = conn_edge_mni(:)+conn_edge_hup(:);
conn_edge_mni(isnan(conn_edge_all)) = [];
conn_edge_hup(isnan(conn_edge_all)) = [];

std_edge_all = std_edge_mni(:)+std_edge_hup(:);
std_edge_mni(isnan(std_edge_all)) = [];
std_edge_hup(isnan(std_edge_all)) = [];

figure(feat);clf;
plot(exp(conn_edge_mni),exp(conn_edge_hup),'ko')
[r,p] = corr(exp(conn_edge_mni)',exp(conn_edge_hup)')
xlabel('MNI delta coherence')
ylabel('HUP delta coherence')
title(sprintf('Correlation between atlas edges, r = %f, p = %f',r,p))
end
%% process bivariate edge features into nodal features

for f = 1:10
    this_feat = [];
    abs_feat = [];
    for s = 1:166
        % take 80th percentile across all edges of each node
        abs_feat = [abs_feat;prctile(abs(bivariate_native(f).subj(s).data),75)'];
    end
    
    % assign into feature matrix
    abs_bivariate_feats(:,f) = abs_feat;
end

% for the single bivariate feature we want the maximum absolute Z score
% across all 10 bivariate features
[single_bivariate_feat, which_feat_bi] = max(abs_bivariate_feats,[],2);

% assign previously calculated univariate z scores into bivariate scores
all_feat_zscores = [univariate_zscores, abs_bivariate_feats];

% incomplete channels are ones which are unlocalized in new atlas
incomplete_channels = find(isnan(sum(all_feat_zscores')));
single_bivariate_feat(incomplete_channels) = []; 

univariate_feats = all_feat_zscores(:,1:10);
univariate_feats(incomplete_channels,:) = [];
[single_univariate_feat, which_feat_uni] = max(abs(univariate_feats),[],2);

%% combine univariate and bivariate
%all_feat_zscores = [univariate_zscores, all_bivariate_feats];
all_feat_zscores(isinf(all_feat_zscores)) = 0;
incomplete_channels = find(isnan(sum(all_feat_zscores')));
all_feat_zscores2 = all_feat_zscores;
all_feat_zscores2(incomplete_channels,:) = [];
EIZ = [zeros(1772,1);[spike_ind - spike_ind.*soz_ch]];
resect_ch_all = [zeros(1772,1);resected_ch];

% 1: patient, 2: abnormal, 3: EIZ, 4: SOZ, 5: resect, 6:outcome
all_possible_labels = [all_pts,[EIZ+[zeros(1772,1);soz_ch]],EIZ,[zeros(1772,1);soz_ch],resect_ch_all,[ones(1772,1);HUP_outcome_all]];
multi_class_all = [~all_possible_labels(:,2)+2.*all_possible_labels(:,3)+3.*all_possible_labels(:,4)];
all_possible_labels(incomplete_channels,:) = [];

multi_class = [~all_possible_labels(:,2)+2.*all_possible_labels(:,3)+3.*all_possible_labels(:,4)];

% get only non EIZ normal versus SOZ (in good outcome patients)
non_normal_non_soz = [zeros(1772,1);[soz_ch==0]];
non_normal_non_soz(incomplete_channels) = [];
eliminate_channels = find(all_possible_labels(:,3)==1);%find(non_normal_non_soz==1);
all_feat_zscores3 = all_feat_zscores2;
all_feat_zscores3(eliminate_channels,:) = [];
new_labels = all_possible_labels;
new_labels(eliminate_channels,:) = [];

all_feats = all_feat_zscores2;
all_labels = new_labels(:,4);
all_pred_1 = zeros(size(all_labels)); 
all_pred_2 = zeros(size(all_labels)); 
all_pred_3 = zeros(size(all_labels)); 

[part] = make_xval_partition(length(all_labels), 20);
for p = 1:10
    p
    X_train = all_feats(find(part~=p),:);
    Y_train = all_labels(find(part~=p));
    X_test = all_feats(find(part==p),:);
    Y_test = all_labels(find(part==p));
    
    Mdl1 = TreeBagger(100,X_train(:,1:10),Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on');

    Mdl2 = TreeBagger(100,X_train(:,11:20),Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on');

    Mdl3 = TreeBagger(100,X_train,Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on');

    [Yfit,pihat1] = predict(Mdl1,X_test(:,1:10));
    [Yfit,pihat2] = predict(Mdl2,X_test(:,11:20));
    [Yfit,pihat3] = predict(Mdl3,X_test);
    
    [~, Y_pred] = max(pihat1,[],2);
    all_pred_1(find(part==p)) = pihat1(:,2);
    all_pred_2(find(part==p)) = pihat2(:,2);
    all_pred_3(find(part==p)) = pihat3(:,2);
    
    xfold_acc(p) = sum(Y_pred==Y_test)./length(Y_test);
    
    Y_pred = round(pihat1(:,2));
    
    
end

    [X1,Y1,T,AUC1] = perfcurve(all_labels,all_pred_1,1)  

    [X2,Y2,T,AUC2] = perfcurve(all_labels,all_pred_2,1)

    [X3,Y3,T,AUC3] = perfcurve(all_labels,all_pred_3,1)
   

%% random forest feature importance
imp = Mdl3.OOBPermutedPredictorDeltaError;

figure(1);clf;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xticks([1:20])

figure(2);clf;
hold on
plot(X1,Y1,'LineWidth',2)%4
plot(X2,Y2,'LineWidth',2)%7
plot(X3,Y3,'LineWidth',2)%10
legend(sprintf('Univariate  - AUC: %f2',AUC1),sprintf('Bivariate  - AUC: %f2',AUC2),sprintf('All - AUC: %f2',AUC3),'Location','SouthEast')%86/85/91
xlabel('False positive rate')
ylabel('True positive rate')
title('All abnormal vs. normal channel classification')

%% need to map back onto all channels and render
all_pred = ones(5203,1);
all_pred(incomplete_channels) = 0;
all_pred(all_pred==1) = 1-xfold_pred3;

%% render abnormalities
soz_ch_all = [zeros(1772,1);soz_ch];
this_pt = 111; %114 res - 45/50 EIZ - 52/64, non-inv - 20,30 % 120 - res - 2., EIZ - 23./64, non-inv - 33/50.
this_pt_abn = all_pred_3(find(all_pts==this_pt));
this_pt_soz = multi_class_all(find(all_pts==this_pt));
this_pt_res = resect_ch_all(find(all_pts==this_pt));
this_pt_coords = all_coords(find(all_pts==this_pt),:);
this_pt_eeg = all_wake_data(:,find(all_pts==this_pt));
%this_pt_delta_coh_abn(isnan(this_pt_delta_coh_abn)) = 0;

this_pt_zscores = all_feat_zscores(find(all_pts==this_pt),:);

figure(1);clf;
hold on
for i = 1:10
    plot([1:12000],[this_pt_eeg(:,4*i)+i*175],'k-')
end
hold off

figure(1);clf;
imagesc(all_feat_zscores(find(all_pts==this_pt),:))
colorbar
caxis([0,5])

figure(2);clf;
subplot(1,3,1)
imagesc(this_pt_zscores(2,:)')
caxis([0,4])
colorbar
subplot(1,3,2)
imagesc(this_pt_zscores(23,:)')
caxis([0,4])
colorbar
subplot(1,3,3)
imagesc(this_pt_zscores(50,:)')
caxis([0,4])
colorbar

figure(3);
subplot(3,1,1)
plot(this_pt_eeg(:,2))
subplot(3,1,2)
plot(this_pt_eeg(:,23))
subplot(3,1,3)
plot(this_pt_eeg(:,50))


figure(4);clf;
final_elec_matrix = [this_pt_coords,this_pt_abn,ones(size(this_pt_coords,1),1)];
%final_elec_matrix = final_elec_matrix([2,23,50],:)
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','abnormality_render.mat')
delete render_elecs.node

%% 
single_bivariate_feat_2 = single_bivariate_feat;
new_roi_2 = new_roi;
new_roi_2(incomplete_channels) = [];

new_roi_2(new_roi_2==17) = 18;

for ch = 1:4702
    if multi_class(ch)==1
        grouping{ch,1} = 'uninvolved';
        class_ind(ch,1) = 1;
    elseif multi_class(ch)==2
        grouping{ch,1} = 'EIZ';
        class_ind(ch,1) = 2;
    else
        grouping{ch,1} = 'SOZ';
        class_ind(ch,1) = 3;
    end
end

univariate_feats = all_feat_zscores(:,1:10);
univariate_feats(incomplete_channels,:) = [];

single_univariate_feat = max(abs(univariate_feats'))';

single_univariate_feat(isinf(single_univariate_feat)) = nanmedian(single_univariate_feat);
single_bivariate_feat_2(isinf(single_bivariate_feat_2)) = nanmedian(single_bivariate_feat_2);

single_univariate_feat = single_univariate_feat(new_roi_2==18);
single_bivariate_feat_2 = single_bivariate_feat_2(new_roi_2==18);
grouping = grouping(new_roi_2==18);
class_ind = class_ind(new_roi_2==18);


figure(1);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat_2,'Group',grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,grouping,'orientation','horizontal',...
     'label',{'','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat_2,grouping,'orientation','horizontal',...
     'label', {'','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

d1 = computeCohen_d(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
p1 = ranksum(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
d2 = computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
p2 = ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
d3 = computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))
p3 = ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))

d4 = computeCohen_d(single_bivariate_feat_2(find(class_ind==2)),single_bivariate_feat_2(find(class_ind==1)))
p4 = ranksum(single_bivariate_feat_2(find(class_ind==2)),single_bivariate_feat_2(find(class_ind==1)))
d5 = computeCohen_d(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==2)))
p5 = ranksum(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==2)))
d6 = computeCohen_d(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==1)))
p6 = ranksum(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==1)))

[d1 d2 d3 d4 d5 d6]
[p1 p2 p3 p4 p5 p6]
%% do node abormality (max univariate+bivariate in plot) for following:
% 1. within/outside soz for lesional / non-lesional

all_lesion = zeros(5203,1);
for i = 1:166
    if i>106
    all_lesion(all_pts==i) = lesion_field(i-106);
    end
end
 
within_soz_lesional = find(all_lesion.*soz_ch_all);
within_soz_nonlesional = find([all_lesion==0].*soz_ch_all);
outside_soz_lesional = find(all_lesion.*[soz_ch_all==0]);
outside_soz_nonlesional = find([all_lesion==0].*[soz_ch_all==0]);

for ch = 1:5203
    if ismember(ch,within_soz_lesional)
        lesional_grouping{ch,1} = 'within_soz_lesional';
    elseif ismember(ch,within_soz_nonlesional)
        lesional_grouping{ch,1} = 'within_soz_nonlesional';
    elseif ismember(ch,outside_soz_lesional)
        lesional_grouping{ch,1} = 'outside_soz_lesional';
    elseif ismember(ch,outside_soz_nonlesional)
        lesional_grouping{ch,1} = 'outside_soz_nonlesional';
    end
end

univariate_feats = all_feat_zscores(:,1:10);
single_univariate_feat = prctile(abs(univariate_feats'),100)';
single_bivariate_feat = prctile(abs_bivariate_feats',100)';

uni_within_soz_lesional = single_univariate_feat(within_soz_lesional);
uni_within_soz_nonlesional = single_univariate_feat(within_soz_nonlesional);
uni_outside_soz_lesional = single_univariate_feat(outside_soz_lesional);
uni_outside_soz_nonlesional = single_univariate_feat(outside_soz_nonlesional);
bi_within_soz_lesional = single_bivariate_feat(within_soz_lesional);
bi_within_soz_nonlesional = single_bivariate_feat(within_soz_nonlesional);
bi_outside_soz_lesional = single_bivariate_feat(outside_soz_lesional);
bi_outside_soz_nonlesional = single_bivariate_feat(outside_soz_nonlesional);

p1 = ranksum(uni_within_soz_lesional,uni_within_soz_nonlesional);
p2 = ranksum(uni_within_soz_lesional,uni_outside_soz_lesional);
p3 = ranksum(uni_within_soz_nonlesional,uni_outside_soz_nonlesional);
p4 = ranksum(uni_outside_soz_lesional,uni_outside_soz_nonlesional);
p5 = ranksum(bi_within_soz_lesional,bi_within_soz_nonlesional);
p6 = ranksum(bi_within_soz_lesional,bi_outside_soz_lesional);
p7 = ranksum(bi_within_soz_nonlesional,bi_outside_soz_nonlesional);
p8 = ranksum(bi_outside_soz_lesional,bi_outside_soz_nonlesional);

single_univariate_feat(incomplete_channels) = [];
single_bivariate_feat(incomplete_channels) = [];
lesional_grouping(incomplete_channels) = [];

figure(2);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat,'Group',lesional_grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,lesional_grouping,'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat,lesional_grouping,'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

[p1 p2 p3 p4 p5 p6 p7 p8]
%% 2. within/outside ablation zone for good / poor outcome
all_outcome = ones(5203,1);
all_abl = NaN(5203,1);
for i = 1:166
    if i>106
    all_outcome(all_pts==i) = late_outcome(i-106);
        if therapy_field(i-106)
            all_abl(all_pts==i) = resect_ch_all(all_pts==i);
        end
    end
end

within_abl_good = find([all_outcome==1].*[all_abl==1]);
outside_abl_good = find([all_outcome==1].*[all_abl==0]);
within_abl_poor = find([all_outcome>1].*[all_abl==1]);
outside_abl_poor = find([all_outcome>1].*[all_abl==0]);

incomplete_channels2 = incomplete_channels;

abl_grouping = cell(5203,1);

for ch = 1:5203
    if ismember(ch,within_abl_good)
        abl_grouping{ch,1} = 'within_abl_good';
    elseif ismember(ch,outside_abl_good)
        abl_grouping{ch,1} = 'outside_abl_good';
    elseif ismember(ch,within_abl_poor)
        abl_grouping{ch,1} = 'within_abl_poor';
    elseif ismember(ch,outside_abl_poor)
        abl_grouping{ch,1} = 'outside_abl_poor';
    else
        incomplete_channels2 = [incomplete_channels2, ch];
    end
end

incomplete_channels2 = unique(incomplete_channels2);

univariate_feats = all_feat_zscores(:,1:10);
single_univariate_feat = prctile(abs(univariate_feats'),100)';
single_bivariate_feat = prctile(abs_bivariate_feats',100)';

uni_within_abl_good = single_univariate_feat(within_abl_good);
uni_outside_abl_good = single_univariate_feat(outside_abl_good);
uni_within_abl_poor = single_univariate_feat(within_abl_poor);
uni_outside_abl_poor = single_univariate_feat(outside_abl_poor);
bi_within_abl_good = single_bivariate_feat(within_abl_good);
bi_outside_abl_good = single_bivariate_feat(outside_abl_good);
bi_within_abl_poor = single_bivariate_feat(within_abl_poor);
bi_outside_abl_poor = single_bivariate_feat(outside_abl_poor);

p1 = ranksum(uni_within_abl_good,uni_outside_abl_good);
p2 = ranksum(uni_within_abl_poor,uni_outside_abl_poor);
p3 = ranksum(uni_within_abl_good,uni_within_abl_poor);
p4 = ranksum(uni_outside_abl_good,uni_outside_abl_poor);
p5 = ranksum(bi_within_abl_good,bi_outside_abl_good);
p6 = ranksum(bi_within_abl_poor,bi_outside_abl_poor);
p7 = ranksum(bi_within_abl_good,bi_within_abl_poor);
p8 = ranksum(bi_outside_abl_good,bi_outside_abl_poor);

single_univariate_feat(incomplete_channels2) = [];
single_bivariate_feat(incomplete_channels2) = [];
abl_grouping(incomplete_channels2) = [];

figure(3);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat,'Group',abl_grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,abl_grouping,'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat,abl_grouping,'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

[p1 p2 p3 p4 p5 p6 p7 p8]

%% Examine correlations between z scores
[r,p] = corr(all_feat_zscores2)
figure(1);clf;
r(p>0.0001) = NaN;
imagesc(r)
colorbar


