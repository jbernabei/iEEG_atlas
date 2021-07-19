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
    late_outcome_raw(s) = engel_scores(end);
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
all_m = {'bandpower'};%,'entropy'}; % metrics

c = 0; % c for column in z score matrix
for m = 1 % m for different metrics
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


%% create atlas data - > edge-wise and patient-wise should be options
% for HUP data, aggregate electrodes to ignore 
% these should be poor-outcome/RZ/SOZ/>1 spike per hour

conn_band = {'delta','theta','alpha','beta','gamma'};
conn_type = {'corr','coh'};

% loop through types of connectivity and frequency bands
for c = 2
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
conn_type = {'coh'};

% which feature
feat = 0;

% loop through connectivity and frequency
for c = 1
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
for feat = 1:5
        a = conn_band{feat};
        aa = 'coherence';
    
    conn_edge_mni = this_feat_conn(feat).mni;
    conn_edge_hup = this_feat_conn(feat).hup;
    
conn_edge_all = conn_edge_mni(:)+conn_edge_hup(:);
conn_edge_mni(isnan(conn_edge_all)) = [];
conn_edge_hup(isnan(conn_edge_all)) = [];

std_edge_all = std_edge_mni(:)+std_edge_hup(:);
std_edge_mni(isnan(std_edge_all)) = [];
std_edge_hup(isnan(std_edge_all)) = [];

figure(1);clf
%subplot(1,5,feat)
hold on
plot((conn_edge_mni),(conn_edge_hup),'k.','Color',[128 128 128]/255)
[r,p1] = corr((conn_edge_mni)',(conn_edge_hup)')
rvals(feat) = r;
xlabel('MNI connectivity')
ylabel('HUP connectivity')
xlim([0,0.15])
ylim([0,0.15])
p = polyfit((conn_edge_mni),(conn_edge_hup),1);
f = polyval(p,(conn_edge_mni));
plot((conn_edge_mni),f,'k-','LineWidth',2)
title(sprintf('%s %s, r = %f, p = %f',a,aa,r,p1))
hold off
end

figure(2);clf
bar(rvals)
%% process bivariate edge features into nodal features

for f = 1:5
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
[single_bivariate_feat] = mean(abs_bivariate_feats')';

% assign previously calculated univariate z scores into bivariate scores
all_feat_zscores = [abs(univariate_zscores), abs_bivariate_feats];

% incomplete channels are ones which are unlocalized in new atlas
incomplete_channels = find(isnan(sum(all_feat_zscores')));
%single_bivariate_feat(incomplete_channels) = []; 

univariate_feats = all_feat_zscores(:,1:5);
univariate_feats(incomplete_channels,:) = [];
%[single_univariate_feat, which_feat_uni] = max(abs(univariate_feats),[],2);

%% combine univariate and bivariate
%all_feat_zscores = [univariate_zscores, all_bivariate_feats];
all_feat_zscores(isinf(all_feat_zscores)) = 0;
incomplete_channels = find(isnan(sum(all_feat_zscores')));
complete_channels = find(~isnan(sum(all_feat_zscores')));
all_feat_zscores2 = all_feat_zscores;
all_feat_zscores2(incomplete_channels,:) = [];
EIZ = [zeros(1772,1);[spike_ind - spike_ind.*soz_ch]];
resect_ch_all = [zeros(1772,1);resected_ch];

% 1: patient, 2: abnormal, 3: EIZ, 4: SOZ, 5: resect, 6:outcome
all_possible_labels = [all_pts,[EIZ+[zeros(1772,1);soz_ch]],EIZ,[zeros(1772,1);soz_ch],resect_ch_all,[ones(1772,1);HUP_outcome_all],[[zeros(1772,1);soz_ch].*resect_ch_all]];
all_possible_labels(incomplete_channels,:) = [];

% get only non EIZ normal versus SOZ (in good outcome patients)
eliminate_channels = [find(all_possible_labels(:,3)==1);find((all_possible_labels(:,4)+all_possible_labels(:,5))==1)];%find(non_normal_non_soz==1);
retain_channels = find(all_possible_labels(:,3)==0);
all_feat_zscores3 = all_feat_zscores2;
all_feat_zscores3(eliminate_channels,:) = [];
new_labels = all_possible_labels;
new_labels(eliminate_channels,:) = [];

all_feats = all_feat_zscores3;
all_labels = new_labels(:,7); % was column 5 before
all_pred_1 = zeros(size(all_labels)); 
all_pred_2 = zeros(size(all_labels)); 
all_pred_3 = zeros(size(all_labels)); 

num_parts = 10;

[pt_part] = make_xval_partition(166, num_parts);

% cross-validate based on patients 
for p = 1:num_parts
    indices = find(pt_part==p)
    for i = 1:length(indices)
    part(find(indices(i)==new_labels(:,1))) = p;
    
    part2(find(indices(i)==all_possible_labels(:,1))) = p;
    end
end

% Define classes, 0 = non seizure, 1 = seizure
ClassNames = [0,1];

% Define costs -> penalty here is 5
cost.ClassNames = ClassNames;
cost.ClassificationCosts = [0 1; 1 0];

for p = 1:num_parts
    p
    X_train = all_feats(find(part~=p),:);
    Y_train = all_labels(find(part~=p));
    X_test = all_feats(find(part==p),:);
    Y_test = all_labels(find(part==p));
    
    X_analysis = all_feat_zscores2(find(part2==p),:);
    
    Mdl1 = TreeBagger(100,X_train(:,1:5),Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on','cost',cost);

    Mdl2 = TreeBagger(100,X_train(:,6:10),Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on','cost',cost);

    Mdl3 = TreeBagger(100,X_train,Y_train,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on','cost',cost);

    [Yfit,pihat1] = predict(Mdl1,X_test(:,1:5));
    [Yfit,pihat2] = predict(Mdl2,X_test(:,6:10));
    [Yfit,pihat3] = predict(Mdl3,X_test);
    [Yfit, pihat_analysis] = predict(Mdl3,X_analysis);
    [Yfit, pihat_analysis2] = predict(Mdl1,X_analysis(:,1:5));
    [Yfit, pihat_analysis3] = predict(Mdl2,X_analysis(:,6:10));
    
    [~, Y_pred] = max(pihat1,[],2);
    all_pred_1(find(part==p)) = pihat1(:,2);
    all_pred_2(find(part==p)) = pihat2(:,2);
    all_pred_3(find(part==p)) = pihat3(:,2);
    all_pred_analysis(find(part2==p),1) = pihat_analysis(:,2);
    all_pred_analysis2(find(part2==p),1) = pihat_analysis2(:,2);
    all_pred_analysis3(find(part2==p),1) = pihat_analysis3(:,2);
    
    xfold_acc(p) = sum(Y_pred==Y_test)./length(Y_test);
    
    Y_pred = round(pihat1(:,2));

    
end

    [X1,Y1,T,AUC1] = perfcurve(all_labels,all_pred_1,1);  
    [X2,Y2,T,AUC2] = perfcurve(all_labels,all_pred_2,1);
    [X3,Y3,T,AUC3] = perfcurve(all_labels,all_pred_3,1);

    AUC1
    AUC2
    AUC3
   

%% random forest feature importance
imp = Mdl3.OOBPermutedPredictorDeltaError;

figure(1);clf;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xticks([1:10])

figure(2);clf;
hold on
plot(X1,Y1,'LineWidth',2,'Color',color4)%4
plot(X2,Y2,'LineWidth',2,'Color',color7)%7
plot(X3,Y3,'LineWidth',2,'Color',color1)%10
legend(sprintf('Univariate  - AUC: %.2f',AUC1),sprintf('Bivariate  - AUC: %.2f',AUC2),sprintf('All - AUC: %.2f',AUC3),'Location','SouthEast')%86/85/91
xlabel('False positive rate')
ylabel('True positive rate')
title('All abnormal vs. normal channel classification')
%% correlation
[r,p] = corr(all_feat_zscores2);
figure(1);clf;
subplot(1,2,1)
imagesc(r-eye(10))
caxis([-1 1])
colormap(color_bar)
colorbar
subplot(1,2,2)
imagesc(p<(0.05./100))
%% need to map back onto all channels and render
all_pred = zeros(5203,1);
all_pred(incomplete_channels) = 0;
all_pred(complete_channels) = all_pred_analysis;

all_pred2 = zeros(5203,1);
all_pred2(incomplete_channels) = 0;
all_pred2(complete_channels) = all_pred_analysis2;

all_pred3 = zeros(5203,1);
all_pred3(incomplete_channels) = 0;
all_pred3(complete_channels) = all_pred_analysis3;

%% render abnormalities
soz_ch_all = [zeros(1772,1);soz_ch];
this_pt = 116;%116/122/154 / 154 %114 res - 45/50 EIZ - 52/64, non-inv - 20,30 % 120 - res - 2., EIZ - 23./64, non-inv - 33/50.
this_pt_abn = all_pred(find(all_pts==this_pt));
this_pt_soz = soz_ch_all(find(all_pts==this_pt));
this_pt_eiz = EIZ(find(all_pts==this_pt));
this_pt_res = resect_ch_all(find(all_pts==this_pt));
this_pt_coords = all_coords(find(all_pts==this_pt),:);
this_pt_eeg = all_wake_data(:,find(all_pts==this_pt));
%this_pt_delta_coh_abn(isnan(this_pt_delta_coh_abn)) = 0;

this_pt_multiclass = 2*this_pt_soz+this_pt_eiz;

this_pt_zscores = all_feat_zscores(find(all_pts==this_pt),:);

this_pt_test = zeros(size(this_pt_res));
this_pt_test([28,9,49])=1;

% figure(1);clf;
% hold on
% for i = 1:10
%     plot([1:2000],[this_pt_eeg(:,4*i)+i*175],'k-')
% end
% hold off

% figure(1);clf;
% imagesc(all_feat_zscores(find(all_pts==this_pt),:))
% colorbar
% caxis([0,5])

figure(2);clf;
subplot(1,3,1)
imagesc(this_pt_zscores(28,:)')
caxis([0,3])
colormap('jet')
colorbar
colorbar
subplot(1,3,2)
imagesc(this_pt_zscores(9,:)')
caxis([0,3]) 
colormap('jet')
colorbar
colorbar
subplot(1,3,3)
imagesc(this_pt_zscores(49,:)')
caxis([0,3])
colormap('jet')
colorbar
colorbar
% 
figure(3);
subplot(3,1,1)
plot(this_pt_eeg(6001:12000,4))
subplot(3,1,2)
plot(this_pt_eeg(6001:12000,23))
subplot(3,1,3)
plot(this_pt_eeg(6001:12000,27))


figure(4);clf;
final_elec_matrix = [this_pt_coords,this_pt_res,ones(size(this_pt_res))];
%final_elec_matrix = final_elec_matrix([2,23,50],:)
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','abnormality_render.mat')
%delete render_elecs.node

%% patient-specific correlation between abnormality features
for s = 107:166
    pt_zscores = all_feat_zscores2(find([all_possible_labels(:,1)==s]),:);
    
    pt_zscores_univar = pt_zscores(:,1:5);
    pt_zscores_bivar = pt_zscores(:,6:10);
    
    corrval2(s-106) = corr(mean(pt_zscores_univar')',mean(pt_zscores_bivar')','Type','Spearman');
    
    for f = 1:5
        corrval((s-106),f) = corr(pt_zscores(:,f),pt_zscores(:,(f+5)),'Type','Spearman');
    end
    
    corrval_outcome((s-106),1) = late_outcome(s-106);
    corrval_lesion((s-106),1) = lesion_field(s-106);
end

p1 = ranksum(corrval(corrval_outcome==1,1),corrval(corrval_outcome>1,1),'tail','right')
p2 = ranksum(corrval(corrval_outcome==1,2),corrval(corrval_outcome>1,2),'tail','right')
p3 = ranksum(corrval(corrval_outcome==1,3),corrval(corrval_outcome>1,3),'tail','right')
p4 = ranksum(corrval(corrval_outcome==1,4),corrval(corrval_outcome>1,4),'tail','right')
p5 = ranksum(corrval(corrval_outcome==1,5),corrval(corrval_outcome>1,5),'tail','right')

p6 = ranksum(corrval2(corrval_outcome==1),corrval2(corrval_outcome>1),'tail','right')

p7 = ranksum(corrval(corrval_lesion==0,1),corrval(corrval_lesion>0,1),'tail','right')
p8 = ranksum(corrval(corrval_lesion==0,2),corrval(corrval_lesion>0,2),'tail','right')
p9 = ranksum(corrval(corrval_lesion==0,3),corrval(corrval_lesion>0,3),'tail','right')
p10 = ranksum(corrval(corrval_lesion==0,4),corrval(corrval_lesion>0,4),'tail','right')
p11 = ranksum(corrval(corrval_lesion==0,5),corrval(corrval_lesion>0,5),'tail','right')

p12 = ranksum(corrval2(corrval_lesion==0),corrval2(corrval_lesion>0),'tail','right')


figure(1);clf;
plot_matrix = padcat(corrval(corrval_outcome==1,2),corrval(corrval_outcome>1,2));
figure(1);clf;
UnivarScatter(plot_matrix)
xticks([1,2])
xticklabels({'Engel 1','Engel 2+'})
ylim([-0.5 0.5])
ylabel('Correlation of beta power & coherence |Z|')

%% compute AUPRC/DRS at a per-patient level
a = 0;
for s = [107:127,128:132,134:139,141:155,157:161,163:166]
    s
    a = a+1;
    this_pt_label = all_possible_labels(find(all_possible_labels(:,1)==s),4);
    this_pt_pred = all_pred_analysis(find(all_possible_labels(:,1)==s));
    this_pt_pred2 = all_pred_analysis2(find(all_possible_labels(:,1)==s));
    this_pt_pred3 = all_pred_analysis3(find(all_possible_labels(:,1)==s));
    [X1,Y1,T,AUC1] = perfcurve(this_pt_label,this_pt_pred,1,'XCrit','TPR','YCrit','PPV');  
    [X2,Y2,T,AUC2] = perfcurve(this_pt_label,this_pt_pred2,1,'XCrit','TPR','YCrit','PPV');
    [X3,Y3,T,AUC3] = perfcurve(this_pt_label,this_pt_pred3,1,'XCrit','TPR','YCrit','PPV');
    
    auc_pt_outcome(a,1) = late_outcome(s-106);
    pt_auc(a,1:3) = [AUC1,AUC2,AUC3];
end

good_inds = find(auc_pt_outcome==1);
poor_inds = find(auc_pt_outcome>1);

plot_matrix = padcat(pt_auc(good_inds,1),pt_auc(poor_inds,1));
figure(1);clf;
UnivarScatter(plot_matrix)

ranksum(plot_matrix(:,1),plot_matrix(:,2))
xticks([1,2])
xticklabels({'Engel 1','Engel 2+'})
ylabel('Area under precision-recall curve')
title('Comparison between AUPRC, rank-sum p = 0.033')

%%
[single_univariate_feat, which_feat_uni] = max(all_feat_zscores(:,1:5),[],2);
[single_bivariate_feat, which_feat_bi] = max(all_feat_zscores(:,6:10),[],2);
figure(1);clf;
plot(abs(single_univariate_feat),single_bivariate_feat,'ko')
%% all-node abnormality score all regions
for ch = 1:4717
    if all_possible_labels(ch,2)==0
        grouping{ch,1} = 'uninvolved';
        class_ind(ch,1) = 1;
    elseif all_possible_labels(ch,3)==1
        grouping{ch,1} = 'EIZ';
        class_ind(ch,1) = 2;
    elseif all_possible_labels(ch,4)==1
        grouping{ch,1} = 'SOZ';
        class_ind(ch,1) = 3;
    end
end

new_roi_2 = new_roi;
new_roi_2(incomplete_channels) = [];

single_bivariate_feat_2 = single_bivariate_feat;
single_bivariate_feat_2(incomplete_channels) = [];

for f = 1:5
    for r = [1:16,17:20]
        
        new_roi_2(new_roi_2==(2*r-1)) = 2*r;
        
        pred_uninvolved = all_feat_zscores2(find((class_ind==1).*(new_roi_2==2*r)),(f+5));
        pred_eiz = all_feat_zscores2(find((class_ind==2).*(new_roi_2==2*r)),(f+5));
        pred_soz = all_feat_zscores2(find((class_ind==3).*(new_roi_2==2*r)),(f+5));
        
        median(pred_uninvolved)
        median(pred_eiz)
        median(pred_soz)
        
        plot_matrix = padcat(pred_uninvolved, pred_eiz,pred_soz);
        
        % figure(1);clf;
        % plot_cell{1,1} = pred_uninvolved';
        % plot_cell{1,2} = pred_eiz';
        % plot_cell{1,3} = pred_soz';
        % violin(plot_cell,'xlabel',{'','',''},'facecolor',[color1;color7;color2],'mc',[],'medc','k');%,'edgecolor','k');
        % xlim([0.5 3.5])
        
        p1_region(f,r) = ranksum(pred_uninvolved,pred_eiz)
        p2_region(f,r) = ranksum(pred_uninvolved,pred_soz)
        p3_region(f,r) = ranksum(pred_eiz,pred_soz)
        
        cohen_d2(f,r) = computeCohen_d(pred_uninvolved,pred_soz)
        cohen_d3(f,r) = computeCohen_d(pred_eiz,pred_soz)
        
    end
end

figure(1);clf;
subplot(3,1,1)
imagesc([p3_region])
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
colorbar
subplot(3,1,2)
imagesc([p3_region]<0.05)
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
subplot(3,1,3)
imagesc([cohen_d3])
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
caxis([0,1])
colorbar

%% patient-level abnormality score all regions

    pt_uninvolved = [];
    pt_eiz = [];
    pt_soz = [];
    
for s = 1:166
    
    pt_uninvolved_inds = find([all_possible_labels(:,2)==0].*[all_possible_labels(:,1)==s].*[new_roi_2==18]);
    pt_eiz_inds = find([all_possible_labels(:,3)==1].*[all_possible_labels(:,1)==s].*[new_roi_2==18]);
    pt_soz_inds = find([all_possible_labels(:,4)==1].*[all_possible_labels(:,1)==s].*[new_roi_2==18]);
    
    pt_uninvolved = [pt_uninvolved; nanmean(all_pred_analysis(pt_uninvolved_inds))];
    pt_eiz = [pt_eiz; nanmean(all_pred_analysis(pt_eiz_inds))];
    pt_soz = [pt_soz; nanmean(all_pred_analysis(pt_soz_inds))];
end

plot_matrix = padcat(pt_uninvolved,pt_eiz,pt_soz);

figure(1);clf;
plot_cell{1,1} = rmmissing(pt_uninvolved);
plot_cell{1,2} = rmmissing(pt_eiz);
plot_cell{1,3} = rmmissing(pt_soz);
violin(plot_cell,'xlabel',{'','',''},'facecolor',[color1;color7;color2],'mc',[],'medc','k');%,'edgecolor','k');
xlim([0.5 3.5])

figure(2);clf;
UnivarScatter(plot_matrix)
xticks([1:3])
xticklabels({'uninvolved','irritative zone','seizure onset zone'})
ylabel('predicted epileptogenicity score')

p1 = ranksum(pt_uninvolved,pt_eiz)
p2 = ranksum(pt_eiz,pt_soz)
p3 = ranksum(pt_uninvolved,pt_soz)

%% patient-level abnormality score res in/out good/poor outcome
good_rz = [];
good_out = [];

poor_rz = [];
poor_out = [];

for s = 107:166
    this_pt_outcome = late_outcome(s-106);
    this_pt_ind_rz = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,5)==1]);
    this_pt_ind_out = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,5)==0]);
    
    if this_pt_outcome==1
        good_rz = [good_rz; nanmedian(all_pred_analysis(this_pt_ind_rz))];
        good_out = [good_out; nanmedian(all_pred_analysis(this_pt_ind_out))];
    else
        poor_rz = [poor_rz; nanmedian(all_pred_analysis(this_pt_ind_rz))];
        poor_out = [poor_out; nanmedian(all_pred_analysis(this_pt_ind_out))];
    end
end

plot_matrix = padcat(good_rz,good_out,poor_rz,poor_out);
figure(2);clf;
UnivarScatter(plot_matrix)
%% patient-level abnormality score soz in/out lesional

pt_in_rz = [];
pt_out_rz = [];

for s = 107:166
    
    w_inds = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,4)==1]);
    o_inds = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,4)==0]);
    
    length_diagnosis = age_surgery - age_onset;
    
    pt_in_rz = [pt_in_rz; nanmedian(single_bivariate_feat_2(w_inds))];
    pt_out_rz = [pt_out_rz; nanmedian(single_bivariate_feat_2(o_inds))];

end

figure(1);clf;
hold on
plot(length_diagnosis, pt_out_rz,'b.','Color',[106 178 180]/255)
p = polyfit(length_diagnosis,pt_out_rz,1);
f = polyval(p,length_diagnosis);
plot(length_diagnosis,f,'b-','LineWidth',2)
[a2,b2] = corr(length_diagnosis, pt_out_rz)
ylim([0.5 2])

%% Examine correlations between z scores
[r,p] = corr(all_feat_zscores2)
figure(1);clf;
r(p>0.0001) = NaN;
imagesc(r)
colorbar

%%
res_roi = new_roi(resect_ch_all==1);
for i = 1:20
    roi_count(i,1) = sum(res_roi==2*i)+sum(res_roi==(2*i-1));
end

