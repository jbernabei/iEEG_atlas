%%
clear all

% set up path - change iEEG_atlas_path to where you download the repository
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas_dev';
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
    early_outcome(s) = floor(engel_scores(1));
    late_outcome(s) = floor(engel_scores(end));
end

therapy_field = metadata{:,6};
implant_field = metadata{:,7};
target_field = metadata{:,8};
laterality_field = metadata{:,9};
lesion_field = metadata{:,10};
gender_field = metadata{:,13};

lesion_field = strcmp(lesion_field,'Lesional')
therapy_field = strcmp(therapy_field,'Ablation')

load('data/MNI_atlas_new.mat')
load('data/HUP_atlas_May11.mat')

all_pts = [Patient; (patient_no+106)];

HUP_outcome_all = zeros(length(soz_ch),1);
for i = 1:60
    HUP_outcome_all(find(patient_no==i)) = late_outcome(i);
end

clear NodesLeft
clear NodesLeftInflated
clear NodesRightInflated
clear NodesRight
clear NodesRegionLeft
clear NodesRegionRight
clear sleep_clip
clear Data_N2
clear Data_N3
clear Data_R

    color1 = [0, 0.4470, 0.7410];
    color2 = [0.6350, 0.0780, 0.1840];
    color3 = [255, 248, 209]./255; 
    color4 = [103 55 155]/255;
    color6 = [78 172 91]/255;
    
%%
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
spike_thresh = 100; % this is empirical

% find indices for which spikes are greater than threshold
spike_ind = [spike_24h>spike_thresh];

% define all abnormal channels
abnormal_ch = find([resected_ch+spike_ind+soz_ch+(HUP_outcome_all>1)]>0)+1772;

% define all seizure onset indices
soz_ch_inds = find(soz_ch)+1772;

% define normal HUP channels
normal_HUP_ch = find([resected_ch+spike_ind+soz_ch+(HUP_outcome_all>1)]==0)+1772;

% define normal MNI channels
normal_MNI_ch = [1:1772]';

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
end

new_counts = [num_ch_roi(:,1),(num_ch_roi(:,2))];

figure(1);clf;
subplot(1,2,1)
barh(roi_count)
set(gca,'TickLabelInterpreter', 'none');
yticks([1:45])
yticklabels(split_roi_name_old)
subplot(1,2,2)
barh(new_counts)
set(gca,'TickLabelInterpreter', 'none');
yticks([1:21])
yticklabels(split_roi_name)
% 106 patients versus 38 patients

%%
all_wake_data = [Data_W,wake_clip];

[pxx,f] = pwelch(all_wake_data,200,100,[0.5:0.5:80],200);
pxx_norm = pxx./sum(pxx);
sum(pxx_norm)

%% new areas to test

figure(1);clf
for i = 1:21
    % extract indices that correspond to a given ROI for HUP and MNI data
    % (combining data from L and R hemispheres)
    roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);

    % median power spectral density across channels
    median_mni = median(pxx_norm(:,roi_mni),2);
    median_hup = median(pxx_norm(:,roi_hup),2);

    % find upper and lower bounds (here 25th and 75th percentiles)
    prct_25_75_mni = prctile(pxx_norm(:,roi_mni)',[25,75]);
    prct_25_75_hup = prctile(pxx_norm(:,roi_hup)',[25,75]);

    x_axis = [0.5:0.5:80];
    
    % 3 x 7 plot for 21 regions (aggregating across hemispheres)
    subplot(3,7,i)
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
    text(2,0.15,txt_d);
    text(5.5,0.15,txt_t);
    text(9.5,0.15,txt_a);
    text(18,0.15,txt_b);
    text(40,0.15,txt_g);
    
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
        [mean_vec,median_vec,std_vec,feat_vals] = create_univariate_atlas(all_wake_data(:,all_normal_ch),[1:42],new_roi(all_normal_ch),sprintf('%s',all_m{m}),sprintf('%s',all_f{f}));
        
        % we call this function again on all of the data, but are only 
        % using it to extract all of the feature values on each channel
        [~,~,~,feat_vals] = create_univariate_atlas(all_wake_data,[1:42],new_roi,sprintf('%s',all_m{m}),sprintf('%s',all_f{f}));
        
        % now we test the univariate atlas to generate z scores
        [zscores, pt_zscores, pt_ind, roi_ind] = test_univariate_atlas(median_vec, std_vec, feat_vals, all_pts, [1:42], new_roi,'patient');
        
        % store the z scores
        univariate_zscores(:,c) = zscores;
        all_pt_zscores(:,c) = pt_zscores;
        
    end
end
%save('univariate_zscores.mat','univariate_zscores')

%% mass univariate testing
% 21 ROI x 5  = 105 tests to correct for
freq_inds = [1,  8; % 0.5 - 4 Hz
             9,  16; % 4 - 8 Hz
             17, 26; % 8 - 13 Hz
             27, 60; % 13 - 30 Hz 
             61, 160]; % 30 - 80 Hz

for i = 1:21
    for j = 1:5
        roi_soz = intersect((find(soz_ch)+1772),[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
        roi_eiz = intersect((EIZ_ch),[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
        roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
        roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
        soz_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_soz))';
        eiz_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_eiz))';
        HUP_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_hup))';
        MNI_feat = median(pxx_norm([freq_inds(j,1):freq_inds(j,2)],roi_mni))';
        [p,h] = ranksum(HUP_feat,MNI_feat); % HUP vs MNI
        pval_matrix_norm(i,j) = p;
        [p1,h] = ranksum([HUP_feat;MNI_feat],soz_feat);
        pval_matrix_soz(i,j) = p1;
        [p2,h] = ranksum([HUP_feat;MNI_feat],eiz_feat);
        pval_matrix_eiz(i,j) = p2;
    end
end

figure(1);clf;
subplot(1,3,1)
imagesc((pval_matrix_norm.*105)<0.05)
title('HUP vs. MNI normal channels')
xticks([1:5])
xticklabels({'Delta','Theta','Alpha','Beta','Gamma'})
yticks([1:21])
set(gca,'TickLabelInterpreter', 'none');
yticklabels(split_roi_name)
subplot(1,3,2)
imagesc((pval_matrix_eiz.*105)<0.05)
title('Irritative zone vs. composite atlas')
xticks([1:5])
xticklabels({'Delta','Theta','Alpha','Beta','Gamma'})
yticks([1:21])
set(gca,'TickLabelInterpreter', 'none');
yticklabels(split_roi_name)
subplot(1,3,3)
imagesc((pval_matrix_soz.*105)<0.05)
title('Seizure onset zone  vs. composite atlas')
xticks([1:5])
xticklabels({'Delta','Theta','Alpha','Beta','Gamma'})
yticks([1:21])
set(gca,'TickLabelInterpreter', 'none');
yticklabels(split_roi_name)

%% wavelet entropy validation across all regions
for i = 1:60
    for j = 1:size(all_wake_data,2)
    start_inds = (i-1)*200+1;
    end_inds = 200*i;
    all_wentropy(i,j) = wentropy(all_wake_data((start_inds:end_inds),j),'shannon');
    end
end

all_mean_wentropy = log(-1*median(all_wentropy));

r1 = all_mean_wentropy([normal_HUP_ch;normal_MNI_ch]);
r2 = all_mean_wentropy(abnormal_ch);
r3 = all_mean_wentropy(EIZ_good_ch);
r4 = all_mean_wentropy(EIZ_poor_ch);
r5 = all_mean_wentropy(soz_good_ch);
r6 = all_mean_wentropy(soz_poor_ch);
r7 = all_mean_wentropy(poor_out_ch);

ranksum(all_mean_wentropy(normal_HUP_ch),all_mean_wentropy(normal_MNI_ch))

%ranksum(r1,r2)
p1 = ranksum(r1,r3)
p2 = ranksum(r1,r4)
p3 = ranksum(r3,r4)
p4 = ranksum(r1,r5)
p5 = ranksum(r1,r6)
p6 = ranksum(r5,r6)
p7 = ranksum(r1,r7)
p8 = ranksum(r3,r5)
p9 = ranksum(r4,r6)

figure(1);clf;
hold on
scatter(ones(length(r1),1),r1,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(r1) median(r1)],'k-','LineWidth',2)
scatter(2*ones(length(r7),1),r7,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(r7) median(r7)],'k-','LineWidth',2)
scatter(3*ones(length(r2),1),r2,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([2.75 3.25],[median(r2) median(r2)],'k-','LineWidth',2)
scatter(4*ones(length(r3),1),r3,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([3.75 4.25],[median(r3) median(r3)],'k-','LineWidth',2)
scatter(5*ones(length(r5),1),r5,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([4.75 5.25],[median(r5) median(r5)],'k-','LineWidth',2)
scatter(6*ones(length(r6),1),r6,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([5.75 6.25],[median(r6) median(r6)],'k-','LineWidth',2)
xlim([0.5, 6.5])
ylim([2,24])
ylabel('Wavelet entropy')
xticks([1:6])
xticklabels({'Engel 1, non-SOZ, non-EIZ','Engel 2+, non-SOZ, non-EIZ', 'Engel 1 EIZ','Engel 2+ EIZ','Engel 1 SOZ','Engel 2+ SOZ'})
hold off


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
        save(sprintf('all_adj_%s_%s.mat',conn_type{c},conn_band{f}),'all_pt_adj')

    end
end

%% assign atlas locations per patient into cell structure
for s = 1:max(all_pts)
    pt_loc{s} = new_roi(all_pts==s)';
end

% indices of all abnormal channels
all_abnormal_ch = zeros(4156,1);
all_abnormal_ch(abnormal_ch) = 1;

for s = 1:166
    
    % check abnormal channels for HUP patients
    if s>106
        
        % abnorma channel indices
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
        
        load(sprintf('all_adj_%s_%s.mat',conn_type{c},conn_band{f}));
        
        if c==2
            for s = 1:166
            all_pt_adj{s} = log(all_pt_adj{s});
            end
        end

        [conn_edge, std_edge, samples_edge, sem_edge, raw_atlas_edge] = create_atlas_by_edge(all_pt_adj, pt_loc, all_abn, [1:42]', 1);
        
        % call for MNI patients
        [~, ~, samples_edge_mni, ~, raw_atlas_edge_mni] = create_atlas_by_edge(all_pt_adj(1:106), pt_loc(1:106), all_abn(1:106), [1:42]', 1);
        
        % call for HUP patients
        [~, ~, samples_edge_hup, ~, raw_atlas_edge_hup] = create_atlas_by_edge(all_pt_adj(107:166), pt_loc(107:166), all_abn(107:166), [1:42]', 1);

        % loop through all patients
        for s = 1:166
            a = a+1;

            % do atlas space
            %[pt_conn, std_pt, num_samples, sem_conn, raw_atlas_pt] = create_atlas({all_pt_adj{s}}, {pt_loc{s}}, {[]}, [1:42]', 1);
            %[score_mat, corr_val, residuals] = test_patient_conn(conn_edge, std_edge, [1:42]', pt_conn, pt_loc{s});
            %bivariate_feats(feat).subj(a).data = score_mat;    

            % do nodal space
            try [native_adj_scores, corr_val] = test_native_adj(all_pt_adj{s}, pt_loc{s}, conn_edge, std_edge, [1:42]');

            bivariate_native(feat).subj(a).data = native_adj_scores;

            catch anyerror
                bivariate_native(feat).subj(a).data = NaN(length(pt_loc{s}));
            end

        end
    end
    
end

%save('bivariate_hup_mni_May13.mat','bivariate_feats','bivariate_native')
%% check bivariate features in order to assess whether they distinguish
all_normal_ch = zeros(4156,1);
all_normal_ch(normal_HUP_ch) = 1;
all_normal_ch(normal_MNI_ch) = 1;
which_feat = 5;

% epileptogenic zone in the amygdala + hippocampus or not
MTL_channels = [new_roi==17]+[new_roi==18];
SOZ_MTL_channels = MTL_channels.*soz_ch_bin;%+MTL_channels.*EIZ;
EIZ_MTL_channels = MTL_channels.*EIZ;
norm_MTL_channels = MTL_channels.*all_normal_ch;

SOZ_MTL_feat = single_bivariate_feat(find(SOZ_MTL_channels))%,which_feat);
EIZ_MTL_feat = single_bivariate_feat(find(EIZ_MTL_channels))%,which_feat);
norm_MTL_feat = single_bivariate_feat(find(norm_MTL_channels))%,which_feat);

figure(1);clf;
hold on
scatter(ones(length(norm_MTL_feat),1),norm_MTL_feat,'jitter','on')
plot([0.75 1.25],[nanmedian(norm_MTL_feat) nanmedian(norm_MTL_feat)],'k-','LineWidth',2)
scatter(2*ones(length(EIZ_MTL_feat),1),EIZ_MTL_feat,'jitter','on')
plot([1.75 2.25],[nanmedian(EIZ_MTL_feat) nanmedian(EIZ_MTL_feat)],'k-','LineWidth',2)
scatter(3*ones(length(SOZ_MTL_feat),1),SOZ_MTL_feat,'jitter','on')
plot([2.75 3.25],[nanmedian(SOZ_MTL_feat) nanmedian(SOZ_MTL_feat)],'k-','LineWidth',2)
hold off

p1 = ranksum(norm_MTL_feat,EIZ_MTL_feat)
p2 = ranksum(EIZ_MTL_feat,SOZ_MTL_feat)
p3 = ranksum(norm_MTL_feat,SOZ_MTL_feat)

%% process bivariate
for f = 1:10
    this_feat = [];
    abs_feat = [];
    for s = 1:166
        this_feat = [this_feat; prctile(bivariate_native(f).subj(s).data,50)'];
        abs_feat = [abs_feat;prctile(abs(bivariate_native(f).subj(s).data),50)'];
    end
    all_bivariate_feats(:,f) = this_feat;
    abs_bivariate_feats(:,f) = abs_feat;
end
single_bivariate_feat = max(abs_bivariate_feats')';
%single_bivariate_feat(incomplete_channels) = [];
%% 
for ch = 1:3656
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
single_bivariate_feat(isinf(single_bivariate_feat)) = nanmedian(single_bivariate_feat);

figure(1);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat,'Group',grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,grouping,'orientation','horizontal',...
     'label',{'','',''},'color',clr);
boxplot(h(3),single_bivariate_feat,grouping,'orientation','horizontal',...
     'label', {'','',''},'color',clr);
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 5])
ylim([0 5])

computeCohen_d(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
ranksum(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))
ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))

computeCohen_d(single_bivariate_feat(find(class_ind==2)),single_bivariate_feat(find(class_ind==1)))
ranksum(single_bivariate_feat(find(class_ind==2)),single_bivariate_feat(find(class_ind==1)))
computeCohen_d(single_bivariate_feat(find(class_ind==3)),single_bivariate_feat(find(class_ind==2)))
ranksum(single_bivariate_feat(find(class_ind==3)),single_bivariate_feat(find(class_ind==2)))
computeCohen_d(single_bivariate_feat(find(class_ind==3)),single_bivariate_feat(find(class_ind==1)))
ranksum(single_bivariate_feat(find(class_ind==3)),single_bivariate_feat(find(class_ind==1)))
%% combine univariate and bivariate
%all_feat_zscores = [univariate_zscores, all_bivariate_feats];
incomplete_channels = find(isnan(sum(all_feat_zscores')));
all_feat_zscores2 = all_feat_zscores;
all_feat_zscores2(incomplete_channels,:) = [];
EIZ = [zeros(1772,1);[spike_ind - spike_ind.*soz_ch]];

% 1: patient, 2: abnormal, 3: EIZ, 4: SOZ, 5: resect, 6:outcome
all_possible_labels = [all_pts,[EIZ+[zeros(1772,1);soz_ch]],EIZ,[zeros(1772,1);soz_ch],resect_ch_all,[ones(1772,1);HUP_outcome_all]];
all_possible_labels(incomplete_channels,:) = [];

multi_class = [~all_possible_labels(:,2)+2.*all_possible_labels(:,3)+3.*all_possible_labels(:,4)];

% get only non EIZ normal versus SOZ
non_normal_non_soz = [zeros(1772,1);[soz_ch==0]];
non_normal_non_soz(incomplete_channels) = [];
eliminate_channels = find(all_possible_labels(:,2)==0);%find(non_normal_non_soz==1);
all_feat_zscores3 = all_feat_zscores2;
all_feat_zscores3(eliminate_channels,:) = [];
new_labels = all_possible_labels;
new_labels(eliminate_channels,:) = [];

all_feats = all_feat_zscores2;
all_labels = all_possible_labels(:,2);%new_labels(:,4);%

[part] = make_xval_partition(length(all_labels), 10);
for p = 1:10
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
    
    xfold_acc(p) = sum(Y_pred==Y_test)./length(Y_test);
    
    Y_pred = round(pihat1(:,2));
    
    [X1,Y1,T,AUC1] = perfcurve(Y_test,pihat1(:,2),1);  
    AUC_all_1(p).X = X1;
    AUC_all_1(p).Y = Y1;
    xfold_auc1(p) = AUC1;

    [X2,Y2,T,AUC2] = perfcurve(Y_test,pihat2(:,2),1);
    AUC_all_2(p).X = X2;
    AUC_all_2(p).Y = Y2;
    xfold_auc2(p) = AUC2;

    [X3,Y3,T,AUC3] = perfcurve(Y_test,pihat3(:,2),1);
    AUC_all_3(p).X = X3;
    AUC_all_3(p).Y = Y3;
    xfold_auc3(p) = AUC3;
    
    
end

%%
mean(xfold_auc1)
mean(xfold_auc2)
mean(xfold_auc3)
signrank(xfold_auc1,xfold_auc2)
signrank(xfold_auc1,xfold_auc3)
signrank(xfold_auc2,xfold_auc3)

% random forest feature importance
imp = Mdl3.OOBPermutedPredictorDeltaError;

figure(1);clf;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xticks([1:20])

figure(2);clf;
hold on
plot(AUC_all_1(1).X,AUC_all_1(1).Y,'LineWidth',2)%4
plot(AUC_all_2(4).X,AUC_all_2(4).Y,'LineWidth',2)%7
plot(AUC_all_3(3).X,AUC_all_3(3).Y,'LineWidth',2)%10
legend('Univariate  - AUC: 0.78','Bivariate  - AUC: 0.81','All - AUC: 0.82','Location','SouthEast')%86/85/91
xlabel('False positive rate')
ylabel('True positive rate')
title('All abnormal vs. normal channel classification')

%% Examine correlations between z scores
[r,p] = corr(all_feat_zscores2)
figure(1);clf;
imagesc(r)
colorbar

%% lesional vs non-lesional
EIZ = [zeros(1772,1);[spike_ind - spike_ind.*soz_ch]];
all_labels = [all_pts,[EIZ+[zeros(1772,1);soz_ch]],EIZ,[zeros(1772,1);soz_ch],resect_ch_all,[ones(1772,1);HUP_outcome_all]];

%% get patient-level SOZ, RZ, spike indicators
all_normal_good = [];
all_normal_poor = [];
all_les_in_in = [];
all_nonles_in_in = [];
all_SOZ_in_in = [];
all_SOZ_in_out = [];
all_spike_in_in = [];
all_spike_in_out = [];

b = 0;
ptile = 75;
for s = 107:166
    b = b+1
    a = s;
    
    pt_roi = pt_loc{s};
    res = resected_ch([patient_no==b]); 
    spike = spike_24h([patient_no==b]);
    soz = soz_ch([patient_no==b]);
    
    roi_res = NaN*zeros(42,1);
    
%     for r = 1:42
%         this_roi = (r);
%         roi_res(r) = mean(res(pt_roi==this_roi));
%         roi_soz(r) = mean(soz(pt_roi==this_roi));
%         roi_spike(r) = mean(spike(pt_roi==this_roi));
%     end
%     
%     HUP_scores(a).res = ceil(roi_res);
%     HUP_scores(a).spike = (roi_spike>24);
%     HUP_scores(a).soz = ceil(roi_soz);
%     
%     res_roi = [HUP_scores(a).res==1];
%     soz_roi = ([HUP_scores(a).soz==1]);
%     spk_roi = ([HUP_scores(a).spike==1]);
%     
%     non_res_roi = ([HUP_scores(a).res==0]);
%     non_soz_roi = ([HUP_scores(a).soz==0]);
%     non_spk_roi = ([HUP_scores(a).spike==0]);
%     
    
    % spared edges
    all_normal = bivariate_native(feat).subj(a).data(find(~res),find(~res));
    HUP_normal(b,1) = nanmean(all_normal(:));
    
    all_normal_good = [all_normal_good;prctile((all_normal(:)),ptile)];
    
    % resected edges 
    %in_out_resect = bivariate_native(feat).subj(a).data(find(res),:);
    in_in_resect = bivariate_native(feat).subj(a).data(find(res),find(~res));
    %HUP_resect(b,1) = nanmedian(in_out_resect(:)); % (in-out)
    %HUP_resect(b,2) = nanmedian(in_in_resect(:)); % in-in
    %HUP_resect(b,3) = nanmedian(test_HUP_scores(:)); % all
    if lesion_field(b)
    all_les_in_in = [all_les_in_in;prctile(in_in_resect(:),ptile)];
    else
    all_nonles_in_in = [all_nonles_in_in;prctile(in_in_resect(:),ptile)];
    end
    %all_RZ_in_out = [all_RZ_in_out;prctile((in_out_resect(:)),ptile)];
    
    % soz edges (in-out)
    in_out_soz = bivariate_native(feat).subj(a).data(find(soz),find(~soz));
    in_in_soz = bivariate_native(feat).subj(a).data(find(soz),find(soz));
    %HUP_soz(b,1) = nanmedian(in_out_soz(:)); % (in-out)
    %HUP_soz(b,2) = nanmedian(in_in_soz(:)); % in-in
    
    all_SOZ_in_in = [all_SOZ_in_in;prctile(in_in_soz(:),ptile)];
    all_SOZ_in_out = [all_SOZ_in_out;prctile((in_out_soz(:)),ptile)];
    
    % spike edges 
    try in_out_spike = bivariate_native(feat).subj(a).data(find(spike),find(~spike));
        all_spike_in_out = [all_spike_in_out;prctile((in_out_spike(:)),ptile)];
    catch anyerror
    end
    try in_in_spike = bivariate_native(feat).subj(a).data(find(spike),find(spike));
        all_spike_in_in = [all_spike_in_in;prctile(in_in_spike(:),ptile)];
    catch anyerror
    end
%     HUP_spike(a,1) = nanmedian(in_out_spike(:)); % (in-out)
%     HUP_spike(a,2) = nanmedian(in_in_spike(:)); % in-in
%     
    
     
%     
end

%
figure(1);clf;
hold on

scatter(ones(length(all_normal_good),1),all_normal_good,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[nanmedian(all_normal_good) nanmedian(all_normal_good)],'k-','LineWidth',2)
scatter(3*ones(length(all_SOZ_in_in),1),all_SOZ_in_in,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([2.75 3.25],[nanmedian(all_SOZ_in_in) nanmedian(all_SOZ_in_in)],'k-','LineWidth',2)
scatter(4*ones(length(all_SOZ_in_out),1),all_SOZ_in_out,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([3.75 4.25],[nanmedian(all_SOZ_in_out) nanmedian(all_SOZ_in_out)],'k-','LineWidth',2)
scatter(5*ones(length(all_spike_in_in),1),all_spike_in_in,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
plot([4.75 5.25],[nanmedian(all_spike_in_in) nanmedian(all_spike_in_in)],'k-','LineWidth',2)
scatter(6*ones(length(all_spike_in_out),1),all_spike_in_out,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
plot([5.75 6.25],[nanmedian(all_spike_in_out) nanmedian(all_spike_in_out)],'k-','LineWidth',2)
xlim([0.5, 6.5])
ylim([0 6])
hold off

ranksum(all_normal_good,all_SOZ_in_in)
ranksum(all_normal_good,all_SOZ_in_out)
ranksum(all_normal_good,all_spike_in_in)
ranksum(all_normal_good,all_spike_in_out)

ranksum(all_SOZ_in_in,all_SOZ_in_out)
ranksum(all_SOZ_in_in,all_spike_in_in)
ranksum(all_spike_in_in,all_spike_in_out)
ranksum(all_SOZ_in_out,all_spike_in_out)
%% boxplots baby

num_good = sum([HUP_outcome==1]);
num_poor = sum([HUP_outcome>1]);

good_out=HUP_normal(HUP_outcome==1);
poor_out=HUP_normal(HUP_outcome>1);

good_in_out_soz=HUP_soz([HUP_outcome==1],1);
good_in_in_soz=HUP_soz([HUP_outcome==1],2);
poor_in_out_soz=HUP_soz([HUP_outcome>1],1);
poor_in_in_soz=HUP_soz([HUP_outcome>1],2);

good_in_out_spike=HUP_spike([HUP_outcome==1],1);
good_in_in_spike=HUP_spike([HUP_outcome==1],2);
poor_in_out_spike=HUP_spike([HUP_outcome>1],1);
poor_in_in_spike=HUP_spike([HUP_outcome>1],2);

all_plot_data = NaN*zeros(57,6);
all_plot_data(1:num_good,1) = good_in_in_soz;
all_plot_data(1:num_poor,2) = poor_in_in_soz;
all_plot_data(1:num_good,3) = good_in_out_soz;
all_plot_data(1:num_poor,4) = poor_in_out_soz;
all_plot_data(1:num_good,5) = good_out;
all_plot_data(1:num_poor,6) = poor_out;

ranksum(good_in_in_soz,poor_in_in_soz)
ranksum(good_in_out_soz,poor_in_out_soz)
ranksum(good_out,poor_out)

signrank(good_in_in_soz,good_in_out_soz)
signrank(good_in_out_soz,good_out)
signrank(good_in_in_soz,good_out)

signrank(poor_in_in_soz,poor_in_out_soz)
signrank(poor_in_out_soz,poor_out)
signrank(poor_in_in_soz,poor_out)

ranksum(good_in_in_spike,poor_in_in_spike)
ranksum(good_in_out_spike,poor_in_out_spike)
ranksum(good_out,poor_out)

signrank(good_in_in_spike,good_in_out_spike)
signrank(good_in_out_spike,good_out)
signrank(good_in_in_spike,good_out)

signrank(poor_in_in_spike,poor_in_out_spike)
signrank(poor_in_out_spike,poor_out)
signrank(poor_in_in_spike,poor_out)

figure(1);clf;
boxplot(all_plot_data)
ylim([-2 2])

%% lesional vs non-lesional
r1 = HUP_soz((lesion_field==0),2);
r2 = HUP_soz(lesion_field,2);

r3 = HUP_normal((lesion_field==0));
r4 = HUP_normal(lesion_field);

ranksum(r1,r2)

figure(1);clf;
hold on
plot(ones(length(r1),1),r1,'ko')
plot(2*ones(length(r2),1),r2,'ko')
xlim([0.5 2.5])
hold off


