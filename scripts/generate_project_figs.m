%% generate project figures
% Set up workspace
clear all

% load pre-processed z scores
load('univariate_zscores.mat')
load('bivariate_zscores.mat')

% set up path - change iEEG_atlas_path to where you download the repository
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas_dev';

% other paths may stay the same
metadata = readtable('data/atlas_metadata_simplified.xlsx');
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
    late_outcome(s) = floor(engel_scores(end)); % late outcome is latest time point up to 24 months. These are used in the project.
end

% extract additional fields of metadata
therapy_field = metadata{:,6};
implant_field = metadata{:,7};
target_field = metadata{:,8};
laterality_field = metadata{:,9};
lesion_field = metadata{:,10};
gender_field = metadata{:,13};
lesion_field = strcmp(lesion_field,'Lesional'); % 1 for lesional
therapy_field = strcmp(therapy_field,'Ablation'); % 1 for ablation

load('data/MNI_atlas_new.mat')
load('data/HUP_atlas_final.mat')

all_pts = [Patient; (patient_no+106)]; % combine patient numbers for MNI & HUP

HUP_outcome_all = zeros(length(soz_ch),1);
for i = 1:60
    HUP_outcome_all(find(patient_no==i)) = late_outcome(i); % we use outcome as latest measured outcome
end

% clear unused data from MNI atlas
clear NodesLeft NodesLeftInflated NodesRightInflated NodesRight NodesRegionLeft NodesRegionRight Data_N2 Data_N3 Data_R

% some colors for plotting, etc (r/g/b)
color1 = [0, 0.4470, 0.7410]; color2 = [0.6350, 0.0780, 0.1840]; color3 = [255, 248, 209]./255; color4 = [103 55 155]/255;
color6 = [78 172 91]/255; color7 = [0.9290, 0.6940, 0.1250]; color8 = [0.8500, 0.3250, 0.0980];

% assign wake clips from MNI and HUP into a final data matrix
all_wake_data = [Data_W,wake_clip];
    
% perform localization based on AAL116 atlas (which has an internal WM mask also)
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

% define normal MNI/HUP and abnormal channels (SOZ/high-spike/RZ)
% set threshold on spike counts
spike_thresh = 24; % this is empirical, 1 spike/hour

% find indices for which spikes are greater than threshold
spike_ind = [spike_24h>spike_thresh];

% define all abnormal channels
abnormal_ch = find([spike_ind+soz_ch+(HUP_outcome_all>1)]>0)+1772;

% define all seizure onset indices
soz_ch_inds = find(soz_ch)+1772;

% define normal HUP channels
normal_HUP_ch = find([spike_ind+soz_ch+(HUP_outcome_all>1)]==0)+1772;

% define normal MNI channels
normal_MNI_ch = [1:1772]';

% define all normal channels
all_normal_ch = [normal_MNI_ch;normal_HUP_ch];

% define exclusive irritative zone
EIZ_ch = find(spike_ind - spike_ind.*soz_ch)+1772;

% binary variables for whether channel is in exclusive irritative zone, resected, or
% seizure onset zone
EIZ_binary = [zeros(1772,1);[spike_ind - spike_ind.*soz_ch]];
resect_binary = [zeros(1772,1);resected_ch];
soz_binary = [zeros(1772,1);soz_ch];

% assign atlas locations per patient into cell structure
for s = 1:max(all_pts)
    pt_loc{s} = new_roi(all_pts==s)';
end

% indices of all abnormal channels
all_abnormal_ch = zeros(5203,1);
all_abnormal_ch(abnormal_ch) = 1;

% make structure of patient-level clinically abnormal channels, useful for
% atlas construction
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

%% *** Figure 4 and Supplemental Figure 3 ***
% calculate power spectrum
all_wake_data = [Data_W,wake_clip];

[pxx,f] = pwelch(all_wake_data,400,200,[0.5:0.5:80],200); % original was 200/100
pxx_norm = pxx./sum(pxx);

% Analyze power spectral density between HUP/MNI, Normal/irritative zone/SOZ 
% change this to be on a per-patient level rather than a per-node level

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
    roi_EIZ = intersect([EIZ_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_soz = intersect([soz_ch_inds],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    
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

% process bivariate edge features into nodal features

for f = 1:5
    this_feat = [];
    abs_feat = [];
    for s = 1:166
        % take 75th percentile across all edges of each node
        abs_feat = [abs_feat;prctile(abs(bivariate_native(f).subj(s).data),75)']; 
    end
    
    % assign into feature matrix
    abs_bivariate_feats(:,f) = abs_feat;
end

% for the single bivariate feature we want the maximum absolute Z score
% across all 10 bivariate features
[single_bivariate_feat] = median(abs_bivariate_feats')';

% assign previously calculated univariate z scores into bivariate scores
all_feat_zscores = [abs(univariate_zscores), abs_bivariate_feats];
all_feat_zscores(isinf(all_feat_zscores)) = 0; % zero out any z-scores that were infinite due to zero variance in edge estimates

% incomplete channels are ones which are unlocalized in new atlas
incomplete_channels = find(isnan(sum(all_feat_zscores')));

% remove incomplete channels from univariate features
univariate_feats = all_feat_zscores(:,1:5);
univariate_feats(incomplete_channels,:) = [];


%% set up machine learning features and labels
% define channels in which sparse atlas coverage yields NaN's 
incomplete_channels = find(isnan(sum(all_feat_zscores')));
complete_channels = find(~isnan(sum(all_feat_zscores')));

% start with channel indices for ML features as all channels
ch_indices = [1:5203]';

% remove channels for which features are incomplete from further analysis
zscores_complete = all_feat_zscores;
zscores_complete(incomplete_channels,:) = []; % remove channels that are incomplete
ch_indices(incomplete_channels) = []; % remove channel indices that are incomplete

% define binary variable for which channels are abnormal
abn_ch_binary = zeros(5203,1);
abn_ch_binary(abnormal_ch) = 1;

% define a matrix containing different labels. The binary variables
% correspond to the following labels
% column 1: patient number
% column 2: whether the channel is normal or abnormal
% column 3: the exclusive irritative zone
% column 4: the seizure onset zone
% column 5: whether or not the channel was resected
% column 6: surgical outcome
% column 7: whether the channel was both resected and seizure onset zone
all_class_labels = [all_pts,abn_ch_binary,EIZ_binary,soz_binary,resect_binary,[ones(1772,1);HUP_outcome_all],[[zeros(1772,1);soz_ch].*resect_binary]];
all_class_labels(incomplete_channels,:) = [];

% we want only resected, soz channels vs normal channels
normal_channels = find(all_class_labels(:,2)==0);
resected_soz_channels = find(all_class_labels(:,7)==1);
abn_ch_list = find(all_class_labels(:,2)==1);

% extract features for these channels from all features and define labels
ML_features = zscores_complete([normal_channels;resected_soz_channels],:);
ML_chs = ch_indices([normal_channels;resected_soz_channels]); % original channel indices for data we are using
ML_labels = all_class_labels([normal_channels;resected_soz_channels],7); % 0 for normal 1 for resected-SOZ
ML_pt_inds = all_class_labels([normal_channels;resected_soz_channels],1); % patient indices for data we are using

%% *** Figure 6A and 6B ***
% random forest and feature importance

clear part all_pred_1 all_pred_2 all_pred_3

num_parts = 10;

[pt_part] = make_xval_partition(166, num_parts);

% cross-validate based on patients 
for p = 1:num_parts
    indices = find(pt_part==p)
    for i = 1:length(indices)
    part(find(indices(i)==ML_pt_inds(:,1))) = p;
    
    part2(find(indices(i)==all_class_labels(:,1))) = p;
    end
end

% Define classes, 0 = non seizure, 1 = seizure
ClassNames = [0,1];

% Define costs -> penalty here is 5
cost.ClassNames = ClassNames;
cost.ClassificationCosts = [0 1; 1 0];

for p = 1:num_parts
    p
    X_train = ML_features(find(part~=p),:);
    Y_train = ML_labels(find(part~=p));
    X_test = ML_features(find(part==p),:);
    Y_test = ML_labels(find(part==p));
    
    X_analysis = zscores_complete(find(part2==p),:);
    
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

    [X1,Y1,T,AUC1] = perfcurve(ML_labels,all_pred_1,1);  
    [X2,Y2,T,AUC2] = perfcurve(ML_labels,all_pred_2,1);
    [X3,Y3,T,AUC3] = perfcurve(ML_labels,all_pred_3,1);

    AUC1
    AUC2
    AUC3
    
    Mdl3 = TreeBagger(100,ML_features,ML_labels,'Method',...
    'classification','PredictorSelection','curvature','OOBPredictorImportance','on','cost',cost);

% random forest feature importance
imp = Mdl3.OOBPermutedPredictorDeltaError;

figure(3);clf;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xticks([1:10])

figure(4);clf;
hold on
plot(X1,Y1,'LineWidth',2,'Color',color4)
plot(X2,Y2,'LineWidth',2,'Color',color7)
plot(X3,Y3,'LineWidth',2,'Color',color1)%
legend(sprintf('Univariate  - AUC: %.2f',AUC1),sprintf('Bivariate  - AUC: %.2f',AUC2),sprintf('All - AUC: %.2f',AUC3),'Location','SouthEast')%86/85/91
xlabel('False positive rate')
ylabel('True positive rate')
title('All abnormal vs. normal channel classification')

%% *** Figure 5A ***
% mean Z-scores for Mesial Temporal Lobe (regions 17-20) in SOZ vs IZ vs uninvolved at nodal level
new_roi_2 = new_roi;
new_roi_2(incomplete_channels) = [];

single_bivariate_feat_2 = single_bivariate_feat;
single_bivariate_feat_2(incomplete_channels) = [];

MTL_uninvolved = [];
MTL_EIZ = [];
MTL_SOZ = [];

for pt = 1:166
    uninv_pt = find([all_class_labels(:,1)==pt].*[all_class_labels(:,2)==0].*[[new_roi_2==17]+[new_roi_2==19]+[new_roi_2==18]+[new_roi_2==20]]);
    EIZ_pt = find([all_class_labels(:,1)==pt].*[all_class_labels(:,3)==1].*[[new_roi_2==17]+[new_roi_2==19]+[new_roi_2==18]+[new_roi_2==20]]);
    SOZ_pt = find([all_class_labels(:,1)==pt].*[all_class_labels(:,4)==1].*[[new_roi_2==17]+[new_roi_2==19]+[new_roi_2==18]+[new_roi_2==20]]);
    
    MTL_uninvolved = [MTL_uninvolved; mean(single_bivariate_feat_2(uninv_pt))];
    MTL_EIZ = [MTL_EIZ; mean(single_bivariate_feat_2(EIZ_pt))];
    MTL_SOZ = [MTL_SOZ; mean(single_bivariate_feat_2(SOZ_pt))];
end

mean_uninv = MTL_uninvolved;
mean_EIZ = MTL_EIZ;
mean_SOZ = MTL_SOZ;

plot_matrix = padcat(mean_uninv,mean_EIZ,mean_SOZ);
figure(5);clf;
UnivarScatter(plot_matrix)
xticks([1:3])
xticklabels({'uninvolved','EIZ','SOZ'})

p1 = ranksum(mean_uninv,mean_EIZ,'tail','left')
cd1 = computeCohen_d(mean_uninv,mean_EIZ)
p2 = ranksum(mean_uninv,mean_SOZ,'tail','left')
cd2 = computeCohen_d(mean_uninv,mean_SOZ)
p3 = ranksum(mean_EIZ,mean_SOZ,'tail','left')
cd3 = computeCohen_d(mean_SOZ,mean_EIZ)

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

%% *** Figure 7 ***
% render abnormalities
soz_ch_all = [zeros(1772,1);soz_ch];
this_pt = 116;%116/122/154 / 154 %114 res - 45/50 EIZ - 52/64, non-inv - 20,30 % 120 - res - 2., EIZ - 23./64, non-inv - 33/50.
this_pt_abn = all_pred(find(all_pts==this_pt));
this_pt_soz = soz_ch_all(find(all_pts==this_pt));
this_pt_eiz = EIZ_binary(find(all_pts==this_pt));
this_pt_res = resect_binary(find(all_pts==this_pt));
this_pt_coords = all_coords(find(all_pts==this_pt),:);
this_pt_eeg = all_wake_data(:,find(all_pts==this_pt));
%this_pt_delta_coh_abn(isnan(this_pt_delta_coh_abn)) = 0;

this_pt_multiclass = 2*this_pt_soz+this_pt_eiz;

this_pt_zscores = all_feat_zscores(find(all_pts==this_pt),:);

this_pt_test = zeros(size(this_pt_res));
this_pt_test([28,9,49])=1;

figure(6);clf;
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
figure(7);
subplot(3,1,1)
plot(this_pt_eeg(6001:12000,4))
subplot(3,1,2)
plot(this_pt_eeg(6001:12000,23))
subplot(3,1,3)
plot(this_pt_eeg(6001:12000,27))


figure(8);clf;
final_elec_matrix = [this_pt_coords,this_pt_abn,ones(size(this_pt_res))];
%final_elec_matrix = final_elec_matrix([2,23,50],:)
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','abnormality_render.mat')
%delete render_elecs.node

%% *** Figure 5C ***
% patient-specific correlation between abnormality features
for s = 107:166
    pt_zscores = zscores_complete(find([all_class_labels(:,1)==s]),:);
    
    pt_zscores_univar = pt_zscores(:,1:5);
    pt_zscores_bivar = pt_zscores(:,6:10);
    
    corrval2(s-106) = corr(pt_zscores_univar(:),pt_zscores_bivar(:),'Type','Pearson');
    
    corrval_outcome((s-106),1) = late_outcome(s-106);
end

p6 = ranksum(corrval2(corrval_outcome==1),corrval2(corrval_outcome>1),'tail','right')
d6 = computeCohen_d(corrval2(corrval_outcome==1),corrval2(corrval_outcome>1))

plot_matrix = padcat(corrval2(corrval_outcome==1)',corrval2(corrval_outcome>1)');
figure(9);clf;
UnivarScatter(plot_matrix);
xticks([1,2])
xticklabels({'Engel 1','Engel 2+'})
ylim([-0.3 0.7])


%% *** Supplemental Figure 8 ***
% Resected vs non-resected epileptogenicity
a = 0;
for s = [107:166]
    a = a+1;
    pred_pt = all_pred3(all_pts==s);
    mean_in(a,1) = mean(pred_pt(find(resect_binary(all_pts==s))));
    mean_out(a,1) = mean(pred_pt(find(~resect_binary(all_pts==s))));
    res_test(a,1) = mean(pred_pt(find(resect_binary(all_pts==s))))>mean(pred_pt(find(~resect_binary(all_pts==s))))
    
end

mean_diff = mean_in-mean_out;

figure(10);clf;
hold on
plot_matrix = padcat(mean_diff(late_outcome==1),mean_diff(late_outcome>1));
UnivarScatter(plot_matrix);
xticks([1:2])
xticklabels({'Engel 1','Engel 2+'})
plot([0.5 2.5],[0 0],'k-.')
ylabel('Mean resected - mean non-resected epileptogenicity')



%% *** Figure 6D ***
% compute AUPRC/DRS at a per-patient level
a = 0;
for s = 107:166
    
    a = a+1;
    this_pt_label = all_class_labels(find(all_class_labels(:,1)==s),4);
    this_pt_pred3 = all_pred_analysis3(find(all_class_labels(:,1)==s));
    try [X3,Y3,T,AUC] = perfcurve(this_pt_label,this_pt_pred3,1,'XCrit','TPR','YCrit','PPV');
    catch
        AUC = NaN;
    end
    
    auc_pt_outcome(a,1) = late_outcome(s-106);
    pt_auc(a,1) = [AUC];
end

good_inds = find(auc_pt_outcome==1);
poor_inds = find(auc_pt_outcome>1);

plot_matrix = padcat(pt_auc(good_inds),pt_auc(poor_inds));
figure(11);clf;
UnivarScatter(plot_matrix)

ranksum(plot_matrix(:,1),plot_matrix(:,2))
computeCohen_d(plot_matrix(:,1),plot_matrix(:,2))
xticks([1,2])
xticklabels({'Engel 1','Engel 2+'})
ylabel('Area under precision-recall curve')

%% *** Figure 6C ***
% patient-level abnormality score all regions

pt_uninvolved = [];
pt_eiz = [];
pt_soz = [];
    
for s = 107:166
    
    pt_uninvolved_inds = find([all_class_labels(:,2)==0].*[all_class_labels(:,1)==s]);
    pt_eiz_inds = find([all_class_labels(:,3)==1].*[all_class_labels(:,1)==s]);
    pt_soz_inds = find([all_class_labels(:,4)==1].*[all_class_labels(:,1)==s]);
    
    pt_uninvolved = [pt_uninvolved; nanmedian(all_pred_analysis(pt_uninvolved_inds))];
    pt_eiz = [pt_eiz; nanmedian(all_pred_analysis(pt_eiz_inds))];
    pt_soz = [pt_soz; nanmedian(all_pred_analysis(pt_soz_inds))];
end

plot_matrix = padcat(pt_uninvolved,pt_eiz,pt_soz);

figure(12);clf;
UnivarScatter(plot_matrix)
xticks([1:3])
xticklabels({'uninvolved','irritative zone','seizure onset zone'})
ylabel('predicted epileptogenicity score')
ylim([0 0.8])

p1 = ranksum(pt_uninvolved,pt_eiz)
d1 = computeCohen_d(pt_uninvolved,pt_eiz)
p2 = ranksum(pt_eiz,pt_soz)
d2 = computeCohen_d(pt_eiz,pt_soz)
p3 = ranksum(pt_uninvolved,pt_soz)
d3 = computeCohen_d(pt_uninvolved,pt_soz)

%% *** Figure 5B ***
% correlation length of diagnosis
pt_in_rz = [];
pt_out_rz = [];

for s = 107:166
    
    w_inds = find([all_class_labels(:,1)==s]);
    o_inds = find([all_class_labels(:,1)==s]);
    
    length_diagnosis = age_surgery - age_onset;
    
    pt_in_rz = [pt_in_rz; nanmedian(single_bivariate_feat(w_inds))];
    pt_out_rz = [pt_out_rz; nanmedian(single_bivariate_feat(o_inds))];
 

end

figure(13);clf;
hold on
plot(length_diagnosis, pt_out_rz,'b.','Color',[106 178 180]/255,'MarkerSize',12)
p = polyfit(length_diagnosis,pt_out_rz,1);
f = polyval(p,length_diagnosis);
plot(length_diagnosis,f,'b-','LineWidth',2)
[R,P] = corr(length_diagnosis, pt_out_rz,'type','Pearson','tail','right')