%% Set up workspace
clear all

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

%% generate test patient Z-scores
% univariate: activity
demo_pt = 110;

all_f = {'delta','theta','alpha','beta','gamma'}; % frequency bands
all_m = {'bandpower'}; % metrics

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
        [zscores, pt_zscores, pt_ind, roi_ind] = test_univariate_atlas(median_vec, std_vec, feat_vals, all_pts, [1:40], new_roi, 'patient');
        
        % store the z scores
        demo_univariate_zscores(:,c) = zscores(all_pts==demo_pt);
        
    end
end

% bivariate: connectivity
conn_band = {'delta','theta','alpha','beta','gamma'};
conn_type = {'coh'};

% which feature
feat = 0;

% loop through connectivity and frequency
for c = 1
    for f = 1:5
        
        feat = feat+1;
        a = 0
        
        load(sprintf('data/adj_matrices/all_adj_%s_%s_May25.mat',conn_type{c},conn_band{f}));
        
        [conn_edge, std_edge, samples_edge, sem_edge, raw_atlas_edge] = create_atlas_by_edge(all_pt_adj, pt_loc, all_abn, [1:40]', 2);
        
        % call for MNI patients
        [conn_edge_mni, std_edge_mni, samples_edge_mni, ~, raw_atlas_edge_mni] = create_atlas_by_edge(all_pt_adj(1:106), pt_loc(1:106), all_abn(1:106), [1:40]', 1);
        
        % call for HUP patients
        [conn_edge_hup, std_edge_hup, samples_edge_hup, ~, raw_atlas_edge_hup] = create_atlas_by_edge(all_pt_adj(107:166), pt_loc(107:166), all_abn(107:166), [1:40]', 1);

        this_feat_conn(feat).mni = conn_edge_mni;
        this_feat_conn(feat).hup = conn_edge_hup;
        
        % loop through all patients
        for s = 110 % change this back
            a = a+1;

            % do nodal space
            try [native_adj_scores, corr_val] = test_native_adj(all_pt_adj{s}, pt_loc{s}, conn_edge, std_edge, [1:40]');

            demo_bivariate_native(feat).data = native_adj_scores;

            catch anyerror
                demo_bivariate_native(feat).data = NaN(length(pt_loc{s}));
            end

        end
    end
    
end

%% process the raw z scores
for f = 1:5
   
    pt_abs_feat = prctile(abs(demo_bivariate_native(f).data),75)'; 
    
    % assign into feature matrix
    demo_abs_bivariate_feats(:,f) = pt_abs_feat;
end

% assign previously calculated univariate z scores into bivariate scores
demo_all_feat_zscores = [abs(demo_univariate_zscores), demo_abs_bivariate_feats];
demo_all_feat_zscores(isinf(all_feat_zscores)) = 0; % zero out any z-scores that were infinite due to zero variance in edge estimates

%% Predict epileptogenicity from random forest
%load RF_model
[Yfit,pihat3] = predict(Mdl3,demo_all_feat_zscores);
this_pt_abn = pihat3(:,2);

%% render results of predicted epileptogenicity
final_elec_matrix = [all_coords(all_pts==110,:),this_pt_abn,ones(size(this_pt_abn))];
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','abnormality_render.mat')

%% compare to resection and irritative zone (this patient is good outcome)
this_pt = 110; % HUP068
this_pt_res = resect_binary(find(all_pts==this_pt));
final_elec_matrix = [all_coords(all_pts==110,:),this_pt_res,ones(size(this_pt_res))];
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','abnormality_render.mat')

