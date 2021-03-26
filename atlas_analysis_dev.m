%% set up workspace
clear all

iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas_dev';
metadata = readtable("data/atlas_metadata_nishant.xlsx");

all_engel_scores = metadata{:,3:5}; % extracted 3 columns of engel scores (1.2 = 1B, 2.1 = 2A, etc)

[all_patients_raw, all_inds, all_locs, coords_field, hasData_field, hasVar_field, id_field,...
    implant_field, outcome_field, target_field,...
    therapy_field, region_list, region_name, lesion_field] = set_up_workspace_dev(iEEG_atlas_path);

% extract outcome
for s = 1:length(all_patients_raw)
    engel_scores = metadata{s,3:5};
    engel_scores = rmmissing(engel_scores);
    early_outcome(s) = floor(engel_scores(1));
    late_outcome(s) = floor(engel_scores(end));
end

remove_patients = ~[all_patients_raw.hasData];

age_onset = metadata{:,11};
age_surgery = metadata{:,12};

all_patients_raw(remove_patients) = [];
early_outcome(remove_patients) = [];
late_outcome(remove_patients) = [];
age_onset(remove_patients) = [];
age_surgery(remove_patients) = [];
therapy_field(remove_patients) = [];
lesion_field(remove_patients) = [];
target_field(remove_patients) = [];

all_engel_scores(remove_patients,:) = [];

metadata(remove_patients,:) = [];

%% get basic atlas info
test_band = 1;
test_threshold = 1;

data_patient_indices = find([all_patients_raw.hasData]);

atlas_patients = all_patients_raw(data_patient_indices);
[mean_conn, std_conn, num_samples, sem_conn] = create_atlas({atlas_patients.conn}, {atlas_patients.roi}, {atlas_patients.resect}, region_list, test_band, test_threshold);

