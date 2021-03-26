function [all_patients, all_inds, all_locs, coords_field, ...
    hasData_field, hasVar_field, id_field, implant_field, outcome_field, target_field,...
    therapy_field, region_list, region_names, ...
    lesion_field] = set_up_workspace_dev(iEEG_atlas_path)

addpath(genpath(iEEG_atlas_path))

% suppress warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% load in data from excel spreadsheet
metadata = readtable("data/atlas_project_metadata.xlsx");

% enable warnings
warning ('on','all')

% place patients in a struct, extracting all relevant metadata
all_patients = struct('patientID',metadata.Patient, ...
'outcome', metadata.Engel_6_mo,'conn',cell(length(metadata.Patient),1), ...
'var',cell(length(metadata.Patient),1), ...
'coords',cell(length(metadata.Patient),1), ...
'roi',cell(length(metadata.Patient),1), ...
'resect',cell(length(metadata.Patient),1), ...
'hasData',cell(length(metadata.Patient),1),...
'therapy',metadata.Therapy,'implant',metadata.Implant,...
'target',metadata.Target,'laterality',metadata.Laterality,...
'lesion_status',metadata.Lesion_status,'age_onset',metadata.Age_onset,...
'age_surgery',metadata.Age_surgery,'gender',metadata.Gender);

% Extract atlas indices and ROIs available from atlas (here AAL116 w/WM)
fileID = fopen('localization/AAL116_WM.txt');
atlas_info = textscan(fileID,'%s %s %d');
all_inds = [double(atlas_info{3})];
all_locs = [atlas_info{2}];

% set up arrays to store data
id_field = {all_patients.patientID};
conn_field = {all_patients.conn};
var_field = {all_patients.var};
coords_field = {all_patients.coords};
roi_field = {all_patients.coords};
in_brain_field = {all_patients.coords};
resect_field = {all_patients.resect};
outcome_field = {all_patients.outcome};
hasData_field = {all_patients.hasData};
therapy_field = {all_patients.therapy};
implant_field = {all_patients.implant};
target_field = {all_patients.target};
lesion_field = {all_patients.lesion_status};
age_onset_field = {all_patients.age_onset};
age_surgery_field = {all_patients.age_surgery};

% if true, the script will automatically move problematic data to another
% directory
move_files = false;

% load in data from all patients
for k = 1:length(metadata.Patient)
    folderpath = sprintf('data/%s',id_field{k});
    datapath = sprintf('%s/patient_data.mat',folderpath);
    if isfile(datapath)
        fprintf('%s: ',datapath)
        d = load(datapath);
        conn_field{k} = d.II_conn;
        try var_field{k} = d.II_var;
            hasVar_field(k) = 1;
        catch ME
            var_field{k} = [];
            hasVar_field(k) = 0;
        end
        if sum(sum(~isnan(d.II_conn(1).data))) == 0
            hasData_field{k} = false;
            fprintf('(connectivity data is all NaNs!)\n')
            if move_files, movefile(folderpath,'data/exclude/no_conn_data'); end
            continue
        else
            hasData_field{k} = true;
        end
        coords_field{k} = d.mni_coords;
        
        try
            resect_field{k} = d.res_elec_inds;
        catch Error
            resect_field{k} = [];
        end
        % convert all electrode coordinates to region names
        try
            [~,electrode_regions,~] = nifti_values(coords_field{1,k},'localization/AAL116_WM.nii');
            roi_field{k} = electrode_regions;
            in_brain_field{k} = d.in_brain;
            fprintf('loaded\n')
        catch ME
            fprintf('failed to load\n')
            warning('Problem converting MNI coordinates to region labels\n(%s)',datapath, ME.identifier)
            hasData_field{k} = false;
            if move_files, movefile(folderpath,'data/exclude/out_of_bound_electrodes'); end
        end
    else
        hasData_field{k} = false;
    end
end

% place data back into main struct
[all_patients.conn] = conn_field{:};
[all_patients.var] = var_field{:};
[all_patients.coords] = coords_field{:};
[all_patients.roi] = roi_field{:};
[all_patients.resect] = resect_field{:};
[all_patients.hasData] = hasData_field{:};
[all_patients.lesion_status] = lesion_field{:};
[all_patients.age_onset] = age_onset_field{:};
[all_patients.age_surgery] = age_surgery_field{:};
[all_patients.in_brain] = in_brain_field{:};

fprintf('\nAll patient data loaded.')

% load in region numbers
region_list = zeros(1,90); % for the 90 AAL regions we will be using
region_names = cell(1,90);
fi = fopen("localization/AAL116_WM.txt");
for j = 1:90
    label = split(fgetl(fi));
    region_list(j) = str2double(label{3});
    region_names{j} = label{2};
end

fprintf('\nRegion list loaded.\n')

end