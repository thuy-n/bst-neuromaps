% compute_parcellation_brainstorm.m
%
% Parcellate surface and volume annotations using the Desikan-Killiany atlas
% and save results in two .csv files:
%
%     - 'parcellation_dk_surface_brainstorm.csv'
%     - 'parcellation_dk_volume_brainstorm.csv'
%
% Surface annotations are in the Brainstorm MNI152 template 15k (both hemispheres)
% Volume annotations are in the MNI152 space
%
% Surface annotations are fetched with the process 'process_nmp_fetch_maps'
% Volume annotations have to be previously fetched with 'preprocess.py', and
% are located at ../../tmp/volume.
%
% Created on 2024

%% Create protocol
ProtocolName = 'bst-neuromaps-parcellation';
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
% Parcellation (for surface and volume)
atlasName  = 'Desikan-Killiany';
% Parcels name order in FreeSurfer (hemisphere less)
% without 'unknown' and without 'corpuscallosum'
dk_fs_labels = {'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus',    ...
                'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal',          ...
                'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', 'lingual', ...
                'medialorbitofrontal', 'middletemporal', 'parahippocampal', 'paracentral', ...
                'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine',   ...
                'postcentral', 'posteriorcingulate', 'precentral', 'precuneus',            ...
                'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal',     ...
                'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole',    ...
                'temporalpole', 'transversetemporal', 'insula'};
% Labels for each hemisphere
dk_fs_labels = [cellfun(@(x) [x ' L'], dk_fs_labels, 'UniformOutput', false), ...
                cellfun(@(x) [x ' R'], dk_fs_labels, 'UniformOutput', false)];


%% Surface files
% Process: Fetch brain annotations from neuromaps
sFiles = bst_process('CallProcess', 'process_nmp_fetch_maps', [], [], ...
    'sspace',        'surface', ...  % Surface
    'brainmaps_srf', {'Acetylcholine: M1_lsn3172176_naganawa2020_N24_Age40', ...
                      'Acetylcholine: VaCht_feobv_aghourian2017_N18_Age67', ...
                      'Acetylcholine: VaCht_feobv_bedard2019_N5_Age68', ...
                      'Acetylcholine: VaCht_feobv_tuominen_N4_Age37', ...
                      'Acetylcholine: a4b2_flubatine_hillmer2016_N30_Age34', ...
                      'Cannabinoid: CB1_omar_normandin2015_N77_Age30', ...
                      'Dopamine: D1_sch23390_kaller2017_N13_Age33', ...
                      'Dopamine: D2_flb457_sandiego2015_N55_Age32', ...
                      'Dopamine: D2_flb457_smith2017_N37_Age48', ...
                      'Dopamine: DAT_fpcit_dukart2018_N174_Age61', ...
                      'GABA: GABAa/bz_flumazenil_dukart2018_N6_Age43', ...
                      'GABA: GABAa/bz_flumazenil_norgaard2021_N16_Age27', ...
                      'Glutamate: mGluR5_abp688_dubois2015_N28_Age33', ...
                      'Glutamate: mGluR5_abp688_rosaneto_N22_Age68', ...
                      'Glutamate: mGluR5_abp688_smart2019_N73_Age20', ...
                      'Histamine: H3_gsk189254_gallezot2017_N8_Age32', ...
                      'Norepinephrine: NET_mrb_ding2010_N77_Age33', ...
                      'Opioid: MOR_carfentanil_kantonen2020_N204_Age32', ...
                      'Serotonin: 5-HT1a_way100635_savli2012_N35_Age26', ...
                      'Serotonin: 5-HT1b_p943_gallezot2010_N23_Age29', ...
                      'Serotonin: 5-HT1b_p943_savli2012_N23_Age29', ...
                      'Serotonin: 5-HT2a_cimbi36_beliveau2017_N29_Age23', ...
                      'Serotonin: 5-HT4_sb207145_beliveau2017_N59_Age26', ...
                      'Serotonin: 5-HT6_gsk215083_radnakrishnan2018_N30_Age37', ...
                      'Serotonin: 5-HTT_dasb_beliveau2017_N100_Age25', ...
                      'Serotonin: 5-HTT_dasb_savli2012_N18_Age31'});
surf_files = {sFiles.FileName};

% Short name: Ntx__subtype_desc_source
tmp = cellfun(@(x) regexprep(x, 'GABAa_bz', 'GABAa---bz'), surf_files, 'UniformOutput', false);
pattern = '^.*results_surface_(.*?__.*?_.*?_.*?)_.*';
surf_shortfiles = cellfun(@(x) regexp(x, pattern, 'tokens', 'once'), tmp, 'UniformOutput', false);
surf_shortfiles = cellfun(@(x) x{1}, surf_shortfiles, 'UniformOutput', false);
% Sort surface annotations filenames
[~, ix] = sort(lower(surf_files));
surf_files = surf_files(ix);
surf_shortfiles = surf_shortfiles(ix);

%% Volume files
% Create Subject with ICBM152 (MNI152) anatomy
SubjectName = 'ICBM152-neuromaps';
% Create Subject that is NOT using Default anatomy NOR Default channel
[~, iSubject] = db_add_subject(SubjectName, [], 0, 0);
% Set template ICBM152 for this Subject
sTemplates = bst_get('AnatomyDefaults');
iTemplate = find(strcmpi({sTemplates.Name}, 'ICBM152'));
db_set_template(iSubject, sTemplates(iTemplate), 0);

% Search '.nii.gz' files in the '../../tmp'
nii_files = dir('../../tmp/**/*.nii.gz');
nii_files = arrayfun(@(nii_file) bst_fullfile(nii_file.folder, nii_file.name), nii_files, 'UniformOutput', false);
n_nii_files = length(nii_files);
vol_files = cell(1, n_nii_files);
vol_shortfiles = cell(1, n_nii_files);
pattern = '^.*volume/(.*?)/subtype-(.*?)_source-(.*?)_desc-(.*?)_N-(.*?)_Age-(.*?)_';
for ix = 1 : n_nii_files
    match = regexp(nii_files{ix}, pattern, 'tokens');
    match = match{1};
    file_comment = sprintf('%s: %s_%s_%s_N%s_Age%s', match{1}, match{2}, match{4}, match{3}, match{5}, match{6});
    vol_shortfiles{ix} = sprintf('%s__%s_%s_%s', match{1}, match{2}, match{4}, match{3});
    % Import annotation files in Brainstorm as MRI volumes
    % Use isAutoAdjust=1, to re-slice 3mm brain annotations to 1mm, otherwise
    % if parcellation is downsampled, there is the risk of removing labels
    vol_files{ix} = import_mri(iSubject, nii_files{ix}, 'ALL-MNI', 0, 1, file_comment);
end
% Sort volume annotations filenames
[~, ix] = sort(lower(nii_files));
vol_files = vol_files(ix);
vol_shortfiles = vol_shortfiles(ix);

%% Surface parcellation
% Process: Scout time series: [68 scouts]
sMatrixFiles = bst_process('CallProcess', 'process_extract_scout', surf_files, [], ...
    'timewindow',     [], ...
    'scouts',         {atlasName, dk_fs_labels}, ...
    'flatten',        0, ...
    'scoutfunc',      'mean', ...  % Mean
    'pcaedit',        struct(...
         'Method',         'pcaa', ...
         'Baseline',       [0, 0.9992], ...
         'DataTimeWindow', [0, 0.9992], ...
         'RemoveDcOffset', 'file'), ...
    'isflip',         0, ...
    'isnorm',         0, ...
    'concatenate',    0, ...
    'save',           1, ...
    'addrowcomment',  1, ...
    'addfilecomment', []);
% Process: Concatenate time
sMatrixFile = bst_process('CallProcess', 'process_concat', sMatrixFiles, []);
% Process: Delete selected files
bst_process('CallProcess', 'process_delete', sMatrixFiles, [], ...
    'target', 1);  % Delete selected files
% Load matrix
sMat = in_bst_matrix(sMatrixFile.FileName);
% Correct scout names (remove "@ file" and preceding spaces)
sMat.Description = cellfun(@(x) regexprep(x, ' *@.*$', ''), sMat.Description, 'UniformOutput', false);
% Resort Scouts to FreeSurfer order
[Lia, Locb] = ismember(dk_fs_labels, sMat.Description);
if sum(Lia) ~= length(dk_fs_labels)
    error('Error with scout labels')
end
parcel_srf_results = sMat.Value(Locb, [1:2:length(surf_files)*2]);
% Save as CSV
T = array2table(parcel_srf_results, 'RowNames', strrep(dk_fs_labels, ' ', '_'), 'VariableNames', surf_shortfiles);
writetable(T, './parcellation_dk_surface_brainstorm.csv', 'WriteRowNames', true);

%% Volume parcellation
% Find Desikan-Killiany anatomical parcellation for ICBM512 anatomy
sSubject = bst_get('Subject', iSubject);
iMriAtlas = find(strcmpi({sSubject.Anatomy.Comment}, atlasName));
if isempty(iMriAtlas)
    % Find atlas with typo (until January 2023)
    iMriAtlas = find(strcmpi({sSubject.Anatomy.Comment}, 'Deskian-Killiany'));
end
sMriAtlas = in_mri_bst(sSubject.Anatomy(iMriAtlas).FileName);
% Keep 68 cortical labels (ids 1001 to 2025, except for x004), already in FreeSurfer order
ix_parcels_del = [sMriAtlas.Labels{:,1}] < 1001 | [sMriAtlas.Labels{:,1}] == 1004 ...
               | [sMriAtlas.Labels{:,1}] > 2035 | [sMriAtlas.Labels{:,1}] == 2004;
sMriAtlas.Labels(ix_parcels_del, :) = [];
% Find cube indices for parcellation
for iParcellation = 1 : size(sMriAtlas.Labels, 1)
    sMriAtlas.Labels{iParcellation, 4} = find(sMriAtlas.Cube == sMriAtlas.Labels{iParcellation, 1});
end
% Values to save
parcel_vol_results = zeros(size(sMriAtlas.Labels, 1), length(n_nii_files));
% Apply parcellation to imported volumes
for ix = 1 : n_nii_files
    % Load MRI
    sMriAnnot = in_mri_bst(vol_files{ix});
    for iParcellation = 1 : size(sMriAtlas.Labels, 1)
        parcel_vol_results(iParcellation, ix) = mean(sMriAnnot.Cube(sMriAtlas.Labels{iParcellation, 4}));
    end
end
% Save as CSV
T = array2table(parcel_vol_results, 'RowNames', strrep(dk_fs_labels, ' ', '_'), 'VariableNames', vol_shortfiles);
writetable(T, './parcellation_dk_volume_brainstorm.csv', 'WriteRowNames', true);
