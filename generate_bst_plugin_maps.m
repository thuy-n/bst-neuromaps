function generate_bst_plugin_maps()
% GENERATE_BST_PLUGIN_MAPS
%
% USAGE: generate_bst_plugin_maps(bstNeuromapsDir)
%

% 'bst-neuromaps' and 'brainstorm3' are expected at the same level
bstNeuromapsDir = bst_fullfile(fileparts(bst_get('BrainstormHomeDir')), 'bst-neuromaps');
if ~exist(bstNeuromapsDir, 'dir')
    return
end

% Annotations fetched from neuromaps
% Surface: FsAverage space with 327,684 vertices (~164k vertices per hemisphere)
% Volume:  MNI152 space 1mm [193x229x193]
nmpSurfMapsDir = bst_fullfile(bstNeuromapsDir, 'tmp', 'surface');
nmpVolMapsDir  = bst_fullfile(bstNeuromapsDir, 'tmp', 'volume');

% Annotations are imported in Brainstorm MNI ICBM152 default anatomy
% https://neuroimage.usc.edu/brainstorm/Tutorials/DefaultAnatomy
% Surface: Cortex (pial) with 15,002 vertices
% Volume:  MNI152 space 1mm [193x229x193]
bstMapsDir = bst_fullfile(bstNeuromapsDir, 'maps');
bstSurfMapsDir = bst_fullfile(bstNeuromapsDir, 'maps', 'surface');
bstVolMapsDir  = bst_fullfile(bstNeuromapsDir, 'maps', 'volume');
if exist(bstMapsDir, 'dir')
    rmdir(bstMapsDir, 's');
end
mkdir(bstMapsDir);
mkdir(bstSurfMapsDir);
mkdir(bstVolMapsDir);

% Surafce: Find all GIfTI (.gii) files, only left hemisphere
nmpSurfMapsFiles = dir(bst_fullfile(nmpSurfMapsDir, '**/*lh.shape.gii'));
nNmpSurfMaps = length(nmpSurfMapsFiles);
% Volume:  Find all NIfTI (.nii.gz) files
nmpVolMapsFiles  = dir(bst_fullfile(nmpVolMapsDir, '**/*.nii.gz'));
nNmpVolMaps = length(nmpVolMapsFiles);
% Check that Brainstorm is running with local database
bstDbDir = bst_get('BrainstormDbDir');
if isempty(regexp(bstDbDir, 'local_db$', 'once'))
    warning('Run this script with brainstorm using local db');
    return
end

%% GENERATE ANNOTATIONS TO IMPORT IN BRAINSTORM
%% PREPARE PROTOCOL AND SUBJECT
% Generate surface and volume annotations to distribute with bst-neuromaps plugin
% Both annotation types rely on the MNI ICBM152 Brainstorm template

ProtocolName = 'bstNeuromaps';
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% Create two subjects:
% 1. Subject using 'FsAverage'   template, to import annotations from neuromaps
% 2. Subject using 'MNI ICBM152' template, to generate annotations for plugin
SubjectNames  = {'Subject_FsA', 'Subject_MNI'};
templateAnats = {'FsAverage_2020',   'ICBM152'};
for iSubjectName = 1 : length(SubjectNames)
    % Create subject
    [~, iSubject] = db_add_subject(SubjectNames{iSubjectName}, [], 0, 0);
    % Set anatomy template
    sTemplates = bst_get('AnatomyDefaults');
    ix = find(strcmpi({sTemplates.Name}, templateAnats{iSubjectName}));
    db_set_template(iSubject, sTemplates(ix(1)), 0);
end

%% SURFACE MAPS
% 1. Import annotations into subject using FsAverage template
% 2. Project from FsAverage cortex (328k) to ICBM152 cortex (15k)
% 3. Save annotations files in the ICBM152 cortex

% 1. Import annotations into subject using FsAverage template
%
% For Subject using FsAverage template, set as default cortical surface the one
% with 327684 vertices (163842 per hemisphere)
% This surface matches the surface for the neuromaps-fetched surface annotations
[sSubjectFsa, iSubjectFsa] = bst_get('Subject', SubjectNames{1});
ix = find(~cellfun(@isempty,(regexpi({sSubjectFsa.Surface.Comment}, '.*cortex_327684V'))));
sSubjectFsa = db_surface_default(iSubjectFsa, 'Cortex', ix, 1);
fsaSurfaceFile = sSubjectFsa.Surface(sSubjectFsa.iCortex).FileName;
fsaSurfMapImportFiles = cell(nNmpSurfMaps, 1);
for iMap = 1 : nNmpSurfMaps
    % Build map Comment from map File info, e.g.:
    %  from: ./Acetylcholine/subtype-VaCht_source-aghourian2017_desc-feobv_N-18_Age-67_space-fsaverage_den-164k_lh.shape.gii
    %  to:     Acetylcholine: VaCht_feobv_aghourian2017_N18_Age67
    [~, mapFolder] = bst_fileparts(nmpSurfMapsFiles(iMap).folder);
    tokens = regexp(nmpSurfMapsFiles(iMap).name, 'subtype-(.*?)_source-(.*?)_desc-(.*?)_N-([0-9]*)_Age-([0-9]*)', 'tokens', 'once');
    mapComment = sprintf('%s: %s_%s_%s_N%s_Age%s', mapFolder, tokens{1}, tokens{3}, tokens{2}, tokens{4:5});
    % Replace '---' with '/', replaced in Python to have a valid filepath
    mapComment = strrep(mapComment, '---', '/');
    % Left and right hemisphere files
    mapFilePathL = bst_fullfile(nmpSurfMapsDir, mapFolder, nmpSurfMapsFiles(iMap).name);
    mapFilePathR = strrep(mapFilePathL, 'lh.shape.gii', 'rh.shape.gii');
    % Condition for map Category
    [~, iStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubjectFsa.Name, mapFolder));
    if isempty(iStudy)
        iStudy = db_add_condition(sSubjectFsa.Name, mapFolder);
    end
    % Import brain map to database
    MapFile = import_sources(iStudy, fsaSurfaceFile, mapFilePathL, mapFilePathR, 'GII', mapComment);
    % Get N (number of subjects averaged to obtained brain map)
    mapN = regexp(mapComment, 'N([0-9]*)', 'tokens', 'once');
    mapN = sscanf(mapN{1}, '%d');
    % Update Leff
    tmp.Leff = mapN;
    bst_save(MapFile, tmp, [], 1);
    fsaSurfMapImportFiles{iMap} = file_short(MapFile);
end

% 2. Project from FsAverage cortex (328k) to ICBM152 cortex (15k)
%
% For Subject using MNI template, set as default cortical surface the one
% Set cortical surface with 15002 vertices MUST be the default one
% This low-resolution is distributed with the FsAverage template in Brainstorm
[sSubjectMni, iSubjectMni] = bst_get('Subject', SubjectNames{2});
ix = find(~cellfun(@isempty,(regexpi({sSubjectMni.Surface.Comment}, '.*cortex_15002V'))));
sSubjectMni = db_surface_default(iSubjectMni, 'Cortex', ix, 1);
bstSurfaceFile = sSubjectMni.Surface(sSubjectMni.iCortex).FileName;
bstSurfMapImportFiles = cell(nNmpSurfMaps, 1);
for iMap = 1 : nNmpSurfMaps
    sNmpMapMat = in_bst_results(fsaSurfMapImportFiles{iMap}, 0, 'Comment');
    sMapProjFile = bst_project_sources(fsaSurfMapImportFiles(iMap), bstSurfaceFile, 0, 0);
    bstSurfMapImportFiles(iMap) = sMapProjFile;
    % Update comment
    tmp.Comment = sNmpMapMat.Comment;
    bst_save(file_fullpath(bstSurfMapImportFiles{iMap}), tmp, [], 1);
end

% 3. Save annotations files in the ICBM152 cortex
%
% Export annotations in Brainstorm default anatomy surface (sources)
% Filenaming follows 'BST' format in the Brainstorm function 'export_result()'
for iMap = 1 : length(bstSurfMapImportFiles)
    [mapFolder, mapFileName, mapExt] = bst_fileparts(bstSurfMapImportFiles{iMap});
    [~, mapFolder] = bst_fileparts(mapFolder);
    if exist(bst_fullfile(bstSurfMapsDir, mapFolder), 'dir') == 0
        mkdir(bst_fullfile(bstSurfMapsDir, mapFolder));
    end
    % Remove string with datetime (_YYMMDD_HHMM) and original SubjectName (_Subject_FsA) from filename
    mapFileName = regexprep(mapFileName, ['_[0-9]{6}_[0-9]{4}_' sSubjectFsa.Name '$'], '');
    % Remove string with 'results_' from filename
    mapFileName = regexprep(mapFileName, '$results_', '');
    % Update extension from '.mat' to '_sources.mat'
    mapFileName = [mapFileName, '_sources', mapExt];
    % Copy to bst-neuromaps
    srcMap = file_fullpath(bstSurfMapImportFiles{iMap});
    dstMap = bst_fullfile(bstSurfMapsDir, mapFolder, mapFileName);
    copyfile(srcMap, dstMap, 'f');
end
% Copy Brainstorm default anatomy surface
bstSurfaceFile = file_fullpath(bstSurfaceFile);
copyfile(bstSurfaceFile, bstSurfMapsDir, 'f');

%% VOLUME MAPS
% Copy annotations fetched with neuromaps as they are also in the ICBM152 volume
nmpVolMapImportFiles = cell(nNmpVolMaps, 1);
for iMap = 1 : length(nmpVolMapsFiles)
    % Build map Comment from map File info, e.g.:
    % from: ./Acetylcholine/subtype-VaCht_source-aghourian2017_desc-feobv_N-18_Age-67_space-mni152_den-1mm.nii.gz
    % to:     Acetylcholine: VaCht_feobv_aghourian2017_N18_Age67
    [~, mapFolder] = bst_fileparts(nmpVolMapsFiles(iMap).folder);
    tokens = regexp(nmpVolMapsFiles(iMap).name, 'subtype-(.*?)_source-(.*?)_desc-(.*?)_N-([0-9]*)_Age-([0-9]*)', 'tokens', 'once');
    mapComment = sprintf('%s: %s_%s_%s_N%s_Age%s', mapFolder, tokens{1}, tokens{3}, tokens{2}, tokens{4:5});
    % Replace '---' with '/', replaced in Python to have a valid filepath
    mapComment = strrep(mapComment, '---', '/');
    % Volume file
    mapFilePath = bst_fullfile(nmpVolMapsDir, mapFolder, nmpVolMapsFiles(iMap).name);
    % TBD for volume files:
    %  1. As MRI files (too large), or
    %  2. As sources in Grid (need to find a way to make deterministic the grid)
end


%% EXIT
% Stop Brainstorm
brainstorm stop
% Quit Matlab
exit

end
