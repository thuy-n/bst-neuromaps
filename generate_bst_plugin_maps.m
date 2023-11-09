function generate_bst_plugin_maps()
% GENERATE_BST_PLUGIN_MAPS
%
% USAGE: generate_bst_plugin_maps()
%
% INPUTS:
%    - bstTestName  : Test to run, usually a script in the ./toolbox/scripts folder
%    - bstUserName  : Cell array of signals {[nSignals1, nSamples1], [nSignals2, nSamples2], ...}
%
% For a given bstTestName, the script does:
%    1. Download tutorial data (if needed)
%    2. Run tutorial script
%    3. Send report by email to bstUserName

%% FIND HIGH-DEFINITION MAPS
% Add path 'bst-tests' with support functions
% 'bst-neuromaps' and 'brainstorm3' are expected at the same level
bstNeuromapsDir = bst_fullfile(fileparts(bst_get('BrainstormHomeDir')), 'bst-neuromaps');
if ~exist(bstNeuromapsDir, 'dir')
    return
end

% High definition maps
% Surface: FsAverage space with 327,684 vertices (~164k vertices per hemisphere)
% Volume:
hdSurfMapsDir = bst_fullfile(bstNeuromapsDir, 'tmp', 'surface');
hdVolMapsDir  = bst_fullfile(bstNeuromapsDir, 'tmp', 'volume');

% Low definition maps
% Surface: FsAverage space with 15,002 vertices (distributed in Brainstorm)
% Volume:
ldMapsDir = bst_fullfile(bstNeuromapsDir, 'maps');
ldSurfMapsDir = bst_fullfile(bstNeuromapsDir, 'maps', 'surface');
ldVolMapsDir  = bst_fullfile(bstNeuromapsDir, 'maps', 'volume');
if exist(ldMapsDir, 'dir')
    rmdir(ldMapsDir, 's');
end
mkdir(ldMapsDir);
mkdir(ldSurfMapsDir);
mkdir(ldVolMapsDir);

% Find all gii files (only left hemisphere)
hdSurfMapsFiles = dir(bst_fullfile(hdSurfMapsDir, '**/*lh.shape.gii'));
nHdSurfMaps = length(hdSurfMapsFiles);
% Find all nii files
hdVolMapsFiles  = dir(bst_fullfile(hdVolMapsDir, '**/*.nii.gz'));
nHdVolMaps = length(hdVolMapsFiles);
% Check that Brainstorm is run with local database
bstDbDir = bst_get('BrainstormDbDir');
if isempty(regexp(bstDbDir, 'local_db$', 'once'))
    warning('Run this script with brainstorm using local db');
    return
end

%% GENERATE LOW DEFINITION MAPS
%% PREPARE PROTOCOL AND SUBJECT
% Generate low resolution maps to distribute with bst-neuromaps plugin
% Surface: FsAverage space with 7.5k vertices per hemisphere
% Volume:

ProtocolName = 'bstNeuromaps';
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% Create Subject
SubjectName = 'Neuromaps';
[~, iSubject] = db_add_subject(SubjectName, [], 0, 0);
% Set anatomy to FsAverage template
sTemplates = bst_get('AnatomyDefaults');
ix = find(~cellfun(@isempty,(regexpi({sTemplates.Name}, 'fsaverage'))));
db_set_template(iSubject, sTemplates(ix), 0);
[sSubject, iSubject] = bst_get('Subject', SubjectName);


%% SURFACE MAPS
% Set cortical surface with 327684 vertices (163842 per hemisphere) MUST be the default one
% This is because the high-definition maps are in the fsaverage 164k space
ix = find(~cellfun(@isempty,(regexpi({sSubject.Surface.Comment}, '.*cortex_327684V'))));
sSubject = db_surface_default(iSubject, 'Cortex', ix, 1);
hdSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;

% Import high-definition Surface maps in Subject01
hdSurfMapImportFiles = cell(nHdSurfMaps, 1);
for iMap = 1 : nHdSurfMaps
    % Build map Comment from map File info
    % e.g. ['./Acetylcholine/source-aghourian2017_desc-feobv_N-18_Age-67_space-fsaverage_den-164k_lh.shape.gii'] > 'Acetylcholine: feobv_aghourian2017_N18_Age67'
    [~, mapFolder] = bst_fileparts(hdSurfMapsFiles(iMap).folder);
    tokens = regexp(hdSurfMapsFiles(iMap).name, 'source-(.*?)_desc-(.*?)_N-([0-9]*)_Age-([0-9]*)', 'tokens', 'once');
    mapComment = sprintf('%s: %s_%s_N%s_Age%s', mapFolder, tokens{2}, tokens{1}, tokens{3:4});
    % Left and right hemisphere files
    mapFilePathL = bst_fullfile(hdSurfMapsDir, mapFolder, hdSurfMapsFiles(iMap).name);
    mapFilePathR = strrep(mapFilePathL, 'lh.shape.gii', 'rh.shape.gii');
    % Condition for map Category
    [~, iStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, mapFolder));
    if isempty(iStudy)
        iStudy = db_add_condition(sSubject.Name, mapFolder);
    end
    % Import brain map to database
    MapFile = import_sources(iStudy, hdSurfaceFile, mapFilePathL, mapFilePathR, 'GII', mapComment);
    % Get N (number of subjects averaged to obtained brain map)
    mapN = regexp(mapComment, 'N([0-9]*)', 'tokens', 'once');
    mapN = sscanf(mapN{1}, '%d');
    % Update MapFile
    sMapMat.Leff = mapN;
    bst_save(MapFile, sMapMat, [], 1);
    hdSurfMapImportFiles{iMap} = file_short(MapFile);
end

% Project imported surface maps to low-definition (15k surface)
% Set cortical surface with 15002 vertices MUST be the default one
% This low-resolution is distributed with the FsAverage template in Brainstorm
ix = find(~cellfun(@isempty,(regexpi({sSubject.Surface.Comment}, '.*cortex_15002V'))));
sSubject = db_surface_default(iSubject, 'Cortex', ix, 1);
ldSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
ldSurfMapImportFiles = cell(nHdSurfMaps, 1);
for iMap = 1 : nHdSurfMaps
    sHdMapMat = in_bst_results(hdSurfMapImportFiles{iMap}, 0, 'Comment');
    sMapProjFile = bst_project_sources(hdSurfMapImportFiles(iMap), ldSurfaceFile, 0, 0);
    ldSurfMapImportFiles(iMap) = sMapProjFile;
    % Update comment
    sMapMat.Comment = sHdMapMat.Comment;
    bst_save(file_fullpath(ldSurfMapImportFiles{iMap}), sMapMat, [], 1);
end

% Export low-definition maps (sources)
% Filenaming follows 'BST' format in the Brainstorm function 'export_result()'
for iMap = 1 : length(ldSurfMapImportFiles)
    [mapFolder, mapFileName, mapExt] = bst_fileparts(ldSurfMapImportFiles{iMap});
    [~, mapFolder] = bst_fileparts(mapFolder);
    if exist(bst_fullfile(ldSurfMapsDir, mapFolder), 'dir') == 0
        mkdir(bst_fullfile(ldSurfMapsDir, mapFolder));
    end
    % Remove string with datetime (_YYMMDD_HHMM) and number of vertices (_nV) from filename
    mapFileName = regexprep(mapFileName, '_[0-9]{6}_[0-9]{4}_15002V$', '');
    % Remove string with 'results_' from filename
    mapFileName = regexprep(mapFileName, '$results_', '');
    % Update extension from '.mat' to '_sources.mat'
    mapFileName = [mapFileName, '_sources', mapExt];
    % Copy to bst-neuromaps
    srcMap = file_fullpath(ldSurfMapImportFiles{iMap});
    dstMap = bst_fullfile(ldSurfMapsDir, mapFolder, mapFileName);
    copyfile(srcMap, dstMap, 'f');
end

% Copy low-definition surface
ldSurfaceFile = file_fullpath(ldSurfaceFile);
copyfile(ldSurfaceFile, ldSurfMapsDir, 'f');

%% VOLUME MAPS
% Import Volume maps in Subject01
hdVolMapImportFiles = cell(nHdVolMaps, 1);
for iMap = 1 : length(hdVolMapsFiles)
    % Build map Comment from map File info
    % e.g. ['./Acetylcholine/source-aghourian2017_desc-feobv_N-18_Age-67_space-mni152_den-1mm.nii.gz'] > 'Acetylcholine: feobv_aghourian2017_N18_Age67'
    [~, mapFolder] = bst_fileparts(hdVolMapsFiles(iMap).folder);
    tokens = regexp(hdVolMapsFiles(iMap).name, 'source-(.*?)_desc-(.*?)_N-([0-9]*)_Age-([0-9]*)', 'tokens', 'once');
    mapComment = sprintf('%s: %s_%s_N%s_Age%s', mapFolder, tokens{2}, tokens{1}, tokens{3:4});
    % Volume file
    mapFilePath = bst_fullfile(hdVolMapsDir, mapFolder, hdVolMapsFiles(iMap).name);
    % TBD for volume files:|
    %  1. As MRI files (too large), or
    %  2. As sources in Grid (need to find a way to make deterministic the grid)
end


%% EXIT
% Stop Brainstorm
brainstorm stop
% Quit Matlab
exit

end
