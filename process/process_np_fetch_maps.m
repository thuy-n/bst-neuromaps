function varargout = process_np_fetch_maps( varargin )
% PROCESS_CORR_MAPS

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Le Thuy Duong Nguyen, Raymundo Cassani 2023

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Get List
    brainmapListSrf = GetBrainMapsList('surface');
    brainmapListVol = GetBrainMapsList('volume');
    % Description the process
    sProcess.Comment     = 'Fetch brain annotations from Neuromaps';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 600;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    % === DESCRIPTION
    sProcess.options.help.Comment = ['This process computes Pearson''s correlation coefficient between <BR>' ...
                                     'a surface source file and a brain map<BR>'];
    sProcess.options.help.Type    = 'label';
    % === SOURCE SPACE
    sProcess.options.sspace.Comment    = {'Surface', 'Volume', '<B>Source space:</B>'; ...
                                          'surface', 'volume', ''};
    sProcess.options.sspace.Type       = 'radio_linelabel';
    sProcess.options.sspace.Value      = 'surface';
    sProcess.options.sspace.Controller = struct('surface', 'surface', 'volume', 'volume');
    % === SURFACE BRAIN MAPS
    sProcess.options.brainmaps_srf.Comment = 'Brain maps from Neuromaps (wiil be a list)';
    sProcess.options.brainmaps_srf.Type    = 'textarea';
    sProcess.options.brainmaps_srf.Value   = strjoin(brainmapListSrf, char(10)); % Add checkbox to select all?
    sProcess.options.brainmaps_srf.Class   = 'surface';
    % === VOLUME BRAIN MAPS
    sProcess.options.brainmaps_vol.Comment = 'Brain maps from Neuromaps (wiil be a list)';
    sProcess.options.brainmaps_vol.Type    = 'textarea';
    sProcess.options.brainmaps_vol.Value   = strjoin(brainmapListSrf, char(10)); % Add checkbox to select all?
    sProcess.options.brainmaps_vol.Class   = 'volume';
    % === METRIC, CORRECTION
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    % Get options
    space = sProcess.options.sspace.Value;
    if strcmpi(space, 'surface')
        brainmapsStr = sProcess.options.brainmaps_srf.Value;
    elseif strcmpi(space, 'volume')
        brainmapsStr = sProcess.options.brainmaps_vol.Value;
    end

    % Get brain maps Comments
    brainmaps = str_split(brainmapsStr, 10);
    tmps = {};
    for iMap = 1 : length(brainmaps)
        tmp = strtrim(brainmaps{iMap});
        if ~isempty(tmp)
            tmps{end+1} = tmp;
        end
    end
    brainmaps = tmps;

    % Get brain maps
    MapFiles = PrepareNeuromap(brainmaps);

    % Output filenames
    OutputFiles = MapFiles;
    % Update whole tree
    panel_protocols('UpdateTree');
end


%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================
function [MapFiles, MapSurfaceFile] = PrepareNeuromap(mapComments)
    MapFiles = {};
    MapSurfaceFile = {};

    % Check existence of requested brain files
    for iMap = 1 : length(mapComments)
        [mapCat, mapFileNameL, mapFileNameR] = mapComment2mapFileName(mapComments{iMap});
        mapFilePathL = bst_fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', 'surface', mapCat, mapFileNameL);
        mapFilePathR = bst_fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', 'surface', mapCat, mapFileNameR);
        if ~exist(mapFilePathL, 'file')
            fprintf('Requested brain map "%s" does not exist', mapFilePathL);
            %return with error
        end
        if ~exist(mapFilePathR, 'file')
            fprintf('Requested brain map "%s" does not exist', mapFilePathR);
            %return with error
        end
    end

    % Get Neuromaps Subject
    SubjectName = 'Neuromaps';
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject)
        % Create Subject
        [~, iSubject] = db_add_subject(SubjectName, [], 0, 0);
        % Set anatomy to FsAverage template
        sTemplates = bst_get('AnatomyDefaults');
        ix = find(~cellfun(@isempty,(regexpi({sTemplates.Name}, 'fsaverage'))));
        if isempty(ix)
            disp('There is not FsAverage template')
            %return with error
        elseif length(ix) > 1
            disp('There are more than one FsAverage templates')
            %return with error
        end
        db_set_template(iSubject, sTemplates(ix), 0);
        [sSubject, iSubject] = bst_get('Subject', SubjectName);
    end
    % Check that it uses FsAverage
    if isempty(regexpi(sSubject.Anatomy(sSubject.iAnatomy).Comment, 'fsaverage', 'match'))
        disp('Subject is not using FsAverage anatomy');
        %return with error
    end
    % Cortical surface with 327684 vertices (163842 per hemisphere) MUST be the default one
    % This is because the maps in the bst-neuromaps are in the Fs 164k space
    ix = find(~cellfun(@isempty,(regexpi({sSubject.Surface.Comment}, '.*cortex.*327684.*'))));
    if ~isequal(sSubject.iCortex, ix)
        sSubject = db_surface_default(iSubject, 'Cortex', ix, 1);
    end
    MapSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;

    % Import brain maps in Neuromaps Subject
    for iMap = 1 : length(mapComments)
        [mapCat, mapFileNameL, mapFileNameR] = mapComment2mapFileName(mapComments{iMap});
        mapFilePathL = bst_fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', 'surface', mapCat, mapFileNameL);
        mapFilePathR = bst_fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', 'surface', mapCat, mapFileNameR);
        % Condition for map Category
        [~, iStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, mapCat));
        if isempty(iStudy)
            iStudy = db_add_condition(sSubject.Name, mapCat);
            if isempty(iStudy)
                error('Study could not be created : "%s".', mapCat);
                %return with error
            end
        end
        % Get brain map file
        sStudy = bst_get('Study', iStudy);
        [~, ix] = intersect({sStudy.Result.Comment}, mapComments{iMap}, 'stable');
        if ~isempty(ix)
            % Brain map already in database
            MapFiles{iMap} = sStudy.Result(ix).FileName;
        else
            % Import brain map to database
            MapFile = import_sources(iStudy, MapSurfaceFile, mapFilePathL, mapFilePathR, 'GII', mapComments{iMap});
            % Get N (number of subjects averaged to obtained brain map)
            mapN = regexp(mapComments{iMap}, 'N([0-9]*)', 'tokens', 'once');
            mapN = sscanf(mapN{1}, '%d');
            % Update MapFile
            sMapMat.Leff = mapN;
            bst_save(MapFile, sMapMat, [], 1);
            MapFiles{iMap} = file_short(MapFile);
        end
    end
end

function mapComments = GetBrainMapsList(space)
    % Brain map Comments from brain FileNames in bst_neuromaps Plugin
    if strcmpi(space, 'surface')
        % Find all gii files (only Lh)
        dir_str = {'surface','**/*lh.shape.gii'};
    elseif strcmpi(space, 'volume')
        % Find all nii.gz files
        dir_str = {'volume','**/*.nii.gz'};
    end
    % Search files
    mapFiles = dir(fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', dir_str{:} ));
    mapComments = {};
    for iMap = 1 : length(mapFiles)
        [~, mapFolder] = bst_fileparts(mapFiles(iMap).folder);
        mapComments{end+1} = mapFileName2mapComment(mapFolder ,mapFiles(iMap).name);
    end
end

function mapComment = mapFileName2mapComment(mapCat, mapFileName)
    % Go from brain map Category and map FileName (either Lh or Rh) to map Comment
    % e.g. ['Acetylcholine', 'source-aghourian2017_desc-feobv_N-18_Age-67_space-fsaverage_den-164k_lh.shape.gii'] > 'Acetylcholine : feobv_aghourian2017_N18_Age67'
    tokens = regexp(mapFileName, 'source-(.*?)_desc-(.*?)_N-([0-9]*)_Age-([0-9]*)', 'tokens', 'once');
    mapComment = sprintf('%s: %s_%s_N%s_Age%s', mapCat, tokens{2}, tokens{1}, tokens{3:4});
end

function [mapCat, mapFileNameL, mapFileNameR] = mapComment2mapFileName(mapComment)
    % Go from brain map Comment to brain map Category and map FileNames (Lh and Rh)
    % e.g., 'Acetylcholine : feobv_aghourian2017_N18_Age67' > ['Acetylcholine'
    %                                                          'source-aghourian2017_desc-feobv_N-18_Age-67_space-fsaverage_den-164k_lh.shape.gii',
    %                                                          'source-aghourian2017_desc-feobv_N-18_Age-67_space-fsaverage_den-164k_rh.shape.gii']
    tokens = regexp(mapComment, '(.*?)\s*:\s*(.*?)_(.*?)_N([0-9]*)_Age([0-9]*)', 'tokens', 'once');
    mapCat = tokens{1};
    mapFileNameL = sprintf('source-%s_desc-%s_N-%s_Age-%s_space-fsaverage_den-164k_lh.shape.gii', tokens{3}, tokens{2}, tokens{4:5});
    mapFileNameR = strrep(mapFileNameL, 'lh.shape.gii', 'rh.shape.gii');
end
