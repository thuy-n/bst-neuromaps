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
    % TODO: Update to GUI Lists
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
    % Load neuromaps plugin if needed
    PlugDesc = bst_plugin('GetDescription', 'neuromaps');
    if ~PlugDesc.isLoaded
        bst_plugin('Load', 'neuromaps');
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
    MapFiles = PrepareNeuromap(space, brainmaps);
    if isempty(MapFiles)
        bst_report('Error', sProcess, [], 'Could not find requested maps.');
        OutputFiles = [];
        return;
    end

    % Output filenames
    OutputFiles = MapFiles;
    % Update whole tree
    panel_protocols('UpdateTree');
end


%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================
function [MapFiles, MapSurfaceFile] = PrepareNeuromap(space, mapComments)
    MapFiles = {};
    MapSurfaceFile = {};

    mapInfos = repmat(struct('Comment',  '', ...
                             'Category', '', ...
                             'FullPath', ''), 0);
    % Check existence of requested brain files
    for iMap = 1 : length(mapComments)
        % Get category
        mapCategory = regexp(mapComments{iMap}, '(.*?):', 'tokens', 'once');
        mapCategory = mapCategory{1};
        % Go from filename to comment
        mapFileName = ['results_', space, '_', mapComments{iMap}, '_sources.mat'];
        mapFileName = regexprep(mapFileName, ': ', '__');
        mapFullPath = bst_fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', space, mapCategory, mapFileName);
        if ~exist(mapFullPath, 'file')
            bst_error(sprintf('Requested brain map "%s" does not exist', mapFullPath));
            return
        end
        mapInfos(end+1).Comment = mapComments{iMap};
        mapInfos(end).Category  = mapCategory;
        mapInfos(end).FullPath  = mapFullPath;
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
            bst_error('FsAverage template could not be found');
            return
        elseif length(ix) > 1
            bst_error('There are more than one FsAverage templates');
            return
        end
        db_set_template(iSubject, sTemplates(ix), 0);
        [sSubject, iSubject] = bst_get('Subject', SubjectName);
    end
    % Check that it uses FsAverage
    if isempty(regexpi(sSubject.Anatomy(sSubject.iAnatomy).Comment, 'fsaverage', 'match'))
        bst_error('Subject is not using FsAverage template anatomy');
        return
    end
    % Cortical surface with 15002 vertices (distributed in Brainstorm) MUST be the default one
    % This is because the maps in the bst-neuromaps were obtained for this surface
    % See the code `generate_bst_plugin_maps` of bst-neuromaps
    ix = find(~cellfun(@isempty,(regexpi({sSubject.Surface.Comment}, '.*cortex_15002V'))));
    if ~isequal(sSubject.iCortex, ix)
        sSubject = db_surface_default(iSubject, 'Cortex', ix, 1);
    end
    MapSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    % Sanity check by comparing with the cortical surface in the plugin
    PluginSurfaceFile = fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', 'surface',  'tess_cortex_pial_low.mat');
    sMapSrfMat = in_tess_bst(MapSurfaceFile);
    sPluSrfMat = in_tess_bst(PluginSurfaceFile);
    if ~isequal(sMapSrfMat.Vertices, sPluSrfMat.Vertices)
        bst_error('The FsAverage templates (neuromaps plugin and Protocol) are different');
        return
    end
    % Import brain maps in Neuromaps Subject
    for iMap = 1 : length(mapInfos)
        % Condition for map Category
        [~, iStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, mapInfos(iMap).Category));
        if isempty(iStudy)
            iStudy = db_add_condition(sSubject.Name, mapInfos(iMap).Category);
            if isempty(iStudy)
                bst_error(sprintf('Study could not be created : "%s".', mapInfos(iMap).Category));
                return
            end
        end
        % Get brain map file
        sStudy = bst_get('Study', iStudy);
        [~, ix] = intersect({sStudy.Result.Comment}, mapInfos(iMap).Comment, 'stable');
        if ~isempty(ix)
            % Brain map already in database
            MapFiles{iMap} = sStudy.Result(ix).FileName;
        else
            % Import brain map to database
            MapFile = import_sources(iStudy, MapSurfaceFile, mapInfos(iMap).FullPath, [], 'BST', mapInfos(iMap).Comment);
            % Get N (number of subjects averaged to obtained brain map)
            mapN = regexp(mapInfos(iMap).Comment, 'N([0-9]*)', 'tokens', 'once');
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
    % Find all Brainstorm source files
    mapFiles = dir(fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'bst-neuromaps-main', 'maps', space,  '**/*_sources.mat'));
    mapComments = {};
    for iMap = 1 : length(mapFiles)
        % Go from filename to comment
        mapComments{end+1} = regexprep(mapFiles(iMap).name, ['^results_', space, '_'], '');
        mapComments{end}   = regexprep(mapComments{end}, '_sources.mat$', '');
        mapComments{end}   = regexprep(mapComments{end}, '__', ': ');
    end
end
