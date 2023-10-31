function varargout = process_corr_maps( varargin )
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
% Authors:  Le Thuy Duong Nguyen, Raymundo Cassani 2023

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Get List
    brainmapList = GetBrainMapsList();
    % Description the process
    sProcess.Comment     = 'Spatial correlation';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 600;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % === DESCRIPTION
    sProcess.options.help.Comment = ['This process computes Pearson''s correlation coefficient between <BR>' ...
                                     'a surface source file and a brain map<BR>'];
    sProcess.options.help.Type    = 'label';
    % === BRAIN MAP List from folders in Plugin
    sProcess.options.brainmaps.Comment = 'Brain maps from Neuromaps (wiil be a list)';
    sProcess.options.brainmaps.Type    = 'textarea';
    sProcess.options.brainmaps.Value   = 'Opioid: carfentanil_kantonen2020_N204_Age32, GABA: flumazenil_norgaard2021_N16_Age27'; % Allow all?


    % === METRIC
    % === CORRECTION
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFile = Run(sProcess, sInput) %#ok<DEFNU>
    % Get options
    brainmapsStr = sProcess.options.brainmaps.Value;
    % Verify that source is surface
    sResultsMat = in_bst_results(sInput.FileName, 0, 'HeadModelType', 'Time');
    if ~strcmpi(sResultsMat.HeadModelType, 'surface')
        disp('Only surface is supported');
        %return with error
    end
    % Check that there is only time point
    if length(sResultsMat.Time) > 2
        disp('Source map MUST be only one point');
        %return with error
    end

    % Get brain maps Comments
    brainmaps = str_split(brainmapsStr, ',');
    brainmaps = cellfun(@strtrim, brainmaps, 'UniformOutput', false);
    % Get Maps and their Surface
    [MapFiles, MapSurfaceFile] = PrepareNeuromap(brainmaps);

    % Project sources
    sResultsProjFileName = bst_project_sources({sInput.FileName}, MapSurfaceFile, 0, 0);
    [~, iStudyProj] = bst_get('ResultsFile', sResultsProjFileName{1});
    sResultsProjMat = in_bst_results(sResultsProjFileName{1}, 1);

    % Perform correlation
    corrVals = zeros(length(MapFiles), 1);
    for iMap = 1 : length(MapFiles)
        % Load map
        sMapMat = in_bst_results(MapFiles{iMap}, 1);
        corrVals(iMap) = bst_corrn(transpose(sMapMat.ImageGridAmp(:,1)), transpose(sResultsProjMat.ImageGridAmp(:,1)));
        fprintf('%s : %f\n', brainmaps{iMap}, corrVals(iMap));
    end

    % Delete projected
    db_delete_studies(iStudyProj);
    % Create matrix file
    sMatrixMat = db_template('matrixmat');
    sMatrixMat.Comment = 'Brain map corr';
    sMatrixMat.Time = [];
    sMatrixMat.Value = corrVals;         % size [nMaps,1]
    sMatrixMat.Description = brainmaps'; % size [nMaps,1]
    % Add history entry
    sMatrixMat = bst_history('add', sMatrixMat, 'process', sprintf('Brain map spatial correlation for %s: ', sInput(1).FileName));

    % === SAVE FILE ===
    % Output filename
    sStudy = bst_get('Study', sInput(1).iStudy);
    OutputFile{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_neuromaps');
    % Save file
    bst_save(OutputFile{1}, sMatrixMat, 'v6');
    % Register in database
    db_add_data(sInput(1).iStudy, OutputFile{1}, sMatrixMat);
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
            mapN = mapN{1};
            % Update MapFile
            sMapMat.nAvg = mapN;
            bst_save(MapFile, sMapMat, [], 1);
            MapFiles{iMap} = file_short(MapFile);
        end
    end
end

function mapComments = GetBrainMapsList()
    % Brain map Comments from brain FileNames in bst_neuromaps Plugin
    % Find all gii files (only Lh)
    mapFiles = dir(fullfile(bst_get('UserPluginsDir'), 'neuromaps', 'neuromaps', 'maps', 'surface','**/*lh.shape.gii'));
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
    mapComment = sprintf('%s : %s_%s_N%s_Age%s', lower(mapCat), tokens{2}, tokens{1}, tokens{3:4});
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
