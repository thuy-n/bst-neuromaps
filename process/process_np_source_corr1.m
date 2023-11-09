function varargout = process_np_source_corr1( varargin )
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
% Authors: Raymundo Cassani 2023

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Get List
    brainmapListSrf = process_np_fetch_maps('GetBrainMapsList', 'surface');
    brainmapListVol = process_np_fetch_maps('GetBrainMapsList', 'volume');
    % Description the process
    sProcess.Comment     = 'Spatial correlation';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 601;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
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
    sProcess.options.brainmaps_vol.Value   = strjoin(brainmapListVol, char(10)); % Add checkbox to select all?
    sProcess.options.brainmaps_vol.Class   = 'volume';
    % === METRIC, CORRECTION
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Get options
    space = sProcess.options.sspace.Value;
    if strcmpi(space, 'surface')
        brainmapsStr = sProcess.options.brainmaps_srf.Value;
    elseif strcmpi(space, 'volume')
        brainmapsStr = sProcess.options.brainmaps_vol.Value;
    end
    % Verify correct source space
    for iInput = 1 : length(sInputs)
        sResultsMat = in_bst_results(sInputs(iInput).FileName, 0, 'HeadModelType');
        if ~strcmpi(sResultsMat.HeadModelType, space)
            disp('Select the proper source space');
            %return with error
        end
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

    % Get Maps and their Surface
    [MapFiles, MapsSurfaceFile] = process_np_fetch_maps('PrepareNeuromap', space, brainmaps);

    for iInput = 1 : length(sInputs)
        % Compute correlations
        sMatrixMat = process_np_source_corr2('CorrelationSurfaceMaps', sInputs(iInput).FileName, MapFiles, MapsSurfaceFile);
        % === SAVE FILE ===
        % Add history entry
        sMatrixMat = bst_history('add', sMatrixMat, 'process', sprintf('Brain map spatial correlation for %s: ', sInputs(iInput).FileName));
        % Output filename
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        OutputFiles{end+1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_neuromaps');
        % Save file
        bst_save(OutputFiles{end}, sMatrixMat, 'v6');
        % Register in database
        db_add_data(sInputs(iInput).iStudy, OutputFiles{end}, sMatrixMat);
    end
    % Update whole tree
    panel_protocols('UpdateTree');
end
