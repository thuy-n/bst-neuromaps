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
    brainmapListSrf = GetBrainMapsList('surface');
    brainmapListVol = GetBrainMapsList('volume');
    % Description the process
    sProcess.Comment     = 'Spatial correlation';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 601;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/Neuromaps';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % === DESCRIPTION
    sProcess.options.help.Comment = ['Spatial comparison between: <BR>' ...
                                     'Each file and selected brain annotations<BR>'];
    sProcess.options.help.Type    = 'label';
    % === SOURCE SPACE
    sProcess.options.sspace.Comment    = {'Surface', 'Volume', '<B>Source space:</B>'; ...
                                          'surface', 'volume', ''};
    sProcess.options.sspace.Type       = 'radio_linelabel';
    sProcess.options.sspace.Value      = 'surface';
    sProcess.options.sspace.Controller = struct('surface', 'surface', 'volume', 'volume');
    sProcess.options.sspace.Hidden     = 1;
    % === SURFACE BRAIN MAPS
    sProcess.options.brainmaps_srf.Comment = [brainmapListSrf, {'Brain annotations from Neuromaps'}];
    sProcess.options.brainmaps_srf.Type    = 'list_vertical';
    sProcess.options.brainmaps_srf.Value   = '';
    sProcess.options.brainmaps_srf.Class   = 'surface';
    % === VOLUME BRAIN MAPS
    sProcess.options.brainmaps_vol.Comment = [brainmapListVol, {'Brain annotations from Neuromaps'}];
    sProcess.options.brainmaps_vol.Type    = 'list_vertical';
    sProcess.options.brainmaps_vol.Value   = '';
    sProcess.options.brainmaps_vol.Class   = 'volume';
    sProcess.options.brainmaps_vol.Hidden  = 1;
    % === METRIC
    sProcess.options.corrmetric.Comment = {'Pearson corr', 'Spearman corr', 'Comparison metric:'; 'Pearson', 'Spearman', ''};
    sProcess.options.corrmetric.Type    = 'radio_linelabel';
    sProcess.options.corrmetric.Value   = 'Pearson';
    % === REMOVE ZEROS
    sProcess.options.removezeros.Comment = 'Ignore zeros when computing correlation (default=True)';
    sProcess.options.removezeros.Type    = 'checkbox';
    sProcess.options.removezeros.Value   = 1;
    % === SPIN TEST
    sProcess.options.nspins.Comment = 'Number of spins for spin test (0 = no spin test): ';
    sProcess.options.nspins.Type    = 'value';
    sProcess.options.nspins.Value   = {10, '', 0};
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
        brainmaps = sProcess.options.brainmaps_srf.Value;
    elseif strcmpi(space, 'volume')
        brainmaps = sProcess.options.brainmaps_vol.Value;
        bst_error('Volume brain annotations are not supported yet.', 'BST-Neuromaps');
        return
    end
    nSpins      = sProcess.options.nspins.Value{1};
    removeZeros = sProcess.options.removezeros.Value;
    processTab  = 1;
    corrMethod  = sProcess.options.corrmetric.Value;


    % Load neuromaps plugin if needed
    PlugDesc = bst_plugin('GetInstalled', 'neuromaps');
    if ~PlugDesc.isLoaded
        bst_plugin('Load', 'neuromaps');
    end

    % Verify correct source space
    for iInput = 1 : length(sInputs)
        sResultsMat = in_bst_results(sInputs(iInput).FileName, 0, 'HeadModelType');
        if ~strcmpi(sResultsMat.HeadModelType, space)
            bst_report('Error', sProcess, sInputs(iInput), 'Input file has a different source space than the requested maps.');
            OutputFiles = [];
            return;
        end
    end

    % Get Maps and their Surface
    [MapFiles, MapsSurfaceFiles] = process_np_fetch_maps('PrepareNeuromap', space, brainmaps);
    if isempty(MapFiles) || isempty(MapsSurfaceFiles)
        bst_report('Error', sProcess, [], 'Could not find requested maps.');
        OutputFiles = [];
        return;
    end

    % Compute and save spatial correlations
    OutputFiles = process_np_source_corr2('CorrelationSurfaceMaps', sInputs, MapFiles, MapsSurfaceFiles, removeZeros, nSpins, corrMethod, processTab, 1);

    % Update whole tree
    panel_protocols('UpdateTree');
end

function mapComments = GetBrainMapsList(space)
    % Brain map Comments from brain FileNames in bst_neuromaps Plugin
    % Find all Brainstorm source files
    plugDesc = bst_plugin('GetInstalled', 'neuromaps');
    mapFiles = dir(fullfile(plugDesc.Path, plugDesc.SubFolder,  'maps', space,  '**/*_sources.mat'));
    mapComments = cell(1, length(mapFiles));
    for iMap = 1 : length(mapFiles)
        % Get comment from file content
        tmp = load(fullfile(mapFiles(iMap).folder, mapFiles(iMap).name), 'Comment');
        mapComments{iMap} = tmp.Comment;
    end
end
