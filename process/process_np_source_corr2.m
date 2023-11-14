function varargout = process_np_source_corr2( varargin )
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
    % Description the process
    sProcess.Comment     = 'Spatial correlation AxB';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 602;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    % === DESCRIPTION
    sProcess.options.help.Comment = ['This process computes Pearson''s correlation coefficient between <BR>' ...
                                     '1 source file and N brain maps<BR>'];
    sProcess.options.help.Type    = 'label';
    % === METRIC, CORRECTION

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
function OutputFiles = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
    OutputFiles = {};
    nSpins = sProcess.options.nspins.Value{1};
    % All files must have surface or volume sources
    sResultsB1 = in_bst_results(sInputsB(1).FileName, 0, 'HeadModelType', 'Time');
    % Verify correct source space
    sInputs = [sInputsA, sInputsB];
    for iInput = 1 : length([sInputsA, sInputsB])
        sResultsMat = in_bst_results(sInputs(iInput).FileName, 0, 'HeadModelType');
        if ~strcmpi(sResultsMat.HeadModelType, sResultsB1.HeadModelType)
            disp('Select the proper source space');
            %return with error
        end
    end
    % Verify time definition
    % FilesB must have the same time axis: TimesB
    for iInputB = 1 : length(sInputsB)
        sResultsMat = in_bst_results(sInputs(iInputB).FileName, 0, 'Time');
        if ~isequal(sResultsMat.Time, sResultsB1.Time)
            disp('Sources in FilesB must have the same time axis');
            %return with error
        end
    end
    % Validate time dimensions for FilesA and FilesB
    % filesA (1 sample)  vs filesB (1 sample)   OK      1 Corr value
    % filesA (N samples) vs filesB (1 sample)   OK      N Corr values
    % filesA (N samples) vs filesB (N samples)  OK      N Corr values
    % filesA (M samples) vs filesB (N samples)  Not OK
    for iInputB = 1 : length(sInputsB)
        sResultsMat = in_bst_results(sInputs(iInputB).FileName, 0, 'Time');
        % FilesA must have the same time axis as TimesB
        if length(sResultsB1.Time) > 1 && ~isequal(sResultsMat.Time, sResultsB1.Time)
            disp('Sources in FilesA must have the same time axis as FilesB');
            %return with error
        end
    end

    % Get maps and their surface from sInputB
    AllMapFiles = {sInputsB.FileName};
    AllMapsSurfaceFiles = cell(1, length(AllMapFiles));
    for iInputB = 1 : length(sInputsB)
        sResultsMat = in_bst_results(sInputs(iInputB).FileName, 0, 'SurfaceFile');
        AllMapsSurfaceFiles{iInputB} = sResultsMat.SurfaceFile;
    end
    % Group by common surface file
    [MapsSurfaceFiles, ia, ic] = unique(AllMapsSurfaceFiles);
    MapFilesGroups = cell(1, length(ia));
    for iia = 1 : length(ia)
        MapFilesGroups{iia} = AllMapFiles(ic == ia(iia));
    end


    % For each FileA compute correlations with all FilesB
    for iInputA = 1 : length(sInputsA)
        sMatrixMat = db_template('matrixmat');
        sMatrixMat.Comment = 'Brain map corr';
        sMatrixMat.Options.SpinTest = [];
        % Accumulator for matrices, one matrix per surface file
        for iMapSurfaceFile = 1 : length(MapsSurfaceFiles)
            % Maps with common Surface file
            MapsSurfaceFile = MapsSurfaceFiles{iMapSurfaceFile};
            MapFiles = MapFilesGroups{iMapSurfaceFile};
            % Compute correlations
            sMatrixTmpMat = CorrelationSurfaceMaps(sInputsA(iInputA).FileName, MapFiles, MapsSurfaceFile, nSpins);
        end
        % Concatename matrices (add results from other maps)
        sMatrixMat.Time = sMatrixTmpMat.Time;
        sMatrixMat.Value =  [sMatrixMat.Value; sMatrixTmpMat.Value];                   % size [nMaps, nTime]
        sMatrixMat.Description = [sMatrixMat.Description; sMatrixTmpMat.Description];  % size [nMaps,1]
        % TODO: Do we need to save this?
        if nSpins > 0
            sMatrixMat.Options.SpinTest = cat(1, sMatrixMat.Options.SpinTest, sMatrixTmpMat.Options.SpinTest); % size [nMaps, nTime, nSpins]
        end
        % Add history entry
        sMatrixMat = bst_history('add', sMatrixMat, 'process', sprintf('Brain map spatial correlation for %s: ', sInputsA(iInputA).FileName));
        % === SAVE FILE ===
        % Output filename
        sStudy = bst_get('Study', sInputsA(1).iStudy);
        OutputFiles{end+1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_neuromaps');
        % Save file
        bst_save(OutputFiles{end}, sMatrixMat, 'v6');
        % Register in database
        db_add_data(sInputsA(iInputA).iStudy, OutputFiles{end}, sMatrixMat);
    end
    % Update whole tree
    panel_protocols('UpdateTree');
end

%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================

function [sMatrixMat] = CorrelationSurfaceMaps(ResultsFile, MapFiles, MapsSurfaceFile, nSpins)
    % No spinning case
    if nargin < 3 || isempty(nSpins) || nSpins < 1
        nSpins = 0;
    end
    % Prepare Surface files
    sResultsMat = in_bst_results(ResultsFile, 0, 'SurfaceFile');
    % Project if needed
    % TODO a better check needs to be done for warped surfaces
    if strcmp(sResultsMat.SurfaceFile, MapsSurfaceFile)
        sResultsProjFileName = '';
        sResultsProjMat = in_bst_results(ResultsFile, 1);
    else
        sResultsProjFileName = bst_project_sources({ResultsFile}, MapsSurfaceFile, 0, 0);
        sResultsProjFileName = sResultsProjFileName{1};
        [~, iStudyProj] = bst_get('ResultsFile', sResultsProjFileName);
        sResultsProjMat = in_bst_results(sResultsProjFileName, 1);
    end
    % Check time dimension for InputA
    TimesA = sResultsProjMat.Time;
    % Correlation values without spins
    corrValues = zeros(length(MapFiles), length(TimesA));
    % Correlation values without spins
    corrSpinValues = zeros(length(MapFiles), length(TimesA), nSpins);

    % Perform correlations between Results and each Map
    mapComments = cell(1, length(MapFiles));
    for iMap = 1 : length(MapFiles)
        % Load Map
        sOrgMapMat = in_bst_results(MapFiles{iMap}, 1);
        mapComments{iMap} = sOrgMapMat.Comment;
        % Comptue correlation for each sample in InputA
        for iTimeA = 1 : length(TimesA)
            if length(sOrgMapMat.Time) == length(TimesA)
                iTimeB = iTimeA;
            else
                iTimeB = 1;
            end
            % Spatial correlation with Map
            corrValues(iMap, iTimeA) = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeA)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeB)));

            % Spatial correlations with spinned Map
            if nSpins > 0
                % Get Subject with Maps' Surface
                [sSubject, iSubject] = bst_get('SurfaceFile', MapsSurfaceFile);
                defSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
                % Load Surface and remove any saved previous tess2tess interpolation
                % Because the registration sphere for the destination surface will change,
                % thus the test2testt_interp needs to be recomputed each time
                sMapSrfMat = in_tess_bst(MapsSurfaceFile);
                tmp.tess2tess_interp = [];
                bst_save(file_fullpath(MapsSurfaceFile), tmp, [], 1);
                % Create a copy to the map surface, this will be the changing sping surface
                rotSrfFileFull = strrep(file_fullpath(MapsSurfaceFile), '.mat', '_spin_test.mat');
                copyfile(file_fullpath(MapsSurfaceFile), rotSrfFileFull);
                spnSrfFile = file_short(rotSrfFileFull);
                iRotSrf = db_add_surface(iSubject, spnSrfFile, [sMapSrfMat.Comment, ' | Spin test']);
                % Spinning...
                for iSpin = 1 : nSpins
                    % Rotate the registration sphere in the original map surface, save it in the Spin test surface
                    sSpnSurfMat = in_tess_bst(MapsSurfaceFile);
                    % =================
                    % Rotate registration sphere (SCS coords) in copy of original surface
                    sSpnSurfMat.Reg.Sphere.Vertices(:,1) =  sSpnSurfMat.Reg.Sphere.Vertices(:,1); % x
                    sSpnSurfMat.Reg.Sphere.Vertices(:,2) =  sSpnSurfMat.Reg.Sphere.Vertices(:,2); % y
                    sSpnSurfMat.Reg.Sphere.Vertices(:,3) = -sSpnSurfMat.Reg.Sphere.Vertices(:,3); % z
                    % =================
                    % Update spin surface file
                    bst_save(file_fullpath(spnSrfFile), sSpnSurfMat);
                    % Project map from original surface to spinned surface
                    WmatSurf = tess_interp_tess2tess(MapsSurfaceFile, spnSrfFile, 0, 0);
                    % Need to clean tess2tess interpolation
                    tmp.tess2tess_interp = [];
                    bst_save(file_fullpath(MapsSurfaceFile), tmp, [], 1);
                    % Interpolation for 1 component
                    spnImageGridAmp = double(WmatSurf * sOrgMapMat.ImageGridAmp);
                    % Spatial correlation with Map
                    corrSpinValues(iMap, iTimeA, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeA)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeB)));
                end
                % Delete Spin test surface
                if (file_delete(rotSrfFileFull, 1) == 1)
                    % Remove from database
                    ProtocolSubjects = bst_get('ProtocolSubjects');
                    if iSubject == 0
                        ProtocolSubjects.DefaultSubject.Surface(iRotSrf) = [];
                    else
                        ProtocolSubjects.Subject(iSubject).Surface(iRotSrf) = [];
                    end
                    bst_set('ProtocolSubjects', ProtocolSubjects);
                    sSubject = bst_get('Subject', iSubject);
                    % Restore default cortex
                    ix = find(~cellfun(@isempty,(regexpi({sSubject.Surface.FileName}, defSurfaceFile))));
                    if ~isequal(sSubject.iCortex, ix)
                        db_surface_default(iSubject, 'Cortex', ix, 1);
                    end
                end
            end
        end
    end
    % Delete projected file and update study
    if ~isempty(sResultsProjFileName)
        file_delete(file_fullpath(sResultsProjFileName), 1);
        db_reload_studies(iStudyProj);
        % TODO Delete study if empty
    end
    % Create matrix structure
    sMatrixMat = db_template('matrixmat');
    sMatrixMat.Comment = 'Brain map corr';
    sMatrixMat.Time = TimesA;
    sMatrixMat.Value = corrValues;         % save [nMaps, nTimes]
    sMatrixMat.Description = mapComments'; % save [nMaps,1]
    % TODO: Do we need to save these values?
    if nSpins > 0
        sMatrixMat.Options.SpinTest = corrSpinValues;   % save [nMaps, nTimes, nSpins]
    end
end
