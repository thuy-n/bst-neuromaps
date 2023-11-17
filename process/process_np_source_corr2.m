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
    % Load neuromaps plugin if needed
    PlugDesc = bst_plugin('GetDescription', 'neuromaps');
    if ~PlugDesc.isLoaded
        bst_plugin('Load', 'neuromaps');
    end

    % All files must have surface or volume sources
    sResultsB1 = in_bst_results(sInputsB(1).FileName, 0, 'HeadModelType', 'Time');
    % Verify correct source space
    sInputs = [sInputsA, sInputsB];
    for iInput = 1 : length([sInputsA, sInputsB])
        sResultsMat = in_bst_results(sInputs(iInput).FileName, 0, 'HeadModelType');
        if ~strcmpi(sResultsMat.HeadModelType, sResultsB1.HeadModelType)
            bst_report('Error', sProcess, sInputs(iInput), 'Input file have different source space.');
            return;
        end
    end
    % Verify time definition
    % FilesB must have the same time axis: TimesB
    for iInputB = 1 : length(sInputsB)
        sResultsMat = in_bst_results(sInputsB(iInputB).FileName, 0, 'Time');
        if ~isequal(sResultsMat.Time, sResultsB1.Time)
            bst_report('Error', sProcess, sInputsB(iInputB), 'Input files B must have the same time definition.');
            return;
        end
    end
    % Validate time dimensions for FilesA and FilesB
    % filesA (1 sample)  vs filesB (1 sample)   OK      1 Corr value
    % filesA (N samples) vs filesB (1 sample)   OK      N Corr values
    % filesA (N samples) vs filesB (N samples)  OK      N Corr values
    % filesA (M samples) vs filesB (N samples)  Not OK
    for iInputA = 1 : length(sInputsA)
        sResultsMat = in_bst_results(sInputsA(iInputA).FileName, 0, 'Time');
        % FilesA must have the same time axis as TimesB
        if length(sResultsB1.Time) > 2 && ~isequal(sResultsMat.Time, sResultsB1.Time)
            bst_report('Error', sProcess, sInputsA(iInputA), 'Input files A must have the same time axis as files B');
            return;
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
        % Concatenate matrices (add results from other maps)
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

function sMatrixMat = CorrelationSurfaceMaps(ResultsFile, MapFiles, MapsSurfaceFile, nSpins)
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
        % Check if Map is one time sample
        isOneTimeSampleMap = length(sOrgMapMat.Time) == 2 && isequal(sOrgMapMat.ImageGridAmp(:,1), sOrgMapMat.ImageGridAmp(:,2));
        % Comptue correlation for each sample in InputA
        for iTimeA = 1 : length(TimesA)
            % Not allow to correlate one sample source map with multiple samples source map
            if isOneTimeSampleA && ~isOneTimeSampleMap
                bst_error(sprintf('Source file %s must be one time sample.', MapFiles{iMap}));
                return
            end
            % Compute correlations for each sample point in InputA and its equivalent in Map
            if length(TimesA) == length(sOrgMapMat.Time)
                iTimeB = iTimeA;
                if isOneTimeSampleA && isOneTimeSampleMap && iTimeA == 2
                    continue
                end
            % Compute correlations for each sample point in InputA and sample 1 in Map
            else
                iTimeB = 1;
            end
            % Spatial correlation with Map
            corrValues(iMap, iTimeA) = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeA)));

            % Spatial correlations with spun Map
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
                % Create a copy to the map surface, this will be the changing spinning surface
                rotSrfFileFull = strrep(file_fullpath(MapsSurfaceFile), '.mat', '_spin_test.mat');
                copyfile(file_fullpath(MapsSurfaceFile), rotSrfFileFull);
                spnSrfFile = file_short(rotSrfFileFull);
                iRotSrf = db_add_surface(iSubject, spnSrfFile, [sMapSrfMat.Comment, ' | Spin test']);

                % Spinning...
                for iSpin = 1 : nSpins
                    % Rotate the registration spheres (L and R) in the original map surface, save it in the Spin test surface
                    sSpnSurfMat = in_tess_bst(MapsSurfaceFile);
                    % Get vertex indices for each hemisphere
                    [ir, il] = tess_hemisplit(sSpnSurfMat);
                    % Get coordinates of sphere center (L and R)
                    offset_coordinatesL= (max(sSpnSurfMat.Reg.Sphere.Vertices(il,:))+min(sSpnSurfMat.Reg.Sphere.Vertices(il,:)))/2;
                    offset_coordinatesR= (max(sSpnSurfMat.Reg.Sphere.Vertices(ir,:))+min(sSpnSurfMat.Reg.Sphere.Vertices(ir,:)))/2;
                    % Alexander-Bloch method applied opposite rotations over the X and Z axes (SCS coords)
                    % https://doi.org/10.1016/j.neuroimage.2018.05.070
                    I1 = eye(3,3);
                    I1(1,1)=-1;
                    % Random uniform sampling procedure, get rotation matrices TL and TR (for each sphere)
                    A = normrnd(0,1,3,3);
                    [TL, temp] = qr(A);
                    TL = TL * diag(sign(diag(temp)));
                    if(det(TL)<0)
                        TL(:,1) = -TL(:,1);
                    end
                    % Reflect across the X-Z plane (SCS coords) for right hemisphere
                    TR = I1 * TL * I1;
                    % Rotate spheres
                    sSpnSurfMat.Reg.Sphere.Vertices(il,:)= sSpnSurfMat.Reg.Sphere.Vertices(il,:) * TL;
                    sSpnSurfMat.Reg.Sphere.Vertices(ir,:)= sSpnSurfMat.Reg.Sphere.Vertices(ir,:) * TR;
                    % Get coordinates for new sphere center (L and R)
                    offset_coordinatesL2= (max(sSpnSurfMat.Reg.Sphere.Vertices(il,:))+min(sSpnSurfMat.Reg.Sphere.Vertices(il,:)))/2;
                    offset_coordinatesR2= (max(sSpnSurfMat.Reg.Sphere.Vertices(ir,:))+min(sSpnSurfMat.Reg.Sphere.Vertices(ir,:)))/2;
                    % Recenter new sphere to old sphere center (L and R)
                    sSpnSurfMat.Reg.Sphere.Vertices(il,:) = sSpnSurfMat.Reg.Sphere.Vertices(il,:)-offset_coordinatesL2 +offset_coordinatesL;
                    sSpnSurfMat.Reg.Sphere.Vertices(ir,:) = sSpnSurfMat.Reg.Sphere.Vertices(ir,:)-offset_coordinatesR2 +offset_coordinatesR;
                    % Update spin surface file
                    bst_save(file_fullpath(spnSrfFile), sSpnSurfMat);
                    % Project map from original surface to spun surface
                    WmatSurf = tess_interp_tess2tess(MapsSurfaceFile, spnSrfFile, 0, 0);
                    % Need to clean tess2tess interpolation
                    tmp.tess2tess_interp = [];
                    bst_save(file_fullpath(MapsSurfaceFile), tmp, [], 1);
                    % Interpolation for 1 component
                    spnImageGridAmp = double(WmatSurf * sOrgMapMat.ImageGridAmp);
                    % Spatial correlation with Map
                    corrSpinValues(iMap, iTimeA, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeA)));
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
    sMatrixMat.Value = corrValues;                     % save [nMaps, nTimes]
    sMatrixMat.Description = mapComments';             % save [nMaps,1]
    if nSpins > 0
        p_spin_values = computeSpinPvalue(corrValues, corrSpinValues, nSpins);
        sMatrixMat.Options.SpinTest = p_spin_values;   % save [nMaps, nTimes, nSpins]
    end
end


function p_spin_values = computeSpinPvalue(corrVal, corrSpinValues, nSpins)
    % Compute p-values for spin test
    p_spin_values= zeros(size(corrVal));
    for iMap = 1:size(corrSpinValues,1)
        for iTime = 1:size(corrSpinValues,2)
            if corrVal(iMap,iTime)> 0
                p_spin_values(iMap,iTime) = (sum(squeeze(corrSpinValues(iMap,iTime,:)) > corrVal(iMap,iTime))+1) ./ (nSpins+1);
            else % corrVal(j,k)<= 0
                p_spin_values(iMap,iTime) = (sum(squeeze(corrSpinValues(iMap,iTime,:)) < corrVal(iMap,iTime))+1) ./ (nSpins+1);
            end
        end
    end
end

