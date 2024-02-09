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
% Authors: Raymundo Cassani 2023-2024

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

    % Load neuromaps plugin if needed
    PlugDesc = bst_plugin('GetInstalled', 'neuromaps');
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
        if length(sResultsMat.Time) ~= length(sResultsB1.Time) || ...
           (length(sResultsMat.Time) > 2 && ( abs(sResultsMat.Time(1) - sResultsB1.Time(1)) > 1e-5) || (abs(sResultsMat.Time(end) - sResultsB1.Time(end)) > 1e-5 ))
            bst_report('Error', sProcess, sInputsB(iInputB), 'Input files B must have the same time definition.');
            return;
        end
    end
    % Validate time dimensions for FilesA and FilesB
    % filesA (1 sample)  vs filesB (1 sample)   OK      1 Corr value
    % filesA (N samples) vs filesB (1 sample)   OK      N Corr values
    % filesA (N samples) vs filesB (N samples)  OK      N Corr values
    % filesA (1 samples) vs filesB (N samples)  Not OK
    % filesA (M samples) vs filesB (N samples)  Not OK
    for iInputA = 1 : length(sInputsA)
        sResultsMat = in_bst_results(sInputsA(iInputA).FileName, 0, 'Time');
        if length(sResultsB1.Time) > 2 && ((abs(sResultsMat.Time(1) - sResultsB1.Time(1)) > 1e-5) || (abs(sResultsMat.Time(end) - sResultsB1.Time(end)) > 1e-5 ))
            bst_report('Error', sProcess, sInputsA(iInputA), 'Input files A must have the same time axis as files B');
            return;
        end
    end

    % Get maps and their surfaces from sInputB
    MapFiles = {sInputsB.FileName};
    MapSurfaceFiles = cell(1, length(MapFiles));
    for iInputB = 1 : length(sInputsB)
        sResultsMat = in_bst_results(sInputsB(iInputB).FileName, 0, 'SurfaceFile');
        MapSurfaceFiles{iInputB} = sResultsMat.SurfaceFile;
    end

    % Compute and save spatial correlations
    OutputFiles = CorrelationSurfaceMaps(sProcess, sInputsA, MapFiles, MapSurfaceFiles, 1);

    % Update whole tree
    panel_protocols('UpdateTree');
end

%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================
function OutputFiles = CorrelationSurfaceMaps(sProcess, sResultsInputs, MapFiles, MapSurfaceFiles, isInteractive)
    % Perform spatial correlations of N Source files and M Maps
    OutputFiles = {};

    % Valida inputs
    if nargin < 5 || isempty(isInteractive)
        isInteractive = 0;
    end
    nSpins = sProcess.options.nspins.Value{1};
    if isempty(nSpins) || nSpins < 1
        nSpins = 0;
    end
    if ischar(MapFiles)
        MapFiles = {MapFiles};
    end
    if ischar(MapSurfaceFiles)
        MapSurfaceFiles = {MapSurfaceFiles};
    end

    nResultsInputs = length(sResultsInputs);
    nMaps = length(MapFiles);

    % Map comments
    mapComments = cell(nMaps, 1);
    % Track progress
    infoStr  = 'Spatial correlations... %s: %%d/%%d, %s: %%d/%%d';
    switch func2str(sProcess.Function)
        case 'process_np_source_corr1'
            infoStr = sprintf(infoStr, 'File', 'Annotation');
        case 'process_np_source_corr2'
            infoStr = sprintf(infoStr, 'FileA', 'FileB');
    end
    if nSpins > 0
        infoStr = [infoStr, ', Spin: %d/%d'];
    end
    if isInteractive
        bst_progress('start', 'Processes', 'Spatial correlations...', 0, 100);
        barStep = 100 ./ (nResultsInputs * nMaps);
        if nSpins > 0
            barStep = barStep ./ nSpins;
        end
    end
    % Compute and save spatial correlations for one source file
    for iResultsInput = 1 : nResultsInputs
        % Input file and input study
        ResultsFile = sResultsInputs(iResultsInput).FileName;
        iStudy = sResultsInputs(iResultsInput).iStudy;
        % Check time dimension for Sources (InputA)
        sResultsMat = in_bst_results(ResultsFile, 0, 'Time', 'SurfaceFile');
        TimesA = sResultsMat.Time;
        % Initialize stat values for MapSurface correlation and their p-values without spins
        r_no_spin = zeros(nMaps, length(TimesA));
        p_no_spin = zeros(nMaps, length(TimesA));
        % Correlation and their p-values for with spin test
        r_spin_test = zeros(nMaps, length(TimesA), nSpins);
        p_spin_test = zeros(nMaps, length(TimesA));
        % Unique MapSurfaceFiles
        UniqueMapSurfaceFiles = unique(MapSurfaceFiles);
        for iUniqueSurfaceFile = 1 : length(UniqueMapSurfaceFiles)
            MapSurfaceFile = UniqueMapSurfaceFiles{iUniqueSurfaceFile};
            if nSpins > 0
                % Backup previous tess2tess interpolation in MapSurface (if any) as Spin process will overwrite it
                sMapSrfMat = in_tess_bst(MapSurfaceFile);
                MapTess2tessBackup = [];
                if isfield(sMapSrfMat, 'tess2tess_interp') && all(isfield(sMapSrfMat.tess2tess_interp, {'Signature', 'Wmat'})) &&  ...
                          ~isempty(sMapSrfMat.tess2tess_interp.Signature) && ~isempty(sMapSrfMat.tess2tess_interp.Wmat)
                    MapTess2tessBackup = sMapSrfMat.tess2tess_interp;
                end
            end
            % Sources and Maps must have same Surface (Project sources if needed)
            if strcmp(sResultsMat.SurfaceFile, MapSurfaceFile)
                sResultsProjFileName = '';
                sResultsProjMat = in_bst_results(ResultsFile, 1);
            else
                sResultsProjFileName = bst_project_sources({ResultsFile}, MapSurfaceFile, 0, 0);
                sResultsProjFileName = sResultsProjFileName{1};
                [~, iStudyProj] = bst_get('ResultsFile', sResultsProjFileName);
                sResultsProjMat = in_bst_results(sResultsProjFileName, 1);
            end
            % Check if is one (time) sample
            isOneSampleA = length(TimesA) == 2 && isequal(sResultsProjMat.ImageGridAmp(:,1), sResultsProjMat.ImageGridAmp(:,2));
            % Create SpinSurface from (MapSurface) if needed
            if nSpins > 0
                [sSubject, iSubjectMapSurface] = bst_get('SurfaceFile', MapSurfaceFile);
                defSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
                % Load Surface and remove any saved previous tess2tess interpolation
                % Because the registration sphere for the destination surface will change,
                % thus the test2test_interp needs to be recomputed each spin
                tmp.tess2tess_interp = [];
                bst_save(file_fullpath(MapSurfaceFile), tmp, [], 1);
                sMapSrfMat = in_tess_bst(MapSurfaceFile);
                % Create a copy to the map surface, this copy will be the changing spinning surface
                rotSrfFileFull = strrep(file_fullpath(MapSurfaceFile), '.mat', '_spin_test.mat');
                copyfile(file_fullpath(MapSurfaceFile), rotSrfFileFull);
                spnSrfFile = file_short(rotSrfFileFull);
                iRotSrf = db_add_surface(iSubjectMapSurface, spnSrfFile, [sMapSrfMat.Comment, ' | Spin test']);
            end
            % Indices of Maps using this Map Surface
            iMaps = find(strcmp(MapSurfaceFiles, MapSurfaceFile));
            for im = 1 : nMaps
                iMap = iMaps(im);
                % Load Map
                sOrgMapMat = in_bst_results(MapFiles{iMap}, 1);
                mapComments{iMap} = sOrgMapMat.Comment;
                % Check if Map is one time sample
                TimesB = sOrgMapMat.Time;
                isOneSampleMap = length(TimesB) == 2 && isequal(sOrgMapMat.ImageGridAmp(:,1), sOrgMapMat.ImageGridAmp(:,2));
                % Not allow to correlate one sample source map with multiple samples source map
                if isOneSampleA && ~isOneSampleMap
                    bst_error(sprintf('Brain annotation file %s must be one time sample.', MapFiles{iMap}));
                    return
                end
                % Spatial correlations, start with 0, as this is for no spin
                for iSpin = 0 : nSpins
                    for iTimeB = 1 : length(TimesB)
                        % Avoid recompute spatial correlation if one-time-sample B
                        if isOneSampleMap && iTimeB == 2
                            continue
                        end

                        % Spatial correlation for witn not spun map
                        if iSpin == 0
                            % Spatial correlation of all time samples in A with the one-time-sample B
                            if isOneSampleMap
                                [r_no_spin(iMap, :), p_no_spin(iMap, :)] = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp));
                            % Spatial correlation of one time sample in A with the same time sample in B
                            else
                                [r_no_spin(iMap, iTimeB), p_no_spin(iMap, iTimeB)] = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeB)));
                            end

                        % Spatial correlations with spun Map
                        else
                            % Spinning...
                            % Rotate the registration spheres (L and R) in the original map surface, save it in the Spin test surface
                            sSpnSurfMat = sMapSrfMat;
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
                            WmatSurf = tess_interp_tess2tess(MapSurfaceFile, spnSrfFile, 0, 0);
                            % Need to clean tess2tess interpolation
                            tmp.tess2tess_interp = [];
                            bst_save(file_fullpath(MapSurfaceFile), tmp, [], 1);
                            % Interpolation for 1 component
                            spnImageGridAmp = double(WmatSurf * sOrgMapMat.ImageGridAmp);
                            % Spatial correlation with Map
                            if isOneSampleMap
                                r_spin_test(iMap, :, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp));
                                %[r_no_spin(iMap, :), p_no_spin(iMap, :)] = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp));
                            else
                                r_spin_test(iMap, iTimeB, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeB)));
                                %[r_no_spin(iMap, iTimeB), p_no_spin(iMap, iTimeB)] = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeB)));
                            end
                            %r_spin_test(iMap, iTimeA, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeA)));
                        end
                    end
                    % Process indicator
                    if isInteractive
                        if nSpins == 0
                            bst_progress('text', sprintf(infoStr, iResultsInput, nResultsInputs, iMap, nMaps));
                            bst_progress('inc', barStep);
                        elseif nSpins > 1 && iSpin > 0
                            bst_progress('text', sprintf(infoStr, iResultsInput, nResultsInputs, iMap, nMaps, iSpin, nSpins));
                            bst_progress('inc', barStep);
                        end
                    else
                        if nSpins == 0
                            fprintf(1, [infoStr, '\n'], iResultsInput, nResultsInputs, iMap, nMaps);
                        elseif nSpins > 1 && iSpin > 0
                            fprintf(1, [infoStr, '\n'], iResultsInput, nResultsInputs, iMap, nMaps, iSpin, nSpins);
                        end
                    end
                end
            end

            % === REMOVE TMP FILES ===
            if isInteractive
                bst_progress('text', 'Spatial correlations... Delete temporary files');
            end
            % Delete Spin test surface from MapSurfaceSubject
            if nSpins > 0 && (file_delete(rotSrfFileFull, 1) == 1)
                % Remove from database
                ProtocolSubjects = bst_get('ProtocolSubjects');
                if iSubjectMapSurface == 0
                    ProtocolSubjects.DefaultSubject.Surface(iRotSrf) = [];
                else
                    ProtocolSubjects.Subject(iSubjectMapSurface).Surface(iRotSrf) = [];
                end
                bst_set('ProtocolSubjects', ProtocolSubjects);
                % Restore default cortex
                sSubject = bst_get('Subject', iSubjectMapSurface);
                ic = find(strcmp({sSubject.Surface.FileName}, defSurfaceFile));
                if ~isequal(sSubject.iCortex, ic)
                    db_surface_default(iSubjectMapSurface, 'Cortex', ic, 1);
                end
            end
            % Restore original tess2tess interpolation in Map surface
            if (nSpins > 0) && ~isempty(MapTess2tessBackup)
                tmp.tess2tess_interp = MapTess2tessBackup;
                bst_save(file_fullpath(MapSurfaceFile), tmp, [], 1);
            end
            % Delete projected file and update study
            if ~isempty(sResultsProjFileName)
                file_delete(file_fullpath(sResultsProjFileName), 1);
                db_reload_studies(iStudyProj);
                % Delete Study if it is empty
                sStudy = bst_get('Study', iStudyProj);
                fieldsStudyCheck = {'Channel', 'Data', 'HeadModel', 'Result', 'Stat', ...
                                    'Image', 'NoiseCov', 'Dipoles', 'Timefreq','Matrix'};
                if all(cellfun(@(x) isempty(sStudy.(x)), fieldsStudyCheck))
                    db_delete_studies(iStudyProj);
                end
            end
        end

        % === SAVE CORRELATION VALUES  ===
        if isInteractive
            bst_progress('text', 'Spatial correlations... Saving correlation values');
        end
        % Compute p-values for spin test
        if nSpins > 0
            p_spin_test = computeSpinPvalue(r_no_spin, r_spin_test, nSpins);
        end
        % Create statmat structure
        sStatMat = db_template('statmat');
        sStatMat.Type    = 'matrix';
        sStatMat.Comment = 'Brain map corr';
        if nSpins > 0
            sStatMat.Comment = [sStatMat.Comment, ' | spintest'];
        end
        sStatMat.Time = TimesA;
        sStatMat.Description = mapComments;             % [nMaps,1]
        sStatMat.tmap        = r_no_spin;               % [nMaps, nTimes]
        if nSpins > 0
            sStatMat.Options.nSpins      = nSpins;
            sStatMat.pmap                = p_spin_test; % [nMaps, nTimes]
            sStatMat.Options.rSpinTest   = r_spin_test; % [nMaps, nTimes, nSpins]
            sStatMat.Options.pNoSpinTest = p_no_spin;   % [nMaps, nTimes]
        else
            sStatMat.pmap    = p_no_spin;               % [nMaps, nTimes]
        end
        sStatMat.Correction   = 'no';
        sStatMat.ChannelFlag  = [];
        sStatMat.ColormapType = 'stat2';
        sStatMat.DisplayUnits = 'correlation';
        % Add history entry
        sStatMat = bst_history('add', sStatMat, 'process', sprintf('Spatial correlation for: %s', ResultsFile));
        % Save file
        sStudy = bst_get('Study', iStudy);
        OutputFiles{end+1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'pmatrix_neuromaps');
        % Save file
        bst_save(OutputFiles{end}, sStatMat, 'v6');
        % Register in database
        db_add_data(iStudy, OutputFiles{end}, sStatMat);
    end
    if isInteractive
        bst_progress('set', 100);
        bst_progress('stop');
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
