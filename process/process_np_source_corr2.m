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
    if isempty(nSpins) || nSpins < 1
        nSpins = 0;
    end
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
        % TODO, better way to detect same time axis
        if length(sResultsB1.Time) > 2 && (max(abs(sResultsMat.Time - sResultsB1.Time)) > diff(sResultsMat.Time(1:2)))
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

    % Progress bar
    bst_progress('start', 'Processes', 'Computing spatial correlation...', 0, 100);
    % For each FileA compute correlations with all FilesB
    for iInputA = 1 : length(sInputsA)
        % Update progress bar
        bst_progress('text', sprintf('Processing file #%d/%d', iInputA, length(sInputsA)));
        % Compute and save spatial correlations
        OutputFiles = CorrelationSurfaceMaps(sInputsA(iInputA).FileName, MapFiles, MapSurfaceFiles, nSpins);
        % Update progress bar
        bst_progress('set', iInputA ./ length(sInputsA));
    end
    % Update whole tree
    panel_protocols('UpdateTree');
end

%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================
function OutputFiles = CorrelationSurfaceMaps(ResultsFiles, MapFiles, MapSurfaceFiles, nSpins)
    % Perform spatial correlations of N Source files and M Maps

    OutputFiles = {};
    % Map comments
    mapComments = cell(length(MapFiles), 1);
    % Compute and save spatial correlations for one source file
    for iResultFile = 1 : length(ResultsFiles)
        ResultsFile = ResultsFiles{iResultFile};
        % Output study
        [~, iStudy] = bst_get('ResultsFile', ResultsFile);
        % Check time dimension for Sources (InputA)
        sResultsMat = in_bst_results(ResultsFile, 0, 'Time', 'SurfaceFile');
        TimesA = sResultsMat.Time;
        % Initialize stat values for MapSurface correlation and their p-values without spins
        r_no_spin = zeros(length(MapFiles), length(TimesA));
        p_no_spin = zeros(length(MapFiles), length(TimesA));
        % Correlation and their p-values for with spin test
        r_spin_test = zeros(length(MapFiles), length(TimesA), nSpins);
        p_spin_test = zeros(length(MapFiles), length(TimesA));

        % Unique MapSurfaceFiles
        UniqueMapSurfaceFiles = unique(MapSurfaceFiles);
        for iUniqueSurfaceFile = 1 : length(UniqueMapSurfaceFiles)
            MapSurfaceFile = UniqueMapSurfaceFiles{iUniqueSurfaceFile};
            % Backup previous tess2tess interpolation in MapSurface (if any) as Spin process will overwrite it
            if nSpins > 0
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
            % Indices of Maps using this Map Surface
            iMaps = find(strcmp(MapSurfaceFiles, MapSurfaceFile));
            for im = 1 : length(iMaps)
                iMap = iMaps(im);
    %             if iMap == 1
    %                 pBarParams = bst_progress('getbarparams');
    %                 strBase = pBarParams.Msg;
    %             end
    %             bst_progress('text', [strBase, sprintf(', map #%d/%d', iMap, length(MapFiles))]);
                % Load Map
                sOrgMapMat = in_bst_results(MapFiles{iMap}, 1);
                mapComments{iMap} = sOrgMapMat.Comment;
                % Check if Map is one time sample
                isOneSampleMap = length(sOrgMapMat.Time) == 2 && isequal(sOrgMapMat.ImageGridAmp(:,1), sOrgMapMat.ImageGridAmp(:,2));
                % Comptue correlation for each sample in InputA
                for iTimeA = 1 : length(TimesA)
                    % Not allow to correlate one sample source map with multiple samples source map
                    if isOneSampleA && ~isOneSampleMap
                        bst_error(sprintf('Source file %s must be one time sample.', MapFiles{iMap}));
                        return
                    end
                    % Compute correlations for each sample point in InputA and its equivalent in Map
                    if length(TimesA) == length(sOrgMapMat.Time)
                        iTimeB = iTimeA;
                        if isOneSampleA && isOneSampleMap && iTimeA == 2
                            continue
                        end
                    % Compute correlations for each sample point in InputA and sample 1 in Map
                    else
                        iTimeB = 1;
                    end
                    % Spatial correlation with Map
                    [r_no_spin(iMap, iTimeA), p_no_spin(iMap, iTimeA)] = bst_corrn(transpose(sOrgMapMat.ImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeA)));
                    % Spatial correlations with spun Map
                    if nSpins > 0
                        % String with backspaces
    %                     bckstr = '';
                        % Get Subject with MapSurface
                        [sSubject, iSubject] = bst_get('SurfaceFile', MapSurfaceFile);
                        defSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
                        % Load Surface and remove any saved previous tess2tess interpolation
                        % Because the registration sphere for the destination surface will change,
                        % thus the test2test_interp needs to be recomputed each spin
                        sMapSrfMat = in_tess_bst(MapSurfaceFile);
                        tmp.tess2tess_interp = [];
                        bst_save(file_fullpath(MapSurfaceFile), tmp, [], 1);
                        % Create a copy to the map surface, this copy will be the changing spinning surface
                        rotSrfFileFull = strrep(file_fullpath(MapSurfaceFile), '.mat', '_spin_test.mat');
                        copyfile(file_fullpath(MapSurfaceFile), rotSrfFileFull);
                        spnSrfFile = file_short(rotSrfFileFull);
                        iRotSrf = db_add_surface(iSubject, spnSrfFile, [sMapSrfMat.Comment, ' | Spin test']);
                        % Spinning...
                        for iSpin = 1 : nSpins
                            % Print indicator of nSpin
    %                         nbytes = fprintf('%sSpin: %d of %d\n', bckstr, iSpin, nSpins);
    %                         bckstr = repmat(char(8), 1, nbytes - length(bckstr));
                            % Rotate the registration spheres (L and R) in the original map surface, save it in the Spin test surface
                            sSpnSurfMat = in_tess_bst(MapSurfaceFile);
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
                            r_spin_test(iMap, iTimeA, iSpin) = bst_corrn(transpose(spnImageGridAmp(:,iTimeB)), transpose(sResultsProjMat.ImageGridAmp(:,iTimeA)));
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
                            ic = find(~cellfun(@isempty,(regexpi({sSubject.Surface.FileName}, defSurfaceFile))));
                            if ~isequal(sSubject.iCortex, ic)
                                db_surface_default(iSubject, 'Cortex', ic, 1);
                            end
                        end
                    end
                end
    %             bst_progress('inc', (iMap ./ length(MapFiles)) * valMapBar);
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
        % Copy correlation values for second sample
        if isOneSampleA && isOneSampleMap
            r_no_spin(:,2)     = r_no_spin(:,1);
            p_no_spin(:,2)     = p_no_spin(:,1);
            r_spin_test(:,2,:) = r_spin_test(:,1,:);
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
