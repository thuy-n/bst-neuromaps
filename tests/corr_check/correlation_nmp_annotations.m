% compute_correlation_brainstorm.m
%
% Computes the correlations among all the annotations in the bst-neuromaps plugin
%
%
%
% Created on 2024

%% Create protocol
ProtocolName = 'bst-neuromaps-correlation';
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% List of maps
brainmapListSrf = process_nmp_fetch_maps('GetBrainMapsList', 'surface');

% Process: Fetch brain annotations from neuromaps
sMapFiles = bst_process('CallProcess', 'process_nmp_fetch_maps', [], [], ...
    'sspace',        'surface', ...  % Surface
    'brainmaps_srf', brainmapListSrf, ...
    'brainmaps_vol', '');

% Process: Spatial correlation
sStatFiles = bst_process('CallProcess', 'process_nmp_source_corr1', sMapFiles, [], ...
    'sspace',        'surface', ...  % Surface
    'brainmaps_srf', brainmapListSrf, ...
    'brainmaps_vol', '', ...
    'corrmetric',    'Pearson', ...  % Pearson corr
    'removezeros',   1, ...
    'nspins',        0);

% Keep all correlation values, statistifcally significant or not
% Process: Apply statistic threshold: alpha=1
sCorrFiles = bst_process('CallProcess', 'process_extract_pthresh', sStatFiles, [], ...
    'pthresh',    1, ...
    'durthresh',  0, ...
    'correction', 1, ...  % Uncorrected
    'control1',   1, ...
    'control2',   1, ...
    'control3',   1);

% Process: Concatenate time
sCorrFile = bst_process('CallProcess', 'process_concat', sCorrFiles, []);

% Process: Delete selected files
sFiles = bst_process('CallProcess', 'process_delete', [sStatFiles, sCorrFiles], [], ...
    'target', 1);  % Delete selected files
% Process: Move files: neuromaps/NewFolder
sCorrFile = bst_process('CallProcess', 'process_movefile', sCorrFile, [], ...
    'subjectname', sCorrFile.SubjectName, ...
    'folder',      'all_vs_all');

% Save file as NxN connectivity
sMatrixMat = in_bst_matrix(sCorrFile.FileName);
% Process: Delete selected files
sFiles = bst_process('CallProcess', 'process_delete', sCorrFile, [], ...
    'target', 1);  % Delete selected files
% Remove duplicated samples
ix = [1 : 2 : length(sMatrixMat.Time)];
R = sMatrixMat.Value(:,ix);
% Create file structure
sCnxMat = db_template('timefreqmat');
sCnxMat.TF = process_compress_sym('Compress', R(:));
sCnxMat.Comment  = 'Corr: maps';
sCnxMat.Time     = [0, 1];
sCnxMat.DataType = 'matrix';
sCnxMat.RowNames = sMatrixMat.Description';
sCnxMat.RefRowNames = sCnxMat.RowNames;
sCnxMat.Method = 'corr';
sCnxMat.Measure = 'other';
% Output filename
NewFile = bst_process('GetNewFilename', bst_fileparts(sCorrFile.FileName), 'timefreq_connectn_corr');
bst_save(NewFile, sCnxMat, 'v6');
% Add file to database structure
db_add_data(sCorrFile.iStudy, NewFile, sCnxMat);
panel_protocols('UpdateNode', 'Study', sCorrFile.iStudy);
panel_protocols('SelectNode', [], NewFile);

