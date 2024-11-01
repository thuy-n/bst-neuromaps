% compare_parcellation_values.m
%
% Read parcellation values (surface and volume) obtained for different annotations
% using the Desikan-Killiany atlas. These values are obtained with
% Brainstorm and neuromaps:
%
%     - 'parcellation_dk_surface_brainstorm.csv'
%     - 'parcellation_dk_volume_brainstorm.csv'
%     - 'parcellation_dk_surface_neuromaps.csv'
%     - 'parcellation_dk_volume_neuromaps.csv'
%
% Created on 2024

t_bst_srf = readtable('./parcellation_dk_surface_brainstorm.csv', 'VariableNamingRule', 'preserve');
t_bst_vol = readtable('./parcellation_dk_volume_brainstorm.csv', 'VariableNamingRule', 'preserve');
t_nmp_srf = readtable('./parcellation_dk_surface_neuromaps.csv', 'VariableNamingRule', 'preserve');
t_nmp_vol = readtable('./parcellation_dk_volume_neuromaps.csv', 'VariableNamingRule', 'preserve');

% Check headers
if (~isequal(t_bst_srf.Properties.VariableDescriptions(2:end), t_bst_vol.Properties.VariableDescriptions(2:end)) || ...
    ~isequal(t_nmp_srf.Properties.VariableDescriptions(2:end), t_nmp_vol.Properties.VariableDescriptions(2:end))  || ...
    ~isequal(t_bst_srf.Properties.VariableDescriptions(2:end), t_nmp_srf.Properties.VariableDescriptions(2:end)))
    error('Headers are not the same')
end

% Check rows
if (~isequal(t_bst_srf.Row, t_bst_vol.Row) || ...
    ~isequal(t_nmp_srf.Row, t_nmp_vol.Row)  || ...
    ~isequal(t_bst_srf.Row, t_nmp_srf.Row))
    error('Row names are not the same')
end

%% Compute correlations
annot_names = t_bst_srf.Properties.VariableNames;
annot_names = setdiff(annot_names, {'Row'});
n_annot = length(annot_names);

bst_srf_x_bst_vol = zeros(n_annot, 1);
nmp_srf_x_nmp_vol = zeros(n_annot, 1);
bst_srf_x_nmp_srf = zeros(n_annot, 1);
bst_vol_x_nmp_vol = zeros(n_annot, 1);

for i_annot = 1 : n_annot
    annot_name = annot_names{i_annot};
    bst_srf = t_bst_srf.(annot_name);
    bst_vol = t_bst_vol.(annot_name);
    nmp_srf = t_nmp_srf.(annot_name);
    nmp_vol = t_nmp_vol.(annot_name);
    % Brainstorm: vol vs srf
    bst_srf_x_bst_vol(i_annot) = corr(bst_srf, bst_vol);
    % neuromaps: vol vs srf
    nmp_srf_x_nmp_vol(i_annot) = corr(nmp_srf, nmp_vol);
    % Surface: Brainstorm vs neuromaps
    bst_srf_x_nmp_srf(i_annot) = corr(bst_srf, nmp_srf);
    % Volume: Brainstorm vs neuromaps
    bst_vol_x_nmp_vol(i_annot) = corr(bst_vol, nmp_vol);
end

corrs = [bst_srf_x_bst_vol, nmp_srf_x_nmp_vol, bst_srf_x_nmp_srf, bst_vol_x_nmp_vol];




