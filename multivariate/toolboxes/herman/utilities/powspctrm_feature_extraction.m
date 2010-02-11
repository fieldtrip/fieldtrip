
function [feature_data] = powspctrm_feature_extraction(freq_data,cfg_sets)

% POWSPCTRM_FEATURE_EXTRACTION supports feature extraction from Fieldtrip
% data through the sequential use of PREPARE_TIMEFREQ_DATA with different
% cfg settings. It allows for example for channel-independent selection of
% frequency bands.
%
%  Use as 
%   [feature_data] = powspctrm_feature_extraction(freq_data,[cfg1 cfg2 ..  cfgn])
%   [feature_data] = powspctrm_feature_extraction(freq_data,{cfg1 cfg2 ..  cfgn})
%
%  INPUT 
%       freq_data    - input data extracted FREQANALYSIS
%       cfg_sets     - set of cfg structures - array or cell array (each cfg
%                      structure is the same as the one used with PREPARE_TIMEFREQ_DATA)
%
%  OUTPUT
%       feature_data - the resulting data set - structure with four fields
%                       .biol   -> two-dimensional feature set
%                       .design -> design matrix (including labels)
%                       .cfg    -> array of accumulated cfgs (cfg1..cfgn)
%                       .ncol   -> array containing the number of variables
%                                  (channels, feature components etc.) 
%                                  extracted through each cfg
%
%  SEE ALSO
%       prepare_timefreq_data.m
%

% Pawel Herman, 2009

feature_data.cfg = [];
feature_data.biol = [];
feature_data.design = [];
feature_data.ncol = [];
Ncol = 0;

for i=1:length(cfg_sets)
    if iscell(cfg_sets),  cfg = cfg_sets{i}; else cfg = cfg_sets(i); end
    [cfg, features_design] = prepare_timefreq_data(cfg, freq_data{:});
    biol = features_design.biol;
    design = features_design.design;
    features_design = rmfield(features_design,'biol');
    features_design = rmfield(features_design,'design');
    features_design.cfg.featurecol = Ncol+1:Ncol+1+size(biol,2);
    features_design.freq = cfg.frequency;
    feature_data.cfg = [feature_data.cfg features_design];
    feature_data.biol = [feature_data.biol biol];
    feature_data.design = design;
    feature_data.ncol = [feature_data.ncol size(biol,2)];
    Ncol = size(feature_data.biol,2);
end



