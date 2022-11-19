function [trl, event] = ft_trialfun_hed(cfg)

% FT_TRIALFUN_HED is a trial function that can be used with HED tags. It demonstrates
% some basic functionality for selecting specific events, but mainly serves as an
% example or template for your own conditial trial definitions. For that you would
% copy this function and giuve it your own name, e.g. FT_TRIALFUN_MYEXPERIMENT.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset           = string with the filename
%   cfg.trialfun          = 'ft_trialfun_hed' % or your own copy
%
% The selection of events and timing of the epochs is specified with
%   cfg.trialdef.regexp     = regular expression that is applied to the HED tags
%   cfg.trialdef.prestim    = number, in seconds
%   cfg.trialdef.poststim   = number, in seconds
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL, FT_TRIALFUN_EXAMPLE1,
% FT_TRIALFUN_EXAMPLE2

% read the header information
hdr = ft_read_header(cfg.dataset);

% determine the number of samples before and after the event
prestim  = -round(cfg.trialdef.prestim  * hdr.Fs);
poststim =  round(cfg.trialdef.poststim * hdr.Fs);

% rather than using FT_READ_EVENT which has a limited representation, we assume that
% the data is formatted according to BIDS and includes HED tags

tsvfile  = bids_sidecar(cfg.dataset, 'events', 'tsv');
jsonfile = bids_sidecar(cfg.dataset, 'events', 'json');

tsv = ft_read_tsv(tsvfile);
if ~isempty(jsonfile)
  json = ft_read_json(jsonfile);
else
  json = [];
end

% this will be returned as the complete set of events without selection
event = tsv;

% make a selection based on the HED tags
if any(strcmp(fieldnames(tsv), 'HED'))
  hedtags = tsv.HED;
elseif ~isempty(jsonfile)
  warning('assembling the HED tags on the fly');
  % this requires functions from an external toolbox, see https://hed-examples.readthedocs.io/en/latest/HedMatlabTools.html
  ft_hastoolbox('hedtools', 1);
  tsv = hed_assemble(tsvfile, jsonfile);
  hedtags = tsv.HED;
else
  error('there are no HED tags')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code uses a regular expression to select events
% you may want to change this and use another mechanism to define the trials/epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the HED tags that match the regular expression
sel = find(~cellfun(@isempty, regexp(hedtags, cfg.trialdef.regexp)));

% these are the first three required columnss
trl = table;
trl.begsample = tsv.sample(sel) + prestim;
trl.endsample = tsv.sample(sel) + poststim - 1;
trl.offset = repmat(prestim, size(trl.begsample));

% concatenate all other details to the table
trl = cat(2, trl, tsv(sel,:));
