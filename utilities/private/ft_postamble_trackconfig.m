% FT_POSTAMBLE_TRACKCONFIG

% otherwise the empty field would end up in the output cfg
global ft_default
ft_default = rmfield(ft_default, 'postamble');

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
if isfield(cfg, 'outputfile')
  cfg.outputfile;
end
if isfield(cfg, 'outputlock')
  cfg.outputlock;
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% the calling ft_postamble expects it to be present
ft_default.postamble = {};
