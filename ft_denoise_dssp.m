function [dataout] = ft_denoise_dssp(cfg, datain)

% FT_DENOISE_DSSP 
%
% Use as 
%   dataout = ft_denoise_dssp(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.grid     = structure, source model with precomputed leadfields (see FT_PREPARE_LEADFIELD)
%   cfg.order    = integer, order of the source space
%
% See also FT_DENOISE_PCA, FT_DENOISE_SYNTHETIC, FT_DENOISE_TSR

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain
ft_preamble trackconfig      

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% get the options
cfg.grid  = ft_getopt(cfg, 'grid');
cfg.order = ft_getopt(cfg, 'order', 42);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do your stuff...
dataout = [];

% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug              
ft_postamble trackconfig        
ft_postamble previous   datain  
ft_postamble provenance dataout 
ft_postamble history    dataout 
ft_postamble savevar    dataout 

