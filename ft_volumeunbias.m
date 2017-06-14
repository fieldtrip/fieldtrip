function mri_unbias = ft_volumeunbias(cfg, mri)

% FT_VOLUMEUNBIAS corrects the image inhomogeneity bias in an anatomical MRI
%
% Use as
%   mri_unbias = ft_volumeunbias(cfg, mri)
%
% See also FT_VOLUMEREALIGN

ft_revision = '$Id$';
ft_nargin = nargin;
ft_nargout = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init
ft_preamble debug
ft_preamble loadvar     mri
ft_preamble provenance  mri
ft_preamble trackconfig

if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', {'volume'});

% set the defaults
cfg.spmversion       = ft_getopt(cfg, 'spmversion', 'spm8');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');

if ~isfield(cfg, 'intermediatename')
  cfg.intermediatename = tempname;
end

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

% check whether the input has an anatomy
if ~isfield(mri, 'anatomy')
  error('no anatomical information available, this is required for bias correction');
end

% do an approximate alignment
mri_acpc = ft_convert_coordsys(mri, 'spm');

% create an spm-compatible file for the anatomical volume data
V = ft_write_mri([cfg.intermediatename '_anatomy.img'], mri_acpc.anatomy, 'transform', mri_acpc.transform, 'spmversion', cfg.spmversion);

flags = struct();
% flags.nbins = 4;
% flags.cutoff = 10;
T = spm_bias_estimate(V, flags);
VO = spm_bias_apply(V,T);

mri_unbias = keepfields(mri, {'coordsys', 'dim', 'transform', 'unit', 'inside'});
mri_unbias.anatomy = VO.dat;

% if strcmp(cfg.keepintermediate, 'no')
% % remove the intermediate files
% for flop=1:length(files)
% [p, f, x] = fileparts(files{flop});
% delete(fullfile(p, [f, '.*']));
% [p, f, x] = fileparts(wfiles{flop});
% delete(fullfile(p, [f, '.*']));
% end
% end

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance mri_unbias
ft_postamble history    mri_unbias
ft_postamble savevar    mri_unbias
