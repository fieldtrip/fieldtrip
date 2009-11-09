function flip = spm_flip_analyze_images
% Do Analyze format images need to be left-right flipped?
%-----------------------------------------------------------------------
% @(#)spm_flip_analyze_images.m	2.4 03/05/07

global defaults
if isempty(defaults) | ~isfield(defaults,'analyze') |...
         ~isfield(defaults.analyze,'flip')
	warning('Cant get default Analyze orientation - assuming flipped');
	flip = 1;
	return;
end;
flip = defaults.analyze.flip;
