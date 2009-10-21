function [norm] = electrodenormalize(cfg);

% ELECTRODENORMALIZE is deprecated, please use ELECTRODEREALIGN 

% Copyright (C) 2005-2006, Robert Oostenveld
%
% $Log: electrodenormalize.m,v $
% Revision 1.11  2006/09/13 07:20:06  roboos
% renamed electrodenormalize to electroderealign, added "deprecated"-warning to the old function
%
% Revision 1.10  2006/09/13 07:09:24  roboos
% Implemented support for cfg.method=interactive, using GUI for specifying and showing transformations. Sofar only for electrodes+headsurface.
%
% Revision 1.9  2006/09/12 15:26:06  roboos
% implemented support for aligning electrodes to the skin surface, extended and improved documentation

warning('ELECTRODENORMALIZE is deprecated, please use ELECTRODEREALIGN');

[norm] = rejectvisual(cfg);

