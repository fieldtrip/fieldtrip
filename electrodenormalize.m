function [norm] = electrodenormalize(cfg);

% ELECTRODENORMALIZE is deprecated, please use ELECTRODEREALIGN 

% Copyright (C) 2005-2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('ELECTRODENORMALIZE is deprecated, please use ELECTRODEREALIGN');

[norm] = electroderealign(cfg);

