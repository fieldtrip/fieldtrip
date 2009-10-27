function [grid] = precompute_leadfield(cfg, data)

% PRECOMPUTE_LEADFIELD is deprecated, please use PREPARE_LEADFIELD

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('PRECOMPUTE_LEADFIELD is deprecated, please use PREPARE_LEADFIELD');

if nargin==1
  [grid] = prepare_leadfield(cfg);
else
  [grid] = prepare_leadfield(cfg, data);
end
