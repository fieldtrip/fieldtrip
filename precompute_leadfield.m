function [grid] = precompute_leadfield(cfg, data)

% PRECOMPUTE_LEADFIELD is deprecated, please use PREPARE_LEADFIELD

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: precompute_leadfield.m,v $
% Revision 1.9  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

warning('PRECOMPUTE_LEADFIELD is deprecated, please use PREPARE_LEADFIELD');

if nargin==1
  [grid] = prepare_leadfield(cfg);
else
  [grid] = prepare_leadfield(cfg, data);
end
