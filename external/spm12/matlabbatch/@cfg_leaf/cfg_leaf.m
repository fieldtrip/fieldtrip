function item = cfg_leaf(varargin)

% This is currently only a "marker" class that should be inherited by all
% leaf classes. It does not add fields and does not have methods.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_leaf.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

item = class(struct('unused',[]), mfilename);