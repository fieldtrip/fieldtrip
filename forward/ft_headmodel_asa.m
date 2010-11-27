function vol = ft_headmodel_asa(filename, varargin)

% FT_HEADMODEL_ASA
% 
% Use as
%   vol = ft_headmodel_asa(filename)

% read the headmodel from file
% this works for ASA version 2.x, perhaps also for ASA 3.x 
vol = ft_read_vol(filename);

