function [headmodel] = ft_fetch_vol(cfg, data)

% FT_FETCH_VOL mimics the behaviour of FT_READ_VOL, but for a FieldTrip
% configuration instead of a file on disk.
%
% Use as
%   [headmodel] = ft_fetch_vol(cfg)
% where you should specify the volume conductor model with
%   cfg.headmodel     = structure with volume conduction model or string with filename
%
% See also FT_READ_VOL, FT_FETCH_DATA

% Copyright (C) 2011, J?rn M. Horschig
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% check input arguments
if nargin > 1 
  ft_error('something is wrong');
end

% get the headmodel definition/volume conduction model
if isfield(cfg, 'headmodel') && ischar(cfg.headmodel)
  fprintf('reading headmodel from file ''%s''\n', cfg.headmodel);
  headmodel = ft_read_vol(cfg.headmodel);
elseif isfield(cfg, 'headmodel') && (isstruct(cfg.headmodel) || isa(cfg.headmodel, 'config'))
  headmodel = cfg.headmodel;
elseif isfield(cfg, 'headmodel') && iscell(cfg.headmodel) 
  % this might represent combined EEG, ECoG and/or MEG
  headmodel = cfg.headmodel;
else
  ft_error('no headmodel specified');
end

% ensure that the headmodel description is up-to-date
headmodel = ft_datatype_headmodel(headmodel);
