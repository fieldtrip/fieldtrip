function [ama] = loadama(filename)

% LOADAMA read an inverted A-matrix and associated geometry information
% from an ama file that was written by Tom Oostendorp's DIPOLI
%
% Use as
%   [ama] = loadama(filename)
%
% See also LOADTRI, LOADMAT

% Copyright (C) 2005, Robert Oostenveld
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

fid = fopen(filename, 'rb', 'ieee-le');

version = fread(fid, 1, 'int');
if version~=10
  ft_error(sprintf('%s is either not an inverted A matrix, or one of an old version', filename));
end

mode = fread(fid, 1, 'int');
ngeo = fread(fid, 1, 'int');

totpnt = 0;
tottri = 0;
nrow   = 0;

% read the boundaries
geo  = [];
for i=1:ngeo
  geo(i).name    = char(fread(fid, [1 80], 'uchar'));
  geo(i).npos    = fread(fid, 1, 'int');
  geo(i).pos     = fread(fid, [3 geo(i).npos], 'float')';
  geo(i).ntri    = fread(fid, 1, 'int');
  geo(i).tri     = fread(fid, [3 geo(i).ntri], 'int')' + 1;  % Matlab indexing starts at 1
  geo(i).sigmam  = fread(fid, 1, 'float');
  geo(i).sigmap  = fread(fid, 1, 'float');
  geo(i).geocon  = fread(fid, ngeo, 'int');
  geo(i).deflat  = fread(fid, ngeo, 'float');
  totpnt = totpnt + geo(i).npos;
  tottri = tottri + geo(i).ntri;
end

% read the electrodes
if mode~=1
  elec.name    = char(fread(fid, [1 80], 'uchar'));
  elec.npos    = fread(fid, 1, 'int');
  for i=1:(elec.npos+1)
    elec.el(i).tri  = fread(fid, 1, 'int') + 1; % Matlab indexing starts at 1
    elec.el(i).la   = fread(fid, 1, 'float');
    elec.el(i).mu   = fread(fid, 1, 'float');
    elec.el(i).name = char(fread(fid, [1 10], 'char'));
    % the ELECTRODE c-structure is padded to word boundaries, i.e. to 4 bytes
    dum = fread(fid, 2, 'char');
  end
  elec.vertex  = fread(fid, 1, 'int');
  elec.surface = fread(fid, 1, 'int');
  nrow = nrow + elec.npos;
else
  elec = [];
end

% read the gradiometers
if mode~=0
  ft_error('gradiometers not yet implemented');
else
  grad = [];
end

% read the inverted A-matrix
bi = fread(fid, [totpnt nrow], 'float')';

% read the isolated source compartment information, if present
iso_sur    = fread(fid, 1, 'int') + 1;  % Matlab indexing starts at 1
inner_only = fread(fid, 1, 'int');
if iso_sur~=0
  iso_totpnt = geo(iso_sur).npos;
  iso_b      = fread(fid, [iso_totpnt iso_totpnt], 'float')';
else
  iso_b = [];
end

fclose(fid);

% put all local variables into a structure, this is a bit unusual programming style
% the output structure is messy, but contains all relevant information
tmp = whos;
ama = [];
for i=1:length(tmp)
  if isempty(strmatch(tmp(i).name, {'tmp', 'fid', 'ans', 'handles'}))
    ama = setfield(ama, tmp(i).name, eval(tmp(i).name));
  end
end

