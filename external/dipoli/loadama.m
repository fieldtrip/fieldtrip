function [ama] = loadama(filename);

% LOADAMA read an inverted A-matrix and associated geometry information
% from an ama file that was written by Tom Oostendorp's DIPOLI
%
% Use as
%   [ama] = loadama(filename)
%
% See also LOADTRI, LOADMAT

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: loadama.m,v $
% Revision 1.1  2008/12/24 10:25:41  roboos
% cleaned up the dipoli wrapper, replaced the binary by a better one and added a copy of the helper functions (from fileio)
%
% Revision 1.2  2008/11/12 17:02:03  roboos
% explicitely specify ieee-le in fopen()
%
% Revision 1.1  2005/11/16 13:50:51  roboos
% new implementation
%

fid = fopen(filename, 'rb', 'ieee-le');

version = fread(fid, 1, 'int');
if version~=10
  error(sprintf('%s is either not an inverted A matrix, or one of an old version', filename));
end

mode = fread(fid, 1, 'int');
ngeo = fread(fid, 1, 'int');

totpnt = 0;
totdhk = 0;
nrow   = 0;

% read the boundaries
geo  = [];
for i=1:ngeo
  geo(i).name    = char(fread(fid, [1 80], 'uchar'));
  geo(i).npnt    = fread(fid, 1, 'int');
  geo(i).pnt     = fread(fid, [3 geo(i).npnt], 'float')';
  geo(i).ndhk    = fread(fid, 1, 'int');
  geo(i).dhk     = fread(fid, [3 geo(i).ndhk], 'int')' + 1;  % Matlab indexing starts at 1
  geo(i).sigmam  = fread(fid, 1, 'float');
  geo(i).sigmap  = fread(fid, 1, 'float');
  geo(i).geocon  = fread(fid, ngeo, 'int');
  geo(i).deflat  = fread(fid, ngeo, 'float');
  totpnt = totpnt + geo(i).npnt;
  totdhk = totdhk + geo(i).ndhk;
end

% read the electrodes
if mode~=1
  elec.name    = char(fread(fid, [1 80], 'uchar'));
  elec.npnt    = fread(fid, 1, 'int');
  for i=1:(elec.npnt+1)
    elec.el(i).dhk  = fread(fid, 1, 'int') + 1; % Matlab indexing starts at 1
    elec.el(i).la   = fread(fid, 1, 'float');
    elec.el(i).mu   = fread(fid, 1, 'float');
    elec.el(i).name = char(fread(fid, [1 10], 'char'));
    % the ELECTRODE c-structure is padded to word boundaries, i.e. to 4 bytes
    dum = fread(fid, 2, 'char');
  end
  elec.vertex  = fread(fid, 1, 'int');
  elec.surface = fread(fid, 1, 'int');
  nrow = nrow + elec.npnt;
else
  elec = [];
end

% read the gradiometers
if mode~=0
  error('gradiometers not yet implemented');
else
  grad = [];
end

% read the inverted A-matrix
bi = fread(fid, [totpnt nrow], 'float')';

% read the isolated source compartment information, if present
iso_sur    = fread(fid, 1, 'int') + 1;  % Matlab indexing starts at 1
inner_only = fread(fid, 1, 'int');
if iso_sur~=0
  iso_totpnt = geo(iso_sur).npnt;
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

