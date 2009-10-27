function [source] = source2full(source);

% SOURCE2FULL recreates the grid locations outside the brain in the source 
% reconstruction, so that the source volume again describes the full grid.
% This undoes the memory savings that can be achieved using SOURC2SPARSE
% and makes it possible again to plot  the source volume and save it to an
% external file.
%
% Use as
%   [source] = source2full(source)
%
% See also SOURCE2SPARSE

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

if ~isfield(source, 'inside')  || ...
   ~isfield(source, 'outside') || ...
   ~isfield(source, 'dim')
  error('one of the required fields is missing in the source structure');
end

if ~isfield(source, 'pos') && (~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ...
                               ~isfield(source, 'zgrid'))
  error('the input data needs at least a ''pos'' field, or ''x/y/zgrid''');
end

if isfield(source, 'xgrid'),
  xgrid = source.xgrid;
  ygrid = source.ygrid;
  zgrid = source.zgrid;
else
  %FIXME this assumes that the voxel data are ordered as if in a regularly spaced 3D grid,
  %but with only the inside voxels present
  xgrid = 1:source.dim(1);
  ygrid = 1:source.dim(2);
  zgrid = 1:source.dim(3);
end
% recreate the positions of the dipole grid
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
pos = [X(:) Y(:) Z(:)];

Nsparse = length(source.inside);
siz     = source.dim;
Nfull   = prod(siz);

% determine the size that each slice takes in memory
sx = 1;
sy = siz(1);
sz = siz(1) * siz(2);

if isfield(source, 'inside') && isfield(source, 'outside') && size(source.pos,1)==Nfull
  % it contains all source positions
  inside = source.inside;
  outside = source.outside;
else
  % it only contains the inside source positions, which are all inside the brain
  % reconstruct the original inside and outside grid locations
  for i=1:Nsparse
    fx = find(source.xgrid==source.pos(i,1));
    fy = find(source.ygrid==source.pos(i,2));
  fz = find(source.zgrid==source.pos(i,3));
      inside(i) = (fx-1)*sx + (fy-1)*sy + (fz-1)*sz + 1;
  end
  outside = setdiff(1:Nfull, inside);
end

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

% first do the non-trial fields
source.dim = [1 length(inside) 1]; %to fool parameterselection
[param]    = parameterselection('all', source);
trlparam   = strmatch('trial', param);
sel        = setdiff(1:length(param), trlparam);
param      = param(sel);

for j = 1:length(param)
  dat = getsubfield(source, param{j});
  if islogical(dat),
    tmp         = false(1,Nfull); 
    tmp(inside) = dat;
  elseif iscell(dat),
    tmp          = cell(1,Nfull);
    tmp(inside)  = dat;
    %tmp(outside) = nan;
  else
    tmp         = zeros(1,Nfull) + nan;
    tmp(inside) = dat;   
  end
  source = setsubfield(source, param{j}, tmp);
end

% then do the trial fields
if     isfield(source, 'trial' ),
  for j = 1:length(source.trial)
    tmpsource     = source.trial(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat = getsubfield(tmpsource, tmpparam{k});
      if strcmp(class(dat), 'logical'),
        tmp         = logical(zeros(1,Nfull)); 
        tmp(inside) = dat;
      elseif strcmp(class(dat), 'cell'),
        tmp          = cell(1,Nfull);
        tmp(inside)  = dat;
        %tmp(outside) = nan;
      else
        tmp         = zeros(1,Nfull) + nan;
        tmp(inside) = dat;   
      end
      tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
    end
    tmpsource       = rmfield(tmpsource, 'dim');
    source.trial(j) = tmpsource;
  end   
elseif isfield(source, 'trialA'),
  for j = 1:length(source.trialA)
    tmpsource     = source.trialA(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat = getsubfield(tmpsource, tmpparam{k});
      if strcmp(class(dat), 'logical'),
        tmp         = logical(zeros(1,Nfull)); 
        tmp(inside) = dat;
      elseif strcmp(class(dat), 'cell'),
        tmp          = cell(1,Nfull);
        tmp(inside)  = dat;
        %tmp(outside) = nan;
      else
        tmp         = zeros(1,Nfull) + nan;
        tmp(inside) = dat;   
      end
      tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
    end
    tmpsource        = rmfield(tmpsource, 'dim');
    source.trialA(j) = tmpsource;   
  end
elseif isfield(source, 'trialB'),
  for j = 1:length(source.trialB)
    tmpsource     = source.trialB(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat = getsubfield(tmpsource, tmpparam{k});
      if strcmp(class(dat), 'logical'),
        tmp         = logical(zeros(1,Nfull)); 
        tmp(inside) = dat;
      elseif strcmp(class(dat), 'cell'),
        tmp          = cell(1,Nfull);
        tmp(inside)  = dat;
        %tmp(outside) = nan;
      else
        tmp         = zeros(1,Nfull) + nan;
        tmp(inside) = dat;   
      end
      tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
    end
    tmpsource        = rmfield(tmpsource, 'dim');
    source.trialB(j) = tmpsource;   
  end
end

% and finally do the coherence-like matrices (size Nvox X Nvox)
fn = fieldnames(source);
for i=1:length(fn)
  d = getfield(source, fn{i});
  m = size(d, 1);
  n = size(d, 2);
  if m==Nsparse && n==Nsparse
    tmp = nan*zeros(Nfull,Nfull);
    tmp(inside,inside) = d;
    source = setfield(source, fn{i}, tmp);
  end
end

% update the inside and outside definitions
source.inside  = inside;
source.outside = outside;
source.pos     = pos;
source.dim     = siz;

cfg = [];
% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: source2full.m,v 1.19 2009/10/12 14:20:19 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;
