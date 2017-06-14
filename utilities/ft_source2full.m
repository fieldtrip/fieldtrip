function [source] = ft_source2full(source)

% FT_SOURCE2FULL recreates the grid locations outside the brain in the source 
% reconstruction, so that the source volume again describes the full grid.
% This undoes the memory savings that can be achieved using FT_SOURCE2SPARSE
% and makes it possible again to plot  the source volume and save it to an
% external file.
%
% Use as
%   [source] = ft_source2full(source)
%
% See also FT_SOURCE2SPARSE

% Copyright (C) 2004, Robert Oostenveld
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

ft_defaults

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
  sparsepos = source.pos;
  
  % recreate the positions of the dipole grid
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  pos = [X(:) Y(:) Z(:)];
else
  %FIXME this assumes that the voxel data are ordered as if in a regularly spaced 3D grid,
  %but with only the inside voxels present
  warning('assuming the voxel data to be ordered as if in a regularly spaced 3D grid');
  xgrid = 1:source.dim(1);
  ygrid = 1:source.dim(2);
  zgrid = 1:source.dim(3);

  %establish a homogeneous transformation matrix from voxels to headspace based on the sparse positions
  sparsepos = source.pos;
  ok  = 0;
  cnt = 0;
  while ok==0,
    cnt  = cnt+1;
    dpos = sparsepos - sparsepos(cnt*ones(size(sparsepos,1),1),:);
    [srt, indx] = sort(sum(dpos.^2,2));
    srt    = dpos(indx,:);
    tmpsrt = abs(srt(2:7,:));
    csrt   = tmpsrt*tmpsrt';
    sel    = find(sum(csrt==0)>=2);
    if numel(sel)>=3, 
      ok = 1;
    end
  end  
  tmppos  = sparsepos(indx([1 sel(:)'+1]),:);
  tmpdpos = dpos(indx([1 sel(:)'+1]),:);
 
  % FIXME the following is a bit experimental and not fully tested yet it works in general case rotation
  M         = pinv(tmpdpos(2:4,:));
  
  % get rotation such that maxima are on diagonal and positive
  m(1) = find(M(1,:)==max(abs(M(1,:))));
  m(2) = find(M(2,:)==max(abs(M(2,:))));
  m(3) = find(M(3,:)==max(abs(M(3,:))));
  [srt, indx] = sort(m);
  M    = M(indx,:);
  M    = M*diag(sign(diag(M)));
  sparsepos = sparsepos*M;
  
  % translation
  T         = -min(sparsepos,[],1)+1;
  sparsepos = sparsepos + T(ones(size(sparsepos,1),1), :);  

  % recreate the positions of the dipole grid
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  pos = [X(:) Y(:) Z(:)];
  pos = ft_warp_apply(inv([M T(:);0 0 0 1]), pos);
end

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
  inside = zeros(Nsparse,1);
  for i=1:Nsparse
    fx = find(xgrid==sparsepos(i,1));
    fy = find(ygrid==sparsepos(i,2));
    fz = find(zgrid==sparsepos(i,3));
      inside(i) = (fx-1)*sx + (fy-1)*sy + (fz-1)*sz + 1;
  end
  outside = setdiff([1:Nfull]', inside);
end

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

% determine whether the source is old or new style
fnames = fieldnames(source);
if any(~cellfun('isempty', strfind(fnames, 'dimord'))),
  stype = 'new';
else
  stype = 'old';
end

if strcmp(stype, 'old')
  % original code
  % first do the non-trial fields
  source.dim = [1 length(inside) 1]; %to fool parameterselection
  [param]    = parameterselection('all', source);
  trlparam   = strmatch('trial', param);
  sel        = setdiff(1:length(param), trlparam);
  ind=find(ismember(param,'inside'));% find the index of 'inside' field
  % because its position varies with isfield('plvspctrm') vs. 'cohspctrm'
  param      = param(sel(ind));
  
  for j = 1:length(param)
    dat = getsubfield(source, param{j});
    if islogical(dat)
      tmp         = false(1,Nfull); 
      tmp(inside) = dat;
    elseif iscell(dat)
      tmp          = cell(1,Nfull);
      tmp(inside)  = dat;
      %tmp(outside) = nan;
    else
      tmp         = nan(1,Nfull);
      tmp(inside) = dat;   
    end
    source = setsubfield(source, param{j}, tmp);
  end
  
  % then do the trial fields
  if     isfield(source, 'trial' )
    for j = 1:length(source.trial)
      tmpsource     = source.trial(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
          tmp(inside) = dat;   
        end
        tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
      end
      tmpsource       = rmfield(tmpsource, 'dim');
      source.trial(j) = tmpsource;
    end   
  elseif isfield(source, 'trialA')
    for j = 1:length(source.trialA)
      tmpsource     = source.trialA(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
          tmp(inside) = dat;   
        end
        tmpsource = setsubfield(tmpsource, tmpparam{k}, tmp);
      end
      tmpsource        = rmfield(tmpsource, 'dim');
      source.trialA(j) = tmpsource;   
    end
  elseif isfield(source, 'trialB')
    for j = 1:length(source.trialB)
      tmpsource     = source.trialB(j);
      tmpsource.dim = source.dim; % to fool parameterselection
      tmpparam      = parameterselection('all', tmpsource);
      for k = 1:length(tmpparam)
        dat = getsubfield(tmpsource, tmpparam{k});
        if islogical(dat)
          tmp         = false(1,Nfull); 
          tmp(inside) = dat;
        elseif iscell(dat)
          tmp          = cell(1,Nfull);
          tmp(inside)  = dat;
          %tmp(outside) = nan;
        else
          tmp         = nan(1,Nfull);
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
      tmp = nan(Nfull,Nfull);
      tmp(inside,inside) = d;
      source = setfield(source, fn{i}, tmp);
    end
  end
  
  % update the inside and outside definitions
  source.inside  = inside;
  source.outside = outside;
  source.pos     = pos;
  source.dim     = siz;
elseif strcmp(stype, 'new')
  % new style conversion
  fn = fieldnames(source);
  for i=1:numel(fn)
    if any(size(source.(fn{i}))==Nsparse)
      if iscell(source.(fn{i}))
        indx = find(size(source.(fn{i}))==Nsparse);
        if all(indx==1)
          tmp            = cell(Nfull,1);
          tmp(inside,1)  = source.(fn{i});
          source.(fn{i}) = tmp;
        elseif all(indx==2)
          tmp            = cell(1,Nfull);
          tmp(1,inside)  = source.(fn{i});
          source.(fn{i}) = tmp;
        else
          warning('sparse to full conversion failed for field %s\n', fn{i});
        end
      else
        indx = find(size(source.(fn{i}))==Nsparse);
        if all(indx==1)
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([Nfull tmpsiz(2:end)]);
          tmp(inside,:,:,:,:) = source.(fn{i});
        elseif all(indx==2)
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([tmpsiz(1) Nfull tmpsiz(3:end)]);
          tmp(:,inside,:,:,:) = source.(fn{i});
        elseif all(indx==[1 2])
          % bivariate matrix
          tmpsiz = [size(source.(fn{i})) 1];
          tmp    = nan([Nfull Nfull tmpsiz(3:end)]);
          tmp(inside,inside,:,:,:) = source.(fn{i});
        else
          warning('sparse to full conversion failed for field %s\n', fn{i});
        end
      end
      % nothing to do
    end
  end
  
  % update the inside and outside definitions and pos
  source.inside  = inside;
  source.outside = outside;
  source.pos     = pos;

end
cfg = [];
% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with MATLAB versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;
