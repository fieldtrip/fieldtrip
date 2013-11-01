function [dipout] = ft_eloreta(dip, grad, vol, dat, Cf, varargin)
%
% Use as
%   [dipout] = ft_eloreta(dipin, grad, vol, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   grad        is the gradiometer definition
%   vol         is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid, Npositions x 3
%   dipin.mom   dipole orientation (optional), 3 x Npositions
%   dipin.filter precomputed filter cell-array of 1 x Npositions, each cell
%               containing nchannels x 3 matrix
%
% Additional options should be specified in key-value pairs and can be
%  'keepfilter'       = remember the beamformer filter,    can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,  can be 'yes' or 'no'
%  'keepmom'          = remember the dipole moment,        can be 'yes' or 'no'
%  'lambda'           = regularisation parameter

% These options influence the forward computation of the leadfield
%  'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'        = normalize the leadfield
%  'normalizeparam'   = parameter for depth normalization (default = 0.5)

% If the dipole definition only specifies the dipole location, a rotating
% dipole (regional source) is assumed on each location. If a dipole moment
% is specified, its orientation will be used and only the strength will
% be fitted to the data.

% Copyright (C) 2013, Marlene Boenstrup, Jan-Mathijs Schoffelen and Guido
% Nolte
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id:$

if mod(nargin-5,2)
    % the first 5 arguments are fixed, the other arguments should come in pairs
    error('invalid number of optional arguments');
end

% these optional settings do not have defaults
keepfilter  = ft_getopt(varargin, 'keepfilter',      'no');
keepmom     = ft_getopt(varargin, 'keepmom',         'no');
keepleadfield = ft_getopt(varargin, 'keepleadfield', 'no');
lambda      = ft_getopt(varargin, 'lambda', 0.05);
reducerank  = ft_getopt(varargin, 'reducerank');
normalize   = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam');

% convert the yes/no arguments to the corresponding logical values
keepfilter     = istrue(keepfilter);
keepmom        = istrue(keepmom);
keepleadfield  = istrue(keepleadfield);

% find the dipole positions that are inside/outside the brain
if ~isfield(dip, 'inside') && ~isfield(dip, 'outside');
    insideLogical  = ft_inside_vol(dip.pos, vol);
    dip.inside     = find(insideLogical);
    dip.outside    = find(~dip.inside);
elseif isfield(dip, 'inside') && ~isfield(dip, 'outside');
    dip.outside    = setdiff(1:size(dip.pos,1), dip.inside);
elseif ~isfield(dip, 'inside') && isfield(dip, 'outside');
    dip.inside     = setdiff(1:size(dip.pos,1), dip.outside);
end

% select only the dipole positions inside the brain for scanning
dip.origpos     = dip.pos;
dip.originside  = dip.inside;
dip.origoutside = dip.outside;

if isfield(dip, 'mom'),
    dip.mom = dip.mom(:,dip.inside);
end
if isfield(dip, 'leadfield'), fprintf('using precomputed leadfields\n');
    dip.leadfield = dip.leadfield(dip.inside);
end
if isfield(dip, 'filter')
    fprintf('using precomputed filters\n');
    dip.filter = dip.filter(dip.inside);
end

dip.pos     = dip.pos(dip.inside, :);
dip.inside  = 1:size(dip.pos,1);
dip.outside = [];

% deal with the regularition in the low level function

% use existing leadfields, or compute them
if ~isfield(dip, 'leadfield')
    % compute the leadfield
    fprintf('computing leadfields\n');
    for i=1:size(dip.pos,1)
        dip.leadfield{i} = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
    end
end

% deal with dip.mom
if isfield(dip, 'mom') && size(dip.leadfield{1},2)==size(dip.mom,1)
    fprintf('projecting the forward solutions using specified dipole moment\n');
    for i=1:size(dip.pos,1)
        dip.leadfield{i} = dip.leadfield{i}*dip.mom(:,i);
    end
end

% deal with reduced rank
% check the rank of the leadfields, and project onto the lower dimensional
% subspace if the number of columns per leadfield > rank: this to avoid
% numerical issues in the filter computation
rank_lf = zeros(size(dip.pos,1));
for i=1:size(dip.pos,1)
    rank_lf(i) = rank(dip.leadfield{i});
end
if ~all(rank_lf==rank_lf(1)),
    error('the forward solutions have a different rank for each location. this is not supported');
end
if rank_lf(1)<size(dip.leadfield{1})
    fprintf('the forward solutions have a rank of %d, but %d orientations\n',rank_lf(1),size(dip.leadfield{1},2));
    fprintf('projecting the forward solutions on the lower dimensional subspace\n');
    for i=1:size(dip.pos,1)
        [u,s,v{i}] = svd(dip.leadfield{i}, 'econ');
        dip.leadfield{i} = dip.leadfield{i}*v{i}(:,1:rank_lf);
    end
end

% convert the leadfield into NchanxNdipxNori
[Nchan, Nori] = size(dip.leadfield{1});
Ndip          = numel(dip.leadfield);
leadfield     = permute(reshape(cat(2,dip.leadfield{:}),Nchan,Nori,Ndip),[1 3 2]);

% use existing filters, or compute them
if ~isfield(dip, 'filter')
    filt = mkfilt_eloreta_v2(leadfield, lambda);
    for i=1:size(dip.pos,1)
        dip.filter{i} = squeeze(filt(:,i,:))';
    end
end

% get the power
dip.pow = zeros(size(dip.pos,1),1);
dip.ori = cell(size(dip.pos,1),1);
for i=1:size(dip.pos,1)
    csd     = dip.filter{i}*Cf*dip.filter{i}';
    [u,s,v] = svd(real(csd));
    dip.pow(i) = s(1);
    dip.ori{i} = u(:,1);
end

% get the dipole moment
if keepmom && ~isempty(dat)
    % remove the dipole moment from the input
    if isfield(dip, 'mom')
        dip = rmfield(dip, 'mom');
    end
    for i=1:size(dip.pos,1)
        dip.mom{i} = dip.filter{i}*dat;
    end
end

if keepfilter,    dipout.filter    = dip.filter;    end
if keepleadfield, dipout.leadfield = dip.leadfield; end
if keepmom && isfield(dip, 'mom'), dipout.mom = dip.mom; end

dipout.pow     = dip.pow;
dipout.inside  = dip.originside;
dipout.outside = dip.origoutside;
dipout.pos     = dip.origpos;
dipout.ori     = dip.ori;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
    dipout.leadfield(dipout.inside)  = dipout.leadfield;
    dipout.leadfield(dipout.outside) = {[]};
end
if isfield(dipout, 'filter')
    dipout.filter(dipout.inside)  = dipout.filter;
    dipout.filter(dipout.outside) = {[]};
end
if isfield(dipout, 'mom')
    dipout.mom(dipout.inside)  = dipout.mom;
    dipout.mom(dipout.outside) = {[]};
end
if isfield(dipout, 'pow') %here pow is cell
    dipout.pow(dipout.inside)  = dipout.pow;
    dipout.pow(dipout.outside) = nan;
end
if isfield(dipout, 'ori') %here pow is cell
    dipout.ori(dipout.inside)  = dipout.ori;
    dipout.ori(dipout.outside) = nan;
end
