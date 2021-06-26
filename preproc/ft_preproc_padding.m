function [dat] = ft_preproc_padding(dat, padtype, prepadlength, postpadlength)

% FT_PREPROC_PADDING performs padding on the data, i.e. adds or removes samples
% to or from the data matrix.
%
% Use as
%   [dat] = ft_preproc_padding(dat, padtype, padlength)
% or as
%   [dat] = ft_preproc_padding(dat, padtype, prepadlength, postpadlength)
% where
%   dat           data matrix (Nchan x Ntime)
%   padtype       'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove'
%   padlength     scalar, number of samples that will be padded
%   prepadlength  scalar, number of samples that will be padded before the data
%   postpadlength scalar, number of samples that will be padded after the data
%
% If padlength is used instead of prepadlength and postpadlength, padding
% will be symmetrical (i.e. padlength = prepadlength = postpadlength)
%
% If the data contains NaNs, these are ignored for the computation, but
% retained in the output. Depending on the type of padding, NaNs may spread
% to the pads.
%
% See also FT_PREPROCESSING

% Copyright (C) 2012, Jorn M. Horschig, Robert Oostenveld, Jan-Mathijs Schoffelen
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

if nargin < 4
  postpadlength = prepadlength;
end

if prepadlength == 0 && postpadlength == 0
  return;
end

[nchans, nsamples] = size(dat);

switch(padtype)
  case 'remove'
    dat = dat(:, prepadlength+1:end-postpadlength);
    
  case 'mirror'
    % create an indexvector to index the data with
    index = (1:(prepadlength+nsamples+postpadlength))-prepadlength;
    while any(index<1|index>nsamples)
      index(index<1)        = -index(index<1) + 1;
      index(index>nsamples) = 2.*nsamples - index(index>nsamples) + 1;
    end
    dat = dat(:, index);
    
  case 'edge'
    dat       = [dat(:,1)*ones(1,prepadlength) dat dat(:,end)*ones(1,postpadlength)];
    
  case 'mean'
    mu        = mean(dat, 2);
    dat       = [mu*ones(1,prepadlength) dat mu*ones(1,postpadlength)];
    
  case 'localmean'
    prepad    = min(prepadlength, floor(size(dat, 2)/2));
    edgeleft  = mean(dat(:, 1:prepad), 2);
    postpad   = min(postpadlength, floor(size(dat, 2)/2));
    edgeright = mean(dat(:, 1+end-postpad:end), 2);
    dat       = [edgeleft*ones(1,prepadlength) dat edgeright*ones(1,postpadlength)];
    
  case 'zero'
    dat       = [zeros(nchans,prepadlength) dat zeros(nchans,postpadlength)];
    
  case 'nan'
    dat       = [nan(nchans,prepadlength) dat nan(nchans,postpadlength)];
    
  otherwise
    ft_error('unknown padding option');
end
