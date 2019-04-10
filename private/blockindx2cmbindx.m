function [cmbindx, n, blocklabel] = blockindx2cmbindx(labelcmb, blockindx, block)

% This is a helper function that is needed for the bookkeeping of the data,
% when requesting (conditional)-blockwise granger causality estimates. Its
% single use is in ft_connectivityanalysis, but in order to keep that code
% clean, it was decided to put this function as a private function.
%
% Use as 
%   [cmbindx, n, blocklabel] = blockindx2cmbindx(labelcmb, blockindx,
%   block)
%
% The purpose is to generate a cell-array (Nx2, same size as input array
% block) of numeric indices, which index into the rows of the Mx2 labelcmb
% array, and which can subsequently be used by lower-level functionality
% (i.e. blockwise_conditionalgranger) to compute the connectivity metric of
% interest. Blockindx is a 1x2 cell-array, which maps the individual
% channels in blockindx{1} to an indexed block in blockindx{2}. Block
% specifies in each row of cells two ordered lists of blocks that are
% needed to compute a conditioned Granger spectrum.

% Copyright (C) 2010-2017 Jan-Mathijs Schoffelen
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

% The variable labelcmb contains the specification of how each
% row of the data has been computed, from chanX to chanY given a multivariate
% decomposition of [chanA,...,chanZ], first generate a cell-array y and z
% that splits the channel pair from the multivariate block
y = cell(size(labelcmb));
z = cell(size(labelcmb,1),1);
cmbindx = cell(size(block));
n       = cell(size(block));
for k = 1:size(y,1)
  x1     = strfind(labelcmb{k,1},'['); 
  x2     = strfind(labelcmb{k,1},']');
  y{k,1} = labelcmb{k,1}(1:x1-1);    % string denoting single channel
  z{k,1} = labelcmb{k,1}(x1+1:x2-1); % concatenated string of all channels in the current blockwise decomp
  x1     = strfind(labelcmb{k,2},'['); 
  y{k,2} = labelcmb{k,2}(1:x1-1);
end
z(:,2) = z(:,1);

uz = unique(z(:,1)); % this is the number of unique blocks
% the unique blocks can be of different size, but should be present in full
% in the labelcmb, i.e. correspond to n^2 rows in the labelcmb

z_idx = zeros(size(z,1),1);
z_block = cell(numel(uz),1);
z_label = cell(numel(uz),1);
for k = 1:numel(uz)
  tmp = strcmp(z(:,1),uz{k});
  tok = tokenize(uz{k}, ',');
  
  z_label{k} = tok(:); % the label of the channels that participate in a block
  
  [dum,i2x] = match_str(y(tmp,1), z_label{k});
  [dum,i2y] = match_str(y(tmp,2), z_label{k});
  
  if sqrt(sum(tmp))~=numel(tok)
    error('incomplete data');
  end
  z_idx(tmp) = k;
  
  tmp = find(tmp);
  list = zeros(numel(tok));
  for kk = 1:numel(i2x)
  list(i2x(kk),i2y(kk)) = tmp(kk);
  end
  z_block{k} = list; % the ordered indices of the combinations in a block, as per the order in z_label
end
z_label_sorted = z_label;
for k = 1:numel(z_label)
  z_label_sorted{k} = sort(z_label{k});
end

% now expand the block vectors to channel labels (ordered)
b_label = cell(size(block));
b_n     = cell(size(block));
for k = 1:numel(block)
  tmplist = cell(0,1);
  tmpn    = zeros(0,1);
  for m = 1:numel(block{k})
    tmplist = cat(1,tmplist,blockindx{1}(blockindx{2}==block{k}(m)));
    tmpn    = cat(1,tmpn,sum(blockindx{2}==block{k}(m)));
  end
  b_label{k} = tmplist;
  b_n{k}     = tmpn;
end

% now identify which of the z_blocks corresponds to the requested blocks,
% and reorder if needed
for k = 1:numel(block)
  indx = [];
  for m = 1:numel(z_label_sorted)
    if ~isequal(sort(b_label{k}), z_label_sorted{m})
      continue;
    else
      indx = m;
      break;
    end
  end
  [dum, i2]    = match_str(b_label{k}, z_label{indx});
  cmbindx{k} = reshape(z_block{indx}(i2,i2),[],1);
  n{k}       = b_n{k}; %(i2);
end

blocklabel = []; % don't know why this is needed, keep for now.
