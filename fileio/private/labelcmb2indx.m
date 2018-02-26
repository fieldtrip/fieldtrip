function [indx, label, blockindx, blocklabel] = labelcmb2indx(labelcmb, label)

% LABELCMB2INDX computes an array with indices, corresponding to the order
% in a list of labels, for an Nx2 list of label combinations
%
% Use as
%   [indx] = labelcmb2indx(labelcmb, label)
% or
%   [indx] = labelcmb2indx(labelcmb)
%
% Labelcmb is an Nx2 cell-array with label combinations, label is an Mx1 
% cell-array with labels. If only one input is provided, the indices are
% with respect to the rows in the labelcmb matrix, where the corresponding
% auto combinations are located. As a consequence, the labelcmb matrix 
% needs to contain rows containing auto-combinations
% 
% Example: 
%  labelcmb = {'a' 'b';'a' 'c';'b' 'c';'a' 'a';'b' 'b';'c' 'c'};
%  label    = {'a';'b';'c'};
%
% indx = labelcmb2indx(labelcmb, label)
%  returns:  [1 2;1 3;2 3;1 1;2 2;3 3]
% 
% indx = labelcmb2indx(labelcmb)
%  returns:  [4 5;4 6;5 6;4 4;5 5;6;6]
%
% This is a helper function to FT_CONNECTIVITYANALYSIS

% Copyright (C) 2009-2010 Donders Institute, Jan-Mathijs Schoffelen
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

% check whether the labelcmb contains any square brackets, indicative of
% blockwise decompositions
if all(~cellfun('isempty',strfind(labelcmb(:),'[')))
  tmp = strfind(labelcmb, '[');
  tmplabelcmb = labelcmb;
  tmpblock    = labelcmb;
  for k = 1:numel(tmplabelcmb)
    tmplabelcmb{k} = labelcmb{k}(         1:(tmp{k}-1));
    tmpblock{k}    = labelcmb{k}((tmp{k}+1):(end-1));
  end
  blocklabel = cell(0,1);
  tmp = tmpblock;
  while ~isempty(tmp)
    blocklabel{end+1, 1} = tmp{1};
    tmp = tmp(~strcmp(tmp(:,1),blocklabel{end}),:);
  end
  [indx,      label]      = labelcmb2indx(tmplabelcmb, unique(tmplabelcmb(:)));
  [blockindx, blocklabel] = labelcmb2indx(tmpblock,    blocklabel);
  blockindx = blockindx(:,1);
  return;
end


if nargin==1,
  label = unique(labelcmb(:));
end

%identify the auto-combinations
ncmb = size(labelcmb,1);
indx = zeros(ncmb,2);

nchan    = numel(label);
for k = 1:nchan
  sel = strcmp(label{k}, labelcmb);
  if nargin==1,
    autoindx = find(sel(:,1) & sel(:,2), 1, 'first');
    if isempty(autoindx), ft_error('the required autocombination is not found in the input'); end
  else
    autoindx = k;
  end
  indx(sel(:,1),1) = autoindx;
  indx(sel(:,2),2) = autoindx;
end
