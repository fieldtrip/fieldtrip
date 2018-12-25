function GridDim = determine_griddim(elec)

% DETERMINE_GRIDDIM uses the labels and positions of electrodes in elec to
% determine the dimensions of each set of electrodes (i.e., electrodes with
% the same string, but different numbers)
%
% use as: 
%   GridDim = determine_griddim(elec)
%   where elec is a structure that contains an elecpos field and a label field
%   and GridDim(1) = number of rows and GridDim(2) = number of columns
%   
% See also FT_ELECTRODEREALIGN

% Copyright (C) 2017-2018, Arjen Stolk, Sandon Griffin
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

% determine electrode range
digits = regexp(elec.label, '\d+', 'match');
elec.maxdigit = 1;
elec.mindigit = Inf;
for l=1:numel(digits)
  ElecStrs{l,1} = regexprep(elec.label{l}, '\d+(?:_(?=\d))?', ''); % without electrode numbers
  elec.ElecLab{l,1} = digits{l}{1}; % use first found digit
  if str2num(digits{l}{1}) > elec.maxdigit
    elec.maxdigit = str2num(digits{l}{1});
  end
  if str2num(digits{l}{1}) < elec.mindigit
    elec.mindigit = str2num(digits{l}{1});
  end
end
elec.ElecStr = cell2mat(unique(ElecStrs));

% determine if any electrodes appear to be misordered or cut out of this grid
pos_ordered = [];
labels_ordered = {};
elec.cutout = []; % index of electrodes that appear to be cut out
dowarn = false;
for e = elec.mindigit:elec.maxdigit
  if ~isempty(match_str(elec.ElecLab, num2str(e))) % in case labels are 1, 2, 3 etc.
    labels_ordered{end+1,1} = [elec.ElecStr num2str(e)];
    pos_ordered(end+1, :) = elec.elecpos(match_str(elec.ElecLab, num2str(e)),:);
  elseif ~isempty(match_str(elec.ElecLab, num2str(e, ['%0' num2str(numel(elec.ElecLab{1})) 'd']))) % in case labels are 001, 002, 003 etc.
    labels_ordered{end+1,1} = [elec.ElecStr num2str(e, ['%0' num2str(numel(elec.ElecLab{1})) 'd'])];
    pos_ordered(end+1, :) = elec.elecpos(match_str(elec.ElecLab, num2str(e, ['%0' num2str(numel(elec.ElecLab{1})) 'd'])),:);
  else
    elec.cutout(end+1) = e;
    labels_ordered{end+1,1} = [num2str(e)];
    pos_ordered(end+1, :) = NaN(1,3); % replace missing numbers with electrodes at position [NaN NaN NaN]
    dowarn = true;
  end
end
if dowarn
  ft_warning('%s appears to be missing electrodes %s or have electrodes that are labeled using an unconventional numbering system', elec.ElecStr, num2str(elec.cutout));
end
elec.label = labels_ordered;
elec.elecpos = pos_ordered;

% find the electrode spacing
d_bwelec = [];
for e = 1:size(elec.elecpos,1)-1
  d_bwelec(e) = sqrt(sum((elec.elecpos(e+1,:)-elec.elecpos(e,:)).^2));
end
d_bwelec(find(d_bwelec == 0)) = NaN;
elec_space = nanmedian(d_bwelec); 

% For each electrode, find the electrodes that are within 1.3*elec_space,
% which should be the adjacent neighbors. Without considering the
% electrodes in the same row (e.g., those that are numbered either +1 or -1 
% to the given electrode), find the numbers of the electrodes in the
% same column. Then, find the difference between the numbers of electrodes
% in the same column. The most common difference should be the row length.
vertical_skip = []; % n x 2 matrix where column 1 refers to a given electrode and column 1 refers to the difference in this electrode's number and the electrodes in the same column
for n = 1:size(elec.elecpos,1)
  adj_elec_nums = []; % numbers of the electrodes in the same row or column
  for e = 1:size(elec.elecpos,1)
    if n~=e && n~=e-1 && n~=e+1 && n<e
      if sqrt(sum((elec.elecpos(e,:)-elec.elecpos(n,:)).^2)) < 1.3*elec_space
        vertical_skip(end+1, 2) = abs(n-e);
        vertical_skip(end, 1) = n;
      end
    end
  end
end

isLshaped = 0;
if ~isempty(vertical_skip) && size(vertical_skip,1)>2 % if this is a two dimensional set of electrodes
  % try to determine if this is a square set of electrodes or L-shaped
  uniq_skips = unique(vertical_skip(:,2));
  count = zeros(numel(uniq_skips),1);
  for k = 1:numel(uniq_skips)
    count(k) = sum(vertical_skip(:,2)==uniq_skips(k));
  end
  mode1 = mode(vertical_skip(:,2));
  n_mode = count(find(uniq_skips == mode1)); % the count of the mode
  modes = uniq_skips(find(count > n_mode*0.75));
  GridDim(2) = mode(vertical_skip(:,2));
  GridDim(1) = ceil((elec.maxdigit-elec.mindigit+1)/GridDim(2));
  
  if numel(modes) == 2 && diff(modes) ~= 1 % this set of electrodes is derived from an L-shaped grid
    LGridDim = NaN(2,2);
    
    mode2 = modes;
    mode2(find(modes == mode(vertical_skip(:,2)))) = [];
    e_kink = vertical_skip(find(vertical_skip(:,2) == mode2, 1),1) + mode(vertical_skip(:,2)) -1; % number of electrode on the interior where the L-shape starts
    LGridDim(1, 2) = mode(vertical_skip(:,2));
    LGridDim(1, 1) = ceil(e_kink/LGridDim(1, 2));
    
    LGridDim(2, 2) = mode2;
    LGridDim(2, 1) = ceil((elec.maxdigit-elec.mindigit+1-e_kink)/LGridDim(2,2));
    ft_warning('%s appears to be a grid that is derived from an L-shaped structure and is best described as a %d x %d grid adjacent to a %d x %d grid', elec.ElecStr, LGridDim(1,1), LGridDim(1,2), LGridDim(2,1), LGridDim(2,2));
    isLshaped = 1;
  end
  
else
  GridDim(2) = elec.maxdigit-elec.mindigit+1;
  GridDim(1) = 1;
end

% Give feedback about the grid dimensions
fprintf('assuming %s has %d x %d grid dimensions\n', elec.ElecStr, GridDim(1), GridDim(2));
