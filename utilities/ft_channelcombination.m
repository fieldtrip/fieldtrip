function [collect] = ft_channelcombination(channelcmb, datachannel, includeauto, dirflag)

% FT_CHANNELCOMBINATION creates a cell-array with combinations of EEG/MEG
% channels for subsequent cross-spectral-density and coherence analysis
%
% You should specify channel combinations as a two-column cell array,
%   cfg.channelcmb = {  'EMG' 'MLF31'
%                       'EMG' 'MLF32'
%                       'EMG' 'MLF33' };
% to compare EMG with these three sensors, or
%   cfg.channelcmb = { 'MEG' 'MEG' };
% to make all MEG combinations, or
%   cfg.channelcmb = { 'EMG' 'MEG' };
% to make all combinations between the EMG and all MEG channels.
%
% For each column, you can specify a mixture of real channel labels
% and of special strings that will be replaced by the corresponding
% channel labels. Channels that are not present in the raw datafile
% are automatically removed from the channel list.
%
% When directional connectivity measures will subsequently be computed, the
% interpretation of each channel-combination is that the direction of the
% interaction is from the first column to the second column.
%
% Note that the default behaviour is to exclude symmetric pairs and
% auto-combinations.
%
% See also FT_CHANNELSELECTION

% Undocumented options:
%   includeauto = 0 (or 1), include auto-combinations
%   dirflag     = 0 (or 1, or 2) specifies the treatment of the order in
%                   the columns of the output. If dirflag = 0, the order is
%                   preserved, if dirflag = 1, the order is reversed,
%                   second input column is put first. If dirflag = 2, both
%                   directions are added to the list of combinations.

% Copyright (C) 2003-2015, Jan-Mathijs Schoffelen & Robert Oostenveld
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

if nargin<3 || isempty(includeauto), includeauto = 0; end
if nargin<4 || isempty(dirflag),     dirflag     = 0; end

% ensure that datachannel is a column vector
datachannel = datachannel(:);

if ischar(channelcmb) && strcmp(channelcmb, 'all')
  % make all possible combinations of all channels
  channelcmb = {'all' 'all'};
end

% it should have a selection of two channels or channelgroups in each row
if size(channelcmb,1)==2 && size(channelcmb,2)~=2
  warning('transposing channelcombination matrix');
  channelcmb = channelcmb';
end

% allow for channelcmb to be a 1x2 cell-array containing cells
if numel(channelcmb)==2 && iscell(channelcmb{1}) && iscell(channelcmb{2})
  channelcmb{1} = ft_channelselection(channelcmb{1}, datachannel);
  channelcmb{2} = ft_channelselection(channelcmb{2}, datachannel);
  n1  = numel(channelcmb{1});
  n2  = numel(channelcmb{2});
  tmp = cell(n1*n2+n1+n2,2);
  for k = 1:n1
    tmp((k-1)*n2+(1:n2), 1) = channelcmb{1}(k);
    tmp((k-1)*n2+(1:n2), 2) = channelcmb{2};
    tmp(n2*k+(1:n1),     1) = channelcmb{1};
    tmp(n2*k+(1:n1),     2) = channelcmb{1};
    tmp(n2*k+n1+(1:n2),  1) = channelcmb{2};
    tmp(n2*k+n1+(1:n2),  2) = channelcmb{2};
  end
  collect = tmp;
  return;
end

if isempty(setdiff(channelcmb(:), datachannel))
  % there is not much to do, since there are no channelgroups with special names
  % each element of the input therefore already contains a proper channel name

  switch dirflag
    case 0
      % nothing to do
    case 1
      % switch the order
      channelcmb = channelcmb(:,[2 1]);
    case 2
      channelcmb = [channelcmb; channelcmb(:,[2 1])];
    otherwise
      error('unknown value for input argument ''dirflag''');
  end
  collect = channelcmb;

  if includeauto
    autochannel = unique(channelcmb(:));
    for ch=1:numel(autochannel)
      collect = cat(1, collect, [autochannel(ch) autochannel(ch)]);
    end
  end

else
  % a combination is made for each row of the input selection after
  % translating the channel group (such as 'all') to the proper channel names
  % and within each set, double occurences and autocombinations are removed

  selmat = false(numel(datachannel));
  for sel=1:size(channelcmb,1)
    % translate both columns and subsequently make all combinations
    channelcmb1 = ft_channelselection(channelcmb(sel,1), datachannel);
    channelcmb2 = ft_channelselection(channelcmb(sel,2), datachannel);

    % translate both columns and subsequently make all combinations
    list1 = match_str(datachannel, channelcmb1);
    list2 = match_str(datachannel, channelcmb2);

    selmat(list1,list2) = true;

    if ~includeauto,
      % exclude the auto-combinations
      selmat = selmat & ~eye(size(selmat));
    else
      % ensure that the appropriate diagonal entries are filled
      autovec = false(numel(datachannel),1);
      autovec([list1(:);list2(:)]) = true;
      selmat = selmat | diag(autovec);
    end

    if dirflag==2,
      % also fill the other 'direction'
      selmat(list2,list1) = true;
    end
  end

  if dirflag<2
    % remove double occurrences
    selmat   = selmat & ~tril(selmat, -1)';
  end
  [i1, i2] = find(selmat);

  switch dirflag
    case 0
      % original behaviour,
      % row-to-column, i.e. outflow according to FT's convention
      indx = [i1 i2];
    case 1
      % column-to-row, i.e. inflow according to FT's convention
      indx = [i2 i1];
    case 2
      % both
      indx = [i1 i2];

    otherwise
      error('unknown value for input argument ''dirflag''');
  end

  collect = [datachannel(indx(:,1)) datachannel(indx(:,2))];
end
