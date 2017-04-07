function test_ft_channelcombination

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_channelcombination

% this function tests the new implementation of ft_channelcombination

load(fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg'),'preproc_ctf151.mat'));
label = data.label; clear data;

x = ft_channelcombination({'all' 'all'}, label);
y = ft_channelcombination_old({'all' 'all'}, label);
assert(isequal(x,y));

x = ft_channelcombination({'MLT' 'all'}, label);
y = ft_channelcombination_old({'MLT' 'all'}, label);
assert(isequal(x,y));

x = ft_channelcombination({'MLT' 'MRC'}, label);
y = ft_channelcombination_old({'MLT' 'MRC'}, label);
assert(isequal(x,y));

x = ft_channelcombination({{'MLC12' 'MLC13'} {'MRO11' 'MRO12' 'MRO21'}}, label);
y = ft_channelcombination_old({{'MLC12' 'MLC13'} {'MRO11' 'MRO12' 'MRO21'}}, label);
assert(isequal(x,y));

% test the new functionality
x = ft_channelcombination({'MLT' 'MRC'}, label, 0, 0);
y = ft_channelcombination({'MLT' 'MRC'}, label, 0, 1);
z = ft_channelcombination({'MLT' 'MRC'}, label, 0, 2);
assert(isequal(x, y(:,[2 1]))); % the columns should be swapped
assert(numel(z)==2*numel(x));


%%%%%%%%%%
% Below is the old code

function [collect] = ft_channelcombination_old(channelcmb, datachannel, includeauto)

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
% Please note that the default behaviour is to exclude symetric
% pairs and auto-combinations.
%
% See also FT_CHANNELSELECTION

% Undocumented local options: optional third input argument includeauto,
% specifies to include the auto-combinations

% Copyright (C) 2003-2011, Robert Oostenveld
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
% $Id: ft_channelcombination.m 10449 2015-06-10 18:34:02Z roboos $

if nargin==2,
  includeauto = 0;
end

if ischar(channelcmb) && strcmp(channelcmb, 'all')
  % make all possible combinations of all channels
  channelcmb = {'all' 'all'};
end

% it should have a selection of two channels or channelgroups in each row
if size(channelcmb,1)==2 && size(channelcmb,2)~=2
  warning('transposing channelcombination matrix');
  channelcmb = channelcmb';
end

% this will hold the output
collect = {};

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
  % there is nothing to do, since there are no channelgroups with special names
  % each element of the input therefore already contains a proper channel name
  collect = channelcmb;
  
  if includeauto
    for ch=1:numel(datachannel)
      collect{end+1,1} = datachannel{ch};
      collect{end,  2} = datachannel{ch};
    end
  end
else
  % a combination is made for each row of the input selection after
  % translating the channel group (such as 'all') to the proper channel names
  % and within each set, double occurences and autocombinations are removed
  
  for sel=1:size(channelcmb,1)
    % translate both columns and subsequently make all combinations
    channelcmb1 = ft_channelselection(channelcmb(sel,1), datachannel);
    channelcmb2 = ft_channelselection(channelcmb(sel,2), datachannel);
    
    % compute indices of channelcmb1 and channelcmb2 relative to datachannel
    [dum,indx,indx1]=intersect(channelcmb1,datachannel);
    [dum,indx,indx2]=intersect(channelcmb2,datachannel);
    
    % remove double occurrences of channels in either set of signals
    indx1   = unique(indx1);
    indx2   = unique(indx2);
    
    % create a matrix in which all possible combinations are set to one
    cmb = zeros(length(datachannel));
    for ch1=1:length(indx1)
      for ch2=1:length(indx2)
        cmb(indx1(ch1),indx2(ch2))=1;
      end
    end
    
    % remove auto-combinations
    cmb = cmb & ~eye(size(cmb));
    
    % remove double occurences
    cmb = cmb & ~tril(cmb, -1)';
    
    [indx1,indx2] = find(cmb);
    
    % extend the previously allocated cell-array to also hold the new
    % channel combinations (this is done to prevent memory allocation and
    % copying in each iteration in the for-loop below)
    num = size(collect,1);               % count the number of existing combinations
    dum = cell(num + length(indx1), 2);  % allocate space for the existing+new combinations
    if num>0
      dum(1:num,:) = collect(:,:);       % copy the exisisting combinations into the new array
    end
    collect = dum;
    clear dum
    
    % convert to channel-names
    for ch=1:length(indx1)
      collect{num+ch,1}=datachannel{indx1(ch)};
      collect{num+ch,2}=datachannel{indx2(ch)};
    end
  end
  
  if includeauto
    cmb           = eye(length(datachannel));
    [indx1,indx2] = find(cmb);
    num           = size(collect,1);
    dum           = cell(num + length(indx1), 2);
    if num>0,
      dum(1:num,:) = collect(:,:);
    end
    collect = dum;
    clear dum
    
    % convert to channel-names for the auto-combinations
    for ch=1:length(indx1)
      collect{num+ch,1} = datachannel{indx1(ch)};
      collect{num+ch,2} = datachannel{indx2(ch)};
    end
  end
end

