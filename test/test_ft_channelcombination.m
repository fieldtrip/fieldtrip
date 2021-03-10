function test_ft_channelcombination

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_channelcombination

% this function tests the new implementation of ft_channelcombination


label = {'STIM'
    'SCLK01'
    'BG1'   
    'BG2'   
    'BG3'   
    'BP1'   
    'BP2'   
    'BP3'   
    'BQ1'   
    'BQ2'   
    'BQ3'   
    'BR1'   
    'BR2'   
    'BR3'   
    'G11'   
    'G12'   
    'G13'   
    'G22'   
    'G23'   
    'P11'   
    'P12'   
    'P13'   
    'P22'   
    'P23'   
    'Q11'   
    'Q12'   
    'Q13'   
    'Q21'   
    'Q22'   
    'Q23'   
    'R11'   
    'R12'   
    'R13'   
    'R22'   
    'R23'   
    'MLC11' 
    'MLC12' 
    'MLC13' 
    'MLC14' 
    'MLC15' 
    'MLC21' 
    'MLC22' 
    'MLC23' 
    'MLC24' 
    'MLC31' 
    'MLC32' 
    'MLC33' 
    'MLC41' 
    'MLC42' 
    'MLC43' 
    'MLF11' 
    'MLF12' 
    'MLF21' 
    'MLF22' 
    'MLF23' 
    'MLF31' 
    'MLF32' 
    'MLF33' 
    'MLF34' 
    'MLF41' 
    'MLF42' 
    'MLF43' 
    'MLF44' 
    'MLF45' 
    'MLF51' 
    'MLF52' 
    'MLO11' 
    'MLO12' 
    'MLO21' 
    'MLO22' 
    'MLO31' 
    'MLO32' 
    'MLO33' 
    'MLO41' 
    'MLO42' 
    'MLO43' 
    'MLP11' 
    'MLP12' 
    'MLP13' 
    'MLP21' 
    'MLP22' 
    'MLP31' 
    'MLP32' 
    'MLP33' 
    'MLP34' 
    'MLT11' 
    'MLT12' 
    'MLT13' 
    'MLT14' 
    'MLT15' 
    'MLT16' 
    'MLT21' 
    'MLT22' 
    'MLT23' 
    'MLT24' 
    'MLT25' 
    'MLT26' 
    'MLT31' 
    'MLT32' 
    'MLT33' 
    'MLT34' 
    'MLT35' 
    'MLT41' 
    'MLT42' 
    'MLT43' 
    'MLT44' 
    'MRC11' 
    'MRC12' 
    'MRC13' 
    'MRC14' 
    'MRC15' 
    'MRC21' 
    'MRC22' 
    'MRC23' 
    'MRC24' 
    'MRC31' 
    'MRC32' 
    'MRC33' 
    'MRC41' 
    'MRC42' 
    'MRC43' 
    'MRF11' 
    'MRF12' 
    'MRF21' 
    'MRF22' 
    'MRF23' 
    'MRF31' 
    'MRF32' 
    'MRF33' 
    'MRF34' 
    'MRF41' 
    'MRF42' 
    'MRF43' 
    'MRF44' 
    'MRF45' 
    'MRF51' 
    'MRF52' 
    'MRO11' 
    'MRO12' 
    'MRO21' 
    'MRO22' 
    'MRO31' 
    'MRO32' 
    'MRO33' 
    'MRO41' 
    'MRO42' 
    'MRO43' 
    'MRP11' 
    'MRP12' 
    'MRP13' 
    'MRP21' 
    'MRP22' 
    'MRP31' 
    'MRP32' 
    'MRP33' 
    'MRP34' 
    'MRT11' 
    'MRT12' 
    'MRT13' 
    'MRT14' 
    'MRT15' 
    'MRT16' 
    'MRT21' 
    'MRT22' 
    'MRT23' 
    'MRT24' 
    'MRT25' 
    'MRT26' 
    'MRT31' 
    'MRT32' 
    'MRT33' 
    'MRT34' 
    'MRT35' 
    'MRT41' 
    'MRT42' 
    'MRT43' 
    'MRT44' 
    'MZC01' 
    'MZC02' 
    'MZF01' 
    'MZF02' 
    'MZF03' 
    'MZO01' 
    'MZO02' 
    'MZP01' 
    'MZP02'
    'EOG'};

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
assert(~isequal(x,y)); % the correct behavior has changed since the update 
% PR 1664 to ft_channelcombination.
% the default behavior should be not to return the auto-combinations, yet
% the old behavior in this case was to actually return the
% auto-combinations, which is incorrect, resulting in a y of 11x2, rather
% then a 6x2

% test the new functionality
x = ft_channelcombination({'MLT' 'MRC'}, label, 0, 0);
y = ft_channelcombination({'MLT' 'MRC'}, label, 0, 1);
z = ft_channelcombination({'MLT' 'MRC'}, label, 0, 2);
assert(isequal(x, y(:,[2 1]))); % the columns should be swapped
assert(numel(z)==2*numel(x));

lab{1,1} = {'a';'b';'c'};
lab{2,1} = {'b';'c';'d'};
lab{3,1} = {'d';'e'};

opts1 = [0 1]; %includeauto
opts2 = [0 1 2]; %dirflag

labx{1} = [1 2 3];
labx{2} = [2 3 4];
labx{3} = [4 5];
cmb     = false(5,5);
dat = zeros(0,5);
for k = 1:numel(lab)
  for m = 1:numel(opts1)
    for p = 1:numel(opts2)
      
      % this implements what I think should be the behavior of
      % ft_channelcombination
      cmb(:) = false;
      cmb(labx{1},labx{k}) = true;
      if opts2(p)<2
        ix = intersect(labx{1},labx{k});
        cmb(ix,ix) = cmb(ix,ix) & tril(ones(numel(ix)),-1)>0;
      end

      if opts2(p)==2
        cmb(labx{k},labx{1}) = true;
      elseif opts2(p)==1
        cmb = cmb';
      end

      if opts1(m)==0
        cmb = cmb & (ones(5)-eye(5))>0;
      elseif opts1(m)==1
        tmp = zeros(5,1);
        tmp(union(labx{1},labx{k})) = 1;
        cmb = cmb | diag(tmp)>0;
      end
      %imagesc(cmb);pause;

      X{k,m,p} = ft_channelcombination({lab{1},lab{k}},union(lab{1},lab{k}), opts1(m), opts2(p));
      
      dat(end+1,:) = [k m p size(X{k,m,p},1) sum(cmb(:))];
    
      if k==1
        % with fully overlapping lists of labels, ft_channelcombination
        % behaves as expected
        assert(isequal(sum(cmb(:)),size(X{k,m,p},1)));
      elseif k==2
        assert(isequal(sum(cmb(:)),size(X{k,m,p},1)));
      elseif k==3
        assert(isequal(sum(cmb(:)),size(X{k,m,p},1)));
      end
      
    end
  end
end

%%%%%%%%%%
% Below is the old code

function [collect] = ft_channelcombination_old(channelcmb, datachannel, includeauto)

% FT_CHANNELCOMBINATION creates a cell-array with combinations of EEG/MEG
% channels for subsequent cross-spectral-density and coherence analysis
%
% You should specify channel combinations as a two-column cell-array,
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
% Please note that the default behavior is to exclude symetric
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

