function [collect] = channelcombination(channelcmb, datachannel, includeauto)

% CHANNELCOMBINATION creates a cell-array with combinations of EEG/MEG
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
% See also CHANNELSELECTION

% Undocumented local options:
% optional third input argument includeauto, specifies to include the 
% auto-combinations

% Copyright (C) 2003-2006, Robert Oostenveld
%
% $Log: channelcombination.m,v $
% Revision 1.16  2009/10/19 13:22:38  jansch
% added undocumented possibility to append auto-combinations to the list,
% which is necessary for functionality in connectivityanalysis
%
% Revision 1.15  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.14  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.13  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.12  2006/03/06 11:52:52  roboos
% prevent memory allocation and copying in a for loop, this was the
% reason for it being tremendously slow for a large number of channel
% combinations
%
% Revision 1.11  2005/09/29 00:46:11  roboos
% renamed the variable sgncmb into channelcmb, both in the help and code
%
% Revision 1.10  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.9  2005/03/18 09:25:50  roboos
% changed handling of multiple rows of the input selection, each row is now processed separately
% removed (defunct) planar channel detection
%
% Revision 1.8  2004/10/19 10:44:43  roboos
% fixed bug in planar-gradient channel detection if labels are shorter than 3 chars
%
% Revision 1.7  2004/08/26 07:16:11  jansch
% re-added to CVS-repository.
%
% Revision 1.5  2004/01/16 08:28:50  roberto
% fixed small bug for 'all'
%
% Revision 1.4  2004/01/06 12:43:47  roberto
% redesign by JM
%
% Revision 1.3  2003/10/29 09:02:31  roberto
% fixed bug in varible names, sgncmb vs channel
%
% Revision 1.2  2003/10/27 16:24:06  roberto
% added check for double channel combinations
%
% Revision 1.1  2003/10/27 16:02:13  roberto
% new implementation after an idea by JM
%

fieldtripdefs

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
end

% this will hold the output
collect = {};

if isempty(setdiff(channelcmb(:), datachannel))
  % there is nothing to do, since there are no channelgroups with special names
  % each element of the input therefore already contains a proper channel name
  collect = channelcmb;

else
  % a combination is made for each row of the input selection after
  % translating the channel group (such as 'all') to the proper channel names
  % and within each set, double occurences and autocombinations are removed

  for sel=1:size(channelcmb,1)
    % translate both columns and subsequently make all combinations
    channelcmb1 = channelselection(channelcmb(sel,1), datachannel);
    channelcmb2 = channelselection(channelcmb(sel,2), datachannel);

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
