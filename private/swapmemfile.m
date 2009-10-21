function [data] = swapmemfile(data);

% SWAPMEMFILE swaps a variable from file into memory and clears it 
% again from the memory on the subsequent call
%
% Use with extreme caution!

% This function depends on PREPARE_TIMEFREQ_DATA which has the following options:
% cfg.avgoverchan
% cfg.avgoverfreq
% cfg.avgovertime
% cfg.channel            
% cfg.channelcmb         
% cfg.datarepresentation
% cfg.frequency
% cfg.latency            
% cfg.precision
% cfg.previous
% cfg.version

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: swapmemfile.m,v $
% Revision 1.4  2006/04/10 16:35:20  ingnie
% updated documentation
%
% Revision 1.3  2005/05/17 17:50:50  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.2  2004/11/08 11:37:36  roboos
% switched file detection from Matlab function "matfinfo" to own function "filetype"
% since matfinfo caused troubles between different matlab versions
%
% Revision 1.1  2004/10/29 11:37:56  roboos
% new implementation, required for prepare_timefreq_data on very large datasets that do not fit into memory simultaneously
%

% this variable will be empty at the first call
persistent file

if nargin<1 || isempty(data)
  % (re)initialize by clearing the persistent variable
  file = [];
  data = [];
elseif ~isempty(file)
  % data is present in memory and originates from a file
  % remove data from memory and replace it by the filename
  fprintf('swapping %s out of memory\n', file);
  data = file;
  file = [];
elseif isempty(file) && (~ischar(data) || ~filetype(data, 'matlab'))
  % data is already present in memory, and does not originate from a file
  % do nothing
elseif isempty(file) && ischar(data) && filetype(data, 'matlab')
  % data is not present in memory but should be read from a Matlab file
  % remember the filename from which the data originates
  file = data;
  fprintf('swapping %s into memory\n', file);
  data = load(file);
  % select the variable of interest, discard the rest
  str = fieldnames(data);
  if length(str)==1
    % select the single variable that is inside the Matlab file
    data = getfield(data, str{1});
  else
    error('Matlab file should contain only one variable');
  end
else
  error('Unknown combination of input variable and file history');
end
