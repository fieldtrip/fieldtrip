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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
