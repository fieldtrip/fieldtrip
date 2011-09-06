function host = gethostname()

% HOSTNAME
%
% Use as
%   str = hostname;

% Copyright 2011, Eelke Spaak

% this is to speed up subsequent calls
persistent previous_argout
if ~isempty(previous_argout)
  host = previous_argout;
  return
end

if (ispc())
  host = getenv('ComputerName');
elseif (isunix())
  host = getenv('HOSTNAME');
  % the HOSTNAME variable is not guaranteed to be set
  if isempty(host)
    [status, host] = system('hostname -s');
    % remove the ^n from the end of the line
    host = deblank(host);
  end
end

host = strtok(host, '.'); % dots in filenames are not allowed by matlab

if (isempty(host))
  host = 'unknownhost';
end

% remember for subsequent calls
previous_argout = host;

