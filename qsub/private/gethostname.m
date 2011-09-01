function host = gethostname()

% HOSTNAME
%
% Use as
%   str = hostname;

% Copyright 2011, Eelke Spaak

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

