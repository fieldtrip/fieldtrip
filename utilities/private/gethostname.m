function host = gethostname()

% this is to speed up subsequent calls
% the hostname will not change over multiple calls
persistent previous_argout
if ~isempty(previous_argout)
  host = previous_argout;
  return;
end

if (ispc())
  host = getenv('COMPUTERNAME');
else
  host = getenv('HOSTNAME');
  if isempty(host)
    [status, host] = system('hostname -s');
  end
end

% remove trailing whitespace and blank lines
host = strtrim(host);

if (isempty(host))
  host = 'unknown';
end

% remember for subsequent calls
previous_argout = host;
