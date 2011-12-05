function host = gethostname()

% this is to speed up subsequent calls
persistent previous_argout
if ~isempty(previous_argout)
  host = previous_argout;
  return;
end

if (ispc())
  host = getenv('COMPUTERNAME');
else
  host = getenv('HOSTNAME');
end

if (isempty(host))
  host = 'unknown';
end

% remember for subsequent calls
previous_argout = host;
