function c = appendevent(a, b)

% APPENDEVENT

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if isempty(a)
  c = b(:);
elseif isempty(b)
  c = a(:);
else
  c = a(:);
  for i=1:numel(b)
    c(end+1).type     = b(i).type;
    c(end  ).value    = b(i).value;
    c(end  ).sample   = b(i).sample;
    if isfield(b, 'timestamp')
      c(end  ).timestamp = b(i).timestamp; % optional
    end
    if isfield(b, 'offset')
      c(end  ).offset    = b(i).offset;    % optional
    end
    if isfield(b, 'duration')
      c(end  ).duration  = b(i).duration;  % optional
    end
  end
end
