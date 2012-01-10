function id = generatesessionid()

% GENERATESESSIONID
%
% See also GENERATEJOBID, GENERATEBATCHID

if nargin~=0
  error('incorrect number of input arguments');
end

id = sprintf('%s_%s_p%d', getusername(), gethostname(), getpid());
