function id = generatebatchid(batch)

% GENERATEBATCHID
%
% See also GENERATEJOBID, GENERATESESSIONID

if nargin~=1
  error('incorrect number of input arguments');
end

id = sprintf('%s_%s_p%d_b%d', getusername(), gethostname(), getpid(), batch);
