function [boolval, pos] = checkpos(varargin)

% last input is always the required string
tol      = varargin{end};
required = varargin{end-1};
varargin = varargin(1:end-2);

Ndata = numel(varargin);
Npos  = zeros(1,Ndata);
pos   = zeros(0,3);
for i=1:Ndata
  Npos(i) = size(varargin{i}.pos,1);
  pos     = [pos;varargin{i}.pos];
end

if strcmp(required, 'unique')
  boolval = size(unique(pos, 'rows'),1)==size(pos,1) && ~all(isnan(pos(:)));
  % the second condition is included when the pos is set to dummy nans,
  % FIXME does this ever happen
elseif strcmp(required, 'identical')
  % the number of positions needs at least to be the same across
  % inputs
  boolval = all(Npos==Npos(1));
  if boolval
    % then check whether the positions are equal up to a certain tolerance,
    % and assuming the same order of the positions
    pos     = reshape(pos, [Npos(1), Ndata 3]);
    boolval = all(all(all(abs(pos - repmat(pos(:,1,:), [1, Ndata, 1]))<tol)==1));
    pos     = reshape(pos(:,1,:), [Npos(1) 3]);
  end
end

