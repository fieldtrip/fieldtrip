function [boolval, taxis] = checktime(varargin)

% last input is always the required string
tol      = varargin{end};
required = varargin{end-1};
varargin = varargin(1:end-2);

Ndata = numel(varargin);
Ntime = zeros(1,Ndata);
taxis = zeros(1,0);
for i=1:Ndata
  Ntime(i) = numel(varargin{i}.time);
  taxis    = [taxis;varargin{i}.time(:)];
end

if strcmp(required, 'unique')
  boolval = numel(unique(taxis))==numel(taxis) && ~all(isnan(taxis));
  % the second condition is included when the time is set to dummy nan
elseif strcmp(required, 'identical')
  % the number of time bins needs at least to be the same across inputs
  boolval = all(Ntime==Ntime(1));
  if boolval
    % then check whether the axes are equal
    taxis   = reshape(taxis, Ntime(1), []);
    boolval = all(all(abs(taxis - repmat(taxis(:,1), 1, Ndata))<tol)==1);
    taxis   = taxis(:,1);
  end
end

