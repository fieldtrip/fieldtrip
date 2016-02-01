function [boolval, faxis] = checkfreq(varargin)

% last input is always the required string
tol      = varargin{end};
required = varargin{end-1};
varargin = varargin(1:end-2);

Ndata = numel(varargin);
Nfreq = zeros(1,Ndata);
faxis = zeros(1,0);
for i=1:Ndata
  Nfreq(i) = numel(varargin{i}.freq);
  faxis    = [faxis;varargin{i}.freq(:)];
end

if strcmp(required, 'unique')
  boolval = numel(unique(faxis))==numel(faxis) && ~all(isnan(faxis));
  % the second condition is included when the freq is set to dummy nan
elseif strcmp(required, 'identical')
  % the number of freq bins needs at least to be the same across
  % inputs
  boolval = all(Nfreq==Nfreq(1));
  if boolval
    % then check whether the axes are equal
    faxis   = reshape(faxis, Nfreq(1), []);
    boolval = all(all(abs(faxis - repmat(faxis(:,1), 1, Ndata))<tol)==1);
    faxis   = faxis(:,1);
  end
end

