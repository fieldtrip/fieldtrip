function argout = cellfun(fname, varargin)

% This is a replacement function for the cellfun in Matlab version 7.0 (R14)
% for which the built-in cellfun has too limited functionality.

UniformOutput = true;

numargin = length(varargin);
njobs    = length(varargin{1});

if UniformOutput
  argout = zeros(size(varargin{1}));
else
  argout = cell(size(varargin{1}));
end

for i=1:njobs
  % redistribute the input arguments
  argin = cell(1, length(varargin));
  for j=1:numargin
    argin{j} = varargin{j}{i};
  end
  % apply the function and collect the results
  if UniformOutput
    argout(i) = feval(fname, argin{:});
  else
    argout{i} = feval(fname, argin{:});
  end
end

