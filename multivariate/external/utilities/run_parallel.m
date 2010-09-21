function res = run_parallel(obj,funname,varargin)
% allow peercellfun to run an arbitrary class method 
  res = obj.(funname)(varargin{:});
end