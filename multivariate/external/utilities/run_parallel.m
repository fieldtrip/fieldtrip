function res = run_parallel(obj,funname,varargin)
% allow peercellfun to run an arbitrary class method or return the input

  try
    res = obj.(funname)(varargin{:});
  catch
    res = cat(2,{class(obj) funname},varargin);
  end
  
end