classdef predictor < mvmethod
%PREDICTOR base class for classifier and regressor
% 
% NOTE:
% optimizer can also act as a predictor but is not of this class
%
% Copyright (c) 2009, Marcel van Gerven
  
    
  methods(Abstract)
    p = predict(obj,X);
  end
  
end
