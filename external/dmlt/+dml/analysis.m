classdef analysis
% ANALYSIS multivariate analysis class.
%   
%   DESCRIPTION
%   A multivariate analysis contains a cell array of multivariate methods 
%   {method1 method2 method3 ...} that are called in this order and where 
%   the output of the previous method acts as input to the next method.
%
%   DEVELOPER 
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

    properties
      
      % the methods that specify the multivariate analysis
      % methods callable using cell array notation; e.g. a{2}
      method
        
    end
    
    methods
      
       function obj = analysis(mvmethods)
       % constructor expects mva methods

       % return the same thing if input is already an analysis
       if isa(mvmethods,'mva.analysis')
          obj = mvmethods;
          return;
        end
        
        % cast to cell array if only one method is specified
        if ~iscell(mvmethods), mvmethods = {mvmethods}; end
        
        obj.method = mvmethods;
       
       end     
       
       function obj = train(obj,X,Y)
         % train just calls the methods' train functions in order to
         % produce an output
         
         if nargin<2, error('input X expected'); end
         if nargin<3, Y = []; end % Y may be absent
      
         for c=1:length(obj.method)

           obj.method{c} = obj.method{c}.train(X,Y);
           
           if c<length(obj.method)
             X = obj.method{c}.test(X);
           end
             
         end
         
       end
       
       function Y = test(obj,X)
         % test just calls the methods' test functions in order to produce a posterior
         
         if nargin<2, error('input X expected'); end
         
         Y = X;
         for c=1:length(obj.method)
           Y = obj.method{c}.test(Y);
         end
         
       end
       
       
       function m = model(obj)
         % get the models belonging to the last method
         
         m = obj.method{end}.model();
         
       end
       
    end
    
end
