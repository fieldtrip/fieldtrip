classdef mulreg < regressor
%MULREG multvariate linear regression method class
%
%   refers to linear regression with multiple responses. Other regressors
%   can be used by overloading the regressors property
%
%   PARAMETERS: 
%   regressors % regressor objects    
%
%   EXAMPLE:
%   
%   mulreg('regressors',{gpregressor() gpregressor()},'prefun',@(x)([sin(x) cos(x)]),'postfun',@(x)(angle(1i*x(:,1)+x(:,2))));
%
%   will take the sin and cos of the design matrix and try to regress on
%   that using gaussian processes; 'postfun' will recombine the multiple responses; here to an
%   angle.
%
%   Copyright (c) 2008, Marcel van Gerven

    properties        
        
        regressors % regressor objects    
      
        ndep % number of dependent variables (relevant n columns of design matrix)
        
        prefun  % optional preprocessing of design matrix to create multiple inputs
        postfun % optional postprocessing of posteriors to create one output
    end
    
    methods
        
        function obj = mulreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function p = estimate(obj,data,design)
                        
            if isempty(obj.ndep)
              p.ndep = design.nfeatures;
            else
              p.ndep = obj.ndep;
            end
            
            if isa(obj.prefun,'function_handle')
                design = obj.prefun(design.X);
            end
            
            if isempty(obj.regressors)
              p.regressors = cell(1,obj.ndep);
              for c=1:obj.ndep
                p.regressors{c} = linreg();
              end
            else
              p.regressors = obj.regressors;
            end
            
            for c=1:obj.ndep
                p.regressors{c} = p.regressors{c}.train(data,dataset(design.X(:,c)));
            end
        end
        
        function res = map(obj,data)                 
          
          res = zeros(data.nsamples,obj.params.ndep);
          for j=1:obj.ndep
            res(:,j) = obj.params.regressors{j}.test(data);
          end
          
          if isa(obj.postfun,'function_handle')
            res = obj.postfun(res);
          end

          res = dataset(res);
          
        end
 
    end
end 
