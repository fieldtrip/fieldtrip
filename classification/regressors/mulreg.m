classdef mulreg < regressor
%MULREG multvariate linear regression method class
%
%   refers to linear regression with multiple responses. Other regressors
%   can be used by overloading the regressors property
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
%
%   $Log: mulreg.m,v $
%
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
        
        function obj = train(obj,data,design)
            
            
            if isempty(obj.ndep)
                obj.ndep = design.nfeatures;
            end                
            
            if isa(obj.prefun,'function_handle')
                design = obj.prefun(design.X);
            end
            
            if isempty(obj.regressors)
              obj.regressors = cell(1,obj.ndep);
              for c=1:obj.ndep
                obj.regressors{c} = linreg();
              end
            end
            
            for c=1:obj.ndep
                obj.regressors{c} = obj.regressors{c}.train(data.X,dataset(design(:,c)));
            end
        end
        
        function res = test(obj,data)                 
          
          data = data.X;
          
          res = zeros(size(data,1),obj.ndep);
          for j=1:obj.ndep
            res(:,j) = obj.regressors{j}.test(data);
          end
          
          if isa(obj.postfun,'function_handle')
            res = obj.postfun(res);
          end

          res = dataset(res);
          
        end
 
    end
end 
