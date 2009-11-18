classdef gp < classifier
%GP gaussian process classifier
%
%   Options:
%   'optimize' : if true tries to optimize hyperparameters
%
%   SEE ALSO:
%   gpml-matlab
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: gp.m,v $
%

    properties

        data;
        targets;
        method = 'laplace';
        optimize = true;
        loghyper = [3.0; 0.0];
        logmarglik;
    end

    methods
       function obj = gp(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)

         [data,design] = obj.check_input(data,design);
            
           if ~exist('gpml-matlab','dir')
               error('this code requires an external toolbox: http://www.gaussianprocess.org/gpml/code/matlab/doc/');
           end
           
           if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
                   
           if obj.nclasses ~= 2, error('GP only accepts binary class problems'); end

           targets = design(:,1);
           targets(design == 1) = -1;
           targets(design == 2) = 1;

           % training mode

           obj.data = data;
           obj.targets = targets;

           if obj.optimize % optimize hyperparameters

               if strcmp(obj.method,'laplace')

                   [obj.loghyper obj.logmarglik] = minimize(obj.loghyper, 'binaryLaplaceGP', -20, 'covSEiso', 'cumGauss', obj.data, obj.targets);

               else % ep

                   [obj.loghyper obj.logmarglik] = minimize(obj.loghyper, 'binaryEPGP', -20, 'covSEiso', obj.data, obj.targets);

               end
           end
                       
       end
       function post = test(obj,data)       

         data = obj.check_input(data);

           if strcmp(obj.method,'laplace')

               probs = binaryLaplaceGP(obj.loghyper, 'covSEiso', 'cumGauss', obj.data, obj.targets, data);

           else % ep

               probs = binaryEPGP(obj.loghyper, 'covSEiso', obj.data, obj.targets, data);

           end

           post = zeros(size(probs,1),2);
           post(:,1) = 1 - probs;
           post(:,2) = probs;

       end

    end
end 
