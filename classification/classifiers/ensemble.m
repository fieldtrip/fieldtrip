classdef ensemble < classifier
%ENSEMBLE combines multiple classifications of one dataset into one
%
%   EXAMPLE:
%
%   % run discriminant analysis and svm on one dataset and combine
%   posteriors using majority voting
%
%   myproc = clfproc({ensemble('procedures',{clfproc({da()}) clfproc(kernelmethod())},'combination','majority')});   
%
%   SEE ALSO:
%   ensemble_example
%   combine_posteriors
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: ensemble.m,v $
%

    properties
       
        procedures = {clfproc({da()})}; % the used classification procedures
        combination = 'voting'; % how to combine classifier output (not how to combine data)
        
    end

    methods
       function obj = ensemble(varargin)
           
             obj = obj@classifier(varargin{:});
             
             if ~iscell(obj.procedures), obj.procedures = { obj.procedures }; end
       end
       function obj = train(obj,data,design)

           [data,design] = obj.check_input(data,design);
           
           % multiple classifiers per dataset
           for k=1:length(obj.procedures)
               obj.procedures{k} = obj.procedures{k}.train(data,design);
           end

       end
       function post = test(obj,data)       
          % follow the structure of the classifiers;
          % retrieve posteriors
          % combine using some combination rule
                      
          data = obj.check_input(data);
          
          % multiple classifiers per dataset
          if iscell(obj.procedures)
              cpost = cell(1,length(obj.procedures));
              for k=1:length(obj.procedures)
                  cpost{k} = obj.procedures{k}.test(data);
              end

          end
          
          post = combine_posteriors(cpost,obj.combination);
           
       end
              
    end
    
end
