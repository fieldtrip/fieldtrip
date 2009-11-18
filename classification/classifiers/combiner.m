classdef combiner < classifier
%COMBINER combines posteriors of multiple datasets into one
%
%   EXAMPLE:
%
%   % run discriminant analysis on an EEG and an MEG dataset
%
%   myproc = clfproc({ combiner('procedure',{da()}) });
%
%   this will combine the posteriors computed from using a combination rule
%
%   SEE ALSO:
%   combine_posteriors.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: combiner.m,v $
%

    properties
       
        procedure = cproc({da()}); % the used classification procedure
        combination = 'product'; % how to combine classifier output (not how to combine data)
        
    end

    methods
       function obj = combiner(varargin)
           
             obj = obj@classifier(varargin{:});
       end
       function obj = train(obj,data,design)

           if ~iscell(data)
               data = {data};
           end    
           assert(~iscell(design));
           
           if ~iscell(obj.procedure)
               procedure = obj.procedure;
               obj.procedure = cell(1,length(data));
               for j=1:length(data)
                   obj.procedure{j} = procedure;
               end
           end

           for j=1:length(data)
               obj.procedure{j} = obj.procedure{j}.train(data{j},design);
           end

           
       end
       function post = test(obj,data)       
          
           if ~iscell(data), data = {data}; end
           
           cpost = cell(1,length(data));
           for j=1:length(data)
               cpost{j} = obj.procedure{j}.test(data{j});
           end
           
           post = combine_posteriors(cpost,obj.combination);           
      
       end

    end
end
