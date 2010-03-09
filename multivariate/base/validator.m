classdef validator
    
  properties
    
    procedure;
    
    post;   % posteriors
    design  % associated class labels
    
    init        = 5;    % initialize RNG
    balanced    = false; % balance classes
    compact     = false; % only retain necessary results
    verbose     = false; % verbose output
    
    
  end
  
  methods
    function obj = validator(varargin)
      % a procedure is a mva
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      
      if ~isempty(obj.procedure)
      
        if ~isa(obj.procedure,'mva')
          % try to create procedure if input is a cell array or a predictor
          obj.procedure = mva(obj.procedure);
        end
      
        if obj.verbose
          fprintf('creating validator for mva %s\n',obj.procedure.name);
        end
      
      end
        
      % initialize random number generator
      if ~isempty(obj.init)
        if obj.verbose
          fprintf('initializing random number generator with seed %d\n',obj.init)
        end
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',obj.init));
      end
      
    end
    
%     function n = nclasses(obj)
%       % return number of classes when known (called by statistics_crossvalidate)
%       
%       n = obj.getpredictor().nclasses;
%       
%     end
    
    function [m,desc] = getmodel(obj)
      % try to return the classifier parameters as a model
          
      if ~iscell(obj.procedure)
        obj.procedure = {obj.procedure};
      end
      
      if isempty(obj.procedure{1})
        
        m = {};
        desc = {};
        return;
      
      end
      
      fm = cellfun(@(x)(x.getmodel()),obj.procedure,'UniformOutput', false);
      
      [tmp,desc] = obj.procedure{1}.getmodel();
      
      m = fm{1};
      for c=2:length(obj.procedure)
        
        for j=1:length(m)
          m{j} = m{j} + fm{c}{j};
        end
      end
      
      % take the mean of the parameters
      if iscell(m)
        for c=1:length(m)
          m{c} = m{c}./length(obj.procedure);
        end
      else
        m = m./length(obj.procedure);
      end
      
    end
    
    
    function result = evaluate(obj,varargin)
      % indirect call to eval to simplify the interface
      
      % create the concatenation of all folds for each of the datasets
     
      post = obj.post;
      design = obj.design;
      
      tpost = cell(1,size(post,2));
      tdesign = cell(1,size(design,2));
      for c=1:length(tpost)
        tpost{c} = cat(1,post{:,c});
        tdesign{c} = cat(1,design{:,c});
      end
      
      result = cell(1,length(tpost));
      
      for c=1:length(tpost)
        result{c} = obj.getpredictor().evaluate(tpost{c},tdesign{c},varargin{:});   
      end
      
      if length(result) == 1
        result = result{1};
      end
      
    end
    
    function p = significance(obj,varargin)
      
      % create the concatenation of all folds for each of the datasets
     
      post = obj.post;
      design = obj.design;
      
      tpost = cell(1,size(post,2));
      tdesign = cell(1,size(design,2));
      for c=1:length(tpost)
        tpost{c} = cat(1,post{:,c});
        tdesign{c} = cat(1,design{:,c});
      end
      
      p = obj.getpredictor().significance(tpost,tdesign,varargin{:});

    end
      
      
    
  end
  
  methods(Access = protected)
 
    
    function p = getpredictor(obj)
       % get the predictor (ending) mvmethod

       if iscell(obj.procedure)
         if iscell(obj.procedure{1}.mvmethods{end})
           p = obj.procedure{1}.mvmethods{end}{1};
         else
           p = obj.procedure{1}.mvmethods{end};
         end
       else
         if iscell(obj.procedure.mvmethods{end})
           p = obj.procedure.mvmethods{end}{1};
         else
           p = obj.procedure.mvmethods{end};
         end
       end
      
    end
    
  end
 
  
end
