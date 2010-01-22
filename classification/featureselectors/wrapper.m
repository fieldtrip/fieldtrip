classdef wrapper < featureselector
%WRAPPER wrapper approach to feature selection
%
%   The wrapper uses a validator to determine the optimal
%   set of features. It creates all the subsets up to maxfeatures. The
%   search strategy determines how we traverse the space of subsets.
%   Currently all strategies are forward selection methods since we assume
%   backward is too expensive for large datasets
%
%   Note:
%   'bestfirst' search currently only works for medium sized datasets
%
%   Todo:
%   rewrite bestfirst strategy
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: wrapper.m,v $
%

    properties
       
      validator = crossvalidator('procedure',clfproc({nb()}),'cvfolds',0.9);
      metric = 'accuracy'; % evaluation metric
      criterion; % keeps track of the evaluation metric
      
      maxfeatures = Inf; % maximum number of used features
      
      search = 'hillclimbing'; % type of search strategy (maxfeatures determines evaluations for hillclimbing)
      expansions = 0; % maximum number of expansions without improvement for bestfirst search
      epsilon = 0.01; % minimum relative improvement for acceptance in bestfirst search

    end

    methods
      
        function obj = wrapper(varargin)
          
          obj = obj@featureselector(varargin{:});
          
        end
        
        function obj = train(obj,data,design)
          
          data = data.collapse();
          design = design.collapse();
          
          switch obj.search
            
            case 'hillclimbing'
              
              cur = [];
              metric = -Inf;
              nfeat = min(obj.maxfeatures,size(data,2));
              
              candidates = 1:size(data,2);
              obj.criterion = zeros(1,nfeat);
              for f=1:nfeat
                
                if obj.verbose
                  fprintf('evaluating %d out of %d features; ',f,nfeat);
                end
                
                oldmetric = metric;
                for j=candidates
                  
                  tmp = [cur j];
                  
                  cv = obj.validator.validate(dataset(data(:,tmp)),dataset(design));
                  m = validator.eval(cv.post,cv.design,'metric',obj.metric);
                  
                  if m > metric
                    subset = tmp;
                    metric = m;
                  end
                end
                
                if obj.verbose
                  fprintf('%.2g\n',metric)
                end
                obj.criterion(f) = metric;
                
                if metric > oldmetric % better candidate found
                  cur = subset;
                else
                  %break; % no improvement so break out of loop
                end
              end
              
              obj.subset = cur;
              
            case 'bestfirst'
              
              nfeat = size(data,2);
              
              % we represent features as binary numbers
              open = 0;
              closed = [];
              
              % performance of the zero subset is assumed infinitely
              % bad
              metric = -Inf;
              openacc = metric;
              
              best = 0;
              v = best;
              lastchange = 0;
              while lastchange <= obj.expansions
                
                if obj.verbose
                  fprintf('evaluating candidate subset; ');
                end
                
                % find children
                candidates = obj.expand(v,nfeat);
                
                vb = (open == v);
                open = open(~vb);
                openacc = openacc(~vb);
                closed = [closed v];
                
                % remove candidates in closed or open
                candidates = setdiff(candidates,open);
                candidates = setdiff(candidates,closed);
                
                acc = zeros(1,length(candidates));
                for j=1:length(candidates)
                  
                  cv = obj.validator.validate(dataset(data(:,boolean(bitget(candidates(j),1:nfeat)))),dataset(design));
                  acc(j) = validator.eval(cv.post,cv.design,'metric',obj.metric);
                end
                
                open = [open candidates];
                openacc = [openacc acc];
                
                [a,b] = max(openacc);
                v = open(b);
                
                if a > (1+obj.epsilon)*metric
                  metric = a;
                  best = v;
                  
                  lastchange = 0;
                else
                  lastchange = lastchange + 1;
                end
                
                if obj.verbose
                  fprintf('%.2g\n',metric)
                end
                
              end
              
              obj.subset = find(bitget(best,1:nfeat));
              
            otherwise
              error('unknown search strategy');
              
          end
          
        end
        
    end
    
    methods(Access = private)
      
      function candidates = expand(obj,v,nfeat)
        
        % create children of v
        b = boolean(bitget(v,1:nfeat));
        
        cnd = find(~b);
        
        candidates = zeros(1,length(cnd));
        for c=1:length(cnd)
          
          % make candidate
          tmp = b;
          tmp(cnd(c)) = true;
          
          % convert binary candidate back to decimal
          cp = fliplr([1 cumprod(2*ones(1,length(tmp)-1))]);
          candidates(c) = sum(cp(boolean(tmp)));
        end
      end
      
    end
end

