classdef graphicalmodel
%GRAPHICALMODEL abstract graphical model class
%
%   A graphicalmodel is a collection of factors. Extra arguments are:
%
%   'include' : which factors to include in parameter learning
%   'ec' : the equivalence class of each factor; this is used to couple
%       expected sufficient statistics.
%
%   SEE ALSO:
%   bayesnet.m
%   markovnet.m
%   dbnet.m
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: graphicalmodel.m,v $
%

   properties
              
       g         % graph as an adjacency matrix       
       factors   % the factors
       include   % factors included in parameter learning       
       emiter = 10; % number of EM iterations       
   end
   
   methods
       function obj = graphicalmodel(factors,varargin)
          
           obj.factors = factors;          
           obj.include = 1:obj.length();
           
           % override defaults
           for i=1:2:length(varargin)
               obj.(varargin{i}) = varargin{i+1};          
           end                      
           
       end
       function len = length(obj)
           % number of factors
           
           len = length(obj.factors);
       end
       function sz = size(obj,i)
           % size of each factor (or factor i if defined)
           
           if nargin > 1
               sz = obj.factors{i}.states; 
           else
               sz = cell2mat(cellfun(@(x)(x.states),obj.factors,'UniformOutput',false));
           end
       end
       function d = discrete(obj,idx)
           % return discrete variables as a logical array

           d = cellfun(@(x)(isa(x,'discrete_cpd') || isa(x,'discrete_pot')),obj.factors);
           if nargin > 1
              d = d(idx); 
           end
       end
       function c = continuous(obj,idx)
           % return continuous variables as a logical array

           c = cellfun(@(x)(isa(x,'continuous_cpd') || isa(x,'continuous_pot')),obj.factors);                      
           if nargin > 1
               c = c(idx);
           end
       end
       function obj = learn_parameters(obj,data)
                                 
           % learn the parameters of each factor and take MAP value

           % learn from complete data
           ucomplete = false(1,obj.length());
           for i=obj.include
              
               complete = data(:,obj.factors{i}.domain);
               complete = complete(~any(isnan(complete),2),:);

               if ~isempty(complete)
                   
                   ucomplete(i) = true;
                   obj.factors{i} = obj.factors{i}.update(complete);
               end
           end

           for i=find(ucomplete), obj.factors{i} = obj.factors{i}.maximize(); end

           % learn using EM in case of incomplete data
           
           incomplete = data(any(isnan(data),2),:);
           if ~isempty(incomplete)

              % number of examples
              szinc = size(incomplete,1);
           
              % create an inference engine
              if any(obj.continuous())
                  engine = canonical_jtree_ie(obj);
              else
                  engine = discrete_jtree_ie(obj);
              end
              
              niter = 1;
              while niter <= obj.emiter

                  fprintf('EM iteration %d of %d\n',niter,obj.emiter);

                  % should or shouldn't we reset ess?
                  % for i=obj.include, obj.factors{i}.ess = obj.factors{i}.essreset(); end
    
                  % iterate over incomplete data
                  for m=1:szinc

                      fprintf('processing case %d of %d\n',m,szinc);

                      % infer parameters
                      engine.enter_evidence(incomplete(m,:));

                      for i=obj.include

                          % only update factors with incomplete data
                          indom = incomplete(m,obj.factors{i}.domain());
                          if any(isnan(indom))

                              % for a BN, the family of node n is the domain of MN factor n
                              marg = engine.marginalize(obj.factors{i}.domain());

                              obj.factors{i} = obj.factors{i}.updateEM(marg,indom);

                          end

                      end

                  end

                  % maximization step
                  for i=obj.include, obj.factors{i} = obj.factors{i}.maximize(); end

                  % recreate inference engine; this could be faster
                  if any(obj.continuous())
                      engine = canonical_jtree_ie(obj);
                  else
                      engine = discrete_jtree_ie(obj);
                  end

                  niter = niter + 1;
              end

          end

       end


   end
end
