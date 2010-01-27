classdef lvq < classifier
%LVQ prototype-based supervised classification method (stems from Kohonen's learning vector
%    quantisation)
%
% It is a wrapper for different variants of LVQ: LVQ1, LVQ3 and DSLVQ. They
% are referred to as 'lvq1', 'lvq2' and 'dslvq', respectively. Options include 
% the type of initialisation algorithm along with settings and LVQ codebook 
% learning parameters. 
% 
%
% Options:
%         method      % lvq method (could be 'lvq1','lvq3' or 'dslvq') - default -> 'lvq3'
%         initmethod  % initialisation method (see init_codebook.m) - default -> 'randinit2'
%         initparam   % initialisation setting (see init_codebook.m)         
%         alpha        % learning rate for basic codebook update (LVQ learning rule) - default -> 0.1   
%         eps          % multiplicative coefficient for codebook vector update (dslvq and lvq3 learning rules) - default -> 0.3
%         beta         % learning rate for fweights (only applicable in dslvq_train.m) - default -> 0.1
%         lambda       % multiplicative update of alpha and beta (i.e. alpha(1+1) = lambda*alpha(i)) - default -> 0.95
%         winsize      % size of the window for updating codebook vectors (LVQ learning rule)- default -> 0.3
%         niter        % number of iterative learning updates - default -> 300
% 
%
% REQUIRES: no external dependencies
% 
% EXAMPLE:
%
% myproc = clfproc({ lvq('method','lvq3','initmethod','fcminit','initparam',[10 1]) });
%
%
% SEE ALSO:
%   lvq1_train.m
%   lvq3_train.m
%   dslvq_train.m
%   lvq_test.m
%   init_codebook.m
%
%  REFERENCES:
%      Somervuo P.; Kohonen T.  Self-Organizing Maps and Learning Vector Quantization for Feature Sequences. 
%                               Neural Processing Letters, Volume 10,
%                               Number 2, October 1999, pp.151-159.
%      Pregenzer, M; Pfurtscheller, G. Frequency component selection for an EEG-based brain to computer interface. 
%                                      IEEE Transactions on Rehabilitation
%                                      Engineering. 1999;7(4):413–419.
% 
% Copyright (c) 2009, Pawel Herman
%
% $Log: lvq.m,v $
%

    properties
        
      method = 'lvq3';    % lvq method (could be 'lvq1','lvq3' or 'dslvq')
      initmethod = 'randinit1'; % initialisation method (see init_codebook.m) -
      % it could also be 'none' assuming that cdb_vectors and cdb_labels have already been initialized/specified
      initparam;          % parameter setup for initialisation (see default settings in init_codebook.m)
      nclasses;
      
      % Learning parameters
      alpha  = 0.1;       % learning rate for basic codebook update (LVQ learning rule)
      eps    = 0.2;       % multiplicative coefficient for codebook vector update (dslvq and lvq3 learning rules)
      beta   = 0.1;       % learning rate for fweights (only applicable in dslvq_train.m)
      lambda = 0.95;      % multiplicative update of alpha and beta (i.e. alpha(1+1) = lambda*alpha(i))
      winsize= 0.3;       % size of the window for updating codebook vectors (LVQ learning rule)
      niter  = 100;       % number of iterative learning updates
      
      
      cdb_vectors=[];    % codebook vectors
      cdb_labels=[];     % codebook labels
      fweights=[]        % optional weighting for feature components
      % predominantly for DSLVQ since the weigths are adjusted in a training process
      
    end
    
    methods
      
      function obj = lvq(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function obj = train(obj,data,design)
        % simply stores input data and design
        
        obj.nclasses = design.nunique;
        
        data = data.X;
        design = design.X;
        
        if strcmp(obj.initmethod,'none') && obj.verbose
          disp('option initmethod is set as none which implies omitting initialisation procedure.');
        else
          
          obj = initlvq(obj,data,design);
          
          if obj.verbose, disp(sprintf('method %s was applied at the initialisation stage',obj.initmethod)); end
          
          if ~isempty(obj.cdb_vectors) && ~isempty(obj.cdb_labels)
            switch obj.method
              case 'lvq1', [obj.cdb_vectors obj.cdb_labels] = ...
                  lvq1_train(data,design,obj.cdb_vectors,obj.cdb_labels,obj.niter,obj.alpha,obj.lambda);
              case 'lvq3', [obj.cdb_vectors obj.cdb_labels] = ...
                  lvq3_train(data,design,obj.cdb_vectors,obj.cdb_labels,obj.niter,obj.alpha,obj.eps,obj.lambda,obj.winsize);
              case 'dslvq',[obj.cdb_vectors obj.cdb_labels obj.fweights] = ...
                  dslvq_train(data,design,obj.cdb_vectors,obj.cdb_labels,obj.niter,obj.alpha,obj.eps,obj.beta,obj.lambda,obj.winsize);
            end
            if obj.verbose, disp(sprintf('method %s was applied at the learning stage',obj.method)); end
          else
            error('there are no codebook vectors available for learning');
          end
        end
        
        
      end
      function post = test(obj,data)
        
        if ~isempty(obj.cdb_vectors) && ~isempty(obj.cdb_labels)
          
          output = lvq_test(data.X,obj.cdb_vectors,obj.cdb_labels,obj.fweights);
          
          % post is just a hard class assignment (either 0 or 1) and does not have a probabilistic interpretation
          post = zeros(size(output,1),obj.nclasses);
          for i=1:obj.nclasses
            post(:,i) = (output==i);
          end
        else
          post = [];
          error('codebook vectors must be initialised and/or trained');
        end
        
        post = dataset(post);
        
      end
      
      function obj = initlvq(obj,data,design)
        
        if iscell(data), error('classifier does not take multiple datasets as input'); end
        if strcmp(obj.initmethod,'none') && obj.verbose
          disp('option initmethod is set as none which implies omitting initialisation procedure.');
        else
          [obj.cdb_vectors obj.cdb_labels] = init_codebook(data,design,obj.initmethod,obj.initparam);
        end
        
      end
      
      
    end
end
