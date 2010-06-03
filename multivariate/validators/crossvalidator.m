classdef crossvalidator < validator
%CROSSVALIDATOR crossvalidation class
%
%   OPTIONS:
%   nfolds: the type of crossvalidation used [10] :
%       - 1                : test on training data
%       - 0 <= nfolds < 1  : stratified proportion of training data
%       - nfolds = inf     : leave-one-out crossvalidation
%       - nfolds = N       : N-fold stratified crossvalidation
%
%   alternatively the user can override this by explicitly
%   speciyfing train and/or testfolds
%   folds are of the form { [i11; .. i1k_1] ... [in1 .. ink_n] }
%   or a cell-array of this form if multiple datasets/designs are given
%   e.g. in case of multitask learning
%
%   multiple datasets are given the same train and test trials if they are
%   qualitatively the same by reinitializing the (non-empty) RNG
%
%   data and design can both be multidimensional arrays
%
%   specification of multiple datasets and one design is allowed
%
%   if the design contains missing labels then it will become part of the 
%   training data by default (useful for testing adaptive methods)
%
%   Copyright (c) 2010, Marcel van Gerven

    properties
      
      nfolds = 10; % crossvalidation procedure
      nsets;       % number of datasets
      
      trainfolds; % partitioning of the training examples 
      testfolds;  % partitioning of the test examples
      
      % if true then the model averaged over folds will be saved here.
      % useful in conjunction with compact=false (classifiers discarded;
      % average model saved)
      model = 0;  
      desc;
      
    end

    methods
      
      function obj = crossvalidator(varargin)
        
        obj = obj@validator(varargin{:});
        
      end
      
      function obj = validate(obj,data,design)
      
        assert(~isempty(obj.procedure));
      
        if ~isa(obj.procedure,'mva')
          % try to create procedure if input is a cell array or a predictor
          obj.procedure = mva(obj.procedure);
        end
        
        if strcmp(obj.verbose,'all')
          for j=1:length(obj.procedure.mvmethods)
            obj.procedure.mvmethods{j}.verbose = true;
          end
          obj.verbose = true;
        end
      
        % move to cell-array representation
        if ~iscell(data)
          data = {data};
        end
        
        if ~iscell(design)
          design = {design};
        end
        
        % replicate designs (useful for transfer learning with the same designs)
        if length(design)== 1 && length(data) > 1
          design = repmat(design,[1 length(data)]);
        end
        
        nsets = length(data);
        obj.nsets = nsets;
        
        % report some properties
        if obj.verbose          
          fprintf('performing crossvalidation for %d dataset(s)\n',nsets);   
          for d=1:nsets
            szin = size(data{d});
            szout = size(design{d});
            fprintf('dataset %d consists of %d examples, %d',d,szin(1),szin(2));
            for dd=3:length(szin)
              fprintf(' x %d',szin(dd));
            end
            if length(szin)>2
              fprintf(' = %d',prod(szin(2:end)));
            end
            fprintf(' input features and %d',szout(2));
            for dd=3:length(szout)
              fprintf(' x %d',szout(dd));
            end
            fprintf(' output features\n');
          end
        end
        
        % create train and test samples
        [obj.trainfolds,obj.testfolds,nfolds] = obj.create_folds(design);
        
        assert(iscell(obj.trainfolds));
        assert(iscell(obj.testfolds));
        assert(size(obj.trainfolds,2) == nsets);
        
        if obj.verbose
          fprintf('using %d folds for %d datasets\n',size(obj.trainfolds,1),size(obj.trainfolds,2));
        end
        
        obj.post = cell(nfolds,nsets);
        obj.design = cell(nfolds,nsets);
        
        % this can happen if the validator has already been validated
        % e.g., in case of the optimizer
        if ~iscell(obj.procedure)
          tproc = obj.procedure;
          obj.procedure = cell(nfolds,1);
          for c=1:numel(obj.procedure)
            obj.procedure{c} = tproc;
          end
        end
        
        for f=1:nfolds % iterate over folds
          
          if obj.verbose
            for j=1:size(obj.trainfolds,2)
              fprintf('dataset %d: validating fold %d of %d using %d training samples and %d test samples\n',...
                j,f,nfolds,numel(obj.trainfolds{f,j}),numel(obj.testfolds{f,j}));
            end
          end
          
          traindata = cell(1,nsets);
          traindesign = cell(1,nsets);
          testdata  = cell(1,nsets);          
          testdesign  = cell(1,nsets);
          
          for d=1:nsets % iterate over datasets            
            
            % construct data and design; reshape to original dimensions
            
            sz = size(data{d}); 
            
            sz(1) = numel(obj.trainfolds{f,d});
            traindata{d} = reshape(data{d}(obj.trainfolds{f,d},:),sz);
            if isempty(traindata{d}), traindata{d} = []; end 
              
            sz(1) = numel(obj.testfolds{f,d});
            testdata{d} = reshape(data{d}(obj.testfolds{f,d},:),sz);
            if isempty(testdata{d}), testdata{d} = []; end 
            
            sz = size(design{d});
            
            sz(1) = numel(obj.trainfolds{f,d});
            traindesign{d} = reshape(design{d}(obj.trainfolds{f,d},:),sz);
            if isempty(traindesign{d}), traindesign{d} = []; end
          
            sz(1) = numel(obj.testfolds{f,d});
            testdesign{d} = reshape(design{d}(obj.testfolds{f,d},:),sz);
            if isempty(testdesign{d}), testdesign{d} = []; end 
          
          end

          if nsets == 1 && ~obj.getpredictor().istransfer()
            tproc = obj.procedure{f}.train(traindata{1},traindesign{1});
            obj.post{f} = tproc.test(testdata{1});
            obj.design{f} = testdesign{1};
          else
            
            tproc = obj.procedure{f}.train(traindata,traindesign);
            
%             % transfer learner is free to return a cell array or not in
%             % case of one input dataset
%             res = tproc.test(testdata);
%             if iscell(res), res = res{1}; end
%             obj.post{f} = res;
            
            obj.post(f,:) = tproc.test(testdata);
            
            obj.design(f,:) = testdesign;  
          
          end
          
          if iscell(obj.model) || obj.model
            if f==1
              [obj.model,obj.desc] = tproc.getmodel();
            else
              try % fails for incompatible models in each fold
                m = tproc.getmodel();
                for mm=1:numel(obj.model)
                  obj.model{mm} = obj.model{mm} + m{mm};
                end
                clear m;
              catch
                obj.model = {};
              end
            end
          end
          
          if ~obj.compact
            obj.procedure{f} = tproc; 
          else            
            
            clear tproc;
            clear traindata;
            clear testdata;
            clear traindesign;
            clear testdesign;
          end
          
        end
       
        if iscell(obj.model) || obj.model
          for mm=1:numel(obj.model)
            obj.model{mm} = obj.model{mm}./nfolds;
          end
        end
        
      end
      
    end
         
   methods
       
       function [trainfolds,testfolds,nfolds] = create_folds(obj,design)
        % workhorse helper function to create the folds in terms of trial indices
         
          trainfolds = obj.trainfolds;
          testfolds = obj.testfolds;
                              
          nsets = length(design);
          
          if ~isempty(obj.trainfolds) && ~isempty(obj.testfolds) % fully specified already
            
            nfolds = size(trainfolds,1); 
            
            if size(trainfolds,2) == 1 && nsets > 1
              trainfolds = repmat(trainfolds,[1 nsets]);
            end
            if size(testfolds,2) == 1 && nsets > 1
              testfolds = repmat(testfolds,[1 nsets]);
            end
            
            return;
            
          elseif isempty(trainfolds) && isempty(testfolds) % samples are not prespecified so create the testfolds
                    
            nfolds = obj.nfolds;
            if isinf(nfolds)
              nfolds = min(cellfun(@(x)(size(x,1)),design));
            end
            
            if (nfolds > 1 && rem(nfolds,1) == 0) % n-fold crossvalidation
              testfolds  = cell(nfolds,nsets);
            else
              testfolds  = cell(1,nsets);
            end
            
            labeled = cellfun(@(x)(mvmethod.labeled(x)),design,'UniformOutput',false);
            unlabeled = cellfun(@(x)(mvmethod.unlabeled(x)),design,'UniformOutput',false);
            
            for d=1:nsets
              
              if obj.verbose && ~isempty(unlabeled{d})
                fprintf('putting %d unlabeled samples into training set\n',length(unlabeled{d}))
              end
              
              % use the same ordering for multiple datasets if possible
              % by reinitializing the random number generator
              if ~isempty(obj.init)                
                 rand('state',obj.init);
                 randn('state',obj.init);
               end
              
              % randomize labeled trials
              nsamples = size(labeled{d},1);
              idxs = labeled{d}(randperm(nsamples));              

              if nfolds == 1 % trainset is testset
                
                if obj.verbose, fprintf('validating using training data\n'); end
          
                trainfolds{1,d} = randperm(nsamples);
                testfolds{1,d}  = idxs;
                              
              elseif nfolds > 0 && nfolds < 1 % percentage
                
                if obj.verbose, fprintf('validating using %.0f%% of the data for testing\n',100*(1-nfolds)); end
          
                % make sure outcomes are evenly represented whenever possible
                [unq,tmp,idx] = unique(design{d},'rows');
                         
                mx = max(idx);
                if mx == nsamples % unique labeled samples
                    testfolds{1,d} = idxs(1:round((1-nfolds)*nsamples));
                else
                
                  % take labeled indices
                  idx = idx(idxs);
                  
                  for j=1:mx                 
                    iidx = find(idx == j);
                    testfolds{1,d} = [testfolds{1,d}; idxs(iidx(1:floor((1-nfolds)*numel(iidx))))];
                  end
                  
                  % add random instances to complete the number
                  n = floor((1-nfolds) * numel(idxs)) - numel(testfolds{1,d});
                  if n>0
                    x = setdiff(idxs,testfolds{1,d});
                    x = x(randperm(numel(x)));
                    testfolds{1,d} = [testfolds{1,d}; x(1:n)];
                  end
                      
                end
                           
              elseif nfolds > 1 && rem(nfolds,1) == 0 % n-fold crossvalidation
                
                if obj.verbose, fprintf('validating using %d-fold cross-validation\n',nfolds); end
          
                % make sure outcomes are evenly represented whenever possible
                [unq,tmp,idx] = unique(design{d}(1:size(design{d}),:),'rows');
                
                mx = max(idx);
                if mx == nsamples % unique samples
                 
                  for f=1:nfolds
                    testfolds{f,d} = idxs((floor((f-1)*(length(idxs)/nfolds))+1):floor(f*(length(idxs)/nfolds)));
                  end
                 
                else
               
                  % take labeled indices
                  idx = idx(idxs);               
               
                  f=1;
                  for j=1:mx
                    iidx = find(idx == j);
                    for k=1:length(iidx)
                      testfolds{f,d} = [testfolds{f,d}; idxs(iidx(k))];
                      f = f+1; if f > nfolds, f=1; end
                    end
                  end
                  
                end
                
              end
              
            end
          
            if nfolds > 0 && nfolds < 1 % percentage
                nfolds = 1;
            end
            
          end
                      
          if ~isempty(trainfolds)
            nfolds = size(trainfolds,1);
            
            % replicate sets for all datasets 
            %(assuming all datasets have the same design matrix)
            if size(trainfolds,2) == 1 && nsets > 1
              trainfolds = repmat(trainfolds,[1 nsets]);
            end
            
          elseif ~isempty(testfolds)
              
            nfolds = size(testfolds,1);
          
            % replicate sets for all datasets
            %(assuming all datasets have the same design matrix)
            if size(testfolds,2) == 1 && nsets > 1 
              testfolds = repmat(testfolds,[1 nsets]);
            end
          else
            nfolds = nan;
          end
          
          if isempty(trainfolds)
            
            trainfolds = cell(nfolds,nsets);            
            for f=1:nfolds
              for d=1:nsets
                
                % use the same ordering for multiple datasets if possible
                % by reinitializing the random number generator
                if ~isempty(obj.init)
                  rand('state',obj.init);
                  randn('state',obj.init);
                end
                
                trainfolds{f,d} = setdiff(1:size(design{d},1),testfolds{f,d})';
                trainfolds{f,d} = trainfolds{f,d}(randperm(size(trainfolds{f,d},1)));
                               
              end
            end
            
          end
          
          if obj.balanced
            % make sure all classes are represented equally in the
            % train samples by sampling with replacement and in the
            % test samples by using the minimum number of trials in a
            % condition. This ensures that training makes most use of
            % the data while testing has no inherent bias in the
            % data. We lose power by removing superfluous trials.
                        
            if obj.verbose
              fprintf('balancing training samples by sampling with replacement\n');
            end
            
            for f=1:nfolds
              for d=1:nsets
            
                [unq,tmp,idx] = unique(design{d}(1:size(design{d}),:),'rows');
                
                idx = idx(trainfolds{f,d});
                
                maxsmp = max(arrayfun(@(x)(sum(idx == x)),unq));
                
                tmp = trainfolds{f,d};
                for j=1:length(unq)
                  
                  iidx = find(idx == unq(j));
                  
                  if obj.verbose
                    fprintf('drawing %d additional samples for label %d\n',maxsmp - length(iidx),unq(j));
                  end
                  
                  % nameclash with fieldtrip_private!
                  tmp = [tmp; randsample(trainfolds{f,d}(iidx),maxsmp - length(iidx),true)];
                  
                end
                
                trainfolds{f,d} = tmp(randperm(length(tmp)));
            
              end
            end
            
          end
          
          if isempty(testfolds)
            
            testfolds = cell(nfolds,nsets);            
            for f=1:nfolds
              for d=1:nsets
                
                % use the same ordering for multiple datasets if possible
                % by reinitializing the random number generator
                if ~isempty(obj.init)
                  rand('state',obj.init);
                  randn('state',obj.init);
                end
                
                testfolds{f,d} = setdiff(1:size(design{d},1),trainfolds{f,d})';
                testfolds{f,d} = testfolds{f,d}(randperm(size(testfolds{f,d},1)));
                
              end
            end
            
          end
          
          if obj.balanced
            % make sure all classes are represented equally in the
            % train samples by sampling with replacement and in the
            % test samples by using the minimum number of trials in a
            % condition. This ensures that training makes most use of
            % the data while testing has no inherent bias in the
            % data. We lose power by removing superfluous trials.
            
            if obj.verbose
              fprintf('balancing test samples by removing superfluous trials\n');
            end
            
            for f=1:nfolds
              for d=1:nsets
                
                [unq,tmp,idx] = unique(design{d}(1:size(design{d}),:),'rows');
                
                idx = idx(testfolds{f,d});
                
                minsmp = min(arrayfun(@(x)(sum(idx == x)),unq));
                
                tmp = [];
                for j=1:length(unq)
                  
                  iidx = find(idx == unq(j));
                  
                  if obj.verbose
                    fprintf('removing %d superfluous samples for label %d\n',length(iidx) - minsmp,unq(j));
                  end
                  
                  tmp = [tmp; testfolds{f,d}(iidx(1:minsmp))];
                  
                end
                
                testfolds{f,d} = tmp(randperm(length(tmp)));
          
              end
            end
            
          end
                 
          
          
       end
   end
end