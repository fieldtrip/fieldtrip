classdef crossvalidator < validator
%CROSSVALIDATOR crossvalidation class
%
%   OPTIONS:
%   nfolds: the type of crossvalidation used [10] :
%       - 1                : test on training data
%       - 0 <= nfolds < 1  : proportion of training data
%       - nfolds = inf     : leave-one-out crossvalidation
%       - nfolds = N       : N-fold crossvalidation
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
      
    end

    methods
      
      function obj = crossvalidator(varargin)
        
        obj = obj@validator(varargin{:});
        
      end
      
      function obj = validate(obj,data,design)

        % move to cell-array representation
        if ~iscell(data)
          data = {data};
        end
        
        if ~iscell(design)
          design = {design};
        end
        
        % cast to datasets if necessary
            
        for c=1:length(data)
          if ~isa(data{c},'dataset')
            data{c} = dataset(data{c});
          end
        end
        
        for c=1:length(design)
          if ~isa(design{c},'dataset')
            design{c} = dataset(design{c});
          end
        end
        
        % replicate designs
        if length(design)== 1 && length(data) > 1
          design = repmat(design,[1 length(data)]);
        end
        
        nsets = length(data);
        obj.nsets = nsets;
        
        % report some properties
        if obj.verbose          
          fprintf('performing crossvalidation for %d dataset(s)\n',nsets);   
          for d=1:nsets
            sz = data{d}.dims;
            fprintf('dataset %d consists of %d examples and %d',d,data{d}.nsamples,sz(1));
            for dd=2:length(sz)
              fprintf(' x %d',sz(dd));
            end
            fprintf(' features\n');
          end
        end
        
        % create train and test samples
        [obj.trainfolds,obj.testfolds,nfolds] = obj.create_folds(design);
                        
        if ~obj.compact, proc = cell(nfolds,1); end
        
        obj.post = cell(nfolds,nsets);
        obj.design = cell(nfolds,nsets);
        
        for f=1:nfolds % iterate over folds
          
          if obj.verbose, fprintf('validating fold %d of %d\n',f,nfolds); end
          
          traindata = cell(1,nsets);
          traindesign = cell(1,nsets);
          testdata  = cell(1,nsets);          
          testdesign  = cell(1,nsets);
          
          for d=1:nsets % iterate over datasets            
            
            % construct data and design
            
            traindata{d} = data{d}.subsample(obj.trainfolds{f,d});
            traindesign{d} = design{d}.subsample(obj.trainfolds{f,d});
            
            testdata{d} = data{d}.subsample(obj.testfolds{f,d});
            testdesign{d} = design{d}.subsample(obj.testfolds{f,d});
          
          end

          if nsets == 1
            tproc = obj.procedure.train(traindata{1},traindesign{1});
            obj.post{f} = tproc.test(testdata{1});
            obj.design{f} = testdesign{1};
          else
            tproc = obj.procedure.train(traindata,traindesign);
            obj.post(f,:) = tproc.test(testdata);
            obj.design(f,:) = testdesign;           
          end
          
          if ~obj.compact, proc{f} = tproc; end
          
        end
       
        if ~obj.compact
          obj.procedure = proc;
        end      
      end
      
    end
         
   methods(Access=protected)
       
       function [trainfolds,testfolds,nfolds] = create_folds(obj,design)
        % workhorse helper function to create the folds in terms of trial indices
         
          trainfolds = obj.trainfolds;
          testfolds = obj.testfolds;
                              
          nsets = length(design);
          
          if isempty(trainfolds) && isempty(testfolds) % samples are not prespecified so create the testfolds
                    
            nfolds = obj.nfolds;
            if isinf(nfolds)
              nfolds = min(cellfun(@(x)(x.nsamples),design));
            end
            
            if (nfolds > 1 && rem(nfolds,1) == 0) % n-fold crossvalidation
              testfolds  = cell(nfolds,nsets);
            else
              testfolds  = cell(1,nsets);
            end
            
            labeled = cellfun(@(x)(x.labeled),design,'UniformOutput',false);
            unlabeled = cellfun(@(x)(x.unlabeled),design,'UniformOutput',false);
            
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
          
                trainfolds{1,d} = randperm(design{d}.nsamples);
                testfolds{1,d}  = idxs;
                              
              elseif nfolds > 0 && nfolds < 1 % percentage
                
                if obj.verbose, fprintf('validating using %.0f%% of the data for testing\n',100*(1-nfolds)); end
          
                % make sure outcomes are evenly represented whenever possible
                [unq,tmp,idx] = unique(design{d}.X,'rows');
                                
                if max(idx) == nsamples % unique samples
                    testfolds{1,d} = idxs(1:(nsamples - floor(nfolds*nsamples)));
                else
                
                  % take labeled indices
                  idx = idx(idxs);
                  
                  for j=1:length(unq)                   
                    iidx = find(idx == unq(j));
                    testfolds{1,d} = [testfolds{1,d}; idxs(iidx(1:(numel(iidx) - floor(nfolds*numel(iidx)))))];
                  end
                  
                end
                           
              elseif nfolds > 1 && rem(nfolds,1) == 0 % n-fold crossvalidation
                
                if obj.verbose, fprintf('validating using %d-fold cross-validation\n',nfolds); end
          
                % make sure outcomes are evenly represented whenever possible
                [unq,tmp,idx] = unique(design{d}.X,'rows');
                
                if max(idx) == nsamples % unique samples
                 
                  for f=1:nfolds
                    testfolds{f,d} = idxs((floor((f-1)*(length(idxs)/nfolds))+1):floor(f*(length(idxs)/nfolds)));
                  end
                 
                else
               
                  % take labeled indices
                  idx = idx(idxs);               
               
                  f=1;
                  for j=1:length(unq)
                    iidx = find(idx == unq(j));
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
                      
          if isempty(trainfolds)
            
            nfolds = size(testfolds,1);
            
            trainfolds = cell(nfolds,nsets);            
            for f=1:nfolds
              for d=1:nsets
                
                % use the same ordering for multiple datasets if possible
                % by reinitializing the random number generator
                if ~isempty(obj.init)
                  rand('state',obj.init);
                  randn('state',obj.init);
                end
                
                trainfolds{f,d} =  setdiff(1:design{d}.nsamples,testfolds{f,d})';
                trainfolds{f,d} = trainfolds{f,d}(randperm(size(trainfolds{f,d},1)));
                
                if obj.balanced
                  % make sure all classes are represented equally in the
                  % train samples by sampling with replacement
                  
                  if obj.verbose
                    fprintf('balancing training samples by sampling with replacement\n');
                  end
                  
                  [unq,tmp,idx] = unique(design{d}.X,'rows');

                  idx = idx(trainfolds{f,d});
                  
                  maxsmp = max(arrayfun(@(x)(sum(idx == x)),unq));
                  
                  tmp = trainfolds{f,d};
                  for j=1:length(unq)                   
                    
                    iidx = find(idx == unq(j));                    
                    
                    if obj.verbose
                      fprintf('drawing %d additional samples for label %d\n',maxsmp - length(iidx),unq(j));
                    end
                    
                    % nameclash with fieldtrip_private
                    tmp = [tmp; randsample(trainfolds{f,d}(iidx),maxsmp - length(iidx),true)];
                    
                  end
                  
                  trainfolds{f,d} = tmp(randperm(length(tmp)));
                  
                end
                               
              end
            end
            
          end
          
          if isempty(testfolds)
            
            nfolds = size(trainfolds,1);
            
            testfolds = cell(nfolds,nsets);            
            for f=1:nfolds
              for d=1:nsets
                
                % use the same ordering for multiple datasets if possible
                % by reinitializing the random number generator
                if ~isempty(obj.init)
                     rand('state',obj.init);
                     randn('state',obj.init);
                end
                
                testfolds{f,d} = setdiff(1:design{d}.nsamples,trainfolds{f,d})';
                testfolds{f,d} = testfolds{f,d}(randperm(size(testfolds{f,d},1)));
              end
            end
            
          end
                    
       end
   end
end