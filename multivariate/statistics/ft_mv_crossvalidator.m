classdef ft_mv_crossvalidator
%FT_MV_CROSSVALIDATOR stratified crossvalidation class
%
%   OPTIONS:
%   nfolds: the type of crossvalidation used [10] :
%       - 1                : test on training data
%       - 0 <= nfolds < 1  : stratified proportion of training data
%       - nfolds = inf     : leave-one-out crossvalidation
%       - nfolds = N       : stratified N-fold crossvalidation
%       - 'online'         : simulated online crossvalidation (sequential)
%
%   alternatively the user can override this by explicitly
%   speciyfing train and/or testfolds
%   folds are of the form { [i11; .. i1k_1] ... [in1 .. ink_n] }
%
%   in case of multiple datasets, folds{d} is a cell-array of folds for the
%   d-th dataset
%
%   multiple datasets are given the same train and test trials if they are
%   qualitatively the same by reinitializing the (non-empty) RNG
%
%   specification of multiple datasets and one design is allowed
%
%   if the design contains missing labels then it will become part of the 
%   training data by default (useful for testing adaptive methods)
%
%   Copyright (c) 2010, Marcel van Gerven

    properties
      
      init        = 1;      % initialize RNG
      compact     = false;  % only retain necessary results
      verbose     = false;  % verbose output
      balance     = false;  % balance data in each fold
      nfolds      = 10;     % crossvalidation procedure
      parallel    = false;  % run in parallel mode 
      
      mva                   % multivariate analysis
      design                % real outputs of size datasets x nr of folds
      post                  % predicted outputs of size datasets x nr of folds
            
      trainfolds;           % partitioning of the training examples 
      testfolds;            % partitioning of the test examples
            
      metric                % performance metric
      sigtest               % significance test 
      
      % the model averaged over folds will be saved here.
      model;  
      description;

    end

    methods
      
      function obj = ft_mv_crossvalidator(varargin)
        
        % parse options
        for i=1:2:length(varargin)
          if ismember(varargin{i},fieldnames(obj))
            obj.(varargin{i}) = varargin{i+1};
          else
            error('unrecognized fieldname %s',varargin{i});
          end
        end
        
      end
      
      function obj = train(obj,X,Y)
      
        % initialize random number generator
        if ~isempty(obj.init)
          if obj.verbose
            fprintf('initializing random number generator with seed %d\n',obj.init)
          end
          %RandStream.setDefaultStream(RandStream('mt19937ar','seed',obj.init));
          rand('state',obj.init);
          randn('state',obj.init);
        end
        
        if isempty(obj.mva), error('mva not specified'); end
      
        if ~isa(obj.mva,'ft_mv_analysis') && ~(iscell(obj.mva) && isa(obj.mva{1},'ft_mv_analysis'))
          % try to create mva if input is a cell array or a predictor
          obj.mva = ft_mv_analysis(obj.mva);
        end
               
        % move to cell-array representation to allow uniform handling
        if ~iscell(X), X = {X}; end
        if ~iscell(Y), Y = {Y}; end

        % create testfolds in sequential order
        if strcmp(obj.nfolds,'online')
          obj.nfolds = [];
          obj.testfolds = 1:size(X{1},1);
        end    
        
        % determine maximum number of folds for leave-one-out
        if isinf(obj.nfolds), obj.nfolds = min(cellfun(@(x)(size(x,1)),Y)); end
    
        % replicate designs in case of one design for multiple datasets
        if length(Y)== 1 && length(X) > 1
          Y = repmat(Y,[1 length(X)]);
        end

        % number of datasets
        nsets = length(X);

        % initialize train and testfolds (if not manually specified)
        if isempty(obj.trainfolds), obj.trainfolds = cell(nsets,1); end
        if isempty(obj.testfolds),  obj.testfolds = cell(nsets,1); end
    
        % change to cell representation in case the folds for just one
        % dataset are specified
        if ~isempty(obj.trainfolds{1}) && ~iscell(obj.trainfolds{1})
          obj.trainfolds = {obj.trainfolds};
        end
        if ~isempty(obj.testfolds{1}) && ~iscell(obj.testfolds{1})
          obj.testfolds = {obj.testfolds};
        end
        
        % report some properties
        if obj.verbose, obj.report(X,Y); end
        
        % create train and test folds
        [obj.trainfolds,obj.testfolds] = obj.create_folds(Y);

        % number of datasets
        nsets = length(obj.trainfolds);
        
        % number of folds
        nf = length(obj.trainfolds{1});

        % initialize mva for each fold
        if ~iscell(obj.mva)
          tproc = obj.mva;
          obj.mva = cell(nf,1);
          for c=1:numel(obj.mva)
            obj.mva{c} = tproc;
          end
        end

        if obj.verbose
          fprintf('using %d folds for %d datasets\n',nf,nsets);
        end
                 
        obj.design = cell(nsets,nf);
        obj.post = cell(nsets,nf);
        if obj.parallel
       
          if obj.verbose
            fprintf('running %s in parallel mode',class(obj));
          end
          train_parallel();
          
        else
          
          train_sequential();
          
        end
   
        for mm=1:numel(obj.model)
          if iscell(obj.model{mm})
            for cc=1:length(obj.model{mm}) % occurs for nested methods
              obj.model{mm}{cc} = obj.model{mm}{cc}./nf;
            end
          else
            obj.model{mm} = obj.model{mm}./nf;
          end
        end
       
        % local function
        function train_sequential()
          
         for f=1:nf % iterate over folds
            
            trainX = cell(nsets,1);
            trainY = cell(nsets,1);
            testX  = cell(nsets,1);
            testY  = cell(nsets,1);
            
            for d=1:nsets % iterate over datasets
              
              if obj.verbose
                fprintf('dataset %d: validating fold %d of %d using %d training samples and %d test samples\n',...
                  d,f,nf,numel(obj.trainfolds{d}{f}),numel(obj.testfolds{d}{f}));
              end
              
              % construct X and Y
              trainX{d} = X{d}(obj.trainfolds{d}{f},:);
              if isempty(trainX{d}), trainX{d} = []; end
              
              testX{d} = X{d}(obj.testfolds{d}{f},:);
              if isempty(testX{d}), testX{d} = []; end
              
              trainY{d} = Y{d}(obj.trainfolds{d}{f},:);
              if isempty(trainY{d}), trainY{d} = []; end
              
              testY{d} = Y{d}(obj.testfolds{d}{f},:);
              if isempty(testY{d}), testY{d} = []; end
              
            end
            
            if nsets == 1
              
              tproc = obj.mva{f}.train(trainX{1},trainY{1});
              obj.post{f} = tproc.test(testX{1});
              obj.design{f} = testY{1};
              
            else
              
              tproc = obj.mva{f}.train(trainX,trainY);
              res = tproc.test(testX);
              if ~iscell(res), res = {res}; end
              obj.post(:,f) = res;
              obj.design(:,f) = testY;
              
            end
            
            if f==1
              [obj.model,obj.description] = tproc.model();
            else
              m = tproc.model();
              for mm=1:numel(obj.model)
                try
                  obj.model{mm} = obj.model{mm} + m{mm};
                catch
                  warning('individual model dimensions do not match');
                end
              end
              clear m;
            end
            
            % attempt to reshape to original dimensions
            sz = size(X); sz = sz(2:end); if numel(sz)==1, sz = [1 sz]; end
            for c=1:length(obj.model)
              if numel(obj.model{c}) == sz
                obj.model{c} = reshape(obj.model{c},sz);
              end
            end
            
            if ~obj.compact, obj.mva{f} = tproc; end            
            
            clear tproc;
            clear trainX;
            clear testX;
            clear trainY;
            clear testY;
            
          end
          
        end
        
        % local function
        function train_parallel()
           
          % we need to construct all trainsets beforehand
          trainX = cell(1,nf);
          trainY = cell(1,nf);
          for f=1:nf 
            
            if nsets == 1
              
              trainX{f} = X{1}(obj.trainfolds{1}{f},:);
              if isempty(trainX{f}), trainX{f} = []; end
              trainY{f} = Y{1}(obj.trainfolds{1}{f},:);
              if isempty(trainY{f}), trainY{f} = []; end
              
            else
              
              trainX{f} = cell(1,nsets);
              trainY{f} = cell(1,nsets);
              for d=1:nsets % iterate over datasets
                
                trainX{f}{d} = X{f}{d}(obj.trainfolds{d}{f},:);
                if isempty(trainX{f}{d}), trainX{f}{d} = []; end
                trainY{f}{d} = Y{f}{d}(obj.trainfolds{d}{f},:);
                if isempty(trainY{f}{d}), trainY{f}{d} = []; end
              
              end
              
            end
          end
          tproc = peercellfun(@run_parallel,obj.mva',repmat({'train'},[1 nf]),trainX,trainY,'UniformOutput',false);
          clear trainX
          clear trainY
            
          % we need to construct all testsets beforehand
          testX  = cell(1,nf);
          testY  = cell(1,nf);
          for f=1:nf
            
            if nsets == 1
              
              testX{f} = X{1}(obj.testfolds{1}{f},:);
              if isempty(testX{f}), testX{f} = []; end
              testY{f} = Y{1}(obj.testfolds{1}{f},:);
              if isempty(testY{f}), testY{f} = []; end
              
            else
              
              testX{f}  = cell(1,nsets);
              testY{f}  = cell(1,nsets);
              for d=1:nsets % iterate over datasets
                
                testX{f}{d} = X{f}{d}(obj.testfolds{d}{f},:);
                if isempty(testX{f}{d}), testX{f}{d} = []; end
                testY{f}{d} = Y{f}{d}(obj.testfolds{d}{f},:);
                if isempty(testY{f}{d}), testY{f}{d} = []; end
                
              end
              
            end
          end
          res = peercellfun(@run_parallel,tproc,repmat({'test'},[nf 1]),testX,'UniformOutput',false);
          if ~iscell(res), res = {res}; end
          obj.post = res;
          obj.design = testY;
          clear testX
          clear testY
          
          if ~obj.compact, obj.mva = tproc; end
   
          for f=1:nf
            if f==1
              [obj.model,obj.description] = tproc{f}.model();
            else
              m = tproc{f}.model();
              for mm=1:numel(obj.model)
                obj.model{mm} = obj.model{mm} + m{mm};
              end
              clear m;
            end
          end
          clear tproc
          
        end
        
      end
   
     function [train,test] = create_folds(obj,Y)
       
       train = obj.trainfolds;
       test = obj.testfolds;
       
       % create train and test folds per dataset
        for c=1:length(Y)
         
          % create testfolds based on set parameters
          if isempty(train{c}) && isempty(test{c})
            test{c} = obj.create_test_folds(Y{c});
          end
            
          if isempty(train{c})
           train{c} = obj.create_complement(test{c},size(Y{c},1));
          end
            
          if isempty(test{c})
            test{c} = obj.create_complement(train{c},size(Y{c},1));
          end
          
          if obj.balance
            train{c} = obj.balance_folds(train{c},Y{c},true);
            test{c}  = obj.balance_folds(test{c},Y{c},false);
          end
          
        end
        
     end
       
     function res = create_complement(obj,folds,nsamples)
       % create complement indices from specified indices
       
       nf = length(folds);
       
       res = cell(nf,1);
       for f=1:nf
         
         % use the same ordering for multiple datasets if possible
         % by reinitializing the random number generator
         if ~isempty(obj.init)
           rand('state',obj.init);
           randn('state',obj.init);
         end
         
         res{f} = setdiff(1:nsamples,folds{f})';
         res{f} = res{f}(randperm(numel(res{f})));
         
       end
              
     end
     
     function res = balance_folds(obj,folds,Y,resample)
       % balance folds
       
        nf = length(folds);
      
        if resample
         % make sure all classes are represented equally in the
         % train samples by sampling with replacement and in the
         % test samples by using the minimum number of trials in a
         % condition. This ensures that training makes most use of
         % the data while testing has no inherent bias in the
         % data. We lose power by removing superfluous trials.
         
         if obj.verbose
           fprintf('balancing training samples by sampling with replacement\n');
         end
         
         for f=1:nf
              
           % get unique rows of Y
           [unq,tmp,idx] = unique(Y,'rows');

           idx = idx(folds{f});
           
           maxsmp = max(arrayfun(@(x)(sum(idx == x)),unq));
           
           tmp = folds{f};
           for j=1:length(unq)
             
             iidx = find(idx == unq(j));
             
             if obj.verbose && maxsmp - length(iidx) ~= 0
               fprintf('drawing %d additional samples for label %d\n',maxsmp - length(iidx),unq(j));
             end
             
             % nameclash with fieldtrip_private!
             tmp = [tmp; randsample(folds{f}(iidx),maxsmp - length(iidx),true)];
             
           end
           
           folds{f} = tmp(randperm(length(tmp)));
           
         end
         
       else
         % make sure all classes are represented equally in the
         % train samples by sampling with replacement and in the
         % test samples by using the minimum number of trials in a
         % condition. This ensures that training makes most use of
         % the data while testing has no inherent bias in the
         % data. We lose power by removing superfluous trials.
         
         if obj.verbose
           fprintf('balancing test samples by removing superfluous trials\n');
         end
         
         for f=1:nf
           
           % get unique rows of Y
           [unq,tmp,idx] = unique(Y,'rows');

           idx = idx(folds{f});
           
           minsmp = min(arrayfun(@(x)(sum(idx == x)),unq));
           
           if minsmp > 0
             
             tmp = [];
             for j=1:length(unq)
               
               iidx = find(idx == unq(j));
               
               if obj.verbose && length(iidx) - minsmp ~=0
                 fprintf('removing %d superfluous samples for label %d\n',length(iidx) - minsmp,unq(j));
               end
               
               tmp = [tmp; folds{f}(iidx(1:minsmp))];
               
             end
             
             folds{f} = tmp(randperm(length(tmp)));
             
           end
         end
        end
       
        res = folds;
        
     end
       
     function res = create_test_folds(obj,Y)
       % create test folds in terms of trial indices
       
       % use the same ordering for multiple datasets by reinitializing the random number generator
       if ~isempty(obj.init)
         rand('state',obj.init);
         randn('state',obj.init);
       end
       
       % determine number of folds
       nf = obj.nfolds;
       if nf<=1, nf = 1; end
       res  = cell(nf,1);
       
       % determine indices of labeled (non-nan) datapoints
       % nan datapoints should end up in the training set
       labeled = find(any(~isnan(Y(1:size(Y,1),:)),2));
       if obj.verbose && numel(labeled)<size(Y,1)
         fprintf('adding unlabeled samples to training set\n');
       end
       
       % randomize labeled trials
       nsamples = size(labeled,1);
       idxs = labeled(randperm(nsamples));
       
       if obj.nfolds > 0 && obj.nfolds < 1 % percentage
         
         if obj.verbose, fprintf('validating using %.0f%% of the data for testing\n',100*(1-obj.nfolds)); end
         
         % make sure outcomes are evenly represented whenever possible
         [unq,tmp,idx] = unique(Y,'rows');
                  
         mx = max(idx);
         if mx == nsamples % unique labeled samples
           res{1} = idxs(1:round((1-obj.nfolds)*nsamples));
         else
           
           % take labeled indices
           idx = idx(idxs);
           for j=1:mx
             iidx = find(idx == j);
             res{1} = [res{1}; idxs(iidx(1:floor((1-obj.nfolds)*numel(iidx))))];
           end
           
           % add random instances to complete the number
           n = floor((1-obj.nfolds) * numel(idxs)) - numel(res{1});
           if n>0
             x = setdiff(idxs,res{1});
             x = x(randperm(numel(x)));
             res{1} = [res{1}; x(1:n)];
           end
           
         end
         
       else % n-fold crossvalidation
         
         if obj.verbose, fprintf('validating using %d-fold cross-validation\n',nf); end
         
         % make sure outcomes are evenly represented whenever possible
         [unq,tmp,idx] = unique(Y(1:size(Y),:),'rows');
         
         mx = max(idx);
         if mx == nsamples % unique samples
           
           for f=1:nf
             res{f} = idxs((floor((f-1)*(length(idxs)/nf))+1):floor(f*(length(idxs)/nf)));
           end
           
         else
           
           % take labeled indices
           idx = idx(idxs);
           
           f=1;
           for j=1:mx
             iidx = find(idx == j);
             for k=1:length(iidx)
               res{f} = [res{f}; idxs(iidx(k))];
               f = f+1; if f > nf, f=1; end
             end
           end
           
         end
         
       end
       
     end
     
     function res = performance(obj,M)

       if exist('M','var')
         obj.metric = M;
       end

       % concatenate all folds
       tpost = cell(1,size(obj.post,1));
       tdesign = cell(1,size(obj.design,1));
       for c=1:length(tpost)
         tpost{c} = cat(1,obj.post{c,:});
         tdesign{c} = cat(1,obj.design{c,:});
       end
       
       % compute metric for each dataset
       res = ft_mv_performance(tdesign,tpost,obj.metric);
       
     end
     
     function res = significance(obj,T)
       % return p-value using some significance test
       
       if exist('T','var')
         obj.sigtest = T;
       end
       
        % concatenate all folds
       tpost = cell(1,size(obj.post,1));
       tdesign = cell(1,size(obj.design,1));
       for c=1:length(tpost)
         tpost{c} = cat(1,obj.post{c,:});
         tdesign{c} = cat(1,obj.design{c,:});
       end
       
       % compute significance test for each dataset
       res = ft_mv_significance(tdesign,tpost,obj.sigtest);
       
     end
     
   end
   
   methods(Static=true)
     
     
    function report(X,Y)
      
      fprintf('validating %d dataset(s)\n',length(X));
      for d=1:length(X)
        szin = size(X{d});
        fprintf('input %d consists of %d examples and %d',d,szin(1),szin(2));
        for dd=3:length(szin)
          fprintf(' x %d',szin(dd));
        end
        if length(szin)>2
          fprintf(' = %d',prod(szin(2:end)));
        end
        fprintf(' features\n');
      end
      
      for d=1:length(Y)
        szin = size(Y{d});
        fprintf('output %d consists of %d examples and %d',d,szin(1),szin(2));
        for dd=3:length(szin)
          fprintf(' x %d',szin(dd));
        end
        if length(szin)>2
          fprintf(' = %d',prod(szin(2:end)));
        end
        fprintf(' features\n');
      end
      
    end
    
   end
end