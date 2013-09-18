classdef crossvalidator
% CROSSVALIDATOR crossvalidation class.
%
%   DESCRIPTION
%   The crossvalidator class performs crossvalidation. It uses a
%   particular crossvalidation on a multivariate analysis specified in the
%   property mva.
%
%   One may use 'type' equal to:
%   - 'nfold' : n-fold cross-validation specified by 'folds' (default: 5)
%   - 'split' : data is split into train and test data; amount of training
%               data is given by 'proportion' (default: 0.75)
%   - 'loo'   : leave-one-out cross-validation; takes 1 trial per fold
%   - 'bloo'  : balanced loo cv; takes 1 trial from each class per fold.
%               Assumes that each class has the same nr of trials.
%  
%   Alternatively one can specify the cell-arrays 'trainfolds' and
%   'testfolds' to contain for each fold the trial indices that need to be
%   used. If one of these is specified then the other one is automatically
%   filled using the complement of the trials. If these cell-arrays are
%   non-empty then the settings above are ignored.
%
%   In order to balance the occurrence of different classes one may set
%   'resample' equal to true (default: false). Resample will upsample less
%   occurring classes during training and downsample often occurring
%   classes during testing.
%
%   In case of memory problems, use 'compact' equal to true (default:
%   false). This will not store the trained multivariate methods but just
%   the classification outcomes and the parameters that are returned by the
%   model function of the multivariate methods.
%
%   EXAMPLE
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]';
%   m = dml.crossvalidator('mva',dml.svm)
%   m = m.train(X,Y);
%   m.statistic('accuracy')
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)
  
    properties
     
      mva % multivariate analysis
      
      result % crossvalidation result
     
      design % the real outputs
      
      model % estimated model per fold
      
      type = 'nfold'; % 'type' : crossvalidation type 'nfold', 'split', 'loo'
      
      folds = 5; % 'folds' : number of folds when using nfold
      
      proportion = 0.75; % 'proportion' : proportion of used training data when using split
      
      resample = false; % 'resample' : whether or not to use resampling to balance classes
            
      % cell array indicating used trials per training fold
      % for multiple datasets this should be come a cell array of cell
      % arrays, e.g. trainfolds = { {f1 ... fn} {f1 ... fn} }
      trainfolds 
      
      % cell array indicating used trials per test fold
      % for multiple datasets this should be come a cell array of cell
      % arrays, e.g. testfolds = { {f1 ... fn} {f1 ... fn} }
      testfolds 
      
      stat = ''; % test statistic used to quantify performance of this cross-validator
      
      verbose = false; % generate output when true

      compact = false; % drop mva details when true 
            
    end

    methods
      
      function obj = crossvalidator(varargin)
        
        % parse options
        for i=1:2:length(varargin)
          if ismember(varargin{i},fieldnames(obj))
            obj.(varargin{i}) = varargin{i+1};
          else
            error('unrecognized fieldname %s',varargin{i});
          end
        end
        
        if ~isa(obj.mva,'dml.analysis'), obj.mva = dml.analysis(obj.mva); end
                
      end
      
      function obj = train(obj,X,Y)
         
        % complete the folds when only train or test is specified
        if isempty(obj.trainfolds) && isempty(obj.testfolds)
          [obj.trainfolds,obj.testfolds] = obj.create_folds(Y);
        elseif isempty(obj.trainfolds)
          obj.trainfolds = obj.complement(Y,obj.testfolds);
        else
          obj.testfolds = obj.complement(Y,obj.trainfolds);
        end
        
        if iscell(Y)
          ndata = length(Y);
        else
          ndata = 1;
        end
        
        if ndata == 1
          nfolds = length(obj.trainfolds);
        else
          nfolds = length(obj.trainfolds{1});
        end
        
        obj.result = cell(nfolds,1);
        obj.design = cell(nfolds,1);
        obj.model = cell(nfolds,1);
        
        if ~iscell(obj.mva)
          obj.mva = repmat({obj.mva},[1 nfolds]);
        end
        
        for f=1:nfolds % iterate over folds
            
          if obj.verbose
            if ndata == 1
              fprintf('validating fold %d of %d for %d datasets\n',f,nfolds,ndata);
            else
              fprintf('validating fold %d of %d for %d datasets\n',f,nfolds,ndata);
            end   
          end

          % construct X and Y for each fold
          if ndata == 1
            trainX = X(obj.trainfolds{f},:);
            testX = X(obj.testfolds{f},:);
            trainY = Y(obj.trainfolds{f},:);
            testY = Y(obj.testfolds{f},:);
          else
            trainX = cell([length(X) 1]);
            testX = cell([length(X) 1]);
            for c=1:length(X)
              trainX{c} = X{c}(obj.trainfolds{c}{f},:);
              testX{c} = X{c}(obj.testfolds{c}{f},:);
            end
            trainY = cell([length(Y) 1]);
            testY = cell([length(Y) 1]);
            for c=1:length(Y)
              trainY{c} = Y{c}(obj.trainfolds{c}{f},:);
              testY{c} = Y{c}(obj.testfolds{c}{f},:);
            end
          end
                    
          tproc = obj.mva{f}.train(trainX,trainY);
          if ~isempty(testY)
            obj.result{f} = tproc.test(testX);
            obj.design{f} = testY;
          end
          obj.model{f} = tproc.model();
              
          if ~obj.compact, obj.mva{f} = tproc; end
            
          clear tproc;
          clear trainX;
          clear testX;
          clear trainY;
          clear testY;
          
        end
        
        % return unique model instead of cell array in case of one fold
        if length(obj.model)==1, obj.model = obj.model{1}; end
        
       
      end
      
      function [train,test] = create_folds(obj,Y)
      % get predefined sequence for train and test folds
      
      if strcmp(obj.type,'loo') || strcmp(obj.type,'bloo')
        
        if iscell(Y)
          obj.folds = min(cellfun(@(x)(size(x,1)),Y));
          if strcmp(obj.type,'bloo')
            nclasses = max(cellfun(@(x)(max(x)),Y));
            obj.folds = obj.folds/nclasses;
          end
        else
          obj.folds = size(Y,1);
          if strcmp(obj.type,'bloo')
            nclasses = max(Y);
            obj.folds = obj.folds/nclasses;
          end
        end
        obj.type = 'nfold';
        
      end
      
      if strcmp(obj.type,'split'), obj.folds = 1; end
      
      if obj.verbose, fprintf('fixing random number generator for reproducibility\n'); end
      
      if ~iscell(Y)
        
        test = create_test_folds(Y);
        
        train = obj.complement(Y,test);
        
        if obj.resample
          train = upsample(train,Y);
          test  = downsample(test,Y);
        end
        
      else
        
        train = cell(length(Y),1);
        test = cell(length(Y),1);
        
        for c=1:length(Y)
          
          test{c} = create_test_folds(Y{c});
          
          train{c} = obj.complement(Y{c},test{c});
          
          if obj.resample
            train{c} = upsample(train{c},Y{c});
            test{c}  = downsample(test{c},Y{c});
          end
          
        end
        
      end
      
        function y = create_test_folds(Y)
          
          % use the same ordering for multiple datasets by reinitializing the random number generator
          try
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
          catch
            rand('twister',1); randn('state',1);
          end
          
          y = cell(obj.folds,1);
          
          % determine indices of labeled (non-nan) datapoints
          % nan datapoints should end up in the training set
          labeled = find(any(~isnan(Y(1:size(Y,1),:)),2));
          if obj.verbose && numel(labeled)<size(Y,1)
            fprintf('adding unlabeled samples to training set\n');
          end
          
          % randomize labeled trials
          nsamples = size(labeled,1);
          idxs = labeled(randperm(nsamples));
          
          % make sure outcomes are evenly represented whenever possible
          [t,t,idx] = unique(Y,'rows');
          mx = max(idx);
          
          switch obj.type
            
            case 'split'
              
              if obj.verbose, fprintf('creating sample indices using %.0f%% of the data for testing\n',100*(1-obj.proportion)); end
              
              if mx == nsamples % unique labeled samples
                
                y{1} = idxs(1:round((1-obj.proportion)*nsamples));
                
              else
                
                % take labeled indices
                idx = idx(idxs);
                for j=1:mx
                  iidx = find(idx == j);
                  y{1} = [y{1}; idxs(iidx(1:floor((1-obj.proportion)*numel(iidx))))];
                end
                
                % add random instances to complete the number
                n = floor((1-obj.proportion) * numel(idxs)) - numel(y{1});
                if n>0
                  x = setdiff(idxs,y{1});
                  x = x(randperm(numel(x)));
                  y{1} = [y{1}; x(1:n)];
                end
                
              end
              
            case 'nfold'
              
              if obj.verbose, fprintf('creating sample indices using %d-fold cross-validation\n',obj.folds); end
              
              if mx == nsamples % unique samples
                
                for f=1:obj.folds
                  y{f} = idxs((floor((f-1)*(length(idxs)/obj.folds))+1):floor(f*(length(idxs)/obj.folds)));
                end
                
              else
                
                % take labeled indices
                idx = idx(idxs);
                
                f=1;
                for j=1:mx
                  iidx = find(idx == j);
                  for k=1:length(iidx)
                    y{f} = [y{f}; idxs(iidx(k))];
                    f = f+1; if f > obj.folds, f=1; end
                  end
                end
                
              end
              
            otherwise
              error('unrecognized option for type');
          end
          
        end
        
        function y = downsample(F,Y)
          % make sure all classes are represented equally in the
          % train samples by sampling with replacement and in the
          % test samples by using the minimum number of trials in a
          % condition. This ensures that training makes most use of
          % the data while testing has no inherent bias in the
          % data. We lose power by removing superfluous trials.
          
          y = cell(size(F));
          
          if obj.verbose
            fprintf('balancing test samples by removing superfluous trials\n');
          end
          
          for f=1:length(F)
            
            % get unique rows of Y
            [unq,t,idx] = unique(Y,'rows');
            
            idx = idx(F{f});
            
            minsmp = min(arrayfun(@(x)(sum(idx == x)),unq));
            
            if minsmp > 0
              
              tmp = [];
              for j=1:length(unq)
                
                iidx = find(idx == unq(j));
                
                if obj.verbose && length(iidx) - minsmp ~=0
                  fprintf('removing %d superfluous samples for label %d\n',length(iidx) - minsmp,unq(j));
                end
                
                tmp = [tmp; F{f}(iidx(1:minsmp))];
                
              end
              
              y{f} = tmp(randperm(length(tmp)));
              
            end
            
          end
          
        end
        
        function y = upsample(F,Y)
          % make sure all classes are represented equally in the
          % train samples by sampling with replacement and in the
          % test samples by using the minimum number of trials in a
          % condition. This ensures that training makes most use of
          % the data while testing has no inherent bias in the
          % data. We lose power by removing superfluous trials.
          
          y = cell(size(F));
          
          if obj.verbose
            fprintf('balancing training samples by sampling with replacement\n');
          end
          
          for f=1:length(F)
            
            % get unique rows of Y
            [unq,t,idx] = unique(Y,'rows');
            
            idx = idx(F{f});
            
            maxsmp = max(arrayfun(@(x)(sum(idx == x)),unq));
            
            tmp = F{f};
            for j=1:length(unq)
              
              iidx = find(idx == unq(j));
              
              if obj.verbose && maxsmp - length(iidx) ~= 0
                fprintf('drawing %d additional samples for label %d\n',maxsmp - length(iidx),unq(j));
              end
              
              % nameclash with fieldtrip_private!
              tmp = [tmp; randsample(F{f}(iidx),maxsmp - length(iidx),true)];
              
            end
            
            y{f} = tmp(randperm(length(tmp)));
            
          end
          
        end
        
      end
      
      
      function s = statistic(obj,stat)
        % report statistics on a trained crossvalidator object; also see
        % dml.statistic
        
        if nargin<2, stat = obj.stat; end
      
        if iscell(obj.design{1}) && iscell(obj.result{1})
          
          s = zeros(length(obj.design),1);
          for c=1:length(obj.design{1}) % iterate over datasets
            D = [];
            P = [];
            for f=1:length(obj.design) % iterate over folds
              D = cat(1,D,obj.design{f}{c});
              P = cat(1,P,obj.result{f}{c});
            end
            s(c) = dml.statistic(stat,D,P);
          end
          
        else
          
          D = cell2mat(obj.design);
          P = cell2mat(obj.result);
          s = dml.statistic(stat,D,P);
        
        end
            
      end
      
    end
    
    methods(Access=private,Static=true)
      
      function y = complement(Y,folds)
        % create complement of trials in each of the folds when only training
        % or test folds are specified
        
        if iscell(Y)
          ndata = length(Y);
        else
          ndata = 1;
        end
        
        y = cell(size(folds));
        for c=1:length(folds)
          if ndata == 1
            y{c} = setdiff(1:size(Y,1),folds{c})';
          else
            for f=1:length(folds{c})
              y{c}{f} = setdiff(1:size(Y{c},1),folds{c}{f})';
            end
          end
        end
        
      end
      
    end
    
    
end