classdef crossvalidator < validator
%CROSSVALIDATOR crossvalidation class
%
%   OPTIONS:
%   cvfolds: the type of crossvalidation used [10] :
%       - 1                 : test on training data
%       - 0 <= cvfolds < 1  : proportion of training data
%       - cvfolds = inf     : leave-one-out crossvalidation
%       - cvfolds = N       : N-fold crossvalidation
%       - cvfolds = { [1:n] [n+1:m] ... } : prespecified partitioning
%
%   if the design contains missing labels then it will become part of the 
%   training data by default
%
%   SEE ALSO:
%   loocrossvalidator.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: crossvalidator.m,v $
%

    properties
        cvfolds = 10; % crossvalidation procedure
    end

   methods
       function obj = crossvalidator(varargin)
           
           obj = obj@validator(varargin{:});
                       
       end
       function obj = validate(obj,data,design)

         % initialize random number generator
         if ~isempty(obj.init)
           if obj.verbose
             fprintf('initializing random number generator with seed %d\n',obj.init)
           end
           RandStream.setDefaultStream(RandStream('mt19937ar','seed',obj.init));
         end
                                             
         % make sure data and design are cell arrays to simplify further
         % processing; we revert at the end of this function
         if ~iscell(data)
           data = {data};
           
           d = cell(1,length(data));
           for c=1:length(data)
             d{c} = design;
           end
           design = d;
           
         end
                
         % make sure data is a (cell array of) matrix
         if ndims(data{1}) > 2, data = obj.collapse(data); end
         
         % make data with missing labels part of training data by default
         % this allows testing of adaptive methods
         [data,design,unlab_data,unlab_design] = obj.get_labeled_unlabeled(data,design);
    
         % report some properties
         if obj.verbose
         
             nex = 0;
             for c=1:length(data)
               nex = nex + size(data{c},1);
             end

             fprintf('performing crossvalidation for %d dataset(s) with %d examples and %d features\n',length(data),nex,size(data{1},2));

         end
         
         if isscalar(obj.cvfolds) && obj.cvfolds==1 % test on training data
           
           if obj.verbose, fprintf('validating using training data for testing\n'); end

           tproc = obj.procedure.train(data,design);
           obj.post = {tproc.test(data)};
           obj.design = {design};
           nfolds = 1;
           
           if ~obj.compact, obj.procedure = {tproc}; end
           
         else % create folds and run cross-validation
       
           % randomize ordering of the data before assigning folds
           if obj.randomize
             [data,design] = obj.shuffle(data,design);
           end
                           
           % folds{v}{d} specifies the test indices for dataset d in fold v
           [folds,nfolds] = obj.create_folds(design);
                 
           if ~obj.compact, proc = cell(1,nfolds); end
           
           obj.post = cell(1,nfolds);
           obj.design = cell(1,nfolds);
           
           for c=1:nfolds % iterate over folds
                      
             if obj.verbose, fprintf('validating fold %d of %d\n',c,nfolds); end
           
             traindata = cell(1,length(data));
             testdata  = cell(1,length(data));
             
             traindesign = cell(1,length(data));
             testdesign  = cell(1,length(data));
             
             for d=1:length(folds{c}) % iterate over datasets

               % construct data and design
               
               testidx = folds{c}{d};
               
               trainidx = ones(1,size(data{d},1)); trainidx(testidx) = 0; trainidx = find(trainidx);
               trainidx = trainidx(randperm(length(trainidx)));
               testidx = testidx(randperm(length(testidx)));
               
               traindata{d} = data{d}(trainidx,:);
               testdata{d}  = data{d}(testidx,:);
               
               traindesign{d} = design{d}(trainidx,:);
               testdesign{d}  = design{d}(testidx,:);
             
             end
             
             if (islogical(obj.balanced) && obj.balanced) || strcmp(obj.balanced,'train')
               if obj.verbose, fprintf('balancing training data\n'); end
               [traindata,traindesign] = obj.balance(traindata,traindesign);
             end
             
             if (islogical(obj.balanced) && obj.balanced) || strcmp(obj.balanced,'test')
               if obj.verbose, fprintf('balancing test data\n'); end
               [testdata,testdesign] = obj.balance(testdata,testdesign);
             end
               
             % this might happen if folds end up empty
             assert(~isempty(traindata{1}) & ~isempty(testdata{1}));
             
             % add possible unlabeled data to train data
             for uc=1:length(design)
               traindata{uc} = [traindata{uc}; unlab_data{uc}];
               traindesign{uc} = [traindesign{uc}; unlab_design{uc}];
             end
             
             % use a separate procedure per fold
             % in order to save all produced results
             tproc = obj.procedure.train(traindata,traindesign);
             
             obj.post{c} = tproc.test(testdata);
             obj.design{c} = testdesign;
             
             if ~obj.compact, proc{c} = tproc; end
             
           end
           
           if ~obj.compact, obj.procedure = proc; end           
                      
         end          
         
         % revert situation for datasets
         for c=1:nfolds
           if length(obj.design{c})==1
             obj.design{c} = obj.design{c}{1};
           end
         end
                  
       end
   end
   
   methods(Access=protected)
       
       function [folds,nfolds] = create_folds(obj,design)
        % helper function to create the folds in terms of trial indices
         
          folds = obj.cvfolds;

          if iscell(folds)
            
            % trial indices are already specified
            
            nfolds = length(folds);
            
            if ~iscell(folds{1})
              for j=1:nfolds
                folds{j} = {folds{j}}; 
              end
            end
            
            return;
            
          elseif folds > 0 && folds < 1 % percentage
            
            if obj.verbose, fprintf('validating using %.0f%% of the data for testing\n',100*(1-folds)); end
            
            folds = cell(1,length(design)); % test indices per dataset
            
            nclasses = max(design{1});
            for c=1:length(design)
              
              % make sure classes are evenly represented
              if strcmp(obj.mode,'classification') && ~isinf(obj.cvfolds)
                testidx = [];
                for k=1:nclasses
                  tmpidx = find(design{c}(:,1) == k);
                  testidx = [testidx; tmpidx(1:(numel(tmpidx) - floor(obj.cvfolds*numel(tmpidx))))];
                end
              else % regression
                tmpidx = randperm(size(design{c},1));
                testidx = tmpidx(1:(numel(tmpidx) - floor(obj.cvfolds*numel(tmpidx))));
              end
              
              folds{c} = testidx;
              
            end
            
            folds = {folds}; % one cv fold at first level
            nfolds = 1;
            
          elseif isinf(folds) || (folds > 1 && rem(folds,1) == 0) % n-fold crossvalidation
            
            nclasses = max(design{1});
            
            if isinf(folds) % leave-one-out; determine # folds              
              nfolds = min(cellfun(@(x)(size(x,1)),design));
            else
              nfolds = obj.cvfolds;
            end
            
            folds = cell(1,nfolds);
            
             for f=1:nfolds
              
              folds{f} = cell(1,length(design));
              for c=1:length(design)
                
                % make sure classes are evenly represented
                if strcmp(obj.mode,'classification') && ~isinf(obj.cvfolds)
                  testidx = [];
                  for k=1:nclasses
                    tmpidx = find(design{c}(:,1) == k);
                    testidx = [testidx; tmpidx((floor((f-1)*(length(tmpidx)/nfolds))+1):floor(f*(length(tmpidx)/nfolds)))];
                  end
                else % regression
                  testidx = (floor((f-1)*(size(design{c},1)/nfolds))+1):floor(f*(size(design{c},1)/nfolds));
                end
                
               folds{f}{c} = testidx;
             end
              
            end
            
          end
       end
          
   end
end
