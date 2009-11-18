classdef neighbourmerger < preprocessor
%NEIGHBOURMERGER neighbour merger
%
%   This preprocessor assumes that each feature has defined neighbours. 
%   It gradually replaces pairs of features by their average.
%   This is useful for instance when the aim is to average
%   over frequencies in order to select optimal frequency bands. 
%
%   Neighbours are specified by a cell array. 
%
%   The criterion to merge is based on comparing a metric between the
%   averaged pair and the nonaveraged pair. Currently the metric is based
%   on a validation procedure. However, we would like to implement
%   (conditional) mutual information
%
% EXAMPLE 
%
%  neighbours = cell(1,46);
%  neighbours{1} = 2;
%  for k=2:45
%     neighbours{k} = [k-1 k+1];
%  end
%  neighbours{end} = 45;
%
%   ncv = crossvalidator('procedure',clfproc({ preprocessor('prefun',@(x)(log10(x))) standardizer() ...
%        da()}),'cvfolds',10,'randomize',true);
%
%   myproc = clfproc({ ...
%     neighbourmerger('neighbours',neighbours,'cv',ncv,'verbose',true,'metric','mi') ...
%     preprocessor('prefun',@(x)(log10(x))) ...
%     standardizer() ...
%     da() ...
%     });
%  
% REMARK
%   This code is untested
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: neighbourmerger.m,v $

    properties            
        value       % mutual information of features 
        candidates  % ordered by mutual information
        neighbours  % cell array of neighbours of features
        mergeord    % features to be merged in chronological order
        merges      % the final merges
        counts      % number of merges per feature
        
        cv          % validation procedure as criterion (mi if empty)
        
        % evaluation metric of the validation procedure
        % currently MI between posteriors and actual class labels
        % as computed from a confusion matrix (assuming balanced classes)
        metric = 'mi';   
    end

    methods
       function obj = neighbourmerger(varargin)
                            
          obj = obj@preprocessor(varargin{:});
          
          if isempty(obj.neighbours)
              error('neighbours are required');
          end
       end
       function obj = train(obj,data,design)
           
           if iscell(data),error('preprocessor does not accept cell array'); end
           
           nfeatures = size(data,2);
           
           % compute function values
           obj.value = zeros(1,nfeatures);
           for j=1:size(data,2)
               
               if isempty(obj.cv)
                   obj.value(j) = mutual_information(data(:,j),design);
               else
                   obj.cv = obj.cv.validate(data(:,j),design);            
                   obj.value(j) = obj.cv.evaluate('metric',obj.metric);
               end
           end
                       
           % get ordering of the features
           [a,obj.candidates] = sort(obj.value,'descend');
           
           obj.counts = ones(1,nfeatures);
      
           while ~isempty(obj.candidates)

               % select candidate
               cnd = obj.candidates(1);
               
               % evaluate neighbours
               nmimax = -inf; nmax = 0;
               for j=1:length(obj.neighbours{cnd})
                   
                   neigh = obj.neighbours{cnd}(j);

                   % compute MI: take counts into account when averaging
                   avgdat = (obj.counts(cnd)*data(:,cnd) + ...
                       obj.counts(neigh)*data(:,neigh)) / sum(obj.counts([cnd neigh]));
                                      
                   if isempty(obj.cv)
                       nmi = mutual_information(avgdat,design);
                   else
                       obj.cv = obj.cv.validate(avgdat,design);
                       nmi = obj.cv.evaluate('metric',obj.metric);
                   end
                   
                   % select the best neighbour
                   if nmi > nmimax
                       nmimax = nmi;
                       nmax = neigh;
                   end                   
               end
               
               % if the averaged variable has higher mi than mi of candidate plus
               % mi of neighbour conditional on the candidate
               if nmax
                   
                   if isempty(obj.cv)
                       paircrit = obj.value(cnd) + conditional_mutual_information(data(:,nmax),design,data(:,cnd));
                   else

                       avgdat = data(:,[cnd nmax]);
                       obj.cv = obj.cv.validate(avgdat,design);
                       paircrit = obj.cv.evaluate('metric',obj.metric);
                   end
               end
                              
               if nmimax > paircrit
                      
                   if obj.verbose
                       fprintf('merging feature %d and %d\n',cnd,nmax);
                   end
                   
                   % merge cnd and best neighbour
                   obj.mergeord = [obj.mergeord; [cnd nmax]];
                   
                   % change data representation
                   data(:,nmax) = (obj.counts(cnd)*data(:,cnd) + obj.counts(nmax)*data(:,nmax)) / sum(obj.counts([cnd nmax]));
                   
                   % adjust mutual information
                   obj.value(nmax) = nmimax;
                   
                   % update counts
                   obj.counts(nmax) = sum(obj.counts([nmax cnd]));
                   
                   % replace cnd as neighbour by nmax
                   for k=[1:(nmax-1) (nmax+1):length(obj.neighbours)]
                       obj.neighbours{k}(obj.neighbours{k} == cnd) = nmax;
                   end
                   
                   % for nmax we take the union minus self
                   obj.neighbours{nmax} = setdiff(union(obj.neighbours{cnd},obj.neighbours{nmax}),[cnd nmax]);
                   
                   % remove data
                   midx = [1:(cnd-1) (cnd+1):size(data,2)];
                   data = data(:,midx);
                   
                   % remove associated mi
                   obj.value = obj.value(midx);
                   
                   % update counts
                   obj.counts = obj.counts(midx);
                   
                   % remove neighbours for cnd
                   obj.neighbours = obj.neighbours(midx);
                   
                   % reduce neighbour indices that are larger than cnd
                   for k=1:length(obj.neighbours)
                       obj.neighbours{k}(obj.neighbours{k} > cnd) = obj.neighbours{k}(obj.neighbours{k} > cnd) - 1;
                       obj.neighbours{k} = unique(obj.neighbours{k});
                   end
                                      
                   % reorder mi
                   [a,obj.candidates] = sort(obj.value,'descend');

               else

                   % no merge found; advance to next candidate
                   obj.candidates = obj.candidates(2:end);
                   
               end
                                                                       
           end
           
           % place mergeord in more readable form
           idx = num2cell(1:nfeatures);
           for j=1:size(obj.mergeord,1)

               merged = idx(obj.mergeord(j,:));

               % update indices
               cnd = obj.mergeord(j,1);
               nmax = obj.mergeord(j,2);

               idx{nmax} = [idx{nmax} idx{cnd}];
               idx = idx([1:(cnd-1) (cnd+1):length(idx)]);

           end         
           
           obj.merges = idx;

       end
       function data = test(obj,tdata)       
                      
           % perform the merges that are implied by the selection
           data = zeros(size(tdata,1),length(obj.merges));
           for j=1:length(obj.merges)
               data(:,j) = mean(tdata(:,obj.merges{j}),2);
           end        
           
       end

    end
end 
