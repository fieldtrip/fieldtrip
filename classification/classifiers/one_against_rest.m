classdef one_against_rest < classifier
%ONE_AGAINST_REST one-against-rest binary classification
%
%   This class evaluates a binary classifier on all a class label and the
%   pooled data for the other class labels
%
%   EXAMPLE:
%   
%   myproc = clfproc({ ...
%        preprocessor('prefun',@(x)(log10(x))) ...
%        standardizer() ... 
%        one_against_rest('procedures',clfproc({gp()})) ...
%        });
%
%   SEE ALSO:
%   one_against_one
%   ensemble.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: one_against_rest.m,v $
%

    properties        
        procedure = clfproc({da()}); % the used classification procedure
    end
    
    methods
        function obj = one_against_rest(varargin)

            obj = obj@classifier(varargin{:});
            
            % exactly one classifier for one-against-one
            assert(~iscell(obj.procedure));
            
            if ~isa(obj.procedure,'clfproc')
              obj.procedure = clfproc(obj.procedure);
            end
            
        end
        
        function obj = train(obj,tdata,tdesign)
          
          [tdata,tdesign] = obj.check_input(tdata,tdesign);
          
          if isnan(obj.nclasses), obj.nclasses = max(tdesign(:,1)); end
          
          % transform the data such that we have a cell element for each
          % class label pair
          nclasses = obj.nclasses;
          
          data = cell(1,nclasses);
          design = cell(size(data));
          
          % create new data representation
          
          for i=1:nclasses
            
            didx = (tdesign == i);
            
            data{i} = tdata;
            design{i} = tdesign;
            design{i}(~didx,:) = 1;
            design{i}(didx,:)  = 2;
            
            %                 % resampling is worse
            %                 didx = (tdesign == i);
            %
            %                 d1 = tdata(~didx,:);
            %                 d2 = tdata(didx,:);
            %
            %                 prm = ceil(size(d2,1)*rand(size(d1,1),1));
            %
            %                 data{i} = cat(1,d1,d2(prm,:));
            %                 design{i} = [ones(size(d1,1),1); 2*ones(size(d1,1),1)];
            
          end
          
          % replicate the classifier
          if ~iscell(obj.procedure)
            procedure = obj.procedure;
            obj.procedure = cell(1,length(data));
            for j=1:length(data)
              obj.procedure{j} = procedure;
            end
          end
          
          for j=1:length(data)
            obj.procedure{j} = obj.procedure{j}.train(data{j},design{j});
          end
          
        end
        
        function post = test(obj,data)
          
          data = obj.check_input(data);
          
          nclasses = obj.nclasses;
          
          % get posteriors for all pairs
          post = zeros(size(data,1),obj.nclasses);
          for i=1:nclasses
            
            tpost = obj.procedure{i}.test(data);
            post(:,i) = tpost(:,2);
          end
          
          % normalize
          post = post ./ repmat(sum(post,2),[1 size(post,2)]);
        end
        
    end
end