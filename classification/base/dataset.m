classdef dataset
% DATASET data object to facilitate handling of input and output data
%
% Copyright (c) 2010, Marcel van Gerven
%

    properties
      
      X; % the data
      
      verbose = false;
      
      dims;      % dimensions of the input data
      
      nsamples; % number of examples
      nfeatures; % number of features
      nunique;   % number of unique trials
    
    end
    
    methods                  
      
      function obj = dataset(D)
                
        obj.X = D;
        
        obj.dims = size(D);        

        obj.nsamples = size(D,1);
        
        obj.nfeatures = prod(obj.dims(2:end));
        
        [tmp,tmp,idx] = unique(D(1:obj.nsamples,:),'rows');
        obj.nunique = max(idx);
                
      end      
      
      function Y = collapse(obj)
        % collapse data to matrix
     
          Y = obj.X(1:obj.nsamples,:);     
      end
      
      function Y = subsample(obj,idx)
        % retrieve a subset of the examples indexed by idx
        
        sel = repmat({':'},[1 ndims(obj.X)]);
        
        sel{1} = idx;
        Y = obj.X(sel{:});
        
      end
      
      function Y = subset(obj,idx)
        % retrieve a subset of the features indexed by idx
        % (collapses the data)
        
        Y = obj.collapse();
        Y = Y(:,idx);
        
      end
      
      function Y = unlabeled(obj)
        % get unlabeled indices
        
        Y = find(all(isnan(obj.collapse),2));
        
      end
      
      function Y = labeled(obj)
        % get labeled indices
        
        Y = find(any(~isnan(obj.collapse),2));        
      end
        
      function image(obj,idx,subdim)
        % plot the samples specified by idx as images
        % works for matrix inputs: data.X(:,:,:)
        % and attempts to converts other dimensions to a square
        % if subdim is specified then we do a subplot with those
        % dimensions 
        
        if length(obj.dims) ~= 3
                    
          D = obj.subsample(idx);
          D = reshape(D,[numel(idx) sqrt(obj.nfeatures) sqrt(obj.nfeatures)]);
          
        end
        
        for j=1:numel(idx)
          if nargin<3 || isempty(subdim)
            figure;
          else
            subplot(subdim(1),subdim(2),j);
          end
          imagesc(squeeze(D(j,:,:)));
          colormap(gray);
          axis off;          
          if size(D,2)==size(D,3)
            axis square;
          end
        end
        
      end
       
      function scatter(obj,labels)
        % scatter plot data if it is in one, two, or three dimensions, possibly with labels 
        
        data = obj.X;
        
        if ~isempty(data)
          
          if ~exist('labels','var') || isempty(labels), labels = ones(size(data,1),1); end
          
          if size(data,2) == 1
          
            scatter(data(:,1),ones(size(data,1),1),50,labels,'filled')            
            xlabel('feature 1');
            ylim([0.8 1.2]);
            set(gca,'YTick',[]);
            title('data distribution');
            
          elseif size(data,2) == 2
            
            scatter(data(:,1),data(:,2),50,labels,'filled')
            xlabel('feature 1');
            ylabel('feature 2');
            title('data distribution');
          
          elseif size(data,2) == 3
            
            scatter3(data(:,1),data(:,2),data(:,3),50,labels,'filled')
            xlabel('feature 1');
            ylabel('feature 2');
            zlabel('feature 3');            
            title('data distribution');
            
          else
            fprintf('too many dimensions to plot\n');
          end
                    
        end
        
      end
      
      
    end
    
end