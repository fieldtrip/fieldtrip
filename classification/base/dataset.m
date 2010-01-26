classdef dataset
% DATASET data object to facilitate handling of input and output data
% 
% data is represented as a matrix but original dimensions are retained
% through data.dims
%
% Copyright (c) 2010, Marcel van Gerven
%

  properties
      
      X; % the data
      
      dims;       % original dimensions of the input data
      
      ndims;      % number of dimensions
      nsamples;   % number of examples
      nfeatures;  % number of features
      
      
      nunique;    % number of unique trials
    
      % standardization parameters
      mu;     % means
      sigma;  % standard deviations
      
      % whitening parameters
      wmat;   % whitening matrix
      uwmat;  % dewhitening matrix
      rdim;   % dimensionality reduction
    
    end
    
    methods                  
      
      function obj = dataset(D)
                        
        obj.dims = size(D);        
        obj.ndims = length(obj.dims);

        obj.nsamples = size(D,1);
        
        obj.nfeatures = prod(obj.dims(2:end));
        
        [tmp,tmp,idx] = unique(D(1:obj.nsamples,:),'rows');
        obj.nunique = max(idx);
                
        if obj.ndims == 2
          obj.X = D;
        else
          obj.X = D(1:obj.nsamples,:);
        end
        
      end      
            
      function Y = subsample(obj,idx)
        % retrieve dataset as a subset of the examples indexed by idx

        Y = dataset(obj.X(idx,:));
        
      end
      
      function Y = subset(obj,idx)
        % retrieve dataset as a subset of the features indexed by linear idx
  
        Y = dataset(obj.X(:,idx));
        
      end
      
      function Y = standardize(obj)
        % return standardized dataset (mean subtracted and SD of one)
        
        if ~isempty(obj.mu)
          error('dataset has already been standardized');
        end
              
        Y = obj.X;
        
        mu = mynanmean(Y);
        
        sigma = mynanstd(Y);
        sigma(sigma==0) = 1; % bug fix
        
        Y = bsxfun(@minus,Y,mu);
        Y = bsxfun(@rdivide,Y,sigma);
        
        Y = dataset(Y);
        
        Y.mu = mu;
        Y.sigma = sigma;
      
      end

      function Y = unstandardize(obj)
        % return dataset while undoing standardization
        
        if isempty(obj.mu)
          error('cannot unstandardize when data is not standardized');
        end
        
        Y = obj.X;
        
        Y = bsxfun(@times,Y,obj.sigma);
        Y = bsxfun(@plus,Y,obj.mu);
        
        Y = dataset(Y);
        
      end
      
      function Y = whiten(obj,rdim)
        % Sort the eigenvalues and select subset, and whiten
        
        if ~isempty(obj.wmat)
          error('dataset has already been whitened');
        end
        
        if nargin < 2
          rdim = obj.nfeatures;
        end
        
        [E, D] = eig(cov(obj.X,1));
        [dummy,order] = sort(diag(-D));
        
        E = E(:,order(1:rdim));
        d = real(diag(D).^(-0.5));
        D = diag(d(order(1:rdim)));
        
        wmat = D*E';
        uwmat = E*D^(-1);
        
        Y = dataset(obj.X*wmat');

        Y.wmat = wmat;
        Y.uwmat = uwmat;
        Y.rdim = rdim;
        
      end
      
      function Y = unwhiten(obj)
        
        if isempty(obj.wmat)
          error('cannot unwhiten when data is not whitened');
        end
        
        Y = dataset(obj.X*obj.uwmat');
        
        Y.wmat = [];
        Y.uwmat = [];
        Y.rdim = [];
        
      end
      
      function Y = unlabeled(obj)
        % get unlabeled indices
        
        Y = find(all(isnan(obj.X),2));
        
      end
      
      function Y = labeled(obj)
        % get labeled indices
        
        Y = find(any(~isnan(obj.X),2));        
      end
        
      function image(obj,idx,subdim)
        % plot the samples specified by idx as images
        % works for matrix inputs: data.X(:,:,:)
        % and attempts to converts other dimensions to a square
        % if subdim is specified then we do a subplot with those
        % dimensions 
        
        if nargin<2 || isempty(idx)
          idx = 1:obj.nsamples;
        end
        
        D = obj.subsample(idx).X;
       
        D = reshape(D,[numel(idx) sqrt(obj.nfeatures) sqrt(obj.nfeatures)]);
        
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
          
          if obj.nfeatures == 1
          
            scatter(data(:,1),ones(size(data,1),1),50,labels,'filled')            
            xlabel('feature 1');
            ylim([0.8 1.2]);
            set(gca,'YTick',[]);
            title('data distribution');
            
          elseif obj.nfeatures == 2
            
            scatter(data(:,1),data(:,2),50,labels,'filled')
            xlabel('feature 1');
            ylabel('feature 2');
            title('data distribution');
          
          elseif obj.nfeatures == 3
            
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