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
    
      % standardization parameters
      mu;     % means
      sigma;  % standard deviations
      
      % whitening parameters
      wmat;   % whitening matrix
      uwmat;  % dewhitening matrix
      rdim;   % dimensionality reduction
    
    end
    
    methods                  
      
      function obj = dataset(D,varargin)
                        
        obj.dims = size(D);        
        obj.ndims = length(obj.dims);

        obj.nsamples = size(D,1);
        
        obj.nfeatures = prod(obj.dims(2:end));
        
        if obj.ndims == 2
          obj.X = D;
        else
          obj.X = D(1:obj.nsamples,:);
        end
        
        % parse options
        for i=1:2:length(varargin)
          if ismember(varargin{i},fieldnames(obj))
            obj.(varargin{i}) = varargin{i+1};
          end
        end
        
      end      
      
      function n = nunique(obj)
        % return the number of unique trials
        
        [tmp,tmp,idx] = unique(obj.X,'rows');
        n = max(idx);
        
      end
      
      function un = unique(obj)
        % return the unique trials
        
        n = unique(obj.X,'rows');
        
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
          warning('dataset has already been standardized');
          Y = obj;
          return
        end
              
        D = obj.X;
        
        mu = mynanmean(D);        
        sigma = mynanstd(D);
        sigma(sigma==0) = 1; % bug fix
        
        D = bsxfun(@minus,D,mu);
        D = bsxfun(@rdivide,D,sigma);
        
        Y = obj;
        Y.X = D;
        Y.dims = size(Y.X);        
        Y.ndims = length(Y.dims);
        Y.nsamples = Y.dims(1);        
        Y.nfeatures = prod(Y.dims(2:end));
        
        Y.mu = mu;
        Y.sigma = sigma;
      
      end

      function Y = unstandardize(obj)
        % return dataset while undoing standardization
        
        if isempty(obj.mu)
          warning('data has already been unstandardized');
          Y = obj;
          return
        end
        
        D = obj.X;        
        D = bsxfun(@times,D,obj.sigma);
        D = bsxfun(@plus,D,obj.mu);
        
        Y = obj;
        Y.X = D;
        Y.dims = size(Y.X);        
        Y.ndims = length(Y.dims);
        Y.nsamples = Y.dims(1);        
        Y.nfeatures = prod(Y.dims(2:end));
      
        Y.mu = [];
        Y.sigma = [];
        
      end
      
      function Y = whiten(obj,rdim)
        % Sort the eigenvalues and select subset, and whiten
        
        % standardization required if not yet performed
        obj = obj.standardize();
        
        if ~isempty(obj.wmat)
          warning('dataset has already been whitened');
          Y = obj;
          return
        end
        
        if nargin < 2
          rdim = obj.nfeatures;
        end
        
        [wmat,uwmat] = dataset.whitening_transform(obj.X,rdim);
        
        Y = obj;
        Y.X = obj.X*wmat';
        Y.dims = size(Y.X);
        Y.ndims = length(Y.dims);
        Y.nsamples = Y.dims(1);
        Y.nfeatures = prod(Y.dims(2:end));
        
        Y.wmat = wmat;
        Y.uwmat = uwmat;
        Y.rdim = rdim;
        
      end
      
      function Y = unwhiten(obj)
        
        if isempty(obj.wmat)
          warning('data has already been unwhitened');
          Y = obj;
          return;
        end
        
        Y = obj;
        Y.X = obj.X*obj.uwmat';
        Y.dims = size(Y.X);        
        Y.ndims = length(Y.dims);
        Y.nsamples = Y.dims(1);        
        Y.nfeatures = prod(Y.dims(2:end));
        
        Y.wmat = [];
        Y.uwmat = [];
        Y.rdim = [];
        
        Y = Y.unstandardize();
        
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
    
    methods(Static=true)
      
       function [wmat,uwmat] = whitening_transform(X,rdim)
        
         [E, D] = eig(cov(X,1));
         
         firstEig = 1;
         
         maxLastEig = sum (diag (D) > 1e-7); % tolerance
       
         if rdim > maxLastEig % tolerance
           error('dimension should be reduced to %d due to the singularity of covariance matrix\n',maxLastEig-firstEig+1);
         end
         
         eigenvalues = sort(diag(D),'descend');
         
         oldDimension = size(X,2);
         
         if rdim < oldDimension
           lowerLimitValue = (eigenvalues(rdim) + eigenvalues(rdim + 1)) / 2;
         else
           lowerLimitValue = eigenvalues(oldDimension) - 1;
         end
         
         lcol = diag(D) > lowerLimitValue;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Drop the larger eigenvalues
         if firstEig > 1
           higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
         else
           higherLimitValue = eigenvalues(1) + 1;
         end
         hcol = diag(D) < higherLimitValue;
         
         % Combine the results from above
         sel = lcol & hcol;
         
         % Select the colums which correspond to the desired range
         % of eigenvalues.
         E = E(:,sel);
         D = D(sel,sel);
         
         wmat = sqrt(D) \ E';
         uwmat = E * sqrt(D);
         
      end
      
    end
    
end