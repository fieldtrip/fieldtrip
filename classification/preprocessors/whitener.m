classdef whitener < preprocessor
%WHITENER whitens the data
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: whitener.m,v $
%

    properties
      
      means; % subtracted means
      
      W; % whitening matrix
      U; % unwhitening matrix 
      
      nc; % number of components
      
    end

    methods
    
        function obj = whitener(varargin)
           
          obj = obj@preprocessor(varargin{:});
          
        end
        
        function obj = train(obj,data,design)

          obj.means = mean(data);
          
          Z = bsxfun(@minus,data, obj.means);
          [E, D] = eig(cov(Z,1));
          
          % Sort the eigenvalues - decending.
          eigenvalues = sort(diag(D),'descend');

          if isempty(obj.nc)
            lastEig = size(data, 2);
          else
            lastEig = obj.nc;
          end
          
          firstEig = 1;
          oldDimension = size (data, 2);

          rankTolerance = 1e-7;
          maxLastEig = sum (diag (D) > rankTolerance);          
          
          if lastEig > maxLastEig
            lastEig = maxLastEig;
            if obj.verbose
              fprintf('Dimension reduced to %d due to the singularity of covariance matrix\n',...
                lastEig-firstEig+1);
            end
          else
            % Reduce the dimensionality of the problem.
            if obj.verbose
              if oldDimension == (lastEig - firstEig + 1)
                fprintf ('Dimension not reduced.\n');
              else
                fprintf ('Reducing dimension...\n');
              end
            end
          end
          
          if lastEig < oldDimension
            lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
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
          
          if obj.verbose, fprintf ('Selected [ %d ] dimensions.\n', sum (sel)); end
          if sum (sel) ~= (lastEig - firstEig + 1), error ('Selected a wrong number of dimensions.'); end
          
          if obj.verbose
            fprintf ('Smallest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(lastEig));
            fprintf ('Largest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(firstEig));
            fprintf ('Sum of removed eigenvalues [ %g ]\n', sum(diag(D) .* (~sel)));
          end
          
          % Select the colums which correspond to the desired range
          % of eigenvalues.
          E = E(:,sel);
          D = D(sel,sel);

          % Some more information
          if obj.verbose
            sumAll=sum(eigenvalues);
            sumUsed=sum(diag(D));
            retained = (sumUsed / sumAll) * 100;
            fprintf('[ %g ] %% of (non-zero) eigenvalues retained.\n', retained);
          end

          obj.W = sqrt(D) \ E';
          obj.U = E * sqrt(D);
                    
        end
        
        function post = test(obj,data)

          post = obj.whiten(data);       
          
        end
        
        function Y = whiten(obj,X)
          
          Z = bsxfun(@minus,X, obj.means);          
          Y = Z * obj.W';  
          
        end
        
        function X = unwhiten(obj,Y)
          
          X = Y * obj.U';
          X = bsxfun(@plus,X, obj.means); 
          
        end

    end
end
