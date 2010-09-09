classdef sparse_pls < reconstructor
% SPARSE_PLS sparse partial least squares reconstructor class
%
% load 69data; X = response; Y = stimuli;
% p = mva({pls_recon('nhidden',3,'verbose',true)});
% p = p.train(X(1:100,:),Y(1:100,:));
% r = p.test(X(101:111,:));
% images(r,1:10,[2 5]);
% figure
% images(Y,101:110,[2 5]);
%
% Copyright (c) 2010, Marcel van Gerven, Tom Heskes


  properties
  
    nhidden = 5; % number of hidden units; estimated using leave-one-out if an array
    L1 = 0.4;
    L2 = 0.01;
    
    prepout; % standardization of images
      
  end

  methods
        
    function obj = sparse_pls(varargin)
      
      obj = obj@reconstructor(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
      
      if isempty(obj.prepout)
        p.prepout = mva({squasher});
      else
        p.prepout = obj.prepout;
      end
      p.prepout = p.prepout.train(Y);
      Y = p.prepout.test(Y);

      if obj.L1 ~=0 || obj.L2 ~= 0
      
        if numel(obj.nhidden) > 1
          
          assert(all(diff(obj.nhidden) > 0));
          
          nsamples = size(X,1);
          score = zeros(numel(obj.nhidden),nsamples);
          basic = sumsqr(Y)/(nsamples-1)^2*nsamples;
          for j=1:nsamples,
            
            XX = X([1:j-1,j+1:nsamples],:)';
            YY = Y([1:j-1,j+1:nsamples],:)';
            RR = YY;
            Rleftout = Y(j,:)';
            Ypredict = zeros(size(Rleftout));
            
            for i=1:length(obj.nhidden)
              
              if i==1
                [AA,BB,CC,YYY] = sparse_bls(XX,RR,obj.nhidden(i),obj.L1,obj.L2,obj.verbose);
              else
                [AA,BB,CC,YYY] = sparse_bls(XX,RR,obj.nhidden(i)-obj.nhidden(i-1),obj.L1,obj.L2,obj.verbose);
              end
              
              Rpredict = AA * bsxfun(@plus,BB'*X(j,:)',CC');
              
              score(i,j) = sumsqr(Rpredict-Rleftout);
              
              % deflation
              RR = RR - YYY;
              Rleftout = Rleftout - Rpredict;
              Ypredict = Ypredict + Rpredict;
              
              if obj.verbose
                fprintf('sample: %d of %d; nhidden: %d; score = %g\n',j,nsamples,obj.nhidden(i),score(i,j)/basic);
              end
              
            end
            
            [d1,d2] = min(score(:,j));
            fprintf('sample: %d of %d; best score = %g (%d hidden)\n',j,nsamples,d1/basic,obj.nhidden(d2));
            
          end
          
          p.score = mean(score,2)';
          [tmp,obj.nhidden] = min(p.score);
        
        end
        
        [p.A,p.B,p.C] = sparse_bls(X',Y',obj.nhidden,obj.L1,obj.L2,obj.verbose);
      
      else
        
        if obj.verbose
          fprintf('using standard PLS\n');
          [XL,YL] = plsregress([X ones(size(X,1),1)],Y,obj.nhidden);
          
          p.A = YL;
          p.B = XL(1:(size(XL,1)-1),:);
          p.C = XL(end,:);
          
        end
        
      end
      
    end
    
    function Y = map(obj,X)
      
      Z = bsxfun(@plus,X * obj.params.B,obj.params.C);
      Y = Z*obj.params.A';
            
      % invert the squashing of the design matrix
      Y = obj.params.prepout.untest(Y);
      
      % reshape to old dimensions
      Y = reshape(Y,[size(X,1) obj.outdims]);
      
    end
    
  end
  
end
