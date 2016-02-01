classdef multiclass_KT < dml.method
% Multiclass classification using Kamitani & Tong approach
%
%   DESCRIPTION
%   Multiclass approach by combining SVM weight vectors. Theoretically
%   doubtful but works well in practice.
%
%   REFERENCE
%   Kamitani Y, Tong F. Decoding seen and attended motion directions from
%   activity in the human visual cortex. Curr. Biol. 2006
%   Jun; 16(11):1096?102.
%
%   Kamitani Y, Tong F. Decoding the visual and subjective contents of the
%   human brain. Nat. Neurosci. 2005;8(5):679?85.
%
%   EXAMPLE
%   X = randn(15,5); Y = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]';
%   m = dml.multiclass_KT
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    svms % SVM objects
    
    nclasses % number of classes
    
  end
  
  methods
    
    function obj = multiclass_KT(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
   
      if isempty(obj.nclasses)  
        C = max(Y);
        obj.nclasses = C;
      else
        C = obj.nclasses;
      end
      
      obj.svms = cell(1,C*(C-1)/2);
      
      % train an SVM for each pair of labels
      idx = 1;
      for i=1:C
        for j=(i+1):C
         
          u = (Y==i) | (Y==j);
          Yt = Y(u); Yt(Yt==i)=1; Yt(Yt==j)=2;
          obj.svms{idx} = svmtrain(X(u,:),Yt,'AutoScale','false');  
          
          idx=idx+1;
          
        end
      end
      
    end
    
    function Y = test(obj,X)
      
      C = obj.nclasses;
      
      N = size(X,1);
      
      Y = zeros(N,C);
      
      % iterate over trials;
      for n=1:N
      
        % compute discriminant function output for each label pair
        
        D = zeros(C,C);
        
        idx=1;
        for i=1:C
          for j=(i+1):C
            
            svm = obj.svms{idx};
            
            for k=1:length(svm.SupportVectorIndices) % loop over support vectors
              D(i,j) = D(i,j) + svm.Alpha(k) * svm.SupportVectors(k,:) * X(n,:)';
            end
            D(i,j) = D(i,j) + svm.Bias;
            
            idx = idx+1;

          end
        end
        
        D = -D + D';
        
        [a,b] = max(sum(D));
        
        Y(n,b) = 1; 
              
     end
        
      
    end

    function m = model(obj)
    
      m=[];
      
    end
  
  end
 
  
end
