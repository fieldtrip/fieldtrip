classdef svmnative < classifier
%SVMNATIVE wrapper for the bioinformatics toolbox SVM
%
% REQUIRES:
% bioinformatics toolbox
%
% Copyright (c) 2010, Marcel van Gerven


    methods
      
      function obj = svmnative(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function p = estimate(obj,X,Y)
        
        p = svmtrain(X,Y);
        
      end
      
      function Y = map(obj,X)
        
        classes = svmclassify(obj.params,X);
             
        Y = zeros(size(classes,1),2);
        Y(:,1) = (classes == 1);
        Y(:,2) = (classes == 2);
          
      end
      
    end
end
