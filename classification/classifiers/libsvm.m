classdef libsvm < classifier
%LIBSVM wrapper for the libsvm toolbox 
%
% example:
% myproc = clfproc({ ...
%    standardizer() ...
%    libsvm('trainparams','-t 0 -c 1') ...
%    });
%
% REQUIRES:
% libsvm toolbox (http://www.csie.ntu.edu.tw/~cjlin/libsvm/)
% matlab interface (http://www.csie.ntu.edu.tw/~cjlin/libsvm/#matlab)
%
% Copyright (c) 2008, Marcel van Gerven
%
% $Log: libsvm.m,v $
%

    properties
         
        model; % the svm object
        trainparams = ''; %  training parameters as a string; e.g., '-c 1 -g 0.07 -b 1'
        testparams = ''; %  test parameters as a string; e.g., '-b 1'
        
    end

    methods
      
      function obj = libsvm(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function obj = train(obj,data,design)
        
        data = data.collapse();
        design = design.collapse();
        
        % regularization parameter
        if isempty(strfind(obj.trainparams,'-c'))
          K = compKernel(data,data,'linear');
          c = .1*(mean(diag(K))-mean(K(:)));
          obj.trainparams = [obj.trainparams ' -c ' num2str(c)];
        end
        
        if isempty(strfind(obj.trainparams,'-b'))
          obj.trainparams = [obj.trainparams ' -b 1'];
        end
        
        obj.model = svmtrain(design, data, obj.trainparams);
        
      end
      
      function post = test(obj,data)
        
        data = data.collapse();
        
        if isempty(strfind(obj.testparams,'-b'))
          obj.testparams = [obj.testparams ' -b 1'];
        end
        
        [a,b,post] = svmpredict(ones(size(data,1),1), data, obj.model, obj.testparams);

        post = dataset(post);
      
      end
      
    end
end
