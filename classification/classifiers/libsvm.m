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
      
      function p = estimate(obj,data,design)
        
        data = data.X;
        design = design.X;
        
        trainparams = obj.trainparams;
        
        % regularization parameter
        if isempty(strfind(trainparams,'-c'))
          K = compKernel(data,data,'linear');
          c = .1*(mean(diag(K))-mean(K(:)));
          trainparams = [trainparams ' -c ' num2str(c)];
        end
        
        if isempty(strfind(trainparams,'-b'))
          trainparams = [trainparams ' -b 1'];
        end
        
        p.model = svmtrain(design, data, trainparams);
        
      end
      
      function post = map(obj,data)
        
        testparams = obj.testparams;
        
       if isempty(strfind(testparams,'-b'))
          testparams = [testparams ' -b 1'];
        end
        
        [a,b,post] = svmpredict(ones(data.nsamples,1), data.X, obj.params.model, testparams);

        post = dataset(post);
      
      end
      
    end
end
