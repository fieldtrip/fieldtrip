classdef indep_recon < reconstructor
% INDEP_RECON independent reconstructor class
%
% assumes outputs can be reconstructed independently from each other
%
% EXAMPLES:
%
% % regression
% load dataset; X = dataset(response); Y = dataset(stimuli);
% p = clfproc({indep_recon('procedure',{ridge},'verbose',true)});
% p = p.train(dataset(X.subsample(1:100)),dataset(Y.subsample(1:100)));
% r = p.test(dataset(X.subsample(101:111)));
% r.image(4);
%
% % probabilistic class labels
% load dataset; X = dataset(response); 
% stimuli = (stimuli > 128) + 1; Y = dataset(stimuli);
% p = clfproc({standardizer indep_recon('procedure',{l2lr},'verbose',true,'prob',true)});
% p = p.train(dataset(X.subsample(1:100)),dataset(Y.subsample(1:100)));
% r = p.test(dataset(X.subsample(101:111)));
% r.image(1:10,[2 5]);
%
% Copyright (c) 2010, Marcel van Gerven


  properties
  
    procedure; % the procedure used for each output 
    
    dims; % dimensions of the output data
    
    prob = false; % give probabilistic output of the maximum class
    
  end

  methods
        
        function obj = indep_recon(varargin)
     
          obj = obj@reconstructor(varargin{:});
          
          assert(~isempty(obj.procedure));

          % cast to clfproc
          obj.procedure = clfproc(obj.procedure);
          
        end        
        
        function obj = train(obj,data,design)
          
          obj.dims = design.dims;
          
          proc = cell(1,design.nfeatures);
          
          for c=1:length(proc)
            if obj.verbose
              fprintf('training outcome %d of %d\n',c,length(proc));
            end
            proc{c} = obj.procedure;
            proc{c} = proc{c}.train(data,dataset(design.subset(c)));
          end
          
          obj.procedure = proc;
          
        end
        
        function out = test(obj,data)
          
          out = zeros(data.nsamples,prod(obj.dims(2:end)));
          
          proc = obj.procedure;
          for c=1:length(proc)
            if obj.verbose
              fprintf('testing outcome %d of %d\n',c,length(proc));
            end
            if obj.prob && isa(proc{c}.clfmethods{end},'classifier')
              % output will be the probability of the last class
              pp = proc{c}.test(data);
              out(:,c) = pp.subset(pp.nfeatures);
            else
              out(:,c) = proc{c}.predict(data);
            end
          end
          
          out = dataset(reshape(out,[data.nsamples obj.dims(2:end)]));
          
        end
        
    end
    
end
