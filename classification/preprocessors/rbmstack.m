classdef rbmstack < preprocessor
%RBMSTACK creates a stack of restricted Boltzmann machines
%
% OPTIONS
%  are passed to the SRBM object
%
% EXAMPLE
%  myproc = { ...    
%    standardizer() ...
%    rbmstack('rbms',mymodel.rbms) ...
%    gslr('maxgroup',10) ...
%    };
%
% SEE ALSO
%     SRBM.m
% 
% REWRITE:
%   add default network for SRBM
%   expectations for hidden units as outputs
%
%   Copyright (c) 2008, Marcel van Gerven

    properties
        
      model
      
      % used to specify which layers are returned (all if empty); note
      % that layer 1 is the visible layer!
      inlayers
      
    end

    methods
        function obj = rbmstack(varargin)
          
          obj = obj@preprocessor(varargin{:});
          
          obj.model = SRBM(varargin{:});
          
        end
        
        function p = estimate(obj,data,design)
          % create hierarchy of rbms
          
          % default architecture
          if isempty(obj.model.rbms)
            
            rbm1 = RBM('nvisible',data.nfeatures,'nhidden',30,'nconditional',1);
            rbm2 = RBM('nvisible',30,'nhidden',30,'nconditional',1);
            p.model = SRBM('rbms',{rbm1 rbm2});
            
          end
          
          p.model = p.model.train(data.X);
          
        end
        
        function M = map(obj,U)            
            % propagate activations and save features as examples
        
            if isempty(obj.inlayers)
              obj.inlayers = 1:(length(obj.parameters.model.rbms)+1);
            end
            
            R = cell(1,length(obj.parameters.model.rbms)+1);
            R{1} = U.X;
            for c=1:length(obj.parameters.model.rbms)
              
              obj.parameters.model.rbms{c}.meanfield = true;
              R{c+1} = obj.parameters.model.rbms{c}.sample_hidden(R{c},ones(size(R{c},1),1));
              
            end
                                
            % concatenation of reconstruction probabilities
            M = dataset(cat(2,R{obj.inlayers}));
            
        end
        
        function U = unmap(obj,M)
          
          % move from feature activations in inlayers to the visible layer
          error('not yet implemented');
          
        
        end
    end
end 
