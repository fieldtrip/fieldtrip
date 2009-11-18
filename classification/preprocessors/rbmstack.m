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
%
%   $Log: rbmstack.m,v $
%

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
                
            if isempty(obj.inlayers)
              obj.inlayers = 1:(length(obj.model.rbms)+1);
            end
            
        end
        
        function obj = train(obj,data,design)
            % create hierarchy of rbms

            % default architecture
            if isempty(obj.model.rbms)
           
              rbm1 = RBM('nvisible',size(data,2),'nhidden',100,'epsilon',0.01);
              rbm2 = RBM('nvisible',100,'nhidden',20,'epsilon',0.01);
              obj.model = SRBM('rbms',{rbm1 rbm2});
            
            end
            
            obj.model = obj.model.train(data);
            
        end
        
        function data = test(obj,data)            
            % propagate activations and save features as examples
        
            R = cell(1,length(obj.model.rbms)+1);
            R{1} = data;
            for c=1:length(obj.model.rbms)
              
              obj.model.rbms{c}.meanfield = true;
              R{c+1} = obj.model.rbms{c}.sample_hidden(R{c},ones(size(R{c},1),1));
              
            end
                                
            % concatenation of reconstruction probabilities
            data = cat(2,R{obj.inlayers});
            
        end
        
    end
end 
