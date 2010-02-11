classdef inference_engine < handle
%INFERENCE_ENGINE abstract inference engine class
%
%   Derives from handle class to prevent copying the whole object when
%   evidence is added.
%
%   SEE ALSO:
%   static_inference_engine.m
%   dynamic_inference_engine.m
%   
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: inference_engine.m,v $
%

   properties
        
        model            % a graphical model
 
        % evidence is a T x K database where NaNs indicate latent variables
        % K equals the number of variables and T the number of samples. 
        % For static data, T=1 and for dynamic data it is as long as the 
        % number of time slices. 
        evidence = [];   
 
   end

   methods
       function enter_evidence(obj,evidence)
           obj.evidence = reshape(evidence',[1 numel(evidence)]);
       end       
   end
   
   methods (Abstract)
       pot = marginalize(obj,query) % marginalize to obtain query variables       
   end
end
