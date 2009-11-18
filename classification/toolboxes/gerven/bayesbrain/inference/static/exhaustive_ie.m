classdef exhaustive_ie < static_inference_engine
%EXHAUSTIVE_IE exhaustive inference engine class
%
%   This engine constructs one big potential and performs marginalization
%   by summing over hidden variables. Is used for debugging purposes and 
%   should be applied to very small models only.
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: exhaustive_ie.m,v $
%
   properties    
        jpd = [];        % joint probability distribution
   end
   
   methods
       function obj = exhaustive_ie(model)
           % constructor
        
           % optionally convert from BN
            if isa(model,'bayesnet') % convert from Bayesian network

                f = model.factors;
                for c=1:length(f)
                    f{c} = cpd2pot(f{c});
                end
                obj.model = markovnet(f);
            else
                obj.model = model;
            end           
           
       end
       function enter_evidence(obj,evidence)
           
           enter_evidence@inference_engine(obj,evidence);
            
           for i=1:obj.model.length

               pot = obj.model.factors{i};
               
               % multiply potentials
               obj.jpd = obj.jpd * pot;

           end
           
           % restrict the observed nodes to their observed values
           obj.jpd = obj.jpd.observe(evidence(obj.jpd.domain));

       end
       function mpot = marginalize(obj,query)
           mpot = obj.jpd.marginalize(query);
       end
   end
end 
