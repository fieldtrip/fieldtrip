classdef hugin_ie < static_inference_engine
%HUGIN_IE wrapper to the hugin inference engine class
%
%   obj = hugin_ie(model,varargin)
%
%   This engine calls the hugin library; license needed!
%
%   A hugin network is created in the default temp folder under the name
%   tmp.net. This network is subsequently used for inference using a
%   library that can be obtained from www.hugin.com
%   
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: hugin_ie.m,v $
%
   
    properties
    
        clusters; % force variables to take part in the same cluster    
    end

   methods
       function obj = hugin_ie(model,varargin)
           % constructor
                   
           assert(isa(model,'bayesnet'));
           
           obj.model = model;
           
           % write model to temporary file
           model.write([tempdir 'tmp'],'hugin');

           % call mex file
           hugin_op('init',[tempdir 'tmp.net']);

       end
       function enter_evidence(obj,evidence)
           
           enter_evidence@inference_engine(obj,evidence); 
           
           hugin_op('enter',evidence);

       end
       function pot = marginalize(obj,query)
           % find a clique which contains our query node

           % for the moment we do not allow multiple query nodes that occur in
           % different cliques. This can always be enforced by creating a combined
           % node in the model

           if length(query) > 1
               error('multiple nodes in marginal not implemented');

               % could be achieved by implementing marginalize to get the
               % full table; however evidence multiplication is not
               % supported by Hugin
           end

           % retrieve single node query by library call
           marray = hugin_op('marg',query);

           if obj.model.size(query) == 1 % continuous

               pot = cpd2pot(gaussian_cpd(query,[],[],marray(1),{},marray(2)));

           else % discrete

               pot = multinomial_pot(query,marray');

           end

       end
   end
end
