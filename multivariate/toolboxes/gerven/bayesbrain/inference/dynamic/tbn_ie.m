classdef tbn_ie < dynamic_inference_engine
%TBN_IE abstract 2tbn inference engine
%
%   constructs two junction trees for the prior and transition slice. 
%
%   NOTE:
%       uses canonical_ie in case there are continuous variables
%
%   SEE ALSO:
%       filtering_ie
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: tbn_ie.m,v $
%

    properties
       outint   % outgoing interface
       ie1      % jtree_ie for slice t = 1
       iet      % jtree_ie for slice t > 1
       
    end

    methods
       
        function obj = tbn_ie(dbn)
            
            nnodes = dbn.length();
            
            % compute outgoing interface (nodes with children in next
            % slice)
            obj.outint = find(any(dbn.g(1:(nnodes/2),((nnodes/2) + 1):nnodes),2))';

            % create junction tree for slice 1
            if any(dbn.continuous)
                obj.ie1 = canonical_jtree_ie(bayesnet(dbn.factors(1:(nnodes/2))),'clusters',{ obj.outint });                
            else
                obj.ie1 = discrete_jtree_ie(bayesnet(dbn.factors(1:(nnodes/2))),'clusters',{ obj.outint });
            end
            
            % create 1.5D model

            % remove nodes in slice 1 that are not in the interface
            factors = dbn.factors;
            for j=1:(nnodes/2)
                factors{j} = factors{j}.disconnect();
            end

            % create junction tree for slice t
            % make sure interface nodes end up in the same cluster

            if any(dbn.continuous)
                  obj.iet = canonical_jtree_ie(bayesnet(factors),'clusters',{ obj.outint obj.outint + (nnodes/2) });               
            else
                obj.iet = discrete_jtree_ie(bayesnet(factors),'clusters',{ obj.outint obj.outint + (nnodes/2) });
            end
        end

    end
    
    methods(Access = protected)
        
        function [marg,iet] = fwd(obj, evidence, outpot)
            % Forward pass for slices >1.
            %
            % returns the engine after evidence instantiation and the marginal for the
            % outgoing clique
            
            % change evidence index to the second slice
            N = obj.iet.model.length();

            evid = nan(1,N);
            evid((N/2+1):end) = evidence;

            % enter evidence and multiply in the old potential
            enter_evidence(obj.iet,evid,'vpots',{outpot});

            iet = obj.iet;
            
            % get the marginal for the outgoing interface
            marg = marginalize(iet,obj.outint + N/2);
            
            marg.cdomain = outpot.cdomain;
            marg.ddomain = outpot.ddomain;
            
        end

        function [marg,ie1] = fwd1(obj, evidence)
            % Forward pass for slice 1.
            %
            % returns the marginal for the outgoing clique and the engine after
            % evidence instantiation (needed during smoothing; can be ignored during
            % filtering)

            % enter evidence
            enter_evidence(obj.ie1,evidence);

            ie1 = obj.ie1;
            
            % get the marginal for the outgoing interface
            marg = marginalize(ie1,obj.outint);
        end
    end
end