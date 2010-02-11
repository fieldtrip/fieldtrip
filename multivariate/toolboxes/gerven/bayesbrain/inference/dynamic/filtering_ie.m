classdef filtering_ie < tbn_ie
%FILTERING_IE 2tbn filtering inference engine
%
%   obj = filtering_ie(dbn)
%
%   constructs a filtering engine for a dbn
%
%   obj.T specifies the current time
%   obj.enter_evidence(evidence) advances obj.T by the number of slices
%   that is implied by the evidence. obj.marginalize(query,dT) computes the
%   distribution of the query variable. If dT is specified then query is
%   evaluated at obj.T + dT
%
%   engine can be reset to zero by setting obj.alpha = [] and obj.T = 0
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: filtering_ie.m,v $
%

    properties
       
        alpha = []; % old interface potential
        T = 0;      % current time
    end

    methods
       
        function obj = filtering_ie(dbn)
            
            obj = obj@tbn_ie(dbn);
        end
        function enter_evidence(obj, evidence)
            % ENTER_EVIDENCE Call the online filter

            obj.evidence = evidence;
            
            % advance the number of slices implied by the evidence

            for u=1:size(evidence,1) % continue filtering until we reach the time of interest

                % note that we take the evidence per time slice up to t
                if u==1 && obj.T == 0
                    [obj.alpha,obj.ie1] = obj.fwd1(evidence(1,:));
                else
                    [obj.alpha,obj.iet] = obj.fwd(evidence(u,:), obj.alpha);
                end
                obj.T = obj.T + 1;
            end
            
        end
        function mpot = marginalize(obj, query, dT)
            % MARGINALIZE Compute the joint distribution on a set of nodes (filter_engine)
            %
            % mpot = marginalize(engine, query, dT)
            %
            
            % if dT > 0 then we first make dT steps into the future

            if nargin == 3
                for u=1:dT % continue filtering until we reach the time of interest

                    % note that we take the evidence per time slice up to t
                    if u==1 && obj.T == 0
                        [obj.alpha,obj.ie1] = obj.fwd1(nan(1,obj.ie1.model.length()));
                    else
                        [obj.alpha,obj.iet] = obj.fwd(nan(1,obj.ie1.model.length()), obj.alpha);
                    end
                    obj.T = obj.T + 1;
                end
            end
            
            % call 2TBN marginalization

            if obj.T == 1 % use prior slice

                mpot = marginalize(obj.ie1, query);

            else % use transition slice

                mpot = marginalize(obj.iet, query + (obj.iet.model.length()/2));
                
            end
        end
    end
end