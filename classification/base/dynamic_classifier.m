classdef dynamic_classifier < classifier
%DYNAMIC_CLASSIFIER class to represent classifiers that explicitly deal with 
% synchronous (finite horizon and trial based) or asynchronous 
% (infinite horizon) temporal data explicitly
%
%   For a finite horizon model the data is assumed to be of the form
%   examples x (variables x slices) whereas for an infinite horizon model
%   the data is assumed to be of the form examples x variables; i.e., here
%   examples represent consecutive time points.
% 
%   obj.horizon < Inf represents a finite horizon model
%   obj.coupled specifies whether variable parameters are coupled over time
%       slices (always true for infinite horizon models)
%
%   A finite horizon model is always converted to a standard Bayesian
%   network.
%
%   The variable of interest is always assumed to be the first variable in
%   each time slice.
%
%   TO DO
%       remove ie/engine double representation
%
%   NOTE
%       this is an abstract class which requires another classifier
%       that specifies the structure of the associated graphical model
%       graphical models are implemented using the Bayesbrain toolbox       
%
%   SEE ALSO
%       hmm
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: dynamic_classifier.m,v $
%

    properties
        ie;
        dbn; % save for info
        numvar; % number of variables per slice
        horizon = inf;  % must be <inf for synchronous (finite horizon) models
        coupled = true; % if false; we decouple all slices
        
        nclasses;
        
    end

    methods
       function obj = dynamic_classifier(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
            
            if iscell(data), error('classifier does not take multiple datasets as input'); end
            
            obj.nclasses = design.nunique;
            
            data = data.collapse();
            design = design.collapse();
                
            if isinf(obj.horizon)
               obj.numvar = 1 + size(data,2); 
            else
                obj.numvar = 1 + (size(data,2)/obj.horizon);
            end
            if mod(obj.numvar,1), error('inconsistent number of variables'); end

            % construct factors (custom code)
            factors = obj.construct_factors();

            % construct DBN
            if obj.verbose, fprintf('creating DBN of length %d\n',length(factors)); end
            dbn = dbnet(factors,'coupled',obj.coupled);
        
            if isinf(obj.horizon) % infinite horizon              
                
                % learn parameters using the standard learner
                dbn = dbn.learn_parameters([design(:,1) data]);
           
                % create inference engine
                obj.ie = filtering_ie(dbn);
                
            else % finite horizon                
                
                if obj.coupled
                    
                    % if coupled then we can learn the parameters for a
                    % two-slice model and unroll later (more efficient)
                    
                    % create two-slice data
                    ncont = size(data,2)/obj.horizon;
                    idx = 1;
                    dbndata = zeros(size(data,1)*(obj.horizon-1),obj.numvar*2);
                    for i=1:size(data,1)
                        for h=1:(obj.horizon-1)

                            % slice 1
                            dbndata(idx,1) = design(i,1);
                            dbndata(idx,2:(1+ncont)) = data(i,(ncont*(h-1)+1):(ncont*h));
                            
                            % slice 2
                            dbndata(idx,2+ncont) = design(i,1);
                            dbndata(idx,(3+ncont):end) = data(i,(ncont*h+1):(ncont*(h+1)));

                            idx = idx + 1;                            
                        
                        end
                    end
                    
                    % learn parameters using the standard learner
                    if obj.verbose, fprintf('learning parameters\n'); end
                    dbn = dbn.learn_parameters(dbndata);
           
                    % unroll the model
                    dbn = dbn.unroll(obj.horizon);

                    % convert to bayesian network
                    dbn = bayesnet(dbn);
                                        
                else

                    % not coupled so each slice must learn its own
                    % parameters
                
                    % reconstruct data to add class variables
                    dbndata = repmat(design(:,1),[1 obj.numvar*obj.horizon]);
                    for h=1:obj.horizon
                        for f=2:obj.numvar
                            dbndata(:,(h-1)*obj.numvar + f) = data(:,(h-1)*(obj.numvar-1) + (f-1));
                        end
                    end
                
                    % unroll to required size
                    dbn = dbn.unroll(obj.horizon);

                    % this unrolled dbn is actually a bn
                    % we don't have time series data as input but rather trial
                    % based data. Therefore, we call the standard BN
                    % parameter learner!
                    dbn = bayesnet(dbn);

                    % learn parameters
                    if obj.verbose, fprintf('learning parameters\n'); end
                    dbn = dbn.learn_parameters(dbndata);
                end
                
                % create inference engine
                if isempty(obj.ie)
                    if obj.verbose, fprintf('using Hugin inference engine\n'); end
                    obj.ie = hugin_ie(dbn);
                else
                    obj.ie = obj.ie(dbn);
                end
                
            end
                      
            obj.dbn = dbn;
            
       end
       function post = test(obj,data)       
           
           if iscell(data), error('classifier does not take multiple datasets as input'); end

           if obj.verbose, fprintf('computing marginals\n'); end
           
           data = data.collapse();
           
           post = zeros([size(data,1) obj.nclasses]);

           if isinf(obj.horizon) % infinite horizon

               for j=1:size(post,1)

                   % add evidence to the inference engine
                   obj.ie.enter_evidence([nan data(j,:)]);

                   % compute marginal for first variable
                   m = normalize(obj.ie.marginalize(1));

                   post(j,:) = m.p';
               end
               
               
           else % finite horizon

               % transform data to incorporate hidden variables
               dbndata = nan([size(data,1) obj.numvar*obj.horizon]);
               for h=1:obj.horizon
                   for f=2:obj.numvar
                       dbndata(:,(h-1)*obj.numvar + f) = data(:,(h-1)*(obj.numvar-1) + (f-1));
                   end
               end

               for j=1:size(post,1)

                   % add evidence to the inference engine
                   obj.ie.enter_evidence(dbndata(j,:));

                   % compute marginal for first variable
                   m = normalize(obj.ie.marginalize(1));

                   post(j,:) = m.p';
               end
           end
                      
       end

       
       function factors = construct_factors(obj)
        % This function should be overloaded. It is currently used 
        % to return the factors of a prespecified DBN.
        
        if isempty(obj.dbn)
          error('unspecified DBN');
        end
        
        factors = obj.dbn.factors;
        
       end
       
    end
end 
