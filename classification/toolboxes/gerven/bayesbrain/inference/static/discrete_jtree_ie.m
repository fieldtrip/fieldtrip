classdef discrete_jtree_ie < jtree_ie
    %DISCRETE_JTREE_IE discrete junction tree inference engine class
    %
    %   discrete_jtree_ie(model,varargin)
    %
    %   is a general inference engine for discrete networks only.
    %
    %   Copyright (C) 2008  Marcel van Gerven
    %
    %   $Log: discrete_jtree_ie.m,v $
    %

    methods
        
        function obj = discrete_jtree_ie(model,varargin)
            % constructs the discrete jtree for our model                        
            
            assert(all(model.discrete));
            
            % optionally convert from BN
            if isa(model,'bayesnet') % convert from Bayesian network

                f = model.factors;
                for c=1:length(f)
                    f{c} = cpd2pot(f{c});
                end
                model = markovnet(f);
            end

            obj@jtree_ie(model,varargin{:});
           
            % add cliques
            obj.cliques = obj.makecliques();
            
            % add potentials
            obj.potentials = obj.makepotentials();
                       
            % construct the graph (junction tree) over the cliques and
            % define sepsets
            fprintf('constructing junction tree\n');

            % build an optimal join tree using Jensen's approach as
            % discussed by Darwiche
            obj.jtree = sparse(length(obj.cliques),length(obj.cliques));
                       
            % compute clique costs
            clqcosts = obj.compute_costs();
            
            % create all possible sepsets
            obj.create_sepsets(clqcosts);
                                    
            % determine neighbours to each clique
            obj.neighbours = arrayfun(@(c)(find(obj.jtree(c,:) | obj.jtree(:,c)')),1:length(obj.cliques),'UniformOutput',false);

            % initialize the messages for collect and distribute evidence

            % we build up the candidate list using a recursive call
            % this determines the messages that need to be sent

            fprintf('computing messages\n');            
            obj.messages = jtree_ie.compute_msgs(1,obj.neighbours,1);

            % assigning fields

            obj.potdom = cell(1,length(obj.potentials));
            for i=1:length(obj.potentials), if ~isempty(obj.potentials{i}), obj.potdom{i} = obj.potentials{i}.domain; end, end

            obj.sepdom = cell(1,length(obj.sepsets));
            for i=1:length(obj.sepsets), if ~isempty(obj.sepsets{i}), obj.sepdom{i} = obj.sepsets{i}.domain; end, end           

            % keep track of to which cliques each node belongs

            obj.evidgraph = sparse(obj.model.length(),length(obj.cliques));
            for i=1:obj.model.length()
                for j=1:length(obj.cliques)
                    if ismembc(i,obj.cliques{j})
                        obj.evidgraph(i,j) = 1;
                    end
                end
            end

        end
        function enter_evidence(obj,evidence,varargin)
            % ENTER_EVIDENCE constructs a large joint potential and takes evidence into
            % account

            % get variable arguments

            vpots = {}; % virtual potentials to multiply in (used by filtering engine)
            for i=1:2:length(varargin)
                switch varargin{i},
                    case 'vpots'
                        vpots = varargin{i+1};
                end
            end

            % initialization
            obj.enter_evidence@inference_engine(evidence);

            potentials = obj.potentials;
            sepsets = obj.sepsets;
            jtree = obj.jtree;

            % find all cliques that contain evidence
            evclq = any(obj.evidgraph(~isnan(obj.evidence),:),1);
            evidclq = find(evclq);
            
            % incorporate the evidence in the potentials

            for i=evidclq

                % restrict the observed nodes to their observed values
                potentials{i} = potentials{i}.observe(obj.evidence(obj.potdom{i}));
            end            

            % multiply in potentials

            if ~isempty(vpots)

                % look for the potential with the same content
                for i=1:length(vpots)

                    for j=1:length(potentials)

                        if all(ismember(vpots{i}.domain,potentials{j}.domain))
                            potentials{j} = potentials{j} * vpots{i};
                            break
                        end

                    end

                end

            end

            
            % send messages for collect evidence

            for i=1:size(obj.messages,2)

                source = obj.messages(1,i);
                sink = obj.messages(2,i);
                idx = jtree(sink,source);

                % compute new separating set potential
                if ~isempty(sepsets{idx}) && ~isempty(potentials{sink}) % && ~isempty(intersect(obj.potdom{sink},obj.sepdom{idx})) % should not happen

                    old_sepset = sepsets{idx};

                    sepsets{idx} = marginalize(potentials{sink},obj.sepdom{idx});

                    % compute new potential after evidence collection
                    potentials{source} = potentials{source} * mrdivide(sepsets{idx},old_sepset);

                end

            end

            % send messages for distribute evidence
            
            for i=size(obj.messages,2):-1:1 % fliplr

                source = obj.messages(1,i); 
                sink = obj.messages(2,i);
                idx = jtree(source,sink);
                
                % compute new separating set potential
                if ~isempty(sepsets{idx}) && ~isempty(potentials{source}) % && ~isempty(intersect(obj.potdom{source},obj.sepdom{idx})) % should not happen

%                     old_sepset = sepsets{idx};
% 
%                     sepsets{idx} = marginalize(potentials{source},obj.sepdom{idx});
% 
%                     % compute new potential after evidence collection
%                     potentials{sink} = potentials{sink} * mrdivide(sepsets{idx}, old_sepset);

                    % compute new potential after evidence collection
                    potentials{sink} = potentials{sink} * mrdivide(marginalize(potentials{source},obj.sepdom{idx}),sepsets{idx});

                end
            end

            obj.posteriors = potentials;

        end 
                
        function potentials = initialize_potentials(obj)
           
            modelsz = obj.model.size();

            potentials = cell(1,length(obj.cliques));
            for j=1:length(obj.cliques)

                sz = modelsz(obj.cliques{j});
                if length(sz)==1, sz = [sz 1]; end
                potentials{j} = multinomial_pot(obj.cliques{j},ones(sz));

            end
            
        end     
        function create_sepsets(obj,clqcosts)
           
            ncliques = length(obj.cliques);
  
            sepsets = cell(1,ncliques*(ncliques-1)); isepset = 1;
            for i=1:ncliques
                
                for j=(i+1):ncliques

                        sepsets{isepset}.domain = obj.cliques{i}(ismembc(obj.cliques{i},obj.cliques{j}));
                        sepsets{isepset}.x = i; sepsets{isepset}.y = j;
                        isepset = isepset + 1;
                end
            end
            
            sepsets = sepsets(1:(isepset-1));
            
            candidates = 1:numel(sepsets);

            % precompute masses
            masses = zeros(1,length(candidates));
            for i=1:length(candidates)
                masses(i) = numel(sepsets{i}.domain);
            end

            % compute costs for each sepset
            costs = zeros(1,length(candidates));
            for i=1:length(candidates)
                costs(i) = clqcosts(sepsets{i}.x) + clqcosts(sepsets{i}.y);
            end

            nsepset = 0;
            esepsets = cell(1,ncliques-1);
            modelsz = obj.model.size();

            % number of candidates is known beforehand
            ncand = ncliques*(ncliques-1)/2;
            candidates = 1:ncand;

            % for hybrid networks, we first connect all discrete cliques
            % this is necessary when we have disconnected components
            connected = false(size(obj.jtree));

            while nsepset < ncliques-1

                selsepset = [];
                if ~isempty(masses)

                    % find the largest mass
                    maxmass = (masses == max(masses));

                    % find the smallest cost
                    mincost = (costs(maxmass) == min(costs(maxmass)));

                    % select the sepset index with largest mass and smallest cost
                    selidx = find(maxmass);
                    selidx = selidx(mincost); selidx = selidx(1);

                    selsepset = candidates(selidx);
                end

                % occurs when we have connected the discrete component for hybrids
                if isempty(selsepset), break; end

                % remove from the list of candidates
                candidates = candidates([1:(selidx-1) (selidx+1):end]);
                masses = masses([1:(selidx-1) (selidx+1):end]);
                costs = costs([1:(selidx-1) (selidx+1):end]);

                % check if the selected sepset can be added to the join tree
                if ~connected(sepsets{selsepset}.x,sepsets{selsepset}.y)

                    nsepset = nsepset + 1;

                    x = sepsets{selsepset}.x; y = sepsets{selsepset}.y;
                    obj.jtree(x,y) = nsepset; obj.jtree(y,x) = nsepset;

                    connected(x,y) = true;
                    connected(y,x) = true;

                    % collect all neighbours
                    neighs = connected(x,:) | connected(y,:);

                    oldneighs = false(size(neighs));
                    while any(xor(oldneighs,neighs))
                        oldneighs = neighs;
                        neighs = any(connected(neighs,:));
                    end

                    connected(neighs,neighs) = true;

                    sepdom = sepsets{selsepset}.domain;

                    if ~isempty(sepdom)
                        sz = modelsz(sepdom);
                        if length(sz)==1, sz = [sz 1]; end

                        esepsets{nsepset} = multinomial_pot(sepdom,ones(sz));
                    end

                end

            end
            
            obj.sepsets = esepsets;
        end
    end
end