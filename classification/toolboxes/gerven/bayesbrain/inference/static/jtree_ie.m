classdef jtree_ie < static_inference_engine
    % JTREE_IE junction tree inference engine class
    %
    %   Should not be called directly
    %
    %   SEE ALSO:
    %       discrete_jtree_ie
    %       canonical_jtree_ie
    %       hugin_ie
    %
    %   Copyright (C) 2008  Marcel van Gerven
    %
    %   $Log: jtree_ie.m,v $
    %
    
    properties
        clusters = {}       % explicit representation of cliques
        posteriors          % storage for posteriors
        neighbours          % neighbours of each clique
        messages            % messages that need to be sent between cliques
        potentials          % the potentials of the junction tree
        cliques             % the cliques of the junction tree
        jtree               % the junction tree graph
        sepsets             % the separating sets
        potdom              % the potential domains
        sepdom              % the separating set domains
        clqasgn             % assignments of each node to its smallest clique
        evidgraph           % keep track of to which cliques each node belongs

    end

    methods
        function obj = jtree_ie(model,varargin)
           
            % model must be a markovnet
            assert(isa(model,'markovnet'));
            
            obj.model = model;
            
            % set explicit clusters
            for i=1:2:length(varargin)
                switch varargin{i},
                    case 'clusters'
                        obj.clusters = varargin{i+1};
                end
            end            
            
        end
        function mpot = marginalize(engine,query)
            % MARGINALIZE computes the marginal for the specified query nodes using the
            % junction tree inference engine

            % find a clique which contains our query node

            % for the moment we do not allow multiple query nodes that occur in
            % different cliques. This can always be enforced by creating a combined
            % node in the model or by passing custom clusters to jtree_ie
            
            if length(query) == 1

                mpot = engine.posteriors{engine.clqasgn(query)}.marginalize(query);

            else

                for c=1:length(engine.cliques)

                    if isempty(setdiff(query,engine.cliques{c}))

                        mpot = engine.posteriors{c}.marginalize(query);

                        return
                    end
                end

                % if the marginal could not be created
                error('query nodes do not belong to the same clique');

            end


        end       
    end

    methods (Access = protected, Static = true)

        function msgs = compute_msgs(node,neighs,marked)
            % compute messages
            
            marked = [marked node];

            nmsgs = neighs{node}(~ismembc(neighs{node},sort(marked)));

            if ~isempty(nmsgs)

                msgs = [cell2mat(arrayfun(@(x)(jtree_ie.compute_msgs(x,neighs,[marked node])),nmsgs,'UniformOutput',false)) [node * ones(1,length(nmsgs)); nmsgs]];
            else
                msgs = [];
            end

        end
    end

    methods(Access = protected)
    
        function cliques = makecliques(obj)

            fprintf('triangulating model\n');

            tgraph = obj.model.g;
                      
            if any(obj.model.continuous)

                cnodes = find(obj.model.continuous());
                if ~isempty(obj.clusters) % force cliques

                    for j=1:length(obj.clusters)
                        tgraph(obj.clusters{j},obj.clusters{j}) = true;
                    end
                    tgraph = triangulate(tgraph,'weights',obj.model.size(),'strong',cnodes);

                else
                    tgraph = triangulate(tgraph,'weights',obj.model.size(),'strong',cnodes);
                end

                % strong mcs
                order = mcs(tgraph,'strong',cnodes);
                if ~order, error('graph is not chordal/decomposable'); end

            else
                if ~isempty(obj.clusters) % force cliques

                    for j=1:length(obj.clusters)
                        tgraph(obj.clusters{j},obj.clusters{j}) = true;
                    end
                    tgraph = triangulate(tgraph,'weights',obj.model.size());

                else
                    tgraph = tgraph.triangulate('weights',obj.model.size());
                end

                % compute node order using maximum cardinality search
                order = tgraph.mcs();
                if ~order, error('graph is not chordal'); end
            end

            % find the cliques of the graph having the running intersection property
            % achieved based on the perfect numbering obtained from MCS (algorithm 4.11
            % pp. 56 Cowell).

            cliques = cell(1,length(obj.model.g)); ncliques = 0;

            % preconstruct Pi_i
            pii = arrayfun(@(x)(length(myintersect(order(1:(x-1)),find(tgraph(order(x),:) | tgraph(:,order(x))')))),1:length(order));

            for i=1:length(order)

                vi = order(i);

                if vi == order(end) || pii(i+1) < 1 + pii(i)
                    ncliques = ncliques + 1;
                    cliques{ncliques} = myunion(vi,myintersect(order(1:(i-1)),find(tgraph(vi,:) | tgraph(:,vi)')));
                end
            end

            cliques = cliques(1:ncliques);

        end        
        function potentials = makepotentials(obj)
           
            % initialize potentials

            fprintf('constructing potentials\n');
                                    
            % initialize potentials
       
            potentials = obj.initialize_potentials();

            % fill potentials

            % find for each variable a clique that contains the variable domain
            for i=1:numel(obj.model.factors)

                f = obj.model.factors{i};
                
                dom = f.domain;

                % initialization and assignment
                for j=1:length(obj.cliques)

                    % if the domain of the potential is contained in a clique
                    if isempty(dom(~ismembc(dom,obj.cliques{j})))

                        potentials{j} = potentials{j} * f;
                        break;
                    end
                end
            end                        
            
        end
        function clqcosts = compute_costs(obj)
            
            % precompute costs of each clique
            modelsz = obj.model.size();
            clqcosts = arrayfun(@(i)(prod(modelsz(obj.cliques{i}))),1:length(obj.cliques));

            % assignments of each node to its smallest clique

            obj.clqasgn = zeros(1,obj.model.length());
            for n=1:obj.model.length()

                % smallest cost

                cost = Inf;
                for c=1:length(obj.cliques)

                    if ismembc(n,obj.cliques{c}) && clqcosts(c) < cost

                        obj.clqasgn(n) = c;
                        cost = clqcosts(c);
                    end
                end
            end
        end
        
    end
    
end