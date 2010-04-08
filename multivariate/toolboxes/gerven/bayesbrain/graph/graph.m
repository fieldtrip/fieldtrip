classdef graph < double
%GRAPH graph class
%
%   The graph is characterized by a sparse adjacency matrix
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: graph.m,v $
%

   methods
       function obj = graph(g)
          % constructor
          
           obj = obj@double(sparse(g));
       end    
       function res = iscomplete(obj,v)
           % ISCOMPLETE returns whether the nodes v are complete in the graph
           
           % convert to double
           d = double(obj);           
           res = all(all(d(v,v) + eye(length(v))));
       end
       function ne = neighbours(obj,v)
           % NEIGHBOURS returns the neigbours of v in the graph
  
           d = double(obj);           
           ne = find(d(v,:) | d(:,v)');
           
           if isempty(ne), ne = []; end
       end
       function pa = parents(obj,v)
           % PARENTS returns parents of node v
           
           d = double(obj);           
           pa = find(d(:,v))';
           
           if isempty(pa), pa = []; end
       end
       function ch = children(obj,v)
           % CHILDREN returns children of node v
           
           d = double(obj);           
           ch = find(d(v,:))';
           
           if isempty(ch), ch = []; end
       end     
       function order = nodeorder(obj)
          % return node order of the graph
          % i.e., parents come before children
          
          d = full(double(obj));

          order = [];          
          while length(order) ~= size(d,2)
          
              for v=setdiff(1:size(d,2),order)
                 
                  if ~any(d(:,v))
                      order = [order v];
                      d(v,:) = 0;
                  end
              end              
          end
                   
       end
       function g = subgraph(obj,nodes)
          
           g = double(obj);           
           g = graph(g(nodes,nodes));
       end
       function [tgraph,tcliques] = triangulate(obj,varargin)
           % TRIANGULATE triangulates a graph
           % using the triangulation strategy of Kjaerulff
           %
           % [tgraph,tcliques] = triangulate(graph,varargin)
           %
           % graph is the to-be-triangulated graph
           % varargin:
           %   - 'weights' defines the node weights (cardinality)
           %   - 'strong' determines the continuous variables for
           %       a strong decomposition
           %
           % tgraph is the triangulated graph
           % tcliques are the obtained cliques
           %

           g = double(obj);

           n = length(g);

           % autoconnection indices
           autoc = (0:n:(n*(n-1))) + (1:n);

           cvars = [];
           for i=1:2:length(varargin)
               switch varargin{i},
                   case 'weights', weights = varargin{i+1};
                   case 'strong', cvars = varargin{i+1};
               end
           end

           % optionally distinguish discrete from continuous rvs for strong
           % triangulation
           if ~isempty(cvars)
               candidates = cvars;
               dvars = setdiff(1:n,cvars);
           else
               candidates = 1:n; dvars = [];
           end

           if ~isvarname('weights'), error('node weights are not specified'); end

           tcliques = cell(1,n); clqi = 0;

           tgraph = g;

           while candidates

               clqi = clqi+1;

               minc = 0;
               minclq = 0;
               minclqval = Inf;
               minclqweight = Inf;

               % select the next candidate
               % by choosing the smallest clique
               for c=candidates

                   % create clique
                   clq = [c find(g(c,:))];

                   % what happens if we use as a criterion the number of
                   % additional arcs added in the original graph?

                   if numel(clq) < minclqval

                       minc = c;
                       minclq = clq;
                       minclqval = length(clq);
                       minclqweight = prod(weights(clq));

                   elseif numel(clq) == minclqval

                       clqweight = prod(weights(clq));
                       if clqweight < minclqweight

                           minc = c;
                           minclq = clq;
                           minclqval = length(clq);
                           minclqweight = clqweight;
                       end
                   end
               end

               tcliques{clqi} = minclq;

               % connect all nodes in the selected clique
               g(minclq,minclq) = true;
               tgraph(minclq,minclq) = true;

               % remove autoconnections
               g(autoc) = false;
               tgraph(autoc) = false;

               % and remove node from the original graph
               g(minc,:) = false;
               g(:,minc) = false;

               % update candidates
               candidates = setdiff(candidates,minc);

               % now address discrete rvs for strong triangulation
               if isempty(candidates), candidates = dvars; dvars = []; end

           end

           tcliques = tcliques(1:clqi);

           tgraph = graph(tgraph);
       end      
       function res = mcs(obj,varargin)
           % MCS performs a maximum cardinality search (pp. 55 and 134 of Cowell)
           %

           g = double(obj);

           res = 0;

           n = length(g);

           cvars = [];
           for i=1:2:length(varargin)
               switch varargin{i},
                   case 'strong', cvars = varargin{i+1};
               end
           end
           dvars = mysetdiff(1:n,cvars);

           L = []; V = 1:n; C = zeros(1,n);

           while ~isequalset(L,V)

               u = V(~myismember(V,L));

               ud = myintersect(u,dvars);
               [x,i] = max(C(ud)); i = ud(i);

               uc = myintersect(u,cvars);
               [y,j] = max(C(uc)); j = uc(j);

               % choose continuous node
               if isempty(i) || (~isempty(y) && y > x)
                   i = j;
               end

               ne = obj.neighbours(i);

               if ~obj.iscomplete(intersect(ne,L))
                   return; % not chordal or decomposable
               end

               w = myintersect(ne,u);

               C(w) = C(w) + 1;

               L = [L i];

           end

           res = L;
       end
       function write(obj,filename,varargin)
           % write graph to graphviz dot format
           
           directed = true;
           for i=1:2:length(varargin)
               switch varargin{i},
                   case 'directed', directed = varargin{i+1};
               end
           end
           
           g = double(obj);
           g = triu(g) | tril(g);
           
           fid = fopen(strcat(filename,'.dot'), 'wt');

           fprintf(fid, 'digraph G{\n');
           %fprintf(fid,'graph [center rankdir=LR]\n');
           
           if ~directed, fprintf(fid,'edge [dir=none]\n'); end
           
           % retrieve arcs
           arcs = ind2subv(size(g),find(g));

           for i=1:size(arcs,1)
               fprintf(fid,'%d -> %d;\n',arcs(i,1),arcs(i,2));
           end

           fprintf(fid, '}\n');

           fclose(fid);

           % try to create and show the obj

           if ispc, shell = 'dos'; else shell = 'unix'; end

           cmdline = strcat([shell '(''dot -Tps ' strcat(filename,'.dot') ' -o ' strcat(filename,'.ps') ''')']); % preserve trailing spaces
           r = eval(cmdline);
           %eval([shell '(''gv ' strcat(filename,'.ps') ''')']); % show ps file
       end
       function sorted = topological_sort(obj)
           % TOPOLOGICAL_SORT generates a well-ordering from parents to children
           %
           % topological_sort(graph)
           %
           % Copyright (C) 2008, Marcel van Gerven

           g = double(obj);
           
           N = size(g,1);

           sorted = zeros(1,N);
           cnd = 1:size(g,1);

           for i=1:N
               for j=cnd
                   if ~any(g(:,j))

                       sorted(i) = j;
                       g(j,:) = false;
                       cnd = mysetdiff(cnd,j);
                       break;
                   end
               end
           end
       end
   end

end
