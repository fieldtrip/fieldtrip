classdef bayesnet < graphicalmodel
%BAYESNET Bayesian network class
%
%   A Bayesian network is constructed using
%
%   bn = bayesnet(factors,varargin)
%
%   where the factors are conditional probability distributions. 
%
%   SEE ALSO:
%   graphicalmodel.m
%   cpd.m
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: bayesnet.m,v $
%

   properties
       ec        % equivalence class          
   end

   methods
       function obj = bayesnet(factors,varargin)
          
           obj = obj@graphicalmodel(factors,varargin{:});
                           
           % construct the graph from the factors
           obj.g = graph(sparse(obj.length(),obj.length()));
           for i=1:obj.length()
               obj.g([factors{i}.parents],factors{i}.child) = true;
           end
               
           % deal with equivalence classes
           if isempty(obj.ec), obj.ec = 1:obj.length(); end
             
           % equate ess parameters with each other
           for i=unique(obj.ec)
               
               eclass = find(obj.ec == i);
               
               ess = obj.factors{eclass(1)}.essclone();               
               for j=eclass                  
                   obj.factors{j}.ess = ess;                   
               end              
           end
           
       end
       function write(obj,filename,type,varargin)
           % write a bn to file (hugin or dot format)
           
           switch type

               case 'hugin'

                   bn2hugin(obj,filename,varargin{:});

               case 'dot'

                   bn2dot(obj, filename,varargin{:});

               otherwise
                   error('unrecognized format');
           end
       end
       function data = simulate(obj,nsamples)
          
           sorted = obj.g.topological_sort();
           M = numel(sorted);

           data = nan(nsamples,M);

           domains = cellfun(@(x)(x.domain),obj.factors,'UniformOutput',false);
           for i=1:nsamples

               bevidence = false(M,1);

               for j=sorted

                   % observe evidence
                   dim = bevidence(domains{j});
                   val = data(i,domains{j}); val = val(dim);

                   % sample from cpd
                   data(i,j) = obj.factors{j}.sample(val);
                   bevidence(j) = true;

               end

           end
           
       end            
       
       function l = loglik(obj,data)
           % compute the log likelihood of the model given the data
           
           l = 0;
           
           % iterate over data samples
           for d=1:size(data,1)

               % iterate over model nodes
               for i=1:obj.length()
               
                   f = obj.factors{i};
                   l = l + f.loglik(data(d,i),data(d,f.parents()));               
               end
           end
           l = l ./ size(data,1);
       end
       
       function d = dim(obj)
           % compute dimensionality (number of free parameters) of the Bayesian network

           d = 0;
           for i=1:length(obj)
               d = d + dim(obj.factors{i});
           end
       end

       function s = bic(obj,data)
           % BIC score of the BN

           s = obj.loglik(data) - dim(obj)*log(length(obj))/2;
       end
       
       function s = aic(obj,data)
           % AIC score of the BN          
           s = obj.loglik(data) - dim(obj);           
       end
       
   end
   
   methods (Access = private)
       function bn2hugin(obj,filename,varargin)
           % write BN to hugin format

           fid = fopen(strcat(filename,'.net'), 'wt');
           fprintf(fid, 'net{}\n\n');

           % iterate over all nodes
           for i=1:obj.length()

               if isa(obj.factors{i},'gaussian_cpd')

                   fprintf(fid, 'continuous node C%d {\n\tposition = (%d %d);\n',i,round(rand*1000),round(rand*500));

                   if obj.factors{i}.name
                       fprintf(fid,'\tlabel = "%s";\n',obj.factors{i}.name);
                   end
                       
                   fprintf(fid,'}\n\n');
                   
               elseif isa(obj.factors{i},'discrete_cpd') % discrete node

                   fprintf(fid, 'node C%d \n{\n\tposition = (%d %d);\n\tstates = ( ',i,round(rand*1000),round(rand*500));
                   if ~isempty(obj.factors{i}.statenames)
                       for j=1:obj.size(i), fprintf(fid, '"%s" ',obj.factors{i}.statenames{j}); end
                   else
                       for j=1:obj.size(i), fprintf(fid, '"%d" ',j); end
                   end
                   fprintf(fid, ');\n');
                   
                   if obj.factors{i}.name
                       fprintf(fid,'\tlabel = "%s";\n',obj.factors{i}.name);
                   end
                   
                   fprintf(fid,'}\n\n');
               else
                   class(obj.factors{i})
                   error('unsupported distribution');
               end
           end

           % write tables
           % iterate over all nodes
           for i=1:obj.length()

               f = obj.factors{i};

               fprintf(fid, 'potential (C%d',i);
               if length(f.domain) > 1, fprintf(fid,' |'); end

               dparents = sort(f.dparents,'descend');
               for j=1:length(f.dparents)
                   fprintf(fid, ' C%d', dparents(j));
               end

               cparents = sort(f.cparents,'descend');
               for j=1:length(f.cparents)
                   fprintf(fid, ' C%d', cparents(j));
               end

               fprintf(fid,')\n{\n\tdata =\n   ');

               szc = cumprod(f.dsize);

               % write data
               if isa(f,'continuous_cpd')

                   mu = f.mu;
                   beta = f.beta;
                   sigma2 = f.sigma2;

                   for j=1:numel(mu)

                       szc = szc(szc > 1);

                       for k=1:length(szc)
                           if mod(j-1,szc(k)) == 0, fprintf(fid,'( '); end
                       end

                       % mean
                       fprintf(fid,'normal ( %f',mu(j));

                       % linear component
                       for c = 1:length(beta{j})

                           if beta{j}(c) >= 0, fprintf(fid,' + '); else fprintf(fid,' -'); end

                           fprintf(fid,'%f * C%d',abs(beta{j}(c)),f.cparents(c));
                       end

                       fprintf(fid,', %f ) ',sigma2(j));

                       for k=1:length(szc)
                           if mod(j,szc(k)) == 0, fprintf(fid,') '); end
                       end

                   end

               else % discrete node

                   data = f.p;

                   for j=1:numel(data)

                       for k=1:length(szc)
                           if mod(j-1,szc(k)) == 0, fprintf(fid,'('); end
                       end

                       fprintf(fid,'%f ',data(j));

                       for k=1:length(szc)
                           if mod(j,szc(k)) == 0, fprintf(fid,')'); end
                       end

                   end

               end

               fprintf(fid,';\n}\n\n');


           end

           fclose(fid);
       end
       function bn2dot(obj,filename,varargin)
           % write BN to GraphViz .dot format
                      
           % file type output
           ext = 'ps';
           
           % get optional parameters
           v = struct(varargin{:});
           if isfield(v,'extension'), ext = v.extension; end 
                          
           fid = fopen(strcat(filename,'.dot'), 'wt');

           fprintf(fid, 'digraph G{\n');
                   
           % write nodes
           for i=1:length(obj.factors)
               
               if ~isempty(obj.factors{i}.name)
                   fprintf(fid,'%d [label="%s"]\n',i,obj.factors{i}.name);
               else
                   fprintf(fid,'%d\n',i);
               end
           end
           
           % retrieve arcs
           for i=1:length(obj.factors)

               if isa(obj.factors{i},'gaussian_cpd')
                   % for a gaussian cpd we check for continuous arcs if
                   % the betas are consistently positive (red), negative
                   % (blue), or mixed (black) given all parent
                   % configurations.

                   b = obj.factors{i}.beta;
                   cparents = obj.factors{i}.cparents;
                   for j=1:length(cparents)

                       signs = zeros(1,numel(b));
                       for k=1:numel(b)
                          signs(k) = sign(b{k}(j)); 
                       end
                       
                       if all(signs == 1)
                           fprintf(fid,'%d -> %d [color=red];\n',cparents(j),i);
                       elseif all(signs == -1)
                           fprintf(fid,'%d -> %d [color=blue];\n',cparents(j),i);
                       else
                           fprintf(fid,'%d -> %d;\n',cparents(j),i);
                       end
                                             
                   end
                   
               else
                   for j=obj.factors{i}.cparents
                    fprintf(fid,'%d -> %d;\n',j,i);
                   end
               end
               
               for j=obj.factors{i}.dparents
                   fprintf(fid,'%d -> %d;\n',j,i);
               end
               
           end
           
           fprintf(fid, '}\n');

           fclose(fid);

           % try to create and show the obj

           if ispc, shell = 'dos'; else shell = 'unix'; end

           cmdline = strcat([shell '(''dot -T' ext ' ' strcat(filename,'.dot') ' -o ' strcat(filename,'.',ext) ''')']); % preserve trailing spaces
           r = eval(cmdline);
       end
            
   end
   
   methods (Static)
       function obj = random(N,varargin)
           % create a random Bayesian network

           sz = 2*ones(1,N);
           continuous = [];
           g = [];
           for i=1:2:length(varargin)
               switch varargin{i},
                   case 'graph', g = varargin{i+1};
                   case 'size', sz = varargin{i+1};
                   case 'continuous', continuous = varargin{i+1};
               end
           end
           
           sz(continuous) = 1;
           if isempty(continuous), continuous = find(sz == 1); end

           % specify structure

           if isempty(g) % random graph

               g = sparse(N,N);

               for i=2:N

                   parents = randperm(i-1);
                   parents = parents(1:min(i-1,1+floor(2*rand))); % min and max number of parents

                   % discrete nodes cannot have continuous parents
                   if ~ismember(i,continuous)
                       parents = setdiff(parents,continuous);
                   end

                   g(parents,i) = 1;
               end
           end

           % create CPDs

           factors = cell(1,N);
           for i=1:N

               parents = find(g(:,i))';
               if isempty(parents)
                   parents = []; 
               end
               
               if ismember(i,continuous) % gaussian rv

                   cparents = intersect(parents,continuous);
                   if isempty(cparents), cparents = []; end
                   dparents = setdiff(parents,continuous);
                   if isempty(dparents)
                       dparents = []; 
                       szd = [1 1];
                   else
                       szd = sz(dparents);
                       if length(szd) == 1, szd = [szd 1]; end            
                   end
                   
                   
                   mu = randn(szd);
                   beta = cell(szd);
                   if ~isempty(cparents)
                       for j=1:prod(szd)
                           beta{j} = randn(numel(cparents),1);
                       end
                   end
                   sigma1 = rand(szd);
                   factors{i} = gaussian_cpd(i,cparents,dparents,mu,beta,sigma1);

               else % discrete rv
    
                   if isempty(parents)
                       szd = [sz(i) 1];
                   else
                       szd = sz([i parents]);
                   end
    
                   factors{i} = multinomial_cpd(i,parents,rand(szd));
               end

           end
           
           obj = bayesnet(factors);

       end
   end
end
