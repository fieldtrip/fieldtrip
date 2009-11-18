classdef dbnet < bayesnet
% DBNET dynamic Bayesian network
%
%   A DBNET is defined as a Bayesian network where nodes belong to a
%   particular timeslice. It is assumed that there are two slices where the 
%   same node in each slice except the first is coupled in terms of 
%   equivalence classes.
%
%   This behavior can be changed using the following options:
%
%   'nslices' : number of slices
%   'coupled' : if false does not assume coupling
%
%   Note that node indices range over timeslices. E.g., 1..n in slice 1 and
%   n+1...2n in slice two.
%
%   If a dbn is defined then it can be unrolled into a larger model consisting
%   of T timeslices using dbn.unroll(T). This copies the last timeslice
%   that is defined in dbn. This is useful to construct finite horizon
%   (hybrid) networks.
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: dbnet.m,v $
%

   properties
       nslices = 2;   % number of slices in the dbnet
       coupled = true; % are timeslices coupled in terms of equivalence classes?
   end

   methods
       function obj = dbnet(factors,varargin)
        
           obj = obj@bayesnet(factors,varargin{:});
        
           % redefine equivalence classes to achieve coupling         
           if obj.coupled

               n = obj.length()/obj.nslices;           
  
               obj.ec = zeros(1,obj.length());
               obj.ec(1:n) = 1:n;
               for j=2:obj.nslices
                   obj.ec((1:n) + (j-1)*n) = ((1:n) + n);
               end
               
               % equate ess parameters with each other
               for i=unique(obj.ec)

                   eclass = find(obj.ec == i);

                   ess = obj.factors{eclass(1)}.essclone();
                   for j=eclass
                       obj.factors{j}.ess = ess;
                   end
               end
               
           end
           
       end
       
       function obj = learn_parameters(obj,data)
          % data is assumed to be a 2-slice representation of a timeseries
          % or converted to it from a timeseries: N times x M features.

          % check if the data is already in proper 2-slice format
           if size(data,2) == obj.length()
               
               obj = obj.learn_parameters@graphicalmodel(data);
               
           else % reshape the N x M timeseries so it can be directly entered as input

               % place time vectors after one another
               tsz = size(data);
               tsz(1) = tsz(1) - obj.nslices + 1;
               tsz(2) = tsz(2)*obj.nslices;

               tdata = zeros(tsz);

               nfeatures = size(data,2);
               for j=1:tsz(1)
                   tdata(j,:) = reshape(data(j:(j+obj.nslices-1),:)',[1 obj.nslices*nfeatures]);
               end

               obj = obj.learn_parameters@graphicalmodel(tdata);
           end
       end
       
       function obj = unroll(obj,nslices,times)
           % UNROLL unrolls a dynamic Bayesian network for a number of time slices
           % node names are also updated to incorporate the slice number
           % and possibly the time when times is given

           % add a number of slices to the dbn
           if nslices > obj.nslices

               % number of nodes per slice
               n = obj.length()/obj.nslices;

               newfactors = cell(1,n*nslices);

               nold = (obj.nslices*n);
               newfactors(1:nold) = obj.factors;
               for j=1:(nslices-obj.nslices)

                   % copy factors
                   newfactors((nold + (j-1)*n + 1):(nold + j*n)) = obj.factors((((obj.nslices-1)*n)+1):nold);
                              
                   % adapt indices
                   for k=1:n

                       idx = nold + (j-1)*n + k;

                       newfactors{idx}.child = newfactors{idx-n}.child + n;

                       newfactors{idx}.cparents = newfactors{idx-n}.cparents;
                       if ~isempty(newfactors{idx}.cparents)
                           newfactors{idx}.cparents = newfactors{idx}.cparents + n;
                       end

                       newfactors{idx}.dparents = newfactors{idx-n}.dparents;
                       if ~isempty(newfactors{idx}.dparents)
                           newfactors{idx}.dparents = newfactors{idx}.dparents + n;
                       end
               
                   end
               end
               
               % if time indices are specified
               if nargin > 2 && ~isempty(times)
               
                   % adapt names
                   for j=1:nslices

                       % adapt indices
                       for k=1:n
                           idx = (j-1)*n + k;
                           if ~newfactors{idx}.name, newfactors{idx}.name = sprintf('%d',k); end
                           newfactors{idx}.name = [newfactors{idx}.name sprintf('(%g)',times(j))];
                       end
                   end
               end
               
               obj = dbnet(newfactors,'nslices',nslices,'coupled',obj.coupled);
           end
           
       end
       
       function bn = bayesnet(obj)
        % convert DBN to standard BN
        
            bn = bayesnet(obj.factors,'ec',obj.ec,'include',obj.include);
                    
       end

       function write(obj,filename,type,varargin)
           % write a dbn to file (dot format)
           
           switch type

               case 'dot'

                   dbn2dot(obj, filename,varargin{:});

               otherwise
                   error('unrecognized format');
           end
       end
   end
   
   methods(Access = private)
      
       
       function dbn2dot(obj,filename,varargin)
           % write DBN to GraphViz .dot format

           % file type output
           ext = 'ps';
           
           % get optional parameters
           v = varargin2struct(varargin);
           if isfield(v,'extension'), ext = v.extension; end
       
           fid = fopen(strcat(filename,'.dot'), 'wt');

           fprintf(fid, 'digraph G{\nrankdir=LR\n');

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

               for j=obj.factors{i}.cparents
                   fprintf(fid,'%d -> %d;\n',j,i);
               end

               for j=obj.factors{i}.dparents
                   fprintf(fid,'%d -> %d;\n',j,i);
               end

           end
           
           % distinguish slices
           n = obj.length / obj.nslices;           
           for s=1:obj.nslices
               
               fprintf(fid,'{ rank=same;');
               for j=(n*(s-1)+1):(n*s)
                    fprintf(fid,' %d',j);
               end
               fprintf(fid,' };\n');
           end
                   
           fprintf(fid, '}\n');

           fclose(fid);

           % try to create and show the obj

           if ispc, shell = 'dos'; else shell = 'unix'; end

           cmdline = strcat([shell '(''dot -T' ext ' ' strcat(filename,'.dot') ' -o ' strcat(filename,'.',ext) ''')']); % preserve trailing spaces
           r = eval(cmdline);
       end

   end
   
   methods(Static)
      
       function obj = random(N,varargin)
           % create a random 2-slice dynamic Bayesian network

           obj = bayesnet.random(N*2,varargin{:});
           
           obj = dbnet(obj.factors);

       end
       
   end
end