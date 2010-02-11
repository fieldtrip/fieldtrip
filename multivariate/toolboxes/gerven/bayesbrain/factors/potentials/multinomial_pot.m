classdef multinomial_pot < discrete_pot
% MULTINOMIAL_POT multinomial potential class 
%
%   obj = multinomial_pot(ddomain,p) 
%
%   ddomain specifies the discrete domain
%   p specifies the potential
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: multinomial_pot.m,v $
%
   properties
       p    % parameters
   end
   
   methods
       function obj = multinomial_pot(ddomain,p)           
           % constructor

           sz = size(p);
           obj = obj@discrete_pot(ddomain,sz(sz > 1));
           obj.p = p;
       end
       function b = eq(o1,o2)
           % check if two potentials are equal
          
           b = (o1.cdomain == o2.cdomain) && ...
               (o1.ddomain == o2.ddomain) && ...
               (o1.dsize == o2.dsize) && ...
               (o1.p == o2.p);          
       end
       function pot = observe(pot, evidence)
       % observe certain values
       
           if any(evidence)

               % get data as indices ':' 1 ':' ...
               idx = num2cell(evidence);
               idx(isnan(evidence)) = {':'};

               % we create a potential that respects the evidence
               eviddata = zeros(size(pot.p));
               eviddata(idx{:}) = 1;

               pot.p = pot.p .* eviddata;
           end
       end
       function pot = mtimes(a,b)
            % multiplication operator

            % deal with other potentials
            if ~isa(b,'multinomial_pot')
                    pot = mtimes(b,a);
                    return
            end
            
           if isempty(a) || isempty(a.ddomain)
               pot = b;
           elseif isempty(b) || isempty(b.ddomain)
               pot = a;
           else

               % create the combined domain
               bidx = ~ismembc(b.ddomain,a.ddomain);
               ddom = [a.ddomain b.ddomain(bidx)];

               % extend both potentials
               if isequal(a.ddomain,ddom)
                   pot = a;
               else
                   sz = b.dsize;
                   pot = extend(a,ddom,[a.dsize sz(bidx)]);
               end

               if isequal(b.ddomain,ddom)
                   pot.p = pot.p .* b.p;
               else
                   sz = b.dsize;
                   bext = extend(b,ddom,[a.dsize sz(bidx)]);
                   pot.p = pot.p .* bext.p;
               end

           end
       end
       function pot = mrdivide(a, b)
           % MRDIVIDE divides two potentials.
           %
           % pot = mrdivide(a, b)
           %
           % extension of potentials is not performed in case noext is true

           if isempty(a) || isempty(a.ddomain)
               pot = b;
           elseif isempty(b) || isempty(b.ddomain)
               pot = a;
           else

               % deal with other potentials
               if ~isa(b,'multinomial_pot')
                   pot = mrdivide(b,a);
                   return
               end

               % create the combined domain
               bidx = ~ismember(b.ddomain,a.ddomain);
               domain = [a.ddomain b.ddomain(bidx)];
               sz = b.dsize;
               nsizes = [a.dsize sz(bidx)];

               % extend both potentials
               if length(a.ddomain) == length(domain) && all(a.ddomain == domain)
                   pot = a;
               else
                   pot = extend(a,domain,nsizes);
               end

               if length(b.ddomain) == length(domain) && all(b.ddomain == domain)
                   bext = b;
               else
                   bext = extend(b,domain,nsizes);
               end

               % create combined potential
               pot.p = pot.p ./ bext.p;

               % deal with divide by zeros
               pot.p(isnan(pot.p)) = 0;

           end
       end
       function pot = extend(pot, ddom, dsize)

           if numel(pot.ddomain) == numel(ddom) % no extension needed

               if all(pot.ddomain == ddom) % the potential is fine

                   return;

               else % we need to swap dimensions

                   % find where the domain of the potential resides in the domain
                   % this code is a fast replacement of [i,i] = ismember(domain,pot.domain);
                   imap(pot.ddomain) = 1:length(pot.ddomain); i = ddom;
                   for k = 1:numel(ddom), i(k) = imap(ddom(k)); end

                   pot.p = permute(pot.p, i);

               end

           else % extension needed
               
               % find where the domain of the potential resides in the
               % domain
               [x,i] = myismember(ddom,pot.ddomain);                              
               
               % permute the potential such that the ordering of the dimension agrees
               % with the full domain
               if length(i(x)) > 1
                   pot.p = permute(pot.p,i(x));
               end

               % reshape parameters such that they agree with the domain
               % repeat the parameters for each singular dimension
               i(x) = pot.dsize; zeroidx = ~i; i(zeroidx) = 1;
               reps = dsize; reps(~zeroidx) = 1;
               pot.p = myrepmat(reshape(pot.p,i),reps);

           end
           
           % change potential domain
           pot.ddomain = ddom;
           pot.dsize = dsize;

       end
       function pot = marginalize(pot,query)
           % return potential over the query nodes
            
           mdom = ismember(pot.ddomain,query);
           
           pot.ddomain = pot.ddomain(mdom);
           pot.dsize = pot.dsize(mdom);
           
           % sum over all variables not in the query
           for i=find(~mdom), pot.p = sum(pot.p,i); end
           
           if length(pot.ddomain) == 1
               pot.p = squeeze(pot.p(:)); % make row vector
           else
               pot.p = squeeze(pot.p);
           end
           
           
          
       end
       function cpd = normalize(pot)
          % normalize 1-D potential to probabilities
            p = pot.p ./ repmat(sum(pot.p),size(pot.p,1),1);
            p(isnan(p)) = 0;
            
            cpd = multinomial_cpd(pot.ddomain(1),[],p);
            
       end
       function pot = canonical_pot(obj)
           % convert to canonical potential

           sz = obj.dsize; if length(sz) == 1, sz = [sz 1]; end
           
           chi = ones(sz); chi(obj.p == 0) = 0;

           nzidx = (obj.p ~= 0);
           g = zeros(sz); g(nzidx) = log(obj.p(nzidx));

           pot = canonical_pot([], obj.ddomain, chi, g, cell(sz), cell(sz));

       end
       function n = states(obj)
            n = size(obj.p,1);
       end
   end
end
