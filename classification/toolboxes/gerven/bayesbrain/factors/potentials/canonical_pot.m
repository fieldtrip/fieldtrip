classdef canonical_pot < continuous_pot
% CANONICAL_POT canonical potential class 
%
%   obj = canonical_pot(cdomain,ddomain,chi,g,h,K)  
%
%   cdomain specifies the continuous domain
%   ddomain specifies the discrete domain
%   chi, g, h, K are the canonical parameters. Please see 
%       Probabilistic Networks and Expert Systems
%       Cowell, R.G., Dawid, A.P., Lauritzen, S.L., Spiegelhalter, D.J. 
%       for additional information.
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: canonical_pot.m,v $
%
   properties
       
       % parameters: canonical characteristics
       chi
       g
       h
       K
   end
   
   methods
       function obj = canonical_pot(cdomain,ddomain,chi,g,h,K)           
           % constructor
           
           sz = size(chi); sz = sz(sz > 1);
           if isempty(sz), sz = []; end
           
           obj = obj@continuous_pot(cdomain,ddomain,sz);

           obj.chi = chi;
           obj.g = g;
           obj.h = h;
           obj.K = K;
       end       
       function b = eq(o1,o2)
          % check if two potentials are equal
          
          if any(o1.chi(:) ~= o2.chi(:)) || any(o1.g(:) ~= o2.g(:))
              b = false;
              return;
          end
          
          for j=1:numel(o1.h)
              
             if any(o1.h{j}(:) ~= o2.h{j}(:)) || any(o1.K{j}(:) ~= o2.K{j}(:))
                 b = false;
                 return;
             end
          end

          b = (o1.cdomain == o2.cdomain) && ...
              (o1.ddomain == o2.ddomain) && ...
              (o1.dsize == o2.dsize);

       end
       function pot = observe(pot,evidence)
           % observe certain values

           % find discrete evidence
           devid = evidence((1+numel(pot.cdomain)):end);

           if any(devid)
 
               ddim = ~isnan(devid); % observed discrete variables

               % we create a potential that respects the evidence
               idx(1:numel(pot.ddomain)) = {':'};        
               idx(ddim) = num2cell(devid(ddim));

                pot.chi = zeros(size(pot.chi));
                pot.chi(idx{:})= 1;

           end

           % find continuous evidence

           cevid = evidence(1:numel(pot.cdomain));

           if any(cevid)

               cdim = ~isnan(cevid); % observed continuous variables
               cval = cevid(cdim); % observed values
               cval2 = (cval.^2)./2; % precomputing

               % iterate over all discrete cells
               for j=1:numel(pot.h)

                   % this is written in such a way that we accommodate all
                   % evidence in one sweep
                   
                   pot.g(j) = pot.g(j) + cval * pot.h{j}(cdim) - cval2 * diag(pot.K{j}(cdim,cdim));

                   h1 = pot.h{j}(~cdim);
                   if isempty(h1)
                       pot.h{j} = [];
                   else

                       % this one seems to be mistaken in Cowell et al.
                       % OR IS IT???????
                       K1g = pot.K{j}(~cdim,cdim); if isempty(K1g), K1g = []; end

                       % is there a more elegant way?
                       pot.h{j} = h1 - sum(repmat(cval,size(K1g,1),1) .* K1g,2);
                   end

                   pot.K{j} = pot.K{j}(~cdim,~cdim); % K11

               end

               % change the continuous domain
               pot.cdomain = pot.cdomain(~cdim);
               if isempty(pot.cdomain), pot.cdomain = []; end
           end

       end
       function pot = mtimes(a,b)
           % multiplication operator

           if isempty(a) || isempty(a.domain)
               pot = b;
           elseif isempty(b) || isempty(b.domain)
               pot = a;
           else

               % convert to canonical potential if possible
               if ~isa(b,'canonical_pot'), b = canonical_pot(b); end
               
               % should we extend?
               if ~isequal(a.cdomain,b.cdomain) || ~isequal(a.ddomain,b.ddomain)
                   
                   % create the combined domain
                   cdomain = myunion(a.cdomain,b.cdomain);

                   % create the combined domain
                   bidx = ~myismember(b.ddomain,a.ddomain);
                   ddomain = [a.ddomain b.ddomain(bidx)];

                   % extend potentials
                   if isempty(ddomain)
                       pot = a.extend(cdomain,[],[]);
                       b = b.extend(cdomain,[],[]);
                   else
                       % determine discrete domain size
                       szb = b.dsize;
                       sz = [a.dsize szb(bidx)];
                       pot = a.extend(cdomain,ddomain,sz);
                       b = b.extend(cdomain,ddomain,sz);
                   end
               else
                   pot = a;
               end
               
               % create combined potential

               pot.chi = pot.chi .* b.chi;

               pot.g = pot.g + b.g;

               if ~isempty(pot.cdomain)

                   for i=1:numel(pot.chi)
                       pot.h{i} = pot.h{i} + b.h{i};
                       pot.K{i} = pot.K{i} + b.K{i};
                   end
               end
           end

       end
       function pot = mrdivide(a,b)
           % division operator

           if isempty(a) || isempty(a.domain)
               pot = b;
           elseif isempty(b) || isempty(b.domain)
               pot = a;
           else

               % convert to canonical potential if possible
               if ~isa(b,'canonical_pot'), b = canonical_pot(b); end

               % should we extend?
               if ~isequal(a.cdomain,b.cdomain) || ~isequal(a.ddomain,b.ddomain)

                   % create the combined domain
                   cdomain = myunion(a.cdomain,b.cdomain);

                   % create the combined domain
                   bidx = ~ismembc(b.ddomain,a.ddomain);
                   ddomain = [a.ddomain b.ddomain(bidx)];

                   if isempty(ddomain)
                       pot = a.extend(cdomain,[],[]);
                       b = b.extend(cdomain,[],[]);
                   else
                       % determine discrete domain size
                       szb = b.dsize;
                       sz = [a.dsize szb(bidx)];
                       pot = a.extend(cdomain,ddomain,sz);
                       b = b.extend(cdomain,ddomain,sz);
                   end
               else
                   pot = a;
               end

               % create combined potential

               pot.chi = pot.chi .* b.chi;
               pot.g = pot.g - b.g;

               if ~isempty(pot.cdomain)
                   for i=1:numel(pot.chi)
                       pot.h{i} = pot.h{i} - b.h{i};
                       pot.K{i} = pot.K{i} - b.K{i};
                   end
               end

           end

       end       
       function pot = extend(pot,cdom,ddom,dsize)
                      
           % first deal with the continuous domain which consists of appending zeros to the parameters
           
           nc = length(cdom);
           if (length(pot.cdomain) < nc) || any(pot.cdomain ~= cdom)

               if nc

                   h = cell(size(pot.h));
                   h(:) = {zeros(nc,1)};

                   K = cell(size(pot.K));
                   K(:) = {zeros(nc)};

                   % find where the domain of the potential resides in the domain
                   [x,i] = myismember(cdom,pot.cdomain); cix = i(x);

                   for j=1:numel(h)
                       h{j}(x) = pot.h{j}(cix);
                       K{j}(x,x) = pot.K{j}(cix,cix);
                   end

                   pot.h = h;
                   pot.K = K;
               end

               pot.cdomain = cdom;
           end
                      
           % now deal with the discrete domain if unequal
           pdom = pot.ddomain;
           psmaller = (length(pdom) < length(ddom));
           if psmaller || any(pdom ~= ddom)
                           
               % resize the discrete domain to incorporate variables
               if psmaller

                   % variables that need to be added
                   dif = ~myismember(ddom,pdom);

                   sz = [ones(1,length(pdom)) dsize(dif)];
                   if length(sz)==1, sz = [sz 1]; end

                   pot.chi = myrepmat(pot.chi,sz);
                   pot.g = myrepmat(pot.g,sz);
                   pot.h = myrepmat(pot.h,sz);
                   pot.K = myrepmat(pot.K,sz);
                    
                   pdom = [pdom ddom(dif)];                   
                   
               end

               pot.ddomain = ddom;
               pot.dsize = dsize;
               
               % permute to make the ordering correct
               [x,i] = myismember(ddom,pdom); ix = i(x);

               if length(ix) > 1

                   pot.chi = permute(pot.chi, ix);
                   pot.g = permute(pot.g, ix);
                   pot.h = permute(pot.h, ix);
                   pot.K = permute(pot.K, ix);
               end
           end
       end
       function pot = marginalize(pot,query)
           % return potential over the query nodes

           % marginalize over continuous variables

           pcdomain = pot.cdomain;
           if ~isempty(pcdomain)

               % indices of continuous variables to keep
               cidx = ismember(pcdomain,query);

               if any(~cidx) % if there are variables to be marginalized over

                   p = numel(find(~cidx));

                   for j=1:numel(pot.chi)

                       % compute integral; valid only if K11 is positive definite

                       K11 = pot.K{j}(~cidx,~cidx);
                       K12 = pot.K{j}(~cidx,cidx);
                       K21 = pot.K{j}(cidx,~cidx);
                       K22 = pot.K{j}(cidx,cidx);

                       h1 = pot.h{j}(~cidx);
                       h2 = pot.h{j}(cidx);

                       % problems with inverse may arise here
                       pot.g(j) = pot.g(j) + ( p * log(2*pi) - log(det(K11)) + h1' * inv(K11) * h1 )/2;
                                       
                       if isempty(h2)
                           pot.h{j} = [];
                       else
                           pot.h{j} = h2 - K21 * inv(K11) * h1;
                       end

                       if isempty(K22)
                           pot.K{j} = [];
                       else
                           pot.K{j} = K22 - K21 * inv(K11) * K12;
                       end

                   end

                   pot.cdomain = pot.cdomain(cidx);
                   if isempty(pot.cdomain), pot.cdomain = []; end
               end
           end

           % marginalize over discrete variables

           pddomain = pot.ddomain;
           if ~isempty(pddomain)

               % inline ismember
               if numel(query) == 1
                   didx = (pddomain ~= query);
               else
                   didx = false(size(pddomain));
                   for i=1:numel(pddomain)
                       didx(i) = ~any(pddomain(i)==query(:));   % ANY returns logical.
                   end
               end

               if any(didx) % if there are variables to be marginalized over

                   % dimensions to marginalize out
                   mdim = find(didx);

                   % first check if we can do a strong marginalization
                   bstrong = true;
                   h = pot.h{1};
                   K = pot.K{1};

                   if ~isempty(h)
                       for j=2:numel(pot.chi)

                           % we allow approximate equality due to roundoff
                           % errors!                           
                           a = abs(pot.h{j} - h);
                           b = abs(pot.K{j} - K);
                           if any(a(:) > 1e-10) || any(b(:) > 1e-10)
                               bstrong = false;
                               break
                           end
                       end
                   end

                   if bstrong % strong marginalization

                       % incorporate evidence
                       pot.g = exp(pot.g) .* pot.chi;

                       % sum over all variables not in the query
                       for i=mdim, pot.g = sum(pot.g,i); end
                       pot.g = squeeze(pot.g);

                       % only take the log of possible entries
                       nzidx = (pot.g ~= 0); pot.g(nzidx) = log(pot.g(nzidx));

                       if size(pot.g,1) == 1, pot.g = pot.g'; end

                       % sum out larger dimensions first
                       smdim = sort(mdim,2,'descend');
                       for i=smdim, pot.chi = squeeze(sum(pot.chi,i)); end
                       if ndims(pot.chi) == 2 && size(pot.chi,1) == 1, pot.chi = pot.chi'; end
                       
                       % bug fix: set chi elements to either one or zero
                       pot.chi = (pot.chi ~= 0);
                      
                       pot.h = cell(size(pot.g));
                       pot.K = cell(size(pot.g));
                       for j=1:numel(pot.h)
                           pot.h{j} = h;
                           pot.K{j} = K;
                       end

                   else % weak marginalization

                       % compute the moments
                       p = zeros(size(pot.chi));
                       ksi = cell(size(pot.chi));
                       sigma = cell(size(pot.chi));

                       % convert from canonical to moment form
                       % THIS CAN BLOW UP

                       for j=1:numel(pot.chi)

                           sigma{j} = inv(pot.K{j});
                           ksi{j} = inv(pot.K{j})*pot.h{j};
                           p(j) = sqrt(det(sigma{j}))*exp(pot.g(j) + pot.h{j}'*sigma{j}*pot.h{j}/2);
                       end

                       % update chi

                       mndim = find(~didx);
                       mchi = pot.chi;
                       for i=mndim, mchi = sum(mchi,i); end
                       sz = pot.dsize; sz(didx) = 1;
                       mchi = repmat(mchi,sz);

                       p = p .* mchi;
                       p = p ./ sum(p(:));

                       % sum out larger dimensions first
                       smdim = sort(mdim,2,'descend');

                       for i=smdim, pot.chi = squeeze(sum(pot.chi,i)); end
                       if ndims(pot.chi) == 2 && size(pot.chi,1) == 1, pot.chi = pot.chi'; end

                       % bug fix: set chi elements to either one or zero
                       pot.chi = (pot.chi ~= 0);
                       
                       % p is an array of numbers
                       mp = p;
                       for i=smdim, mp = sum(mp,i); end

                       % ksi is a cell array of column vectors
                       mksi = ksi;
                       for j=1:numel(mksi)
                           mksi{j} = mksi{j} * p(j);
                       end

                       for i=mdim, mksi = cellsum(mksi,i); end

                       for j=1:numel(mksi), mksi{j} = mksi{j} ./ mp(j); end

                       % msigma is a cell array of matrices
                       msz = size(sigma); msz(~didx) = 1;
                       mksi2 = repmat(mksi,msz);

                       msigma = sigma;
                       for j=1:numel(msigma)
                           msigma{j} = (msigma{j} + (ksi{j} - mksi2{j}) * (ksi{j} - mksi2{j})') * p(j);
                       end

                       for i=mdim, msigma = cellsum(msigma,i); end
                       if ndims(msigma) == 2 && size(msigma,1) == 1, msigma = msigma'; end

                       for j=1:numel(msigma)
                           msigma{j} = msigma{j} ./ mp(j);
                       end

                       % squeeze everything
                       mp = squeeze(mp);
                       mksi = squeeze(mksi);
                       msigma = squeeze(msigma);

                       % convert back to canonical form
                       pot.K = cellfun(@inv,msigma,'UniformOutput',false);
                       pot.h = cell(size(pot.K));
                       for j=1:numel(pot.K)
                           pot.h{j} = pot.K{j}*mksi{j};
                       end
                       pot.g = zeros(size(pot.K));
                       for j=1:numel(pot.K)
                           pot.g(j) = log(mp(j)) + (log(det(pot.K{j})) - numel(pot.cdomain) * log(2*pi) - mksi{j}' * pot.K{j} * mksi{j})/2;
                       end

                   end

                   pot.ddomain = pot.ddomain(~didx); if isempty(pot.ddomain), pot.ddomain = []; end
                   pot.dsize = pot.dsize(~didx); if isempty(pot.dsize), pot.dsize = []; end

               end

           end
          
       end
       function cpd = normalize(pot)
           % normalize potential to probabilities

           if isempty(pot.h{1}) % convert to tabular cpd

               P = exp(pot.g);
               cpd = multinomial_cpd(pot.ddomain,[],pot.chi .* P);

           else

               invK = inv(cell2mat(pot.K));
               h = cell2mat(pot.h);
               cpd = gaussian_cpd(pot.cdomain(1),[],[],invK*h,{[]},invK);

           end
       end
   end
   
   methods(Static)
      
       function pot = empty(cvars,dvars,dsize)
         % return empty potential where dsize specifies the size of the
         % discrete variables
         
         h = cell(dsize);
         K = cell(dsize);

         nc = length(cvars);
         if nc
             for i=1:numel(h)
                 h{i} = zeros(nc,1);
                 K{i} = zeros(nc);
             end
         end

         % create empty potential
         pot = canonical_pot(cvars,dvars,ones(dsize),zeros(dsize),h,K);
           
       end
   end
end
