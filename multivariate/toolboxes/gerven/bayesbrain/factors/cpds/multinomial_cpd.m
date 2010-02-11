classdef multinomial_cpd < discrete_cpd
%MULTINOMIAL_CPD multinomial conditional probability distribution class
%   
%   obj = multinomial_cpd(child,dparents,p)
%
%   child is the child node index
%   dparents is a vector of discrete parent node indices
%   p is a multidimensional array where the first dimension specifies the
%   child states and subsequent dimensions specify parent states.
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: multinomial_cpd.m,v $
%

   properties
       p;       % parameters       
   end
   
   methods
       function obj = multinomial_cpd(child,dparents,p)           
           % constructor
                      
           obj = obj@discrete_cpd(child,dparents);
           
           % normalize table
           p = p ./ repmat(sum(p),size(p,1),1);
           p(isnan(p)) = 0;
           
           obj.p = p;
           obj.ess = obj.essclone();
                      
       end
       function ess = essclone(obj)
           % remove handle
           ess.p = param(zeros(size(obj.p)));           
       end
       function ess = essreset(obj)
           % reset but keep handle
           
           ess = obj.ess;           
           ess.p.value = zeros(size(obj.p));           
       end
       function pot = cpd2pot(obj)
                          
           ddomain = [obj.child obj.dparents]; 
           p = obj.p;
           pot = multinomial_pot(ddomain,p);

       end
       function d = dim(obj)
          % number of free parameters

          d = (size(obj.p,1)-1)*size(obj.p,2);           
       end
       function sz = dsize(obj)
          
           if isempty(obj.dparents)
               sz = size(obj.p,1);
           else
               sz = size(obj.p);             
           end
       end
       function n = states(obj)
           n = size(obj.p,1);
       end
       function obj = update(obj,data)
        % update expected sufficient statistics

        sz = size(obj.p);  childsz = sz(1);
        if sz(2) == 1, sz = childsz; end

        % compute counts for the family
        indices = subv2ind(sz, data);
        pfam = hist(indices, 1:prod(sz));

        % note that we add to the ESS; needs an explicit reset
        obj.ess.p.value = obj.ess.p.value + reshape(pfam, size(obj.p));

       end
       function obj = maximize(obj)
           % MAXIMIZE sets the parameters of a MULTINOMIAL_CPD to their ML/MAP values.

           obj.p = obj.ess.p.value;

           % normalize table
           obj.p = obj.p ./ repmat(sum(obj.p),size(obj.p,1),1);
           obj.p(isnan(obj.p)) = 0;
                      
       end
       function obj = disconnect(obj)
          % disconnect from other nodes
          obj = multinomial_cpd(obj.child,[],ones(obj.states(),1));
       end
       function obj = updateEM(obj,pot,indom)
           % UPDATE_ESS updates the Expected Sufficient Statistics of a MULTINOMIAL_CPD.
           % indom is not used here          
           
           % used to permute domain
           [x,i] = myismember([obj.child obj.dparents],pot.ddomain);                              
               
           if isa(pot,'multinomial_pot') % add marginal data (expectation step)

               % normalize and order domain
               p = permute(pot.p,i(x)) ./ sum(pot.p(:));
               obj.ess.p.value = obj.ess.p.value + p;

           elseif isa(pot,'canonical_pot') % convert from potential

               % normalize and order domain
               p = pot.chi .* exp(pot.g);
               if length(i) > 1, p = permute(p,i(x)) ./ sum(p(:)); end
               obj.ess.p.value = obj.ess.p.value + p;
           end
       end
       function state = sample(obj,val)
           % SAMPLE takes one sample from a MULTINOMIAL_CPD
           %
           % Val are the observed values for the parents
           %
           % Copyright (C) 2008, Marcel van Gerven
           %

           % we create a potential that respects the evidence
           %eviddata = zeros(size(cpd.P));

           % incorporate observations
           index = cell(1,numel(obj.domain)); index{1} = ':';
           for i=2:length(obj.domain), index{i} = val(i-1); end

           %eviddata(index{:}) = 1;
           %P = cpd.P .* eviddata;

           % take a sample
           state = find(cumsum(obj.p(index{:})) >= rand);
           state = state(1);


       end
       
       function l = loglik(obj,childval,parentval)
           % LOGLIK computes the log likelihood of the parameters given the
           % data sample
           %
           % Copyright (C) 2008, Marcel van Gerven
           %

           % default behaviour for unobserved cases
           if any(isnan([childval parentval]))
               l = 0;
               return;
           end
           
           % incorporate observations
           index = cell(1,numel(obj.domain)); index{1} = ':';
           for i=2:length(obj.domain), index{i} = parentval(i-1); end

           l = obj.p(index{:});
           l = log(l(childval));

       end
       function plot(obj,dparentvals)
           % PLOT plots the multinomial prior distribution. 
           % dparentvals must specify a parent configuration.
          
           if nargin < 2, dparentvals = []; end
           
           assert(length(obj.dparents) == length(dparentvals));
           index = cell(1,numel(obj.dparents));
           for i=1:length(obj.dparents), index{i} = dparentvals(i); end

           p = obj.p(index{:});
           
           bar(p,'k');
           ylim([0 1]);
           
           if obj.name, title(obj.name); else title(obj.child); end
        
           if ~isempty(obj.statenames), set(gca,'XTickLabel',obj.statenames); end
           
       end
           
%        function cpd = extend(cpd, dpar, nsizes)
%             % extend discrete parent domain
%             
%             ddomain = [cpd.child cpd.dparents];
% 
%             ddom = [cpd.child dpar];
%             
%             if numel(ddomain) == numel(ddom) % no extension needed
% 
%                 if all(ddomain == ddom) % the potential is fine
% 
%                     return;
% 
%                 else % we need to swap dimensions
% 
%                     % find where the domain of the potential resides in the domain
%                     % this code is a fast replacement of [i,i] = ismember(domain,pot.domain);
%                     imap(ddomain) = 1:length(ddomain); i = ddom;
%                     for k = 1:numel(ddom), i(k) = imap(ddom(k)); end
% 
%                     cpd.p = permute(cpd.p, i);
% 
%                     cpd.dparents = dpar;
% 
%                 end
% 
%             else % extension needed
% 
%                 % find where the domain of the potential resides in the
%                 % domain
%                 [x,i] = ismember(ddom,ddomain);
% 
%                 % permute the potential such that the ordering of the dimension agrees
%                 % with the full domain
%                 if length(i(x)) > 1
%                     cpd.p = permute(cpd.p,i(x));
%                 end
% 
%                 % reshape parameters such that they agree with the domain
%                 % repeat the parameters for each singular dimension
%                 i(x) = cpd.dsize; zeroidx = ~i; i(zeroidx) = 1;
%                 reps = nsizes; reps(~zeroidx) = 1;
%                 cpd.p = myrepmat(reshape(cpd.p,i),reps);
% 
%                 % change potential domain
%                 cpd.dparents = dpar;
% 
%             end
%             
%            
%        end
   end
end 
