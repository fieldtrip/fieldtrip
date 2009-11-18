classdef gaussian_cpd < continuous_cpd
%GAUSSIAN_CPD conditional linear gaussian probability distribution class
%   
%   obj = gaussian_cpd(child,cparents,dparents,mu,beta,sigma2)
%
%   child is the child node index
%   cparents is a vector of continuous parent node indices
%   dparents is a vector of discrete parent node indices
%   mu specifies the mean (for each discrete parent configuration)
%   beta specifies the linear contributions of each continuous parent (for
%   each discrete parent configuration)
%   sigma2 specifies the variance (for each discrete parent configuration)
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: gaussian_cpd.m,v $
%

   properties
       mu      % means per discrete parent configuration
       beta    % linear terms per discrete parent configuration
       sigma2   % variance per discrete parent configuration
   end
   
   methods
       function obj = gaussian_cpd(child,cparents,dparents,mu,beta,sigma2)           
           % constructor
                    
           obj = obj@continuous_cpd(child,cparents,dparents);
           
           obj.mu = mu;
           obj.beta = beta;
           obj.sigma2 = sigma2;                 

           % expected sufficient statistics
           obj.ess = obj.essclone();
           
       end    
       function ess = essclone(obj)
        % return a new uncoupled ess

           ess = obj.essreset();

           % remove handle 
           ess.phi = param(ess.phi.value); % the IG scale parameter beta
           ess.mu  = param(ess.mu.value); % mean of the mean
           ess.tau = param(ess.tau.value); % partially determines variance of the mean
           ess.rho = param(ess.rho.value); % the IG shape parameter alpha
        
       end
       function ess = essreset(obj)
           % reset ESS parameters while maintaining the reference
           % this should allow for an easier way of initializing ESS
           % parameters
           
           cdomain = [obj.child obj.cparents];
           
           tau = cell(size(obj.mu));
           for j=1:numel(tau), tau{j} = zeros(length(cdomain),length(cdomain)); end

           mu  = cell(size(obj.mu));
           for j=1:numel(mu), mu{j} = zeros(length(cdomain),1); end

           rho = num2cell(zeros(size(obj.mu))); % rho = 2*a (shape parameter)
           phi = num2cell(zeros(size(obj.mu))); % phi = 2*b (scale parameter)
           
           ess = obj.ess;
           
           ess.tau.value = tau;
           ess.mu.value =  mu;
           ess.rho.value = rho;
           ess.phi.value = phi;
           
       end
       function pot = cpd2pot(obj) 
           % transforms cpd to canonical potential
                      
           % convert to canonical form
           sz = size(obj.mu);
           
           chi = ones(sz);
           g = zeros(sz);
           h = cell(sz);
           K = cell(sz);

           for i=1:numel(obj.mu)

               g(i) = (- (obj.mu(i)^2/(2*obj.sigma2(i))) - (log(2*pi*obj.sigma2(i))/2));

               if isempty(obj.beta)
                   h{i} = (obj.mu(i)/obj.sigma2(i));
                   K{i} = (1/obj.sigma2(i));
               else
                   h{i} = ((obj.mu(i)/obj.sigma2(i)) * [1; -obj.beta{i}]);
                   K{i} = ((1/obj.sigma2(i)) * [1 -obj.beta{i}'; -obj.beta{i} obj.beta{i}*obj.beta{i}']);
               end

           end

           pot = canonical_pot([obj.child obj.cparents],obj.dparents,chi,g,h,K);
          
       end
       function d = dim(obj)
          % number of free parameters
          
          d = numel(obj.mu) + numel(obj.sigma2) + numel(obj.beta{1})*length(obj.beta);           
       end
       function sz = dsize(obj)
           % size of the discrete parents
           
           if isempty(obj.dparents)
               sz = 1;
           else
               sz = size(obj.mu);
               sz = sz(sz > 1);
           end
       end
       function obj = disconnect(obj)
          % disconnect from other nodes
          obj = gaussian_cpd(obj.child,[],[],[],[],[]);
       end
       function obj = update(obj,data)
        % update expected sufficient statistics

        cdomain = [obj.child obj.cparents];        
        
        % iterate over discrete parent configurations        
        sz = size(obj.mu); if sz(2) == 1, sz = sz(1); end
        if sz(1) == 1
%            indices = ones(1,size(data,2)); % BUG FIX
             indices = ones(size(data,1),1);
        else
            indices = subv2ind(sz, data(:,(numel(cdomain)+1):end));
        end
      
        for j=1:numel(obj.mu);

            % select all continuous data under this discrete configuration
            bindices = (indices == j);
            y = data(bindices,1); % self
            z = [ones(length(y),1) data(bindices,2:numel(cdomain))]; % parents

            % update inverse gamma parameters
            tauprime = obj.ess.tau.value{j} + z' * z;

             if rank(tauprime) == size(tauprime,1)
                muprime = inv(tauprime) * ((obj.ess.tau.value{j} * obj.ess.mu.value{j}) + z' * y);
            else
                % treating singular matrices
                % CHECK VALIDITY OF THIS ASSUMPTION!
                muprime = pinv(tauprime) * ((obj.ess.tau.value{j} * obj.ess.mu.value{j}) + z' * y);                
            end

            % second component may become negative due to no variance in
            % the data!
            obj.ess.phi.value{j} = obj.ess.phi.value{j} + (y - z * muprime)' * y + ...
                (obj.ess.mu.value{j} - muprime)' * obj.ess.tau.value{j} * obj.ess.mu.value{j};

            obj.ess.mu.value{j} = muprime;
            obj.ess.tau.value{j} = tauprime;

            obj.ess.rho.value{j} = obj.ess.rho.value{j} + size(y,1);

        end
        
       end
       function obj = maximize(obj)
           % MAXIMIZE sets the parameters of a MULTINOMIAL_CPD to their ML/MAP values.

           for j=1:numel(obj.mu)

               % set variance equal to the mode of the inverse
               % gamma distribution

              obj.sigma2(j) = ( (obj.ess.phi.value{j}/2) / (obj.ess.rho.value{j}/2 + 1) );

               % mean and linear parents are normally distributed with the following mean
               obj.mu(j) = obj.ess.mu.value{j}(1);

               if length(obj.ess.mu.value{j}) > 1
                   obj.beta{j} = obj.ess.mu.value{j}(2:end);
               end

           end
       end
       
       function factor = updateEM(factor,tpot,upot)
           % UPDATEEM updates the Expected Sufficient Statistics of a GAUSSIAN_CPD.
           %
           % factor = update_ess(factor, tpot, upot)
           %
           % - upot denotes which variables are unobserved
           %
           % REMARK:
           % problems with singular matrices may arise from bad
           % initializations of the distributions. Variance should be
           % chosen large enough to prevent zero probability events.
           % Handling of zero probability events should be handled in the
           % future.
           %
           % NOTE:
           % this code may contain bugs...
           % this code currently depends on potentials in canonical form!
           %
           % Copyright (c) 2008, Marcel van Gerven
                    
           cdomain = [factor.child factor.cparents];
           nc = length(cdomain);
           
           % first order the domains;
           pot = tpot.extend(cdomain,factor.dparents,factor.dsize);

           % find the cdomain that is actually defined (unobserved)
           cdom = pot.h{1} ~= 0;

           % find (un)observed continuous variables
           oc = upot(1:nc);
           unobserved = isnan(oc); 
           observed = ~isnan(oc);

           % normalize chi
           chi = pot.chi ./ repmat(sum(pot.chi),size(pot.chi,1),1);
           chi(isnan(chi)) = 0;

           % get probability of discrete states
           P = exp(pot.g);
           P = P ./ sum(P(:));

           % iterate over discrete configurations
           for j=1:numel(factor.mu)

               % discrete weight
               w = chi(j) * P(j);

               if chi(j) % ignore zero probability events
                   
                   % set IG parameters
                   % compute <z' * z> , <z' * y> and <(y - z * muprime)' * y> for the posteriors

                   % convert to covariance matrix and means
                   m = (inv(pot.K{j}(cdom,cdom))*pot.h{j}(cdom));

                   fm = zeros(1,nc);
                   if any(unobserved)
                       fm(unobserved)  = m;
                   end

                   if any(observed)
                       fm(observed) = oc(observed);
                   end

                   if nc > 1 % if we have continuous parents

                       Sigma = inv(pot.K{j}(cdom,cdom));

                       % create full covariance matrix
                       fSigma = zeros(nc);

                       if any(unobserved)
                           fSigma(unobserved,unobserved) = Sigma;
                       end

                       % compute expectations

                       % compute <z' * z>
                       % Eq. 299 of the matrix cookbook
                       zzm = fm; zzm(1) = 1;
                       zz = trace(fSigma) + zzm'  * zzm;
                       zz = w * zz; % incorporate discrete

                       tauprime = factor.ess.tau.value{j} + zz;

                       % compute <z' * y>
                       zy = fSigma(:,1) + zzm' * fm(1);
                       zy = w * zy; % incorporate discrete

                       muprime = inv(tauprime) * (factor.ess.tau.value{j} * factor.ess.mu.value{j} + zy);

                       % compute <(y - z * muprime)' * y> =
                       % <y'y> - muprime' <z'y>
                       yy = fSigma(1,1) + fm(1) * fm(1);
                       zym = yy - muprime' * zy;

                   else % no continuous parents

                       zz = w; % incorporate discrete
                       if isnan(upot(1)), zy = fm(1); else zy = upot(1); end
                       zy = w * zy;

                       tauprime = factor.ess.tau.value{j} + zz;

                       muprime = inv(tauprime) * (factor.ess.tau.value{j} * factor.ess.mu.value{j} + zy);                                        
                                              
                       if isnan(upot(1))
                           zym = (fm(1) - muprime)' * fm(1);
                       else
                           zym = (upot(1) - muprime)' * upot(1);
                       end

                   end

                   % incorporate discrete
                   zym = w * zym;

                   factor.ess.phi.value{j} = factor.ess.phi.value{j} + zym + ...
                       (factor.ess.mu.value{j} - muprime)' * factor.ess.tau.value{j} * factor.ess.mu.value{j};

                   factor.ess.mu.value{j} = muprime;
                   factor.ess.tau.value{j} = tauprime;

                   factor.ess.rho.value{j} = factor.ess.rho.value{j} + w;

               end
           end
           
       end
       function state = sample(obj,val)
           % SAMPLE takes one sample from a GAUSSIAN_CPD
           %
           % Val are the observed values for the parents
           %
           % Copyright (C) 2008, Marcel van Gerven
           %
           
           % discrete evidence

           dval = val((length(obj.cparents)+1):end);

           % incorporate observations
           index = cell(1,numel(obj.dparents));
           for i=1:length(obj.dparents), index{i} = dval(i); end

           % continuous evidence

           cval = val(1:length(obj.cparents));

           mu = obj.mu(index{:});
           beta = obj.beta(index{:});
           sigma2 = obj.sigma2(index{:});

           state = sigma2(1) * randn + mu(1) ;
           if ~isempty(beta{1})
               state = state + cval * beta{1};
           end
       end
       function l = loglik(obj,childval,parentval)
          % loglik computes the log likelihood of the parameters given the
          % data sample
           
           % default behaviour for unobserved cases
           if any(isnan([childval parentval]))
               l = 0;
               return;
           end
          
           dval = parentval((length(obj.cparents)+1):end);

           % incorporate observations
           index = cell(1,numel(obj.dparents));
           for i=1:length(obj.dparents), index{i} = dval(i); end

           % continuous evidence

           cval = parentval(1:length(obj.cparents));

           mu = obj.mu(index{:});
           beta = obj.beta(index{:});
           sigma = sqrt(obj.sigma2(index{:}));

           if ~isempty(beta{1})
               mu = mu + cval * beta{1};
           end

           z = (childval - mu) ./ sigma;

           l = -.5.*z.*z - log(sqrt(2.*pi) .* sigma);
                                
       end
       function plot(obj,cparentvals,dparentvals)
           % PLOT plots the gaussian distribution. If cparentvals is
           % unspecified then it will assume a zero vector (linear terms are
           % then ignored. If dparentvals is unspecified then all gaussians
           % will be drawn
           
           if nargin < 2 || isempty(cparentvals)
               cparentvals = zeros(1,length(obj.beta{1}));
           end
           
           if nargin < 3 || isempty(dparentvals)
               
               lim = [min(obj.mu(:)) - 3*sqrt(max(abs(obj.sigma2(:)))) max(obj.mu) + 3*sqrt(max(abs(obj.sigma2(:))))];
               colors = {'k' 'r' 'g' 'b' 'c' 'm' 'y'}; colidx = 1;
               for j=1:numel(obj.mu)                   

                   if ~isempty(cparentvals)
                       fplot(@(x)(normpdf(x,obj.mu(j)+cparentvals*obj.beta{j},sqrt(obj.sigma2(j)))),lim,colors{colidx});
                   else
                       fplot(@(x)(normpdf(x,obj.mu(j),sqrt(obj.sigma2(j)))),lim,colors{colidx});
                   end
                   hold on;
                   colidx = colidx + 1; if colidx > 7, colidx = 1; end;
               end
               
           else
              
               assert(length(dparentvals)==length(obj.mu));
               index = cell(1,numel(obj.dparents));
               for i=1:length(obj.dparents), index{i} = dparentvals(i); end

               beta = obj.beta(index{:});
               mu = obj.mu(index{:});
               if ~isempty(cparentvals)
                   mu = mu + cparentvals*beta{1};
               end               
               sigma = sqrt(obj.sigma2(index{:}));
               
               lim = [mu - 3 * sigma mu + 3 * sigma];
               fplot(@(x)(normpdf(x,mu,sigma)),lim,'k');               
           end

           if obj.name, title(obj.name); end
           
       end

   end
end 
