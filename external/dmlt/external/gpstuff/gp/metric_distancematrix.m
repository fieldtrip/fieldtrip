function metric = metric_distancematrix(varargin)
%METRIC_distancematrix An Euclidean distance metric for Gaussian process models.
%
%	Description
%
%	METRIC = METRIC_distancematrix(NIN, COMPONENTS) Constructs an data
%       structure for using a distance matrix input in a GP model. The GP
%       input can be any precomputed symmetric distance matrix. Note that
%       certain distances can only be used with certain covariance
%       functions, for more info, see e.g., Curriero (2006): On the use of Non-Euclidean
%       distance metrics in Geostatistic. Mathematical Geology
%       38(8):907-926.
%
%	The fields and (default values) in METRIC_distancematrix are:
%	  type        = 'metric_distancematrix'
%	  components  = Cell array of vectors specifying which inputs are grouped together
%                       with a same scaling parameter. 
%                       For example, the component specification {[1 2] [3]} means that
%                       distance between 3 dimensional vectors x and z is computed as
%                       r = sqrt(((x_1-z_1)^2+(x_2-z_2)^2)/l_1 + (x_3-z_3)/l_2),
%                       where l_1 and l_2 are lengthscales for corresponding component sets.
%    lengthScales      = Hyperparameters of the metric, which in this case are 
%                       lengthscales for each input component set.
%    Kstarstar         = complete precalculated distance matrix for training&test
%                        data, from which smaller parts are chosen via indices.
%                        NEEDS to be preset!
%    lengthScales_prior  = Prior structure for metric parameters. 
%                       (e.g. p.params.) [prior_unif]
%
%	METRIC = METRIC_distancematrix(METRIC, 'FIELD1', VALUE1, 'FIELD2', VALUE2, ...)
%       Set the values of fields FIELD1... to the values VALUE1... in METRIC.
%
    
% Copyright (c) 2009-2010 Heikki Peura

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

    ip=inputParser;
    ip.FunctionName = 'METRIC_DISTANCEMATRIX';
    ip.addOptional('metric', [], @isstruct);
%     ip.addParamValue('nin', [], @(x) isreal(x));
    ip.addParamValue('components',[], @(x) isempty(x) || iscell(x));
    ip.addParamValue('lengthScales',[] , @(x) isvector(x));
    ip.addParamValue('Kstarstar',[], @(x) ismatrix(x));
    ip.addParamValue('lengthScales_prior',prior_unif, ...
                   @(x) isstruct(x) || isempty(x));
%     ip.addParamValue('Kff', [], @(x) ismatrix(x));
%     ip.addParamValue('Kfu', [], @(x) ismatrix(x));
%     ip.addParamValue('Kffstar', [], @(x) ismatrix(x));
%     ip.addParamValue('Kfstaru', [], @(x) ismatrix(x));
%     ip.addParamValue('X_u', [], @(x) ismatrix(x));
    ip.parse(varargin{:});
    metric=ip.Results.metric;
    
    if isempty(metric)
        % Initialize a Gaussian process
        init=true;
    else
        % Modify a Gaussian process
        if ~isfield(metric,'type') && isequal(metric.type,'metric_distancematrix')
         error('First argument does not seem to be a metric structure')
        end
        init=false;
    end
    
    if init
        % Type
        metric.type = 'metric_distancematrix';
    end
    if init || ~ismember('components', ip.UsingDefaults)
        metric.components = ip.Results.components;
    end
    if init || ~ismember('Kstarstar', ip.UsingDefaults)
        metric.Kstarstar=ip.Results.Kstarstar;
    end
    
%     if init || ~ismember('nin', ip.UsingDefaults)
%         metric.nin = ip.Results.nin;
%     end      
%     if init || ~ismember('Kff', ip.UsingDefaults)
%         metric.Kff = ip.Results.Kff;
%     end
%     if init || ~ismember('Kffstar', ip.UsingDefaults)
%         metric.Kffstar=ip.Results.Kffstar;
%     end
%     if init || ~ismember('X_u', ip.UsingDefaults)
%         metric.X_u=ip.Results.X_u;
%     end
%     if init || ~ismember('Kfu', ip.UsingDefaults)
%         metric.Kfu=ip.Results.Kfu;
%     end
%     if init || ~ismember('Kfstaru', ip.UsingDefaults)
%         metric.Kfstaru=ip.Results.Kfstaru;
%     end
    
    if init || ~ismember('lengthScales', ip.UsingDefaults)
%         if size(metric.lengthScales) ~= size(ip.Results.lengthScales)
%             error('Incorrect number of parameters given')
%         end
        metric.lengthScales=ip.Results.lengthScales;
        if isempty(metric.lengthScales)
          metric.lengthScales = repmat(1,1,length(metric.components));
        end
    end
        
    if init || ~ismember('lengthScales_prior', ip.UsingDefaults)
        % Initialize prior structure
        metric.p=[];
        metric.p.lengthScales=ip.Results.lengthScales_prior;
    end
        
    if init
        % Set the function handles to the nested functions
        metric.fh.pak        = @metric_distancematrix_pak;
        metric.fh.unpak      = @metric_distancematrix_unpak;
        metric.fh.lp         = @metric_distancematrix_lp;
        metric.fh.lpg        = @metric_distancematrix_lpg;
        metric.fh.dist       = @metric_distancematrix_distance;        
        metric.fh.distg      = @metric_distancematrix_ghyper;
        metric.fh.ginput     = @metric_distancematrix_ginput;
        metric.fh.recappend  = @metric_distancematrix_recappend;
        %metric.matrix     = @metric_distancematrix_matrix;
        %metric.initmatrices = @metric_distancematrix_initmatrices;
    end
    end

    
    function w = metric_distancematrix_pak(metric)
    %METRIC_distancematrix_PAK	 Combine the metric parameters into one vector.
    %
    %	Description
    %	W = METRIC_distancematrix_PAK(METRIC, W) takes a metric data structure METRIC and
    %	combines the parameters into a single row vector W.
    %
    %	See also
    %	METRIC_distancematrix_UNPAK
    
    
%         mp=metric.p;
% 
%         i1=0;i2=1;
%         if ~isempty(w)
%             i1 = length(w);
%         end
%         i2=i1+length(metric.params);
%         i1=i1+1;
%         w(i1:i2)=metric.params;
%         i1=i2;
%         
%         % Parameters
%         if isfield(mp.params, 'p') && ~isempty(mp.params.p)
%             i1=i1+1;
%             w(i1)=mp.params.a.s;
%             if any(strcmp(fieldnames(mp.params.p),'nu'))
%                 i1=i1+1;
%                 w(i1)=mp.params.a.nu;
%             end
%         end

        if ~isempty(metric.p.lengthScales)
            w = log(metric.lengthScales);
            
            % Hyperparameters of lengthScale
            w = [w metric.p.lengthScales.fh.pak(metric.p.lengthScales)];
        end

    end
    
    function [metric, w] = metric_distancematrix_unpak(metric, w)
    %METRIC_distancematrix_UNPAK  Separate metric parameter vector into components.
    %
    %	Description
    %	[METRIC, W] = METRIC_distancematrix_UNPAK(METRIC, W) takes a metric data structure
    %   METRIC parameter vector W, and returns a metric data structure  identical to the
    %   input, except that the parameters has been set to the values in W. Deletes the values
    %   set to METRIC from W and returns the modified W. 
    %
    %	See also
    %	METRIC_distancematrix_PAK
    %
        
        if ~isempty(metric.p.lengthScales)
            i2=length(metric.lengthScales);
            i1=1;
            metric.lengthScales = exp(w(i1:i2));
            w = w(i2+1:end);
        
            % Hyperparameters of lengthScale
            [p, w] = metric.p.lengthScales.fh.unpak(metric.p.lengthScales, w);
            metric.p.lengthScales = p;
        end
        
%         mp=metric.p;
%         i1=0;i2=1;
%         i2=i1+length(metric.params);
%         i1=i1+1;
%         metric.params=w(i1:i2);
%         i1=i2;
%         % Parameters
%         if isfield(mp.params, 'p') && ~isempty(mp.params.p)
%             i1=i1+1;
%             metric.p.params.a.s=w(i1);
%             if any(strcmp(fieldnames(mp.params.p),'nu'))
%                 i1=i1+1;
%                 metric.p.params.a.nu=w(i1);
%             end
%         end        
%         w = w(i1+1:end);
    end

    function eprior = metric_distancematrix_lp(metric)
    %METRIC_DISTANCEMATRIX_LP  Evaluate the log prior of metric parameters
    %
    %  Description
    %    LP = METRIC_EUCLIDEAN_LP(METRIC) takes a metric structure
    %    METRIC and returns log(p(th)), where th collects the
    %    parameters.
    %
    %  See also
    %    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN_G, GP_E
    %
%         [n, m] = size(x);

        % Evaluate the prior contribution to the error. The parameters that
        % are sampled are from space W = log(w) where w is all the "real" samples.
        % On the other hand errors are evaluated in the W-space so we need take
        % into account also the  Jakobian of transformation W -> w = exp(W).
        % See Gelman et al., 2004, Bayesian data Analysis, second edition, p24.

        eprior = -metric.p.lengthScales.fh.lp(metric.lengthScales, metric.p.lengthScales) + sum(log(metric.lengthScales));
        
%         eprior = 0;
%         mp=metric.p;
% 
%         if isfield(mp.params, 'p') && ~isempty(mp.params.p)
%             eprior=eprior...
%                    +feval(mp.params.p.s.fe, ...
%                           gpp.params.a.s, mp.params.p.s.a)...
%                    -log(mp.params.a.s);
%             if any(strcmp(fieldnames(mp.params.p),'nu'))
%                 eprior=eprior...
%                        +feval(mp.p.params.nu.fe, ...
%                               mp.params.a.nu, mp.params.p.nu.a)...
%                        -log(mp.params.a.nu);
%             end
%         end
%         eprior=eprior...
%                +feval(mp.params.fe, ...
%                       metric.params, mp.params.a)...
%                -sum(log(metric.params));

    end
    
    function lpg = metric_distancematrix_lpg(metric) 
    %METRIC_DISTANCEMATRIX_LPG  d log(prior)/dth of the metric parameters th
    %
    %  Description
    %    LPG = METRIC_DISTANCEMATRIX_LPG(METRIC) takes a likelihood
    %    structure METRIC and returns d log(p(th))/dth, where th
    %    collects the parameters.
    %
    %  See also
    %    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN, GP_E
    %

        % Evaluate the prior contribution of gradient with respect to lengthScale
          if ~isempty(metric.p.lengthScales)
            i1=1; 
            lll = length(metric.lengthScales);
            lpgs = metric.p.lengthScales.fh.lpg(metric.lengthScales, metric.p.lengthScales);
            lpg(i1:i1-1+lll) = lpgs(1:lll).*metric.lengthScales + 1;
            lpg = [lpg lpgs(lll+1:end)];
          end
        end
   
    function [gdist, gprior]  = metric_distancematrix_ghyper(metric, x1, x2, mask) 
    %METRIC_distancematrix_GHYPER     Evaluate the gradient of the metric function and hyperprior 
    %                            w.r.t to it's hyperparameters.
    %
    %	Description
    %	[GDIST, GPRIOR_DIST] = METRIC_distancematrix_GHYPER(METRIC, X) takes a
    %   metric data structure METRIC together with a matrix X of input vectors and 
    %   return the gradient matrices GDIST and GPRIOR_DIST for each hyperparameter.
    %
    %	[GDIST, GPRIOR_DIST] = METRIC_distancematrix_GHYPER(METRIC, X, X2) forms the gradient
    %   matrices between two input vectors X and X2.
    %     
    %	[GDIST, GPRIOR_DIST] = METRIC_distancematrix_GHYPER(METRIC, X, X2, MASK) forms
    %   the gradients for masked covariances matrices used in sparse approximations.
    %
    %	See also
    %	METRIC_distancematrix_PAK, METRIC_distancematrix_UNPAK, METRIC_distancematrix, GP_E
    %

        mp=metric.p;
        components = metric.components;
        
        %n = size(x1,1);
        %m = length(components);
        m=1;
        i1=0;i2=1;
        %[n1,m1]=size(x1);
        
%         if nargin < 3
%             n2=n1;
%         end

        if ~isempty(metric.p.lengthScales)
        if nargin <= 3
            
            x1=round(x1);
            s = 1./metric.lengthScales;
%             if nargin < 3
%                 dist=s.*metric.x1x2matrix;
%             elseif n1==size(metric.x1,1)
%                 dist=s.*metric.x1matrix;
%                 
%             else
%                 dist=s.*metric.x2matrix;
%             end
            n1=max(size(x1));
            if nargin<3
            
%             if size(metric.Kff,1)==n1
%                 for i=1:m 
%                     gdist{i}=-s.*sqrt(squeeze(metric.Kff));
%                 end
%             elseif length(metric.X_u)==n1
%                 for i=1:m 
%                     gdist{i}=-s.*sqrt(squeeze(metric.Kff(x1matrix,x1matrix)));
%                 end
%             else
                 for i=1:m 
                     gdist{i}=-s.*(squeeze(metric.Kstarstar(x1,x1)));
                 end
%                 
%            end
            else
%                 n2=max(size(x2matrix));
%                 if size(metric.Kff,1)==n1 && length(metric.X_u)==n2
%                     for i=1:m 
%                     gdist{i}=-s.*sqrt(squeeze(metric.Kff(:,x2matrix)));
%                     end
%                 elseif size(metric.Kff,1)==n2 && length(metric.X_u)==n1
%                     for i=1:m 
%                     gdist{i}=-s.*sqrt(squeeze(metric.Kff(x1matrix,:)))';
%                     end
%                 elseif size(metric.Kstarstar,1)==n1 && length(metric.X_u)==n2
%                     for i=1:m 
%                         gdist{i}=-s.*sqrt(squeeze(metric.Kstarstar(:,x2matrix)));
%                     end
%                 elseif size(metric.Kstarstar,1)==n2 && length(metric.X_u)==n1
%                     for i=1:m
%                         gdist{i}=-s.*sqrt(squeeze(metric.Kstarstar(x1matrix,:)))';
%                     end
%                 elseif  size(metric.Kff,1)==n1 && size(metric.Kstarstar,1)==n2
%                     for i=1:m
%                         gdist{i}=-s.*sqrt(squeeze(metric.Kffstar));
%                     end
%                 elseif size(metric.Kff,1)==n2 && size(metric.Kstarstar,1)==n1
                     x2=round(x2);
                     for i=1:m
                         gdist{i}=-s.*(squeeze(metric.Kstarstar(x1,x2)));
                     end
%                     
%                 end
            end
% 
%             for i=1:m 
%                 gdist{i}=-s.*sqrt(squeeze(x1matrix));
%             end
%             
            
%             ii1=0;            
% 
%             dist  =  0;
%             distc = cell(1,m);
%             % Compute the distances for each component set
%             for i=1:m
%                 s = 1./metric.params(i).^2;
%                 distc{i} = 0;
%                 for j = 1:length(components{i})
%                     distc{i} = distc{i} + bsxfun(@minus,x(:,components{i}(j)),x2(:,components{i}(j))').^2;
%                 end
%                 distc{i} = distc{i}.*s;
%                 % Accumulate to the total distance
%                 dist = dist + distc{i};
%             end
%             dist = sqrt(dist);
%             % Loop through component sets 
%             for i=1:m
%                 D = -distc{i};
%                 D(dist~=0) = D(dist~=0)./dist(dist~=0);
%                 ii1 = ii1+1;
%                 gdist{ii1} = D;
%             end

% $$$         elseif nargin == 3
% $$$             if size(x,2) ~= size(x2,2)
% $$$                 error('metric_distancematrix -> _ghyper: The number of columns in x and x2 has to be the same. ')
% $$$             end
        elseif nargin == 4
            gdist = cell(1,length(metric.lengthScales));
        end
        
        % Evaluate the prior contribution of gradient with respect to lengthScale
%         for i2=1:m
%             i1=i1+1;
%             gprior(i1)=feval(mp.params.fg, ...
%                              metric.params(i2), ...
%                              mp.params.a, 'x').*metric.params(i2) - 1;
%         end
% 
%         % Evaluate the prior contribution of gradient with respect to lengthScale.p.s (and lengthScale.p.nu)

        if ~isempty(metric.p.lengthScales)
            i1=1; 
            lll = length(metric.lengthScales);
            gg = -metric.p.lengthScales.fh.lpg(metric.lengthScales, metric.p.lengthScales);
            gprior(i1:i1-1+lll) = gg(1:lll).*metric.lengthScales - 1;
            gprior = [gprior gg(lll+1:end)];
        end
        
%         if isfield(mp.params, 'p') && ~isempty(mp.params.p)
%             i1=i1+1;
%             gprior_dist(i1)=...
%                 feval(mp.params.p.s.fg, ...
%                       mp.params.a.s,...
%                       mp.params.p.s.a, 'x').*mp.params.a.s - 1 ...
%                 +feval(mp.params.fg, ...
%                        metric.params, ...
%                        mp.params.a, 's').*mp.params.a.s;
%             if any(strcmp(fieldnames(mp.params.p),'nu'))
%                 i1=i1+1;
%                 gprior_dist(i1)=...
%                     feval(mp.params.p.nu.fg, ...
%                           mp.params.a.nu,...
%                           mp.params.p.nu.a, 'x').*mp.params.a.nu -1 ...
%                     +feval(mp.params.fg, ...
%                            metric.params, ...
%                            mp.params.a, 'nu').*mp.params.a.nu;
%             end
%         end
        end
    end
    
    

    function dist  = metric_distancematrix_distance(metric, x1, x2)         
    %METRIC_distancematrix_DISTANCE   Compute the distancematrix distence between
    %                            one or two matrices.
    %
    %	Description
    %	[DIST] = METRIC_distancematrix_DISTANCE(METRIC, X) takes a metric data
    %   structure METRIC together with a matrix X of input vectors and 
    %   calculates the distancematrix distance matrix DIST.
    %
    %	[DIST] = METRIC_distancematrix_DISTANCE(METRIC, X1, X2) takes a metric data
    %   structure METRIC together with a matrices X1 and X2 of input vectors and 
    %   calculates the distancematrix distance matrix DIST.
    %
    %	See also
    %	METRIC_distancematrix_PAK, METRIC_distancematrix_UNPAK, METRIC_distancematrix, GP_E
    %
    
    
%         if nargin == 2 || isempty(n2)
%             n2=n1;
%         end
        
%         [n1,m1]=size(x1);
%         [n2,m2]=size(x2);
        
%         if m1~=m2
%             error('the number of columns of X1 and X2 has to be same')
%         end

  
%         [m1,n1]=size(x1);
%         m1=max(m1,n1);
%         dist=0;
        s = 1./metric.lengthScales;
        x1=round(x1);
        if (nargin == 3)
%             [m2,n2]=size(x2);
%             m2=max(m2,n2);
            x2=round(x2);
%             if m2~=n2
%                 if n2==1
%                    dist=s.*sqrt(squeeze(x1matrix(:,x2matrix)))'; 
%                 else
%                 dist=s.*sqrt(squeeze(x2matrix))';
%                 end
%             else
%                 dist=s.*sqrt(squeeze(x1matrix))';
%             end

            dist=s.*(squeeze(metric.Kstarstar(x1,x2)));
%              
%              if size(metric.,1)==m1 && length(metric.X_u)==m2
%                  dist=s.*sqrt(squeeze(metric.Kff(:,x2matrix)));
%              elseif size(metric.Kff,1)==m2 && length(metric.X_u)==m1
%                  dist=s.*sqrt(squeeze(metric.Kff(x1matrix,:)));
%              elseif size(metric.Kff,1)==m1 && size(metric.Kstarstar,1)==m2
%                  dist=s.*sqrt(squeeze(metric.Kffstar));
%              elseif size(metric.Kff,1)==m2 && size(metric.Kstarstar,1)==m1
%                  dist=s.*sqrt(squeeze(metric.Kffstar))';
%              elseif size(metric.Kstarstar,1)==m1 && length(metric.X_u)==m2
%                  dist=s.*sqrt(squeeze(metric.Kstarstar(:,x2matrix)));
%              elseif size(metric.Kstarstar,1)==m2 && length(metric.X_u)==m1
%                  dist=s.*sqrt(squeeze(metric.Kstarstar(x1matrix,:)))';
%              else
%                  dist=s.*sqrt(squeeze(x1matrix));
%              end
%          elseif n1~=metric.n1
%              dist=s.*squeeze(metric.x2matrix);
        else
%             if length(metric.X_u)==m1
%                 dist=s.*sqrt(squeeze(metric.Kff(x1matrix,x1matrix)));
%             elseif size(metric.Kff,1)==m1
%                 dist=s.*sqrt(squeeze(metric.Kff));
%             else
                dist=s.*(squeeze(metric.Kstarstar(x1,x1)));
%             end
        end

%         components = metric.components;
%         m = length(components);
%         dist  =  0;        
%         
%         if isempty(metric.dmatrix)
%             for i=1:m
%                 s = 1./metric.params(i).^2;
%                 for j = 1:length(components{i})
%                     dist = dist + s.*bsxfun(@minus,x1(:,components{i}(j)),x2(:,components{i}(j))').^2;
%                 end
%             end
%         dist = sqrt(dist);
%         end
        
        %dist=metric.dmatrix;
        

    end
    
    function [ginput, gprior_input]  = metric_distancematrix_ginput(metric, x1, x2)         
    %METRIC_distancematrix_GINPUT   Compute the gradient of the distancematrix distance
    %                          function with respect to input.
    %[n, m] =size(x);
        ii1 = 0;
        components = metric.components;
%                
%         if nargin == 2 || isempty(x2)
%             x2=x1;
%         end
%         
%         [n1,m1]=size(x1);
%         [n2,m2]=size(x2);
%         
%         if m1~=m2
%             error('the number of columns of X1 and X2 has to be same')
%         end

%         if n1~=n2
%             dist=metric.x1x2matrix;
%         elseif n1==length(metric.x1)
%             dist=metric.x1matrix;
%         else
%             dist=metric.x2matrix;
%         end
 
        if nargin == 2 || isempty(x2)
            x2=x1;
            dist=metric.fh.dist(metric,x1);          
        else
            dist=metric.fh.dist(metric,x1,x2);
        end
        
        
        
%         s = 1./metric.params.^2;
%         dist = 0;
%         for i=1:length(components)
%             for j = 1:length(components{i})
%                 dist = dist + s(i).*bsxfun(@minus,x1(:,components{i}(j)),x2(:,components{i}(j))').^2;
%             end
%         end
%         dist = sqrt(dist);
        x1=round(x1);
        x2=round(x2);
        s=1./metric.lengthScales.^2;
        n1=length(x1);
        n2=length(x2);
        if n1==n2
%             if n1==length(metric.X_u)
%                 xx1=metric.Kff(x1,x2);
%             elseif n1==size(metric.Kff,1)
%                 xx1=metric.Kff;
%             else
                xx1=metric.Kstarstar(x1,x2);
%            end
        else
%             if n1==length(metric.X_u)
%                 if n2==size(metric.Kff,1)
%                     xx1=metric.Kff(x1,:);
%                 elseif n2==size(metric.Kstarstar,1)
%                     xx1=metric.Kstarstar(x1,:);
%                 end
%             elseif n1==size(metric.Kff,1)
%                 if n2==length(metric.X_u)
%                     xx1=metric.Kff(:,x2);
%                 elseif n2==size(metric.Kstarstar,1)
%                     xx1=metric.Kffstar;
%                 end
%             else
%                if n2==length(metric.X_u)
%                     xx1=metric.Kfstaru;
%                elseif n2==size(metric.Kff,1)
                     xx1=metric.Kstarstar(x1,x2);
%                end
%             end
        end

            for j = 1:n1
                DK = zeros(n1,n2);                
                for k = 1:length(components)
                    DK(j,:) = DK(j,:)+s(k).*xx1(j,:);
                end
                if nargin == 2
                    DK = DK + DK';
                end
                DK(dist~=0) = DK(dist~=0)./dist(dist~=0);
                DK=DK./2;                        
                ii1 = ii1 + 1;
                ginput{ii1} = DK;
                gprior_input(ii1) = 0; 
            end
        
    end
    
    
    function recmetric = metric_distancematrix_recappend(recmetric, ri, metric)
    % RECAPPEND - Record append
    %          Description
    %          RECMETRIC = METRIC_distancematrix_RECAPPEND(RECMETRIC, RI, METRIC) takes old covariance
    %          function record RECMETRIC, record index RI and covariance function structure. 
    %          Appends the parameters of METRIC to the RECMETRIC in the ri'th place.
    %
    %          RECAPPEND returns a structure RECMETRIC containing following record fields:
    %          lengthHyper    
    %          lengthHyperNu  
    %          lengthScale    
    %
    %          See also
    %          GP_MC and GP_MC -> RECAPPEND

    % Initialize record
        if nargin == 2
            recmetric.type = 'metric_distancematrix';
%             recmetric.nin = ri;
            metric.components = recmetric.components;
            
            % Initialize parameters
            recmetric.lengthScales = [];

            % Set the function handles
            recmetric.fh.pak       = @metric_distancematrix_pak;
            recmetric.fh.unpak     = @metric_distancematrix_unpak;
            recmetric.fh.lp        = @metric_distancematrix_lp;
            recmetric.fh.lpg       = @metric_distancematrix_lpg;
            recmetric.fh.dist      = @metric_distancematrix_distance;
            recmetric.fh.distg     = @metric_distancematrix_ghyper;
            recmetric.fh.ginput    = @metric_distancematrix_ginput;            
            recmetric.fh.recappend = @metric_distancematrix_recappend;
            %recmetric.matrix    = @metric_distancematrix_matrix;
            return
        end
        mp = metric.p;

        % record parameters
         if ~isempty(metric.lengthScales)
            recmetric.lengthScales(ri,:)=metric.lengthScales;
            recmetric.p.lengthScales = metric.p.lengthScales.fh.recappend(recmetric.p.lengthScales, ri, metric.p.lengthScales);
        elseif ri==1
            recmetric.lengthScales=[];
        end
        
%         if ~isempty(metric.params)
%             if ~isempty(mp.params)
%                 recmetric.lengthHyper(ri,:)=mp.params.a.s;
%                 if isfield(mp.params,'p')
%                     if isfield(mp.params.p,'nu')
%                         recmetric.lengthHyperNu(ri,:)=mp.params.a.nu;
%                     end
%                 end
%             elseif ri==1
%                 recmetric.lengthHyper=[];
%             end
%             recmetric.params(ri,:)=metric.params;
%         elseif ri==1
%             recmetric.params=[];
%         end
    end
