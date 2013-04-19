function metric = metric_euclidean(varargin)
%METRIC_EUCLIDEAN An euclidean metric function
%
%  Description
%    METRIC = METRIC_EUCLIDEAN('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a an euclidean metric function structure in which the
%    named parameters have the specified values. Either
%    'components' or 'deltadist' has to be specified. Any
%    unspecified parameters are set to default values.
%   
%    METRIC = METRIC_EUCLIDEAN(METRIC,'PARAM1',VALUE1,'PARAM2,VALUE2,...)
%    modify a metric function structure with the named parameters
%    altered with the specified values.
%
%    Parameters for Euclidean metric function [default]
%       components        - cell array of vectors specifying which 
%                           inputs are grouped together with a same
%                           scaling parameter. For example, the
%                           component specification {[1 2] [3]}
%                           means that distance between 3
%                           dimensional vectors computed as 
%                           r = (r_1^2 + r_2^2 )/l_1 + r_3^2/l_2,
%                           where r_i are distance along component
%                           i, and l_1 and l_2 are lengthscales for
%                           corresponding component sets. If
%                           'components' is not specified, but
%                           'deltadist' is specified, then default
%                           is {1 ... length(deltadist)}
%       deltadist         - indicator vector telling which component sets
%                           are handled using the delta distance 
%                           (0 if x=x', and 1 otherwise). Default is
%                           false for all component sets.
%       lengthScale       - lengthscales for each input component set
%                           Default is 1 for each set
%       lengthScale_prior - prior for lengthScales [prior_unif]
%
%  See also
%    GP_SET, GPCF_SEXP
  
% Copyright (c) 2008 Jouni Hartikainen 
% Copyright (c) 2008 Jarno Vanhatalo     
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'METRIC_EUCLIDEAN';
  ip.addOptional('metric', [], @isstruct);
  ip.addParamValue('components',[], @(x) isempty(x) || iscell(x));
  ip.addParamValue('deltadist',[], @(x) isvector(x));
  ip.addParamValue('lengthScale',[], @(x) isvector(x) && all(x>0));
  ip.addParamValue('lengthScale_prior',prior_unif, ...
                   @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  metric=ip.Results.metric;

  if isempty(metric)
    % Initialize a Gaussian process
    init=true;
  else
    % Modify a Gaussian process
    if ~isfield(metric,'type') && isequal(metric.type,'metric_euclidean')
      error('First argument does not seem to be a metric structure')
    end
    init=false;
  end

  if init
    % Type
    metric.type = 'metric_euclidean';
  end
  
  % Components
  if init || ~ismember('components',ip.UsingDefaults)
    metric.components = ip.Results.components;
  end
  % Deltadist
  if init || ~ismember('deltadist',ip.UsingDefaults)
    metric.deltadist = ip.Results.deltadist;
  end
  % Components+Deltadist check and defaults
  if isempty(metric.components) && isempty(metric.deltadist)
    error('Either ''components'' or ''deltadist'' has to be specified')
  elseif isempty(metric.components)
    metric.components=num2cell(1:length(metric.components));
  elseif isempty(metric.deltadist)
    metric.deltadist = false(1,length(metric.components));
  end
  % Lengthscale
  if init || ~ismember('lengthScale',ip.UsingDefaults)
    metric.lengthScale = ip.Results.lengthScale;
    if isempty(metric.lengthScale)
      metric.lengthScale = repmat(1,1,length(metric.components));
    end
  end
  % Prior for lengthscale
  if init || ~ismember('lengthScale_prior',ip.UsingDefaults)
    metric.p=[];
    metric.p.lengthScale = ip.Results.lengthScale_prior;
  end
  
  if init
    % Set the function handles to the subfunctions
    metric.fh.pak       = @metric_euclidean_pak;
    metric.fh.unpak     = @metric_euclidean_unpak;
    metric.fh.lp        = @metric_euclidean_lp;
    metric.fh.lpg       = @metric_euclidean_lpg;
    metric.fh.dist      = @metric_euclidean_dist;
    metric.fh.distg     = @metric_euclidean_distg;
    metric.fh.ginput    = @metric_euclidean_ginput;
    metric.fh.recappend = @metric_euclidean_recappend;
  end

end

function [w s] = metric_euclidean_pak(metric)
%METRIC_EUCLIDEAN_PAK  Combine GP covariance function
%                      parameters into one vector.
%
%  Description
%    W = METRIC_EUCLIDEAN_PAK(GPCF) takes a covariance function
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W and takes a logarithm of the covariance function
%    parameters.
%
%       w = [ log(metric.lengthScale(:))
%             (hyperparameters of metric.lengthScale)]'
%       
%  See also
%    METRIC_EUCLIDEAN_UNPAK
  
  w = []; s = {};
  if ~isempty(metric.p.lengthScale)
    w = log(metric.lengthScale);
    if numel(metric.lengthScale)>1
      s = [s; sprintf('log(metric.lengthScale x %d)',numel(metric.lengthScale))];
    else
      s = [s; 'log(metric.lengthScale)'];
    end
    % Hyperparameters of lengthScale
    [wh sh] = metric.p.lengthScale.fh.pak(metric.p.lengthScale);
    w = [w wh];
    s = [s; sh];
  end
  
end

function [metric, w] = metric_euclidean_unpak(metric, w)
%METRIC_EUCLIDEAN_UNPAK  Separate metric parameter vector into components
%
%  Description
%    METRIC, W] = METRIC_EUCLIDEAN_UNPAK(METRIC, W) takes a
%    metric structure GPCF and a parameter vector W, and returns
%    a covariance function structure identical to the input,
%    except that the covariance parameters have been set to the
%    values in W. Deletes the values set to GPCF from W and
%    returns the modified W.
%
%    The covariance function parameters are transformed via exp
%    before setting them into the structure.
%
%  See also
%    METRIC_EUCLIDEAN_PAK
%
  
  if ~isempty(metric.p.lengthScale)
    i2=length(metric.lengthScale);
    i1=1;
    metric.lengthScale = exp(w(i1:i2));
    w = w(i2+1:end);
    
    % Hyperparameters of lengthScale
    [p, w] = metric.p.lengthScale.fh.unpak(metric.p.lengthScale, w);
    metric.p.lengthScale = p;
  end
end

function lp = metric_euclidean_lp(metric)
%METRIC_EUCLIDEAN_LP  Evaluate the log prior of metric parameters
%
%  Description
%    LP = METRIC_EUCLIDEAN_LP(METRIC) takes a metric structure
%    METRIC and returns log(p(th)), where th collects the
%    parameters.
%
%  See also
%    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN_G, GP_E
%
  
% Evaluate the prior contribution to the error. The parameters that
% are sampled are from space W = log(w) where w is all the "real" samples.
% On the other hand errors are evaluated in the W-space so we need take
% into account also the  Jakobian of transformation W -> w = exp(W).
% See Gelman et.all., 2004, Bayesian data Analysis, second edition, p24.
  if ~isempty(metric.p.lengthScale)
    lp = metric.p.lengthScale.fh.lp(metric.lengthScale, metric.p.lengthScale) + sum(log(metric.lengthScale));
  else
    lp=0;
  end
  
end

function lpg = metric_euclidean_lpg(metric) 
%METRIC_EUCLIDEAN_LPG  d log(prior)/dth of the metric parameters th
%
%  Description
%    LPG = METRIC_EUCLIDEAN_LPG(METRIC) takes a likelihood
%    structure METRIC and returns d log(p(th))/dth, where th
%    collects the parameters.
%
%  See also
%    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN, GP_E
%

% Evaluate the prior contribution of gradient with respect to lengthScale
  if ~isempty(metric.p.lengthScale)
    i1=1; 
    lll = length(metric.lengthScale);
    lpgs = metric.p.lengthScale.fh.lpg(metric.lengthScale, metric.p.lengthScale);
    lpg(i1:i1-1+lll) = lpgs(1:lll).*metric.lengthScale + 1;
    lpg = [lpg lpgs(lll+1:end)];
  end
end

function gdist  = metric_euclidean_distg(metric, x, x2, mask) 
%METRIC_EUCLIDEAN_DISTG  Evaluate the gradient of the metric function
%
%  Description
%    DISTG = METRIC_EUCLIDEAN_DISTG(METRIC, X) takes a metric
%    structure METRIC together with a matrix X of input
%    vectors and return the gradient matrices GDIST and
%    GPRIOR_DIST for each parameter.
%
%    DISTG = METRIC_EUCLIDEAN_DISTG(METRIC, X, X2) forms the
%    gradient matrices between two input vectors X and X2.
%     
%    DISTG = METRIC_EUCLIDEAN_DISTG(METRIC, X, X2, MASK) forms
%    the gradients for masked covariances matrices used in sparse
%    approximations.
%
%  See also
%    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN, GP_E
%

  gdist=[];
  components = metric.components;
  
  n = size(x,1);
  m = length(components);
  i1=0;i2=1;

  % NOTE! Here we have already taken into account that the parameters
  % are transformed through log() and thus dK/dlog(p) = p * dK/dp
  
  if ~isempty(metric.p.lengthScale)
    if nargin <= 3
      if nargin == 2
        x2 = x;
      end
      ii1=0;            

      dist  =  0;
      distc = cell(1,m);
      % Compute the distances for each component set
      for i=1:m
        if length(metric.lengthScale)==1
          s=1./metric.lengthScale.^2;
        else
          s=1./metric.lengthScale(i).^2;
        end
        distc{i} = 0;
        for j = 1:length(components{i})
          if metric.deltadist(i)
            distc{i} = distc{i} + double(bsxfun(@ne,x(:,components{i}(j)),x2(:,components{i}(j))'));
          else
            distc{i} = distc{i} + bsxfun(@minus,x(:,components{i}(j)),x2(:,components{i}(j))').^2;
          end
        end
        distc{i} = distc{i}.*s;
        % Accumulate to the total distance
        dist = dist + distc{i};
      end
      dist = sqrt(dist);
      % Loop through component sets
      if length(metric.lengthScale)==1
        D = -distc{1};
        D(dist~=0) = D(dist~=0)./dist(dist~=0);
        ii1 = ii1+1;
        gdist{ii1} = D;
      else
        for i=1:m
          D = -distc{i};
          ind = dist~=0;
          D(ind) = D(ind)./dist(ind);
          ii1 = ii1+1;
          gdist{ii1} = D;
        end
      end
% $$$         elseif nargin == 3
% $$$             if size(x,2) ~= size(x2,2)
% $$$                 error('metric_euclidean -> _ghyper: The number of columns in x and x2 has to be the same. ')
% $$$             end
    elseif nargin == 4
      gdist = cell(1,length(metric.lengthScale));
    end

    % Evaluate the prior contribution of gradient with respect to lengthScale
    if ~isempty(metric.p.lengthScale)
      i1=1; 
      lll = length(metric.lengthScale);
      gg = -metric.p.lengthScale.fh.lpg(metric.lengthScale, metric.p.lengthScale);
      gprior(i1:i1-1+lll) = gg(1:lll).*metric.lengthScale - 1;
      gprior = [gprior gg(lll+1:end)];
    end
  end
end

function dist = metric_euclidean_dist(metric, x1, x2)         
%METRIC_EUCLIDEAN_DIST  Compute the euclidean distence between
%                       one or two matrices.
%
%  Description
%    DIST = METRIC_EUCLIDEAN_DIST(METRIC, X) takes a metric
%    structure METRIC together with a matrix X of input
%    vectors and calculates the euclidean distance matrix DIST.
%
%    DIST = METRIC_EUCLIDEAN_DIST(METRIC, X1, X2) takes a
%    metric structure METRIC together with a matrices X1 and
%    X2 of input vectors and calculates the euclidean distance
%    matrix DIST.
%
%  See also
%    METRIC_EUCLIDEAN_PAK, METRIC_EUCLIDEAN_UNPAK, METRIC_EUCLIDEAN, GP_E
%
  if (nargin == 2 || isempty(x2))
    % use fast c-code for self-distance
    x2=x1;
    % force deltadist to be logical for simplified c-code
    metric.deltadist=logical(metric.deltadist);
    dist = dist_euclidean(metric,x1);
    if ~any(isnan(dist))
      % if c-code was available, result is not NaN
      return
    end        
  end
    
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);
  
  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end
  
  components = metric.components;
  m = length(components);
  dist  =  0;        
  
  s=1./metric.lengthScale.^2;
  if m>numel(s)
    s=repmat(s,1,m);
  end
  for i=1:m
    for j = 1:length(components{i})
      if metric.deltadist(i)
        dist = dist + s(i).*double(bsxfun(@ne,x1(:,components{i}(j)),x2(:,components{i}(j))'));
      else
        dist = dist + s(i).*bsxfun(@minus,x1(:,components{i}(j)),x2(:,components{i}(j))').^2;
      end
    end
  end
  dist=sqrt(dist); % euclidean distance
  
end

function [ginput, gprior_input]  = metric_euclidean_ginput(metric, x1, x2)
%METRIC_EUCLIDEAN_GINPUT  Compute the gradient of the
%                         euclidean distance function with
%                         respect to input. [n, m]=size(x);
  ii1 = 0;
  components = metric.components;
  
  if nargin == 2 || isempty(x2)
    x2=x1;
  end
  
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);
  
  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end
  
  s = 1./metric.lengthScale.^2;
  dist = 0;
  for i=1:length(components)
    for j = 1:length(components{i})
      if metric.deltadist(i)
        dist = dist + s(i).*double(bsxfun(@ne,x1(:,components{i}(j)),x2(:,components{i}(j))'));
      else
        dist = dist + s(i).*bsxfun(@minus,x1(:,components{i}(j)),x2(:,components{i}(j))').^2;
      end
    end
  end
  dist = sqrt(dist);
  
  for i=1:m1
    for j = 1:n1
      DK = zeros(n1,n2);                
      for k = 1:length(components)
        if ismember(i,components{k})
          if metric.deltadist(i)
            DK(j,:) = DK(j,:)+s(k).*double(bsxfun(@ne,x1(j,i),x2(:,i)'));
          else
            DK(j,:) = DK(j,:)+s(k).*bsxfun(@minus,x1(j,i),x2(:,i)');
          end
        end
      end
      if nargin == 2
        DK = DK + DK';
      end
      DK(dist~=0) = DK(dist~=0)./dist(dist~=0);
      
      ii1 = ii1 + 1;
      ginput{ii1} = DK;
      gprior_input(ii1) = 0; 
    end
  end
  %size(ginput)
  %ginput
  
end


function recmetric = metric_euclidean_recappend(recmetric, ri, metric)
%RECAPPEND  Record append
%
%  Description
%    RECMETRIC = METRIC_EUCLIDEAN_RECAPPEND(RECMETRIC, RI,
%    METRIC) takes old metric function record RECMETRIC, record
%    index RI and metric function structure. Appends the
%    parameters of METRIC to the RECMETRIC in the ri'th place.
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND

% Initialize record
  if nargin == 2
    recmetric.type = 'metric_euclidean';
    metric.components = recmetric.components;
    
    % Initialize parameters
    recmetric.lengthScale = [];

    % Set the function handles
    recmetric.fh.pak       = @metric_euclidean_pak;
    recmetric.fh.unpak     = @metric_euclidean_unpak;
    recmetric.fh.lp        = @metric_euclidean_lp;
    recmetric.fh.lpg       = @metric_euclidean_lpg;
    recmetric.fh.dist      = @metric_euclidean_dist;
    recmetric.fh.distg     = @metric_euclidean_distg;
    recmetric.fh.ginput    = @metric_euclidean_ginput;            
    recmetric.fh.recappend = @metric_euclidean_recappend;
    return
  end
  mp = metric.p;

  % record parameters
  if ~isempty(metric.lengthScale)
    recmetric.lengthScale(ri,:)=metric.lengthScale;
    recmetric.p.lengthScale = metric.p.lengthScale.fh.recappend(recmetric.p.lengthScale, ri, metric.p.lengthScale);
  elseif ri==1
    recmetric.lengthScale=[];
  end

end
