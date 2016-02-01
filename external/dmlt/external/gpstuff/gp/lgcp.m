function [l,lq,xt,gp] = lgcp(x,varargin)
% LGCP - Log Gaussian Cox Process intensity estimate for 1D and 2D data
%   
%   LGCP(X)
%   [P,PQ,XT,GP] = LGCP(X,XT,OPTIONS)
% 
%     X  is 1D or 2D point data
%     XT is optional test points
%     OPTIONS are optional parameter-value pairs
%       'gridn' is optional number of grid points used in each axis direction
%         default is 100 for 1D, 15 for grid 2D
%       'range' tells the estimation range, default is data range
%         for 1D [XMIN XMAX]
%         for 2D [XMIN XMAX YMIN YMAX]
%       'gpcf' is optional function handle of a GPstuff covariance function 
%         (default is @gpcf_sexp)
%       'latent_method' is optional 'EP' (default) or 'Laplace'
%       'int_method' is optional 'mode' (default), 'CCD' or 'grid'
% 
%     P  is the estimated intensity  
%     PQ is the 5% and 95% percentiles of the intensity estimate
%     XT contains the used test points
%     GP is the Gaussian process formed. As the grid is scaled to
%        unit range or unit square, additional field 'scale' is
%        included which includes the range for the grid in the
%        original x space.
  
% Copyright (c) 2009-2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LGCP';
  ip.addRequired('x', @(x) isnumeric(x) && size(x,2)==1 || size(x,2)==2);
  ip.addOptional('xt',NaN, @(x) isnumeric(x) && size(x,2)==1 || size(x,2)==2);
  ip.addParamValue('gridn',[], @(x) isscalar(x) && x>0 && mod(x,1)==0);
  ip.addParamValue('range',[], @(x) isreal(x)&&(length(x)==2||length(x)==4));
  ip.addParamValue('gpcf',@gpcf_sexp,@(x) ischar(x) || isa(x,'function_handle'));
  ip.addParamValue('latent_method','EP', @(x) ismember(x,{'EP' 'Laplace'}))
  ip.addParamValue('int_method','mode', @(x) ismember(x,{'mode' 'CCD', 'grid'}))
  ip.addParamValue('normalize',false, @islogical);
  
  ip.parse(x,varargin{:});
  x=ip.Results.x;
  xt=ip.Results.xt;
  gridn=ip.Results.gridn;
  xrange=ip.Results.range;
  gpcf=ip.Results.gpcf;
  latent_method=ip.Results.latent_method;
  int_method=ip.Results.int_method;
  normalize=ip.Results.normalize;
  
  [n,m]=size(x);
  
  switch m
    case 1 % 1D
      % Parameters for a grid
      if isempty(gridn)
        % number of points
        gridn=100;
      end
      xmin=min(x);xmax=max(x);
      if ~isempty(xrange)
        xmin=min(xmin,xrange(1));
        xmax=max(xmax,xrange(2));
      end
      % Discretize the data
      xx=linspace(xmin,xmax,gridn)';
      yy=hist(x,xx)';
      ye=ones(gridn,1)./gridn.*n;
      % Test points
      if isnan(xt)
        xt=linspace(xmin,xmax,max(gridn,200))';
      end
      % normalise to unit range, so that same prior is ok for different scales
      xxn=(xx-min(xx))./range(xx)-0.5;
      xtn=(xt-min(xx))./range(xx)-0.5;
      % smooth...
      [Ef,Varf,gp]=gpsmooth(xxn,yy,ye,xtn,gpcf,latent_method,int_method);
      gp.scale=range(xx);
      
      % compute mean and quantiles
      A=range(xx);
      lm=exp(Ef+Varf/2)./A.*n;
      lq5=exp(Ef-sqrt(Varf)*1.645)./A*n;
      lq95=exp(Ef+sqrt(Varf)*1.645)./A*n;
      lq=[lq5 lq95];

      if nargout<1
        % no output, do the plot thing
        newplot
        hp=patch([xt; xt(end:-1:1)],[lq(:,1); lq(end:-1:1,2)],[.9 .9 .9]);
        set(hp,'edgecolor',[.9 .9 .9])
        xlim([xmin xmax])
        line(xt,lm,'linewidth',2);
      else
        l=lm;
      end
      
    case 2 % 2D
      
      % Find unique points
      [xu,I,J]=unique(x,'rows');
      % and count number of repeated x's
      counts=crosstab(J); 
      nu=length(xu);
  
      % Parameters for a grid
      if isempty(gridn)
        % number of points in direction
        gridn=15;
      end
      x1min=min(x(:,1));x1max=max(x(:,1));
      x2min=min(x(:,2));x2max=max(x(:,2));
      if ~isempty(xrange)
        % range extension
        x1min=min(x1min,xrange(1));
        x1max=max(x1max,xrange(2));
        x2min=min(x2min,xrange(3));
        x2max=max(x2max,xrange(4));
      end
      % Form regular grid to discretize the data
      zz1=linspace(x1min,x1max,gridn)';
      zz2=linspace(x2min,x2max,gridn)';
      [z1,z2]=meshgrid(zz1,zz2);
      z=[z1(:),z2(:)];
      nz=length(z);
      % form data for GP (xx,yy,ye)
      xx=z;
      yy=zeros(nz,1);
      zi=interp2(z1,z2,reshape(1:nz,gridn,gridn),xu(:,1),xu(:,2),'nearest');
      for i1=1:nu
        yy(zi(i1),1)=yy(zi(i1),1)+counts(i1);
      end
      ye=ones(nz,1)./nz.*n;

      % Default test points
      if isnan(xt)
        [xt1,xt2]=meshgrid(linspace(x1min,x1max,max(100,gridn)),...
                           linspace(x2min,x2max,max(100,gridn)));
        xt=[xt1(:) xt2(:)];
      end
      % normalise to unit square, so that same prior is ok for different scales
      xxn=bsxfun(@rdivide,bsxfun(@minus,xx,min(xx,[],1)),range(xx,1))-.5;
      xtn=bsxfun(@rdivide,bsxfun(@minus,xt,min(xx,[],1)),range(xx,1))-.5;
      % smooth...
      [Ef,Varf,gp]=gpsmooth(xxn,yy,ye,xtn,gpcf,latent_method,int_method);
      gp.scale=[range(xx(:,1)) range(xx(:,2))];
      
      % compute mean
      A = range(xx(:,1)).*range(xx(:,2));
      lm=exp(Ef+Varf/2)./A.*n;
      lq5=exp(Ef-sqrt(Varf)*1.645)./A.*n;
      lq95=exp(Ef+sqrt(Varf)*1.645)./A.*n;
      lq=[lq5 lq95];

      if nargout<1
        % no output, do the plot thing
        G=zeros(size(xt1));
        G(:)=lm;
        pcolor(xt1,xt2,G);
        shading flat
        colormap('jet')
        cx=caxis;
        cx(1)=0;
        caxis(cx);
        colorbar
      else
        l=lm;
      end
      
    otherwise
      error('X has to be Nx1 or Nx2')
  end

end

function [Ef,Varf,gp] = gpsmooth(xx,yy,ye,xt,gpcf,latent_method,int_method)
% Make inference with log Gaussian process and EP or Laplace approximation

  nin = size(xx,2);
  % init gp
  if strfind(func2str(gpcf),'ppcs')
    % ppcs still have nin parameter...
    gpcf1 = gpcf('nin',nin);
  else
    gpcf1 = gpcf();
  end
  % default vague prior
  pm = prior_sqrtt('s2', 1^2, 'nu', 4);
  pl = prior_t('s2', 2^2, 'nu', 4);
  %pm = prior_logunif();
  %pl = prior_logunif();
  pa = prior_t('s2', 10^2, 'nu', 4);
  % different covariance functions have different parameters
  if isfield(gpcf1,'magnSigma2')
     gpcf1 = gpcf(gpcf1, 'magnSigma2', .1, 'magnSigma2_prior', pm);
  end
  if isfield(gpcf1,'lengthScale')
     gpcf1 = gpcf(gpcf1, 'lengthScale', .1, 'lengthScale_prior', pl);
  end
  if isfield(gpcf1,'alpha')
    gpcf1 = gpcf(gpcf1, 'alpha', 20, 'alpha_prior', pa);
  end
  if isfield(gpcf1,'weightSigma2')
    gpcf1 = gpcf(gpcf1, 'weightSigma2_prior', prior_logunif(), 'biasSigma2_prior', prior_logunif());
  end
  
  % Create the GP structure
  gp = gp_set('lik', lik_poisson(), 'cf', {gpcf1}, 'jitterSigma2', 1e-4);

  % Set the approximate inference method
  gp = gp_set(gp, 'latent_method', latent_method);
  
  % Optimize hyperparameters
  opt=optimset('TolX', 1e-3, 'Display', 'off');
  if exist('fminunc')
    gp = gp_optim(gp, xx, yy, 'z', ye, 'optimf', @fminunc, 'opt', opt);
  else
    gp = gp_optim(gp, xx, yy, 'z', ye, 'optimf', @fminlbfgs, 'opt', opt);
  end

  % Make prediction for the test points
  if strcmpi(int_method,'mode')
    % point estimate for the hyperparameters
    [Ef,Varf] = gp_pred(gp, xx, yy, xt, 'z', ye);
  else
    % integrate over the hyperparameters
    %[~, ~, ~, Ef, Varf] = gp_ia(opt, gp, xx, yy, xt, param);
    [notused, notused, notused, Ef, Varf]=...
        gp_ia(gp, xx, yy, xt, 'z', ye, 'int_method', int_method);
  end
  
end
