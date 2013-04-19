function [Ef, Varf, xtnn] = gp_cpred(gp,x,y,xt, ind,varargin)
%GP_CPRED Conditional predictions using specific covariates
%
%  Description
%    GP_CPRED(GP,X,Y,XT,IND,OPTIONS) does predictions using only
%    covariates specified in vector IND. Other covariates are fixed to
%    either mean, median or values chosen by user. Returns predictions for
%    latent values, variance and corresponding inputs. If IND=0, only time
%    is used as a covariate for coxph model.
%
%   OPTIONS is optional parameter-value pair
%      method - which value to fix the not used covariates, 'mean'
%               (default) or 'median'
%      var    - vector specifying optional values for not used covariates,
%               elements corresponding to mean/median values should 
%               be set to NaN. 
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%      plot   - Option for plotting, 'off' (default) or 'on'
%      tr     - Euclidean distance treshold for not using grid points when
%               doing predictions with 2 covariates, default 0.25


ip=inputParser;
ip.FunctionName = 'GP_CPRED';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('xt',  @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('ind', @(x) ~isempty(x) && isvector(x))
ip.addParamValue('var',  [], @(x) isreal(x))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                 isvector(x) && isreal(x) && all(isfinite(x)&x>0))
ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                 (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
ip.addParamValue('method', 'mean', @(x)  ismember(x, {'median', 'mean'}))
ip.addParamValue('plot', 'off', @(x)  ismember(x, {'on', 'off'}))
ip.addParamValue('tr', 0.25, @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(gp, x, y, xt, ind, varargin{:});
options=struct();
options.predcf=ip.Results.predcf;
options.tstind=ip.Results.tstind;
method = ip.Results.method;
vars = ip.Results.var;
plot_results = ip.Results.plot;
tr = ip.Results.tr;
z=ip.Results.z;
if ~isempty(z)
  options.zt=z;
  options.z=z;
end

[tmp, nin] = size(x);

if ~isempty(vars) && (~isvector(vars) || length(vars) ~= nin)
  error('Vector defining fixed variable values must be same length as number of covariates')
end

xto=xt; [n,tmp]=size(xto);

if length(ind)==1
  
  if ind~=0
    [xtnn, iu] = unique(xt(:,ind));
  else
    xtnn = xt(1,:);
    iu = 1;
  end
  if ~isempty(z)
    options.zt = options.zt(iu);
  end
  meanxt=mean(xt);
  if isequal(method, 'mean')
    xt = repmat(meanxt, size(xtnn,1), 1);
  else
    xt = repmat(median(xt), size(xtnn,1), 1);
  end
  if ~isempty(vars)
    xt(:,~isnan(vars)) = repmat(vars(~isnan(vars)), length(xtnn), 1);
  end
  
  if ind>0
    xt(:,ind) = xtnn;
  end
  if ~strcmp(gp.lik.type, 'Coxph')
    [Ef, Varf] = gp_pred(gp, x, y, xt, options);
  else
    [Ef1,Ef2,Covf] = pred_coxph(gp,x,y,xt, options);
    if ind>0
      Ef = Ef2; Varf = diag(Covf(size(Ef1,1)+1:end,size(Ef1,1)+1:end));
    else
      Ef = Ef1; Varf = diag(Covf(1:size(Ef1,1), 1:size(Ef1,1)));
      xtnn = gp.lik.xtime;
    end
  end
  if isequal(plot_results, 'on')
    if ind>0
      plot(xtnn, Ef, 'ob', xtnn, Ef, '-k', xtnn, Ef-1.96*sqrt(Varf), '--b', xtnn, Ef+1.96*sqrt(Varf), '--b')
    else
      % use stairs for piecewise constant baseline hazard
      xtnn = gp.lik.stime;
      [xx,yy]=stairs(xtnn, [Ef;Ef(end)]);
      [xx,yyl]=stairs(xtnn, [Ef-sqrt(Varf);Ef(end)-sqrt(Varf(end))]);
      [xx,yyu]=stairs(xtnn, [Ef+sqrt(Varf);Ef(end)+sqrt(Varf(end))]);
      plot(xx, yy, '-k', xx, yyl, '--b', xx, yyu, '--b')
    end
  end
  
elseif length(ind)==2
  
  uu1=unique(xt(:,ind(1)));
  uu2=unique(xt(:,ind(2)));
  nu1=numel(uu1);
  nu2=numel(uu2);
  if nu1==2 || nu2==2
    % First or second covariate binary
    
    if nu1>2 && nu2==2
      % switch indeces, so that binary covariate is first
      tmp=ind(1);ind(1)=ind(2);ind(2)=tmp;
      tmp=uu1;uu1=uu2;uu2=tmp;
    end
    
    xt1=xt(xt(:,ind(1))==uu1(1),:);
    xt2=xt(xt(:,ind(1))==uu1(2),:);
    [xtnn1, iu1] = unique(xt1(:,ind(2)));
    [xtnn2, iu2] = unique(xt2(:,ind(2)));
    
    options1=options;
    options2=options;
    if ~isempty(z)
      options1.zt = options.zt(iu1);
      options2.zt = options.zt(iu2);
    end

    if isequal(method, 'mean')
      xt1 = repmat(mean(xt1), length(xtnn1), 1);
      xt2 = repmat(mean(xt2), length(xtnn2), 1);
    else
      xt1 = repmat(median(xt1), length(xtnn1), 1);
      xt2 = repmat(median(xt2), length(xtnn2), 1);
    end
    if ~isempty(vars)
      xt1(:,~isnan(vars)) = repmat(vars(~isnan(vars)), length(xtnn1), 1);
      xt2(:,~isnan(vars)) = repmat(vars(~isnan(vars)), length(xtnn2), 1);
    end
    xt1(:,ind(1)) = uu1(1); xt1(:,ind(2)) = xtnn1;
    xt2(:,ind(1)) = uu1(2); xt2(:,ind(2)) = xtnn2;
    
    if ~strcmp(gp.lik.type, 'Coxph')
      [Ef1, Varf1] = gp_pred(gp, x, y, xt1, options1);
      [Ef2, Varf2] = gp_pred(gp, x, y, xt2, options2);
    else
      [Ef11,Ef12,Covf] = pred_coxph(gp,x,y,xt1, options1);
      Ef1 = Ef12; Varf1 = diag(Covf(size(Ef11,1)+1:end,size(Ef11,1)+1:end));
      [Ef21,Ef22,Covf] = pred_coxph(gp,x,y,xt2, options2);
      Ef2 = Ef22; Varf2 = diag(Covf(size(Ef21,1)+1:end,size(Ef21,1)+1:end));
    end
    
    if isequal(plot_results, 'on')
      if nu1>2 && nu2==2
        plot(xtnn1, Ef1, 'ob', xtnn1, Ef1, '-r', xtnn1, Ef1-sqrt(Varf1), '--r', xtnn1, Ef1+sqrt(Varf1), '--r'); hold on;
        plot(xtnn2, Ef2, 'ob', xtnn2, Ef2, '-k', xtnn2, Ef2-sqrt(Varf2), '--k', xtnn2, Ef2+sqrt(Varf2), '--k');
      else
        plot(xtnn1, Ef1, 'ob', xtnn1, Ef1, '-k', xtnn1, Ef1-sqrt(Varf1), '--k', xtnn1, Ef1+sqrt(Varf1), '--k'); hold on;
        plot(xtnn2, Ef2, 'ob', xtnn2, Ef2, '-r', xtnn2, Ef2-sqrt(Varf2), '--r', xtnn2, Ef2+sqrt(Varf2), '--r');
      end
    end
    Ef = [Ef1; Ef2]; Varf = [Varf1; Varf2]; xtnn=[xtnn1;xtnn2];
    
  else
    % first or second covariate is not binary
    xtnn1 = linspace(min(xt(:,ind(1))), max(xt(:,ind(1))), 20);
    xtnn2 = linspace(min(xt(:,ind(2))), max(xt(:,ind(2))), 20);
    [XT1, XT2] = meshgrid(xtnn1, xtnn2); XT1=XT1(:); XT2=XT2(:);
    if ~isempty(z)
      options.zt = repmat(options.zt(1), 400, 1 );
    end
    if isequal(method, 'mean')
      xt = repmat(mean(xt), length(XT1), 1);
    else
      xt = repmat(median(xt), length(XT1), 1);
    end
    if ~isempty(vars)
      xt(:,~isnan(vars)) = repmat(vars(~isnan(vars)), length(XT1), 1);
    end
    xt(:,ind) = [XT1 XT2];
    if ~strcmp(gp.lik.type, 'Coxph')
      [Ef, Varf] = gp_pred(gp, x, y, xt, options);
    else
      [Ef1,Ef2,Covf] = pred_coxph(gp,x,y,xt, options);
      Ef = Ef2; Varf = diag(Covf(size(Ef1,1)+1:end,size(Ef1,1)+1:end));
    end
    
    indd = zeros(size(Ef));
    
    for i2=1:n
      for i3=1:400
        if sqrt(sum((xto(i2,ind)-xt(i3,ind)).^2)) < tr
          indd(i3) = 1;
        end
      end
    end
    
    XT1(indd==0) = NaN; XT2(indd==0) = NaN; Ef(indd==0) = NaN; Varf(indd==0) = NaN;
    
    if isequal(plot_results, 'on')
      contourf(reshape(XT1,20,20), reshape(XT2,20,20), reshape(Ef,20,20))
    end
    
    xtnn = [XT1(indd==1), XT2(indd==1)]; Ef = Ef(indd==1); Varf = Varf(indd==1);
  end
  
else
  error('Only 1 or 2 covariates can be defined for predicting')
end


end

