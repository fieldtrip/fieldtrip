function [stat, cfg, dat] = ft_statfun_gcmi(cfg, dat, design)

% FT_STATFUN_GCMI computes mutual information between the dependent variable
% and a discrete-valued design vector.
%
% You can specify the following configuration options:
%   cfg.preconditionflag = 0 (default), or 1, performs Gaussian copula transform
%                          Preconditioning is computationally efficient, because for given data it needs to be done only once.
%   cfg.gcmi.method      = ['cc', 'cd_model' 'cd_mixture'], type of calculation
%   cfg.gcmi.complex     = ['abs' 'real' 'imag' 'complex' 'angle' ], how to treat complex data
%   cfg.gcmi.tra         = matrix which specifies multivariate structure

if ~strcmp(cfg.precondition,'before')
  error('ft_statfun_gcmi: must use precondition=before')
end
if ~isfield(cfg,'preconditionflag')
  error('ft_statfun_gcmi: must specify preconditionflag')
end
cfg.gcmi             = ft_getopt(cfg, 'gcmi', []);
if ~isfield(cfg.gcmi, 'method')
  error('ft_statfun_gcmi: must specify gcmi method')
end
cfg.gcmi.complex     = ft_getopt(cfg.gcmi, 'complex', false);
% get this one as a local variable only
tra                  = ft_getopt(cfg.gcmi, 'tra', speye(size(dat,1)));

% check the validity of the design
if strcmp(cfg.gcmi.method,'cd_model') || strcmp(cfg.gcmi.method,'cd_mixture')
  % discrete valued ivar, 1-based
  Y = design(cfg.ivar, :)';
  % number of classes
  Ym = max(Y);
  if ~all(ismember(Y, 1:Ym))
    error('ft_statfun_gcmi: the design vector is ill-specified');
  end
  Y = Y-1; % 0-based
elseif strcmp(cfg.gcmi.method,'cc')
  % continuous values (univariate)
  Y = design(cfg.ivar, :)';
  cY = copnorm(Y);
else
  error(sprintf('ft_statfun_gcmi: unknown method \"%s\"', cfg.gcmi.method));
end

datT = dat.';
if cfg.preconditionflag
  
  if ischar(cfg.gcmi.complex) && isreal(dat)
    error('ft_statfun_gcmi: data not complex')
  end
  if ~ischar(cfg.gcmi.complex) && ~isreal(dat)
    error('ft_statfun_gcmi: complex option not specified')
  end
  
  if cfg.gcmi.complex
    tra = speye(size(dat,1)); % for now, this overrules the user-specification
    switch cfg.gcmi.complex
      case 'complex'
        % tease apart the real/imag parts, treat as 2D-variable
        datT = cat(2, real(datT), imag(datT));
        tra = cat(1, tra, tra);
      case 'abs'
        % take the amplitude
        datT = abs(datT);
      case 'angle'
        % tease apart the real/imag parts, after amplitude normalization,
        datT = datT ./ abs(datT);
        datT = cat(2, real(datT), imag(datT));
        tra = cat(1, tra, tra);
      case 'real'
        % just real part
        datT = real(datT);
      case 'imag'
        % just imag part
        datT = imag(datT);
      otherwise
        error(sprintf('ft_statfun_gcmi: unsupported value \"%s\" for gcmi.complex',complex));
    end
  else
    % no multivariate or other transforms, use the user-specified 'tra',
    % which defaults to identity
  end
  fprintf('performing the copula-transform\n');
  % FIXME here we should deal with NaNs in the data, note that these can be different across rows (althought that is rare), which prevents the possibility of doing the transform in a single call.
  % TODO - unconditional copnorm for cd_mixture?
  dat  = copnorm(datT)';
  % save tra for permutation non-precondition runs
  cfg.gcmi.tra = tra;
  
  % build multivariate response once
  % Q: will this work with bootstrap?
  Nvar     = full(sum(tra,1));
  % number of multivariate dimensionalities
  uNvar    = unique(Nvar);
  datT = dat';
  % create multivariate data matrix from tra
  % loop over dimensionality k
  datcel = cell(numel(uNvar),1);
  uNvarN = zeros(numel(uNvar),1);
  selcol = cell(numel(uNvar),1);
  for k = 1:numel(uNvar)
    % multivariate output
    sel_col  = Nvar==uNvar(k);
    % input index
    sel_row  = full(sum(tra(:,sel_col),2)>0);
    
    sel_dat = datT(:,sel_row);
    sel_tra = tra(sel_row,sel_col);
    [ix, iy] = find(sel_tra);
    % multivariate data matrix
    datmc = zeros(size(sel_dat,1), length(sel_col), k);
    for vi1=1:uNvar(k)
      % all voxels for this dimension
      indx1 = ix(vi1:uNvar(k):numel(ix));
      datmc(:,:,vi1) = sel_dat(:,indx1);
    end
    % only do info calculation on variables
    % with no nans, in any multivariate dimension
    idx = all(all(isfinite(datmc),1),3);
    datmcclean = datmc(:,idx,:);
    datcel{k} = datmcclean(:,:)';
    uNvarN(k) = size(datcel{k},1);
    selcol{k} = sel_col & idx;
  end

  if strcmp(cfg.gcmi.method,'cc')
    % check for repeated values (slow, remove?)
    if numel(unique(Y))./numel(Y) < 0.9
      warning('Ivar has more than 10% repeated values.')
    end
  end
  
  cfg.gcmi.uNvar = uNvar;
  cfg.gcmi.uNvarN = uNvarN;
  cfg.gcmi.selcol = selcol;
  cfg.gcmi.Nvar = size(tra,2);
  dat = cell2mat(datcel);
end

stat = [];

% for each size of multivariate calculation unpack the required data
uNvar = cfg.gcmi.uNvar;
uNvarN = cfg.gcmi.uNvarN;
selcol = cfg.gcmi.selcol;

startidx = 1;
datT = dat';
Ntrl = size(dat,2);
stat.stat = NaN(cfg.gcmi.Nvar,1);
for k=1:numel(uNvar)
  endidx = startidx + uNvarN(k) - 1;
  datmc = datT(:,startidx:endidx);
  datmc = reshape(datmc, [Ntrl uNvarN(k)./uNvar(k) uNvar(k)]);
  if strcmp(cfg.gcmi.method,'cd_model')
    %mi = mi_model_gd_vec(datmc, Y, Ym, true, true);
    for m = 1:size(datmc,2)
      mi(m) = mi_model_gd(datmc(:,m), Y, Ym, true, true);
    end

  elseif strcmp(cfg.gcmi.method,'cd_mixture')
    mi = mi_mixture_gd_vec(datmc, Y, Ym);
  elseif strcmp(cfg.gcmi.method,'cc')
    mi = mi_gg_vec(datmc, cY, true, true);
%     for m = 1:size(datmc,2)
%       mi(m) = mi_gg(datmc(:,m), cY, true, true);
%     end
  end
  stat.stat(selcol{k},1) = mi;
end
