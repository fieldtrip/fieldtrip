function [data] = ft_regressconfound(cfg, datain)

% FT_REGRESSCONFOUND estimates the regression weight of a set of confounds
% using a General Linear Model (GLM) and removes the estimated contribution
% from the single-trial data.
%
% Use as
%   timelock = ft_regressconfound(cfg, timelock)
% or as
%   freq     = ft_regressconfound(cfg, freq)
% or as
%   source   = ft_regressconfound(cfg, source)
%
% where timelock, freq, or, source come from FT_TIMELOCKANALYSIS,
% FT_FREQANALYSIS, or FT_SOURCEANALYSIS respectively, with keeptrials = 'yes'
%
% The cfg argument is a structure that should contain
%   cfg.confound    = matrix, [Ntrials X Nconfounds], may not contain NaNs
%
% The following configuration options are supported:
%   cfg.reject      = vector, [1 X Nconfounds], listing the confounds that
%                     are to be rejected (default = 'all')
%   cfg.normalize   = string, 'yes' or 'no', normalization to
%                     make the confounds orthogonal (default = 'yes')
%   cfg.output      = 'residual' (default), 'beta', or 'model'.
%                     If 'residual' is specified, the output is a data
%                     structure containing the residuals after regressing
%                     out the in cfg.reject listed confounds. If 'beta' or 'model'
%                     is specified, the output is a data structure containing
%                     the regression weights or the model, respectively.
%
% This method is described by Stolk et al., Online and offline tools for head
% movement compensation in MEG (Neuroimage, 2013)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REJECTCOMPONENT, FT_REJECTARTIFACT

% Undocumented local options:
%   cfg.ftest       = string array, {N X Nconfounds}, to F-test whether
%                     the full model explains more variance than reduced models
%                     (e.g. {'1 2'; '3 4'; '5'} where iteratively the added value of
%                     regressors 1 and 2, and then 3 and 4, etc., are tested)
%   cfg.statistics  = string, 'yes' or 'no', whether to add the statistics
%                     on the regression weights to the output (default = 'no',
%                     applies only when cfg.output = 'beta')

% Copyright (C) 2011-2017, Arjen Stolk, Robert Oostenveld, Lennart Verhagen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
datain = ft_checkdata(datain, 'datatype', {'timelock', 'freq', 'source'}, 'feedback', 'yes');

if isfield(cfg, 'beta') || isfield(cfg, 'model')
 ft_error('The options cfg.beta and cfg.model have been removed as of Aug 2017, please use cfg.output instead');
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'confound'}, 'renamed', {'Ftest','ftest'}, 'forbidden', {'beta','model'});

% specify the defaults
cfg.confound   = ft_getopt(cfg, 'confound');
cfg.reject     = ft_getopt(cfg, 'reject', 'all');
cfg.normalize  = ft_getopt(cfg, 'normalize', 'yes');
cfg.output     = ft_getopt(cfg, 'output', 'residual');
cfg.statistics = ft_getopt(cfg, 'statistics', 'no');
cfg.ftest      = ft_getopt(cfg, 'ftest');
cfg.parameter  = ft_getopt(cfg, 'parameter'); % the default is handled further down

regr = cfg.confound;
if any(isnan(regr(:)))
  ft_error('the confounds may not contain NaNs');
end
nconf     = size(regr,2);
conflist  = 1:nconf;
if strcmp(cfg.reject, 'all')
  cfg.reject = conflist(1:end); % to be removed
else
  cfg.reject = intersect(conflist, cfg.reject); % to be removed
end

fprintf('removing confound %s \n', num2str(cfg.reject));
kprs = setdiff(conflist, cfg.reject); % to be kept
fprintf('keeping confound %s \n', num2str(kprs));

% confound normalization for orthogonality
if strcmp(cfg.normalize, 'yes')
  fprintf('normalizing the confounds, except the constant \n');
  for c = 1:nconf
    AVG = mean(regr(:,c));
    STD = std(regr(:,c),0,1);
    if abs(STD/AVG)<10*eps
      fprintf('confound %s is a constant \n', num2str(c));
    else
      regr(:,c) = (regr(:,c) - AVG) / STD;
    end
    clear AVG STD;
  end
else
  fprintf('skipping normalization procedure \n');
end

switch ft_datatype(datain)
  case 'freq'
    cfg.parameter = ft_getopt(cfg, 'parameter', 'powspctrm');
  case 'timelock'
    cfg.parameter = ft_getopt(cfg, 'parameter', 'trial');
  case 'source'
    cfg.parameter = ft_getopt(cfg, 'parameter', 'pow');
end

dimord = getdimord(datain, cfg.parameter);
dimsiz = getdimsiz(datain, cfg.parameter);
dimtok = tokenize(dimord, '_');
rptdim = find(strcmp(dimtok, 'rpt'));
datdim = setdiff(1:length(dimtok), rptdim);

nrpt = dimsiz(rptdim);

dat = datain.(cfg.parameter);
if strcmp(dimtok{1}, '{pos}')
  indx = find(datain.inside);
  npos = length(indx);
  tmp = nan([nrpt npos datdim(2:end)]); % only positions inside the brain
  for i=indx'
    tmp(:,i,:,:,:) = dat{i}(:,:,:,:);
  end
  dat = tmp;
else
  dat = permute(dat, [rptdim datdim]);
end

dat = reshape(dat, nrpt, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM MODEL
%   Y = X * B + err, where Y is data, X is the model, and B are beta's
% which means
%   Best = X\Y ('matrix division', which is similar to B = inv(X)*Y)
% or when presented differently
%   Yest = X * Best
%   Yest = X * X\Y
%   Yclean = Y - Yest (the true 'clean' data is the recorded data 'Y' -
%   the data containing confounds 'Yest')
%   Yclean = Y - X * X\Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate and remove the confounds
fprintf('estimating the regression weights and removing the confounds \n');
if isempty(find(isnan(dat))) % if there are no NaNs, process all at once 
  beta = regr\dat;                                                        % B = X\Y 
else % otherwise process per colum set as defined by the nan distribution  
  [u,i,j] = unique(~isnan(dat)','rows','first'); % find unique rows
  uniquecolumns = u'; % unique column types
  Nuniques = numel(i); % number of unique types
  beta_temp = NaN(Nuniques, nconf, size(dat,2)); % declare empty variable
  for n = 1:Nuniques % for each unique type
    rowidx = find(uniquecolumns(:,n)==1); % row indices for unique type
    colidx = find(j==n); % column indices for unique type
    if any(uniquecolumns(:,n)) % if vector contains a nonzero number
      beta_temp(n,:,colidx) = regr(rowidx,:)\dat(rowidx,colidx);         % B = X\Y
    end
  end
  beta = reshape(nansum(beta_temp,1),[nconf size(dat,2)]); % sum the betas
  clear beta_temp
end

model = regr(:, cfg.reject) * beta(cfg.reject, :);                        % model = confounds * weights = X * X\Y
Yc = dat - model;                                                         % Yclean = Y - X * X\Y

% reduced models analyses
if ~isempty(cfg.ftest)  
  dfe        = nrpt - nconf;                                              % degrees of freedom
  err        = dat - regr * beta;                                         % err = Y - X * B
  tmse       = sum((err).^2)/dfe;                                         % mean squared error
  for iter = 1:numel(cfg.ftest)    
    % regressors to test if they explain additional variance
    r          = str2num(cfg.ftest{iter});
    fprintf('F-testing explained additional variance of regressors %s \n', num2str(r));
    % regressors in reduced design (that is the original design)
    ri         = ~ismember(1:size(regr,2),r);
    rX         = regr(:,ri);               % reduced design
    rnr        = size(rX,2);               % number of regressors in reduced design
    % estimate reduced model betas
    rXcov      = pinv(rX'*rX);             % inverse design covariance matrix
    rb         = rXcov*rX'*dat;          	 % beta estimates using pinv
    % calculate mean squared error of reduced model
    rdfe       = size(dat,1) - size(rX,2); % degrees of freedom of the error
    rerr       = dat-rX*rb;                % residual error
    rmse       = sum(rerr'.^2,2)./rdfe;	   % mean squared error
    % F-test
    F(iter,:)  = ((rmse'-tmse)./(nconf-rnr)) ./ (tmse./(dfe-2));
    % Rik Henson defined F-test
    % F = ( ( rerr'*rerr - err'*err ) / ( nconf-rnr ) ) / ( err'*err/ ( nrpt-nconf ) );
    % convert F-value to p-value
    idx_pos    = F(iter,:) >= 0;
    idx_neg    = ~idx_pos;
    p(iter,:)     = nan(1,size(F(iter,:),2));
    p(iter,idx_pos) = (1-fcdf(F(iter,idx_pos),rnr,rdfe));
    p(iter,idx_neg) = fcdf(-F(iter,idx_neg),rnr,rdfe);
    clear rerr rmse
    % FIXME: drop in replace tcdf from the statfun/private dir   
  end
  clear dfe err tmse
end

% organize the output
dataout = keepfields(datain, {'label', 'time', 'freq', 'pos', 'dim', 'transform', 'inside', 'outside', 'trialinfo', 'sampleinfo', 'dimord'});
switch cfg.output
  case 'residual'
    dataout.(cfg.parameter) = reshape(Yc, [nrpt dimsiz(datdim)]); % either powspctrm, trial, or pow
    clear Yc   
  case 'beta'
    dataout.beta = reshape(beta, [nconf, dimsiz(datdim)]);
    if strcmp(cfg.statistics, 'yes') % beta statistics
      fprintf('performing statistics on the regression weights \n');
      dfe        = nrpt - nconf;                                              % degrees of freedom
      err        = dat - regr * beta;                                         % err = Y - X * B
      mse        = sum((err).^2)/dfe;                                         % mean squared error
      covar      = diag(regr'*regr)';                                         % regressor covariance
      bvar       = repmat(mse',1,size(covar,2))./repmat(covar,size(mse,2),1); % beta variance
      tval       = (beta'./sqrt(bvar))';                                      % betas -> t-values
      prob       = (1-tcdf(tval,dfe))*2;                                      % p-values
      clear err dfe mse bvar
      % FIXME: drop in replace tcdf from the statfun/private dir
      dataout.stat = reshape(tval, [nconf dimsiz(datdim)]);
      dataout.prob = reshape(prob, [nconf dimsiz(datdim)]);
      clear tval prob
    end    
  case 'model'
    dataout.model = keepfields(datain, {'label', 'time', 'freq', 'pos', 'dim', 'transform', 'inside', 'outside', 'trialinfo', 'sampleinfo', 'dimord'});
    dataout.model.(cfg.parameter) = reshape(model, [nrpt, dimsiz(datdim)]);
  otherwise
    error('output ''%s'' is not supported', cfg.output);    
end

% reduced models analyses
if ~isempty(cfg.ftest)
  dataout.fvar   = reshape(F, [numel(cfg.ftest) dimsiz(datdim)]);
  dataout.pvar   = reshape(p, [numel(cfg.ftest) dimsiz(datdim)]);
  clear F p
end

% discard the gradiometer information because the weightings have been changed
if isfield(dataout, 'grad')
  ft_warning('discarding gradiometer information because the weightings have been changed');
  dataout = rmfield(dataout, 'grad');
end

% discard the electrode information because the weightings have been changed
if isfield(dataout, 'elec')
  ft_warning('discarding electrode information because the weightings have been changed');
  dataout = rmfield(dataout, 'elec');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous datain

% rename the output variable to accomodate the savevar postamble
data = dataout;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
