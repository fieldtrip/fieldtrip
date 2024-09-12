function [dataout] = ft_denoise_sss(cfg, datain)

% FT_DENOISE_SSS implements an spherical harmonics based
% projection algorithm to suppress interference outside an sphere
% spanned by an MEG array. It is based on: REFERENCE.
%
% Use as
%   dataout = ft_denoise_sss(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.pertrial         = 'no', or 'yes', compute the temporal projection per trial (default = 'no')
%   cfg.demean           = 'yes', or 'no', demean the data per epoch (default = 'yes')
%   cfg.updatesens       = 'yes', or 'no', update the sensor array with the spatial projector
%   cfg.sss              = structure with parameters that determine the behavior of the algorithm
%   cfg.sss.order_in     = scalar. Order of the spherical harmonics basis that spans the in space (default = 8) 
%   cfg.sss.order_out    = scalar. Order of the spherical harmonics basis that spans the out space (default = 3) 
%
% The implementation is based on Tim Tierney's code written for spm
%
% See also FT_DENOISE_PCA, FT_DENOISE_SYNTHETIC, FT_DENOISE_TSR, FT_DENOISE_DSSP, FT_DENOISE_HFC

% Copyright (C) 2024, Jan-Mathijs Schoffelen
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
ft_preamble loadvar    datain
ft_preamble provenance datain

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check the input data
datain = ft_checkdata(datain, 'datatype', {'raw'}); % FIXME how about timelock and freq?

% ensure the external cellfunction toolbox is on the path
ft_hastoolbox('cellfunction', 1);

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.trials            = ft_getopt(cfg, 'trials',   'all', 1);
cfg.channel           = ft_getopt(cfg, 'channel',  'all');
cfg.pertrial          = ft_getopt(cfg, 'pertrial', 'yes');
cfg.demean            = ft_getopt(cfg, 'demean',   'yes');
cfg.updatesens        = ft_getopt(cfg, 'updatesens', 'yes');
cfg.sss               = ft_getopt(cfg, 'sss');         % sub-structure to hold the parameters
cfg.sss.order_in      = ft_getopt(cfg.sss, 'order_in',  8); 
cfg.sss.order_out     = ft_getopt(cfg.sss, 'order_out', 3);
cfg.sss.thr           = ft_getopt(cfg.sss, 'thr', 0.95); % threshold value for removal of correlated components  

pertrial = istrue(cfg.pertrial);
if ~pertrial
  ft_error('a custom chunksize for the time dependent filtering has not yet been implemented, the temporal filtering is applied per trial');
  
  cfg.sss.chunksize = ft_getopt(cfg.sss, 'chunksize', 10); 
else
  cfg.sss.chunksize = ft_getopt(cfg.sss, 'chunksize', inf);
end

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);


if istrue(cfg.demean)
  ft_info('demeaning the time series');
  tmpcfg = keepfields(cfg, {'demean', 'updatesens', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  datain = ft_preprocessing(tmpcfg, datain);
  % restore the provenance information
  [cfg, datain] = rollback_provenance(cfg, datain);
end

ft_info('Computing the spatial subspace projector\n');
options = keepfields(cfg.sss, {'order_in' 'order_out' 'regularize' 'channel' 'bad' 'thr' 'st_only'});
S       = sss_spatial(datain, options);

% compute the temporal subspace projector and the clean the data
ft_info('Computing the subspace projector based on signal correlations\n');
options = keepfields(cfg.sss, {'chunksize', 'thr'});
options.SSS = S;
datain = sss_temporal(datain, options);

% apply the spatial projector to the sensors
if istrue(cfg.updatesens)
  montage     = [];
  montage.tra = S.Pin;
  montage.labelold = S.labelold;
  montage.labelnew = S.labelnew;
  datain.grad = ft_apply_montage(datain.grad, montage, 'keepunused', 'yes', 'balancename', 'amm');
end

% keep some additional information in the subspace struct
subspace.S     = S;

% put some diagnostic information in the output cfg.
cfg.sss.subspace = subspace;

% create the output argument
dataout = keepfields(datain, {'label', 'time', 'trial', 'fsample', 'trialinfo', 'sampleinfo', 'grad'}); 

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrices
% adjusted from the SPM imlpementation by Jan-Mathijs Schoffelen
function [varargout] = sss_spatial(data, options)

% SSS_SPATIAL computes a collection of spatial projectors based on spherical
% harmonics. The spherical harmonics computation is done by code that has been
% adjusted from Tim Tierney's implementation (adjusted for esthetics). The 
% scaling of the spatial harmonics is different than in MNE-Python. The original
% intention was to get as close as possible to MNE-Python's implementation, but
% that has been difficult in practice, specifically because the heuristic that
% is used for regularisation relies on a certain relative scaling, which I have
% not been able to figure out.
%
% Use as
%
%  [sss] = sss_spatial(data, options), or
%  [sss, datain, dataout] = sss_spatial(data, options)
%
% The input argument data is a FieldTrip-style raw data structure containing a 
% correct grad-structure, and the input arguments options specifies the behavior
% of the algorithm.
%
%  options.order_in = scalar (default: 8) order of the in-compartment spherical harmonics
%  options.order_out = scalar (default: 3) order of the out-compartment spherical harmonics
%  options.origin
%  options.bad
%  options.regularize

if nargin<2
  options = [];
end

% ft_hastoolbox('opm', 1);

options.order_in  = ft_getopt(options, 'order_in',  8);
options.order_out = ft_getopt(options, 'order_out', 3); 
options.bad       = ft_getopt(options, 'bad', []);
options.regularize = ft_getopt(options, 'regularize', 'no');
options.origin    = ft_getopt(options, 'origin', [0 0 0]);
options.channel   = ft_getopt(options, 'channel', 'all');

grad = ft_convert_units(data.grad, 'm');
grad = ft_datatype_sens(grad);

ismag = strcmp(grad.chantype, 'mag')|strcmp(grad.chantype, 'megmag');
extended_remove = []; % placeholder

% for now only support unbalanced grad structures, it's the user's responsibility to unbalance
assert(isfield(grad, 'balance') && strcmp(grad.balance.current, 'none'));

% select the list of channels that is required for the output
label   = ft_channelselection(options.channel, grad.label);
selchan = match_str(grad.label, label);
label   = grad.label(selchan);

% coil to channel transformation matrix
tra = grad.tra(selchan, :);

% check whether there are any bad channels defined
if ~isempty(options.bad)
  options.bad = ft_channelselection(options.bad, label);
  badchan     = match_str(label, options.bad);
  goodchan    = setdiff(1:numel(label), badchan(:)');
else
  badchan     = [];
  goodchan    = 1:numel(label);
end

% FIXME build in optional origin of sphere + realignment possibility

% regular spherical harmonics basis functions for the outside field
opt     = [];
opt.li  = options.order_out;
opt.v   = grad.coilpos;
opt.o   = grad.coilori;
opt.or  = options.origin;
opt.reg = 1;
Qout    = tra * spm_opm_vslm(opt);
nout    = size(Qout,2);

% % in comparison to the MNE-python implementation, the harmonics are
% % differently scaled, apply some ad hoc scaling here; this brings the
% % harmonics into the same range across both implementations, which allows
% % for a better probability that the copied regularisation heuristic works 
% r   = mean(grad.coilpos - options.origin);
% 
% [d, o] = get_degrees_orders(options.order_out);
% scl    = pi.*4e-7.*mean(r).^((1:options.order_out)-1);
% scl    = scl(d); % this is not fully OK (yet);
% Qout   = Qout*diag(scl);

% irregular spherical harmonics basis functions for the inside field
opt.li  = options.order_in;
opt.reg = 0;
Qin     = tra * spm_opm_vslm(opt);
nin     = size(Qin,2);

% [d, o] = get_degrees_orders(options.order_in);
% scl    = pi.*4e-7.*mean(r).^-((1:options.order_in)+2);
% scl    = scl(d); % this is not fully OK (yet);
% Qin    = Qin*diag(scl);


% inverse matrix to map from data to in/out space combined
Q       = [Qin Qout];
% if cond(Q) > options.condition_threshold && istrue(options.regularize)
%   [Q, sss_indices, nin] = basis_condition_adjustment(Q, nin, options.condition_threshold);
% end
if options.regularize==1
  % this is based on a heuristic that I got from the MNE-python
  % implementation, and is based on a snr estimate per harmonic basis
  % function. Some pruning is done to exclude the basis functions with the
  % lowest snr. It requires the basis functions to be scaled differently
  % with respect to one another. So far I (JM) have only been able to get
  % this scaling by trial and error approximately right.
  [d, o] = get_degrees_orders(options.order_in);
  r      = grad.coilpos - options.origin;
  r      = sqrt(sum(r.^2,2));
  q      = 3.63859533511; % only tested on a single test case....
  scl    = max(r).^-((1:max(d))-q);
  scl    = scl(d);
  thisQ  = Q;
  thisQ(:,1:nin) = Q(:,1:nin).*scl;

  [in_remove, out_remove] = regularize_in(options.order_in, options.order_out, thisQ, ismag, extended_remove);
  nin  = nin - numel(in_remove);
  Q(:, [in_remove out_remove]) = [];
elseif options.regularize==2
  kappa = ft_getopt(options, 'kappa', []);
  if isempty(kappa)
    ft_error('kappa should be specified if options.regulariz==2');
  end

  % use (implicit) kappa truncated version of the Qin
  [U,S,V] = svd(Qin(goodchan,:));
  S((kappa+1):end,(kappa+1):end) = 0;
  Qin = U*S*V';
  
  Q = [Qin Qout];
end

if options.regularize~=2
  kappa   = size(Q,2);
else
  kappa   = kappa+size(Qout,2);
end

[U,S,V] = svd(Q(goodchan,:));
S       = diag(1./diag(S(1:kappa,1:kappa)));
iQ      = V(:,1:kappa)*S*U(:,1:kappa)';

SSS.P  = Q * iQ;
SSS.Q  = Q;
SSS.iQ = iQ;
SSS.n  = size(Q,2);
SSS.nin = nin;
SSS.nout = SSS.n - SSS.nin;

% this is when the inverses are computed separately
%SSS.iQin  = iQin;
%SSS.iQout = iQout;

% this is how it seems to be done in practice, which makes the in and out projectors to interact which each other (upon the inversion step)
SSS.Qin   = Q(:, 1:SSS.nin);
SSS.iQin  = iQ(1:SSS.nin, :);
SSS.Pin   = SSS.Qin * SSS.iQin;
SSS.Qout  = Q(:, (SSS.nin+1):end);
SSS.iQout = iQ((SSS.nin+1):end, :);
SSS.Pout  = SSS.Qout * SSS.iQout;

SSS.labelold = label(goodchan);
SSS.labelnew = label;

SSS.labelin = cell(size(SSS.Qin,2),1);
for k = 1:numel(SSS.labelin)
  SSS.labelin{k} = sprintf('spharm%03din',k);
end
SSS.labelout = cell(size(SSS.Qout,2),1);
for k = 1:numel(SSS.labelout)
  SSS.labelout{k} = sprintf('spharm%03dout',k);
end

varargout{1} = SSS;

if nargout>1
  % Make montage for the next step
  montage     = [];
  montage.tra = SSS.Pin;
  montage.labelold = SSS.labelold;
  montage.labelnew = SSS.labelnew;

  % FIXME think of the mixing of different channel types
%   montage.chantypeold = data.grad.chantype(i2);
%   montage.chantypenew = data.grad.chantype(i2);
%   montage.chanunitold = data.grad.chanunit(i2);
%   montage.chanunitnew = data.grad.chanunit(i2);
  varargout{2} = ft_apply_montage(data, montage, 'keepunused', 'no');
  montage.tra = SSS.Pout;
  varargout{3} = ft_apply_montage(data, montage, 'keepunused', 'no');
end

function [Qnew, sss_indices, ninnew] = basis_condition_adjustment(Q, nin, thr)

n  = size(Q, 2);
cQ = cond(Q);
sss_indices_in = 1:nin;
Qout = Q(:, (nin+1):end);

while cQ > thr
  for j = 1:length(sss_indices_in)
    Q2 = [Q(:,setdiff(sss_indices_in, sss_indices_in(j))) Qout];
    c(j) = cond(Q2);
  end
  [cQ,drop] = min(c);
  sss_indices_in = setdiff(sss_indices_in,sss_indices_in(drop));
end
Qnew = [Q(:,sss_indices_in) Qout];
sss_indices = [sss_indices_in (nin+1):n];
ninnew = length(sss_indices_in);


function [dataclean,vv,ss] = sss_temporal(data, options)

% SSS_TEMPORAL implements the temporal projection step of the
% tSSS algorithm, and follows a call to the companion function
% SSS_SPATIAL. 
% 
% Use as:
%   [dataclean] = sss_temporal(data, options)
%
% Where data is a Fieldtrip-style data structure, and options
% is a structure that at least contains a field called SSS,
% which contains the spatial projectors, as computed by sss_spatial.
% 
% Other options are
%  options.thr = scalar (default 0.98), correlation threshold for 
%                rejection of a 'temporal' component
%  options.bad = cell-array (default []) of channel labels marked as bad
%  options.st_only = only apply the spatial projection for cleaning
%
% The temporal projectors are computed per trial, i.e. the length of the
% trials determines the temporal support.


options.thr = ft_getopt(options, 'thr', 0.98);
options.bad = ft_getopt(options, 'bad', []);
options.st_only = ft_getopt(options, 'st_only', false);
options.chunksize = ft_getopt(options, 'chunksize', 10); % in seconds


% check whether there are any bad channels defined for the temporal
% projection. In principle this shouldn't be needed, because the spatial
% projection matrices take care of that. it could be that the numerical
% differences that JM observed while comparing this implementation with the
% MNE-python implementation are caused by not excluding the bad channels at
% this stage.
if ~isempty(options.bad)
  options.bad = ft_channelselection(options.bad, data.grad.label);
end

SSS  = options.SSS;

% it could be that there are fewer channels in the actual data than in the sensors description
[i1, i2] = match_str(data.label, SSS.labelnew);
%assert(numel(SSS.labelnew)==numel(data.label));
%assert(isequal(sort(i1), (1:numel(SSS.labelnew))'));
%assert(numel(data.trial)==1); % this may change in the future

% project the data into spherical harmonic space, use ft_apply_montage to
% ensure correct matching of the order of the channels
montage     = [];
montage.tra = SSS.iQ;
montage.labelold = SSS.labelold;
montage.labelnew = [SSS.labelin;SSS.labelout];
dataQ       = ft_apply_montage(data, montage, 'keepunused', 'no');

montage.tra = SSS.Qin;
montage.labelold = SSS.labelin;
montage.labelnew = SSS.labelnew;
datain  = ft_apply_montage(dataQ, montage, 'keepunused', 'no');

montage.tra = SSS.Qout;
montage.labelold = SSS.labelout;
montage.labelnew = SSS.labelnew;
dataout = ft_apply_montage(dataQ, montage, 'keepunused', 'no');

cfg         = [];
cfg.channel = dataQ.label(contains(dataQ.label, 'in'));
dataQ       = ft_selectdata(cfg, dataQ);

datas = data;
datas.label = datas.label(i1);
for k = 1:numel(data.trial)
  datas.trial{k} = data.trial{k}(i1,:) - datain.trial{k}(i2,:) - dataout.trial{k}(i2,:);
end
dataclean = data;
clear dataout data

% note the above is memory inefficient

% the below is inteded to mimick MNE-python, which currently estimates the
% temporal projector on the spatially in and out projected data with
% omission of the bad channels
if ~isempty(options.bad)
  cfg = [];
  cfg.channel = setdiff(datas.label, options.bad);
  datas  = ft_selectdata(cfg, datas);
  datain = ft_selectdata(cfg, datain);
end

% 

for k = 1:numel(datain.trial)
  % norm normalise
  Bin  = datain.trial{k}./norm(datain.trial{k}, 'fro');
  Bres = datas.trial{k}./norm(datas.trial{k}, 'fro');
  
  %Uin  = orth(Bin');
  %Ures = orth(Bres');

  % MNE-Python obtains the orthonormal basis with an svd
  [Uin, Sin, Vin] = svd(Bin','econ');
  tol = max(size(Bin)) * Sin(1) * eps;
  sel = diag(Sin)>tol;
  Uin = Uin(:, sel);

  [Ures, Sres, Vres] = svd(Bres','econ');
  tol = max(size(Bres)) * Sres(1) * eps;
  sel = diag(Sres)>tol;
  Ures = Ures(:, sel);

  [qin,  rin]  = qr(Uin,0);
  [qres, rres] = qr(Ures,0);
  [U, S, V] = svd(qin'*qres); clear U

  V     = qres*V;
  diagS = diag(S);
  nint  = find(diagS>options.thr,1,'last');
  V     = V(:,1:nint); % temporal basis functions

  % remove the correlated signals from the data
  %dataclean = datain;
  %dataclean.trial{1} = datain.trial{1} - (datain.trial{1}*U)*U';

  % in the SPM implementation the temporal projection is applied to the data
  % in spherical harmonic space, the 'in' channels have been selected above
  cfg = [];
  cfg.trials = k;
  if options.st_only
    tmp = dataclean.trial{k};
    tmp = tmp - (tmp*V)*V';
    dataclean.trial{k} = tmp;
  else
    if k==1
      datacleanQ          = ft_selectdata(cfg, dataQ);
    else
      datacleanQ.trial    = dataQ.trial(k);
      datacleanQ.time     = dataQ.time(k);
    end
    datacleanQ.trial{1} = datacleanQ.trial{1} - (datacleanQ.trial{1}*V)*V';
    
    montage          = [];
    montage.tra      = SSS.Qin;
    montage.labelold = SSS.labelin;
    montage.labelnew = SSS.labelnew;
    tmp = ft_apply_montage(datacleanQ, montage, 'keepunused', 'no');
    dataclean.trial{k} = tmp.trial{1};
    if k==1
      dataclean.label = tmp.label;
    end
  end
  
  vv{k}               = V;
  ss{k}               = diagS(1:nint);
  
end
