function [dataout] = ft_denoise_sss(cfg, datain)

% FT_DENOISE_sss implements an spherical harmonics based
% projection algorithm to suppress interference outside an sphere
% spanned by an MEG array. It is based on: REFERENCE.
%
% Use as
%   dataout = ft_denoise_sss(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'MEG'), see FT_CHANNELSELECTION for details
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
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
cfg.trials             = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel            = ft_getopt(cfg, 'channel',    'all');
cfg.demean             = ft_getopt(cfg, 'demean',     'yes');
cfg.updatesens         = ft_getopt(cfg, 'updatesens', 'yes');
cfg.updateheadposition = ft_getopt(cfg, 'updateheadposition', 'no');
cfg.sss                = ft_getopt(cfg, 'sss');         % sub-structure to hold the parameters
cfg.sss.order_in       = ft_getopt(cfg.sss, 'order_in',  8); 
cfg.sss.order_out      = ft_getopt(cfg.sss, 'order_out', 3);
cfg.sss.thr            = ft_getopt(cfg.sss, 'thr',       0.95); % threshold value for removal of correlated components  
cfg.sss.chunksize      = ft_getopt(cfg.sss, 'chunksize', 10);

if ~isequal(cfg.sss.chunksize, 'none')
  tmpcfg             = [];
  tmpcfg.length      = cfg.sss.chunksize;
  tmpcfg.keeppartial = 'yes';
  datain = ft_redefinetrial(tmpcfg, datain);
end

if istrue(cfg.updateheadposition)
  % currently only possible for input CTF MEG data, which has HLC channels,
  % and for which the grad structure is in dewar coordinates
  assert(startsWith(ft_senstype(datain.grad), 'ctf'));
  assert(strcmp(datain.grad.coordsys, 'dewar'));
  assert(sum(strncmp(datain.label, 'HLC', 3))==9);
  
  tmpcfg          = [];
  tmpcfg.feedback = 'no';
  tmpcfg.computecircumcenter = 'no';
  tmpcfg.method   = 'pertrial';
  [dum, gradin]   = ft_headmovement(tmpcfg, datain);
  clear dum;
 
  gradout = ft_average_sens(gradin);
  % FIXME consider to allow a gradout and gradin to be provided in the
  % input cfg (so that the function behaves a bit like ft_megrealign)

  datain.grad = gradout;

  % remove the HLC channels from the data
  tmpcfg = [];
  tmpcfg.channel = 'MEG';
  datain = ft_selectdata(tmpcfg, datain);
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
options = cfg.sss;
if ~istrue(cfg.updateheadposition)
  S           = sss_spatial(datain.grad, options);
  options.sss = S;
else
  [Sout, S]      = sss_spatial(datain.grad, gradin, options);
  options.sss    = S;
  options.sssout = Sout;
end

% compute the temporal subspace projector and the clean the data
ft_info('Computing the subspace projector based on signal correlations\n');
[datain] = sss_temporal(datain, options);

% apply the spatial projector to the sensors
if istrue(cfg.updatesens) && isscalar(S)
  montage     = [];
  montage.tra = S.Qin*S.iQin;
  montage.labelold = S.labelold;
  montage.labelnew = S.labelnew;
  datain.grad = ft_apply_montage(datain.grad, montage, 'keepunused', 'yes', 'balancename', 'sss');
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
% adjusted from the SPM implementation by Jan-Mathijs Schoffelen
function [varargout] = sss_spatial(grad, varargin)

% sss_SPATIAL computes a collection of spatial projectors based on spherical
% harmonics. The spherical harmonics computation is done by code that has been
% adjusted from Tim Tierney's implementation (adjusted for esthetics). The 
% scaling of the spatial harmonics is different than in MNE-Python. The original
% intention was to get as close as possible to MNE-Python's implementation, but
% that has been difficult in practice, specifically because the heuristic that
% is used for regularisation relies on a certain relative scaling, which I have
% not been able to figure out.
%
% Use as
%  [sss] = sss_spatial(gradout, options)
% or as
%  [sssout, sssin] = sss_spatial(gradout, gradin, options)
%
% The input argument gradout is a FieldTrip-style grad structure and the input argument(s) options specifies the behavior
% of the algorithm. The optional input gradin (as a struct(-array)) remaps the coil positions from gradin to gradout
%
%  options.order_in = scalar (default: 8) order of the in-compartment spherical harmonics
%  options.order_out = scalar (default: 3) order of the out-compartment spherical harmonics
%  options.origin
%  options.bad
%  options.regularize

if numel(varargin)>=1
  % it could be that the second argument is a grad struct-array
  if isfield(varargin{1}, 'coilpos')
    if numel(varargin)>1
      options = varargin{2};
    else
      options = [];
    end

    % we can recurse into the function and combine the output later
    sssout = sss_spatial(grad, options);
    options.origin = sssout.origin;

    for i = 1:numel(varargin{1})
      sssin(i) = sss_spatial(varargin{1}(i), options);
    end
    varargout{2} = sssin;
    varargout{1} = sssout;
    return;
  else
    options = varargin{1};
  end
end

% ft_hastoolbox('opm', 1);

options.order_in  = ft_getopt(options, 'order_in',  8);
options.order_out = ft_getopt(options, 'order_out', 3); 
options.bad       = ft_getopt(options, 'bad', []);
options.regularize = ft_getopt(options, 'regularize', 'no');
options.origin    = ft_getopt(options, 'origin', []);
options.channel   = ft_getopt(options, 'channel', 'all');

grad = ft_datatype_sens(grad);
grad = ft_convert_units(grad, 'm');

if isempty(options.origin)
  ft_error('the origin of the spherical harmonics expansion needs to be provided');
end

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
  ft_warning('using regularization on the basis functions is experimental code, and not thoroughly tested, use at your own risk');
  % this is based on a heuristic that I got from the MNE-python implementation, and is based 
  % on an snr estimate per harmonic basis function. Some pruning is done to exclude the basis
  % functions with the lowest snr. It requires the basis functions to be scaled differently
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
    ft_error('kappa should be specified if options.regularize==2');
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

sss.Q  = Q;
sss.iQ = iQ;
sss.n  = size(Q,2);
sss.nin = nin;
sss.nout = sss.n - sss.nin;

% this is how it seems to be done in practice, which makes the in and out projectors to interact which each other (upon the inversion step)
sss.Qin   = Q(:, 1:sss.nin);
sss.iQin  = iQ(1:sss.nin, :);
sss.Qout  = Q(:, (sss.nin+1):end);
sss.iQout = iQ((sss.nin+1):end, :);

sss.labelold = label(goodchan);
sss.labelnew = label;

sss.labelin = cell(size(sss.Qin,2),1);
for k = 1:numel(sss.labelin)
  sss.labelin{k} = sprintf('sphharm%03din',k);
end
sss.labelout = cell(size(sss.Qout,2),1);
for k = 1:numel(sss.labelout)
  sss.labelout{k} = sprintf('sphharm%03dout',k);
end
sss.origin = options.origin;

varargout{1} = sss;

% if nargout>1
%   % Make montage for the next step
%   montage     = [];
%   montage.tra = sss.Pin;
%   montage.labelold = sss.labelold;
%   montage.labelnew = sss.labelnew;
% 
%   % FIXME think of the mixing of different channel types
% %   montage.chantypeold = data.grad.chantype(i2);
% %   montage.chantypenew = data.grad.chantype(i2);
% %   montage.chanunitold = data.grad.chanunit(i2);
% %   montage.chanunitnew = data.grad.chanunit(i2);
%   varargout{2} = ft_apply_montage(data, montage, 'keepunused', 'no');
%   montage.tra = sss.Pout;
%   varargout{3} = ft_apply_montage(data, montage, 'keepunused', 'no');
% end

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


function [dataclean,vv,ss,n] = sss_temporal(data, options)

% sss_TEMPORAL implements the temporal projection step of the
% tsss algorithm, and follows a call to the companion function
% sss_SPATIAL. 
% 
% Use as:
%   [dataclean] = sss_temporal(data, options)
%
% Where data is a Fieldtrip-style data structure, and options
% is a structure that at least contains a field called sss,
% which contains the spatial projectors, as computed by sss_spatial.
% 
% Other options are
%  options.thr     = scalar (default 0.98), correlation threshold for 
%                    rejection of a 'temporal' component
%  options.bad     = cell-array (default []) of channel labels marked as bad
%  options.st_only = only apply the spatial projection for cleaning
%
% The temporal projectors are computed per trial, i.e. the length of the
% trials determines the temporal support.


options.thr = ft_getopt(options, 'thr', 0.98);
options.bad = ft_getopt(options, 'bad', []);
options.st_only = ft_getopt(options, 'st_only', false);

% check whether there are any bad channels defined for the temporal
% projection. In principle this shouldn't be needed, because the spatial
% projection matrices take care of that. it could be that the numerical
% differences that JM observed while comparing this implementation with the
% MNE-python implementation are caused by not excluding the bad channels at
% this stage.
if ~isempty(options.bad)
  options.bad = ft_channelselection(options.bad, data.grad.label);
end

sss  = options.sss;
if numel(sss)>1
  assert(numel(data.trial)==numel(sss));
end

% it could be that there are fewer channels in the actual data than in the sensors description
[i1, i2] = match_str(data.label, sss(1).labelnew);

dataclean = data;
for i = 1:numel(data.trial)
  cfg        = [];
  cfg.trials = i;
  thisdata   = ft_selectdata(cfg, data);

  if numel(sss)>1
    thissss = sss(i);
  else
    thissss = sss(1);
  end

  if isfield(options, 'sssout')
    sssout = options.sssout;
  else
    sssout = thissss;
  end
  
  % project the data into spherical harmonic space, use ft_apply_montage to ensure correct matching of the order of the channels
  montage          = [];
  montage.tra      = thissss.iQ;
  montage.labelold = thissss.labelold;
  montage.labelnew = [thissss.labelin;thissss.labelout];
  dataQ            = ft_apply_montage(thisdata, montage, 'keepunused', 'no', 'feedback' ,'none');

  % forward project, in only
  montage.tra      = thissss.Qin;
  montage.labelold = thissss.labelin;
  montage.labelnew = thissss.labelnew;
  datain           = ft_apply_montage(dataQ, montage, 'keepunused', 'no', 'feedback', 'none');

  % forward project, out only
  montage.tra      = thissss.Qout;
  montage.labelold = thissss.labelout;
  montage.labelnew = thissss.labelnew;
  dataout          = ft_apply_montage(dataQ, montage, 'keepunused', 'no', 'feedback', 'none');

  cfg              = [];
  cfg.channel      = dataQ.label(contains(dataQ.label, 'in'));
  dataQ            = ft_selectdata(cfg, dataQ);

  datas            = thisdata;
  datas.label      = datas.label(i1);
  datas.trial{1}   = thisdata.trial{1}(i1,:) - datain.trial{1}(i2,:) - dataout.trial{1}(i2,:);
  
  % the below is inteded to mimick MNE-python, which currently estimates the
  % temporal projector on the spatially in and out projected data with
  % omission of the bad channels
  if ~isempty(options.bad)
    cfg = [];
    cfg.channel = setdiff(datas.label, options.bad);
    datas  = ft_selectdata(cfg, datas);
    datain = ft_selectdata(cfg, datain);
  end

  % norm normalise
  Bin  = datain.trial{1}./norm(datain.trial{1}, 'fro');
  Bres = datas.trial{1}./norm(datas.trial{1},   'fro');
  
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
  [U, S, V] = svd(qin'*qres);

  V     = qres*V;
  diagS = diag(S);
  nint  = find(diagS>options.thr,1,'last');
  V     = V(:,1:nint); % temporal basis functions
 
  % in the SPM implementation the temporal projection is applied to the data
  % in spherical harmonic space, the 'in' channels have been selected above
  
  if options.st_only
    tmp = thisdata.trial{1};
    tmp = tmp - (tmp*V)*V';
    dataclean.trial{i} = tmp;
  else
    datacleanQ = dataQ;
    datacleanQ.trial{1} = datacleanQ.trial{1} - (datacleanQ.trial{1}*V)*V';
    
    montage          = [];
    montage.tra      = sssout.Qin;
    montage.labelold = sssout.labelin;
    montage.labelnew = sssout.labelnew;
    tmp = ft_apply_montage(datacleanQ, montage, 'keepunused', 'no', 'feedback', 'none');
    dataclean.trial{i} = tmp.trial{1};
    if i==1
      dataclean.label = tmp.label;
    end
  end
  
  vv{i}               = V;
  ss{i}               = diagS;
  n(1,i)              = nint;
  
end
