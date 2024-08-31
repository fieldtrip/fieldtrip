function [dataout] = ft_denoise_amm(cfg, datain)

% FT_DENOISE_AMM implements an adaptive multipole modelling based
% projection algorithm to suppress interference outside an ellipsoid
% spanned by an MEG array. It is based on: REFERENCE.
%
% Use as
%   dataout = ft_denoise_amm(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'MEG'), see FT_CHANNELSELECTION for details
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.pertrial         = 'no', or 'yes', compute the temporal projection per trial (default = 'no')
%   cfg.demean           = 'yes', or 'no', demean the data per epoch (default = 'yes')
%   cfg.updatesens       = 'yes', or 'no', update the sensor array with the spatial projector
%   cfg.amm              = structure with parameters that determine the behavior of the algorithm
%   cfg.amm.order_in     = scalar. Order of the spheroidal harmonics basis that spans the in space (default = 9) 
%   cfg.amm.order_out    = scalar. Order of the spheroidal harmonics basis that spans the out space (default = 2) 
%   cfg.amm.reducerank
%   cfg.amm.thr 
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
cfg.trials            = ft_getopt(cfg, 'trials',  'all', 1);
cfg.channel           = ft_getopt(cfg, 'channel', 'MEG');
cfg.pertrial          = ft_getopt(cfg, 'pertrial', 'yes');
cfg.demean            = ft_getopt(cfg, 'demean', 'yes');
cfg.updatesens        = ft_getopt(cfg, 'updatesens', 'yes');
cfg.amm               = ft_getopt(cfg, 'amm');         % sub-structure to hold the parameters
cfg.amm.order_in      = ft_getopt(cfg.amm, 'order_in',  9); % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.amm.order_out     = ft_getopt(cfg.amm, 'order_out', 2); % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.amm.thr           = ft_getopt(cfg.amm, 'thr', 1); % threshold value for removal of correlated components  

pertrial = istrue(cfg.pertrial);
if ~pertrial
  cfg.amm.chunksize = ft_getopt(cfg.amm, 'chunksize', 10); 
else
  cfg.amm.chunksize = ft_getopt(cfg.amm, 'chunksize', inf);
end

% ensure spm12 on the path
ft_hastoolbox('spm12',1);

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

% compute the spatial projectors, don't apply them yet, will be done below
ft_info('Computing the spatial subspace projector\n');
options = keepfields(cfg.amm, {'order_in' 'order_out' 'reducerank' 'channel' 'bad' 'updatesens'});
options.channel = cfg.channel; % the channel argument needs to be dealt with below, not in ft_selectdata above
S       = amm_spatial(datain, options);

% compute the temporal subspace projector and the clean the data
ft_info('Computing the subspace projector based on signal correlations\n');
options = keepfields(cfg.amm, {'chunksize', 'thr'});
options.AMM = S;
datain = amm_temporal(datain, options);

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
cfg.amm.subspace = subspace;

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
function [varargout] = amm_spatial(data, options)

% AMM_SPATIAL computes a set of spatial projection matrices based on spheroidal harmonics
%
% use as
%
% [amm, datain, dataout] = amm_spatial(data, options)
%
% options.order_in
% options.order_out
% options.bad 
% options.reducerank percentage of total spatial variance (in compartment) explained
% options.channel
% options.updatesens

%function [mfD,Yinds] = spm_opm_amm(S)
% models brain signal and interference as a set of geometrically adaptive
% multipole moments
% FORMAT D = spm_opm_amm(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.li             -  internal harmonic order   - Default: 9
%   S.le             -  external harmonic order   - Default: 2
%   S.window        - temporal window size (s)        - 10
%   S.prefix        - prefix to filename          - Default 'm'
%   S.corrLim       - correlation limit          - Default 1
%   S.plotSpheroid  - flag to plot spheroid      - Default 1
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright  Tim Tierney

if nargin<2
  options = [];
end

options.order_in   = ft_getopt(options, 'order_in',  9);
options.order_out  = ft_getopt(options, 'order_out', 2); 
options.bad        = ft_getopt(options, 'bad', []);
options.reducerank = ft_getopt(options, 'reducerank'); % should be a number
options.channel    = ft_getopt(options, 'channel', 'all');
options.updatesens = ft_getopt(options, 'updatesens', 'yes');

grad = ft_convert_units(data.grad, 'm');
grad = ft_datatype_sens(grad);
grad = ft_convert_coordsys(grad, 'ras'); % to ensure that the second axis is the longest

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
pos = grad.chanpos(selchan, :); % FIXME: consider making this the same as SSS, where the harmonics are computed on the (integration points of the) coils 
ori = grad.chanori(selchan, :);

% check whether there are any bad channels defined
if ~isempty(options.bad)
  options.bad = ft_channelselection(options.bad, label);
  badchan     = match_str(label, options.bad);
  goodchan    = setdiff(1:numel(label), badchan(:)');
else
  badchan     = [];
  goodchan    = 1:numel(label);
end
    
%-fit the ellipsoid -> FIXME for axial gradiometer coil systems this makes more sense to do only on the bottom coils
%--------------------------------------------------------------------------
v = pos;
n = ori;
vrange   = abs((max(v)-min(v)));
[tmp, ind] = max(vrange);
[o, r]   = spheroid_fit(v, ind);

if (ind~=2)
  % FIXME this assumes RAS(-like) coordinates, so for ALS the coordinates need to be temporarily adjusted, but I don't see back below why this would be relevant
  inderror('Y is not longest axis.... fix please')
end

inside = (v(:,1)-o(1)).^2/r(1)^2 + (v(:,2)-o(2)).^2/r(2)^2 + (v(:,3)-o(3)).^2/r(3)^2;
c      = sum(inside>1);
stepsize = max(r*.005);

while c~=size(v,1)
  rt = r-stepsize;
  inside = (v(:,1)-o(1)).^2/rt(1)^2 + (v(:,2)-o(2)).^2/rt(2)^2 + (v(:,3)-o(3)).^2/rt(3)^2;
  cc = sum(inside>1);
  if(cc>=c)
    r = r-stepsize;  
    c = cc;
  end 
end

plotSpheroid = false;
if plotSpheroid
  figure()
  plot3(v(:,1),v(:,2),v(:,3),'.k')
  daspect([1,1,1])
  hold on
  [X,Y,Z]=ellipsoid(o(1),o(2),o(3),r(1),r(2),r(3),10);
  plot3(X(:),Y(:),Z(:),'.')
  daspect([1,1,1])
end

%-construct the projectors
%--------------------------------------------------------------------------
a = max(r);
b = min(r);

vtest     = double(bsxfun(@minus,v,o'));
external  = spm_epharm(vtest, n, a, b, options.order_out);
inelipse  = spm_ipharm(vtest, n, a, b, options.order_in);

if ~isempty(options.reducerank)
  [Q,s] = svd(inelipse,'econ'); 
  Ve    = cumsum(diag(s))/sum(diag(s));
  inelipse = Q(:, Ve < options.reducerank);
end

Pout = external*pinv(external);
M    = eye(size(external,1))-Pout; % hmmm, this is slightly different from the spatial in/out projectors for SSS
Pin  = M*inelipse*pinv(M*inelipse)*M; % this projector is different from how it's computed for SSS, corresponds almost to eq.18 in Tim's biorxiv manuscript (apart from the post multiplication with M).

AMM.Pout  = Pout;
AMM.Pin   = Pin;
AMM.Qout  = external;
AMM.iQout = pinv(external);
AMM.Qin   = M*inelipse;

AMM.labelold = label(goodchan);
AMM.labelnew = label;

AMM.labelin = cell(size(AMM.Pin,2),1);
for k = 1:numel(AMM.labelin)
  AMM.labelin{k} = sprintf('amm%03din',k);
end
AMM.labelout = cell(size(AMM.Pout,2),1);
for k = 1:numel(AMM.labelout)
  AMM.labelout{k} = sprintf('amm%03dout',k);
end

varargout{1} = AMM;

if nargout>1
  % Make montage for the next step
  montage     = [];
  montage.tra = AMM.Pin;
  montage.labelold = AMM.labelold;
  montage.labelnew = AMM.labelnew;

  % FIXME think of the mixing of different channel types
  varargout{2} = ft_apply_montage(data, montage, 'keepunused', 'no');
  if istrue(options.updatesens)
    varargout{2}.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'yes', 'balancename', 'amm');
  end
  montage.tra = AMM.Pout;
  varargout{3} = ft_apply_montage(data, montage, 'keepunused', 'no');
  if istrue(options.updatesens)
    varargout{3}.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'yes', 'balancename', 'amm');
  end
end


function [ o, r] = spheroid_fit( X, ax )
%  [o, r] = ellipsoid_fit( X, ax);
%
% Parameters:
%  X   - Coordinates  n x 3 matrix
%  ax  - numeric indicating longer axis
%
% Output:
% o   - origin
% r   - radii


x =X(:,1);
y =X(:,2);
z =X(:,3);
on = ones(size(x,1),1);
b = x.^2 + y.^2 + z.^2;
if ax==1
    A = [ y.^2 + z.^2 - 2*x.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b;
    v(1) = -2 * beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==2
    A = [ x.^2 + z.^2 - 2*y.^2, 2*x,2*y,2*z,on];
    beta = pinv(A)*b;
     v(1) = beta(1)-1;
    v(2) = -2*beta(1)-1;
    v(3) = beta(1)-1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==3
    A = [ x.^2 + y.^2 - 2*z.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b; 
    v = beta;
    v(1) = beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = -2*beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

A = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
    

o = -A( 1:3, 1:3 ) \ v( 7:9 )';
T = eye( 4 );
T( 4, 1:3 ) = o';
R = T * A * T';
[ vec, s ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
r = sqrt( 1 ./ diag( abs( s ) ) );
sgns = sign( diag( s ) );
r = r .* sgns;
r = vec*r;


function dataclean = amm_temporal(data, options)

% AMM_TEMPORAL removes the correlated components (between the 'in' and
% 'inter' space from the 'in' space signal, 


% DOCSTRING of the function from which the relevant code was harvested
% function [mfD,Yinds] = spm_opm_amm(S)
% models brain signal and interference as a set of geometrically adaptive
% multipole moments
% FORMAT D = spm_opm_amm(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.li             -  internal harmonic order   - Default: 9
%   S.le             -  external harmonic order   - Default: 2
%   S.window        - temporal window size (s)        - 10
%   S.prefix        - prefix to filename          - Default 'm'
%   S.corrLim       - correlation limit          - Default 1
%   S.plotSpheroid  - flag to plot spheroid      - Default 1
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright  Tim Tierney

options.thr       = ft_getopt(options, 'thr', 1);
options.chunksize = ft_getopt(options, 'chunksize', 10); % in seconds

AMM = options.AMM;

%-Update forward modelling information
%--------------------------------------------------------------------------
%fprintf('%-40s: %30s\n','Updating Sensor Information',spm('time'));
grad = ft_convert_units(data.grad, 'm');
grad = ft_datatype_sens(grad);
grad = ft_convert_coordsys(grad, 'ras');

dataclean = data; % FIXME this needs to address bad channels etc

%-canonical correlations
%--------------------------------------------------------------------------
Yinds = 1:numel(data.label); % FIXME this is intended to address bad channels I think
for i = 1:numel(data.trial)
  Ytemp = data.trial{i};
  Y     = Ytemp(Yinds,:);
  inner = AMM.Pin*Y;
  Ytemp(Yinds,:) = inner;
  
  if options.thr<1
    % otherwise only the spatial projection will be applied  
    nsmp = round(data.fsample.*options.chunksize);
    if size(Ytemp,2) > nsmp
      smpinds = (1:nsmp:(size(Ytemp,2)-1));
      smpinds(2,:) = nsmp.*(1:size(smpinds,2));
      smpinds = smpinds'; 
    else
      smpinds = [1 size(Ytemp,2)];
    end 

    outer = AMM.Pout*Y;
    inter = Y-inner-outer;
    for k = 1:size(smpinds,1)
      inner_ = inner(:,smpinds(k,1):smpinds(k,2));
      inter_ = inter(:,smpinds(k,1):smpinds(k,2));

      Oinner = orth(inner_');
      Ointer = orth(inter_');
      C = Oinner'*Ointer;
      [tmp,Sc,Z] = svd(C);
      noise = Ointer*Z;
      s     = diag(Sc);
      noisevec = noise(:,1:sum(s>options.thr));
      Beta = noisevec'*inner_';
      mod  = noisevec*Beta;
      binnew = inner_-mod';
      Ytemp(Yinds,smpinds(k,1):smpinds(k,2))=binnew;
    end
  end
  dataclean.trial{i}= Ytemp;
end

