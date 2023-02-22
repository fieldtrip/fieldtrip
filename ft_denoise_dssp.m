function [dataout] = ft_denoise_dssp(cfg, datain)

% FT_DENOISE_DSSP implements a dual signal subspace projection algorithm
% to suppress interference outside a predefined source region of
% interest. It is based on: Sekihara et al. J. Neural Eng. 2016 13(3), and
% Sekihara et al. J. Neural Eng. 2018 15(3).
%
% Use as
%   dataout = ft_denoise_dssp(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials           = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.pertrial         = 'no', or 'yes', compute the temporal projection per trial (default = 'no')
%   cfg.sourcemodel      = structure, source model with precomputed leadfields, see FT_PREPARE_LEADFIELD
%   cfg.demean           = 'yes', or 'no', demean the data per epoch (default = 'yes')
%   cfg.dssp             = structure with parameters that determine the behavior of the algorithm
%   cfg.dssp.n_space     = 'all', or scalar. Number of dimensions for the
%                          initial spatial projection.
%   cfg.dssp.n_in        = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field inside the ROI.
%   cfg.dssp.n_out       = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field outside the ROI.
%   cfg.dssp.n_intersect = scalar (default = 0.9). Number of dimensions (if
%                          value is an integer>=1), or threshold for the
%                          included eigenvalues (if value<1), determining
%                          the dimensionality of the intersection.
%
% See also FT_DENOISE_PCA, FT_DENOISE_SYNTHETIC, FT_DENOISE_TSR

% Copyright (C) 2018-2023, Jan-Mathijs Schoffelen
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
cfg.channel           = ft_getopt(cfg, 'channel', 'all');
cfg.pertrial          = ft_getopt(cfg, 'pertrial', 'no');
cfg.sourcemodel       = ft_getopt(cfg, 'sourcemodel');
cfg.demean            = ft_getopt(cfg, 'demean', 'yes');
cfg.dssp              = ft_getopt(cfg, 'dssp');         % sub-structure to hold the parameters
cfg.dssp.n_space      = ft_getopt(cfg.dssp, 'n_space', 'interactive'); % number of spatial components to retain from the Gram matrix
cfg.dssp.n_in         = ft_getopt(cfg.dssp, 'n_in',    'interactive'); % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.dssp.n_out        = ft_getopt(cfg.dssp, 'n_out',   'interactive'); % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.dssp.n_intersect  = ft_getopt(cfg.dssp, 'n_intersect', 'interactive'); % dimensionality of the intersection
cfg.output            = ft_getopt(cfg, 'output', 'original');

pertrial = istrue(cfg.pertrial);

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);


if istrue(cfg.demean)
  ft_info('demeaning the time series');
  tmpcfg = [];
  tmpcfg.demean = 'yes';
  datain = ft_preprocessing(tmpcfg, datain);
  % restore the provenance information
  [cfg, datain] = rollback_provenance(cfg, datain);
end

% match the input data's channels with the labels in the leadfield
sourcemodel = cfg.sourcemodel;
if ~isfield(sourcemodel, 'leadfield')
  ft_error('cfg.sourcemodel needs to contain leadfields');
end
[indx1, indx2] = match_str(datain.label, sourcemodel.label);
if ~isequal(indx1(:),(1:numel(datain.label))')
  ft_error('unsupported mismatch between data channels and leadfields');
end
if islogical(sourcemodel.inside)
  inside = find(sourcemodel.inside);
else
  inside = sourcemodel.inside;
end
for k = inside(:)'
  sourcemodel.leadfield{k} = sourcemodel.leadfield{k}(indx2,:);
end

% compute the Gram-matrix of the supplied forward model
lf = cat(2, sourcemodel.leadfield{:});
G  = lf*lf';

%dat     = cat(2,datain.trial{:});
[Bclean, subspace] = dssp(datain.trial, G, cfg.dssp.n_in, cfg.dssp.n_out, cfg.dssp.n_space, cfg.dssp.n_intersect, pertrial);
% datAe   = datain.trial*cellfun(@transpose, Ae, 'UniformOutput', false); % the projection is a right multiplication
% with a matrix (eye(size(Ae,1))-Ae*Ae'), since Ae*Ae' can become quite
% sizeable, it's computed slightly differently here.

% put some diagnostic information in the output cfg.
cfg.dssp.subspace = subspace;

% replace the input cfg values
cfg.dssp.n_space = subspace.S(1).n;
cfg.dssp.n_in    = subspace.Sin(1).n;
cfg.dssp.n_out   = subspace.Sout(1).n;
cfg.dssp.n_intersect = subspace.T(1).n;

% compute the cleaned data and put in a cell-array
switch cfg.output
  case 'original'
    trial = Bclean;
  case 'complement'
    trial = datain.trial-Bclean;
  otherwise
    ft_error(sprintf('cfg.output = ''%s'' is not implemented',cfg.output));
end

% create the output argument
dataout       = keepfields(datain, {'label', 'time', 'fsample', 'trialinfo', 'sampleinfo', 'grad', 'elec', 'opto'}); % grad can be kept and does not need to be balanced, since the cleaned data is a mixture over time, not space.
dataout.trial = trial;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrix
% kindly provided by Kensuke, and adjusted a bit by Jan-Mathijs
function [Bclean, subspace] = dssp(B, G, Nin, Nout, Nspace, Nintersect, pertrial)

% Nc: number of sensors
% Nt: number of time points
% inputs
% B(Nc,Nt):  interference overlapped sensor data
% G(Nc,Nc): Gram matrix of voxel lead field
% Nout and Nin: dimensions of the two row spaces
% recom_Nspace: recommended value for the dimension of the pseudo-signal subspace
% outputs
% Bclean(Nc,Nt): cleaned sensor data
% Nintersect: dimension of the intersection
% Nspace: dimension of the pseudo-signal subspace
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%
% The code below is modified by Jan-Mathijs, no functional changes
% merely cosmetics, added the possibility to run the temporal subspace per
% trial

% eigen decomposition of the Gram matrix, matrix describing the spatial components of the defined 'in' compartment
fprintf('Computing the spatial subspace projection\n');
fprintf('Eigenvalue decomposition of the Gram matrix\n');
[Uspace,S] = eig(G);
Sspace     = abs(diag(S));

[Sspace, iorder] = sort(-Sspace);
Sspace           = -Sspace;
Uspace(:,:)      = Uspace(:,iorder);
Nspace           = getN(Nspace, Sspace, 'spatial');

% spatial subspace projector
Us   = Uspace(:,1:Nspace);

%USUS = Us*Us';
% % Bin and Bout creation
%Bin  =                  USUS  * B;
%Bout = (eye(size(USUS))-USUS) * B;

% computationally more efficient than the above
fprintf('Applying the spatial subspace projector\n');
Bin  = Us*(Us'*B);
Bout = B - Bin; 

fprintf('Computing the subspace projector based on signal correlations\n');
[Ae, subspace] = CSP01(Bin, Bout, Nin, Nout, Nintersect, pertrial);

% add the first spatial subspace projection information as well
subspace.S.U = Uspace;
subspace.S.S = Sspace;
subspace.S.n = Nspace;
subspace.trial = Ae;

fprintf('Applying the subspace projector\n');
%Bclean    = B - (B*Ae)*Ae'; % avoid computation of Ae*Ae'
Bclean = B - (B*cellfun(@transpose, Ae, 'UniformOutput', false))*Ae;

function [Ae, subspace] = CSP01(Bin, Bout, Nin, Nout, Nintersect, pertrial)
%
% interference rejection by removing the common temporal subspace of the two subspaces
% K. Sekihara,  March 28, 2012
% Golub and Van Loan, Matrix computations, The Johns Hopkins University Press, 1996
%
%  Nc: number of channels
%  Nt: number of time points
% inputs
%  Bout(1:Nc,1:Nt): interference data
%  Bin(1:Nc,1:Nt): signal plus interference data
%  Nout: dimension of the interference subspace
%  Nin: dimension of the signal plus interference subspace
%  Nintersect: dimension of the intersection of the two subspaces
% outputs
% Ae = matrix from which the projector onto the intersection can
%      be obtained:
% subspace: struct containing information about the different subspace
%      projections
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%

if ~pertrial
  % compute the projection across trials
  trllist = 1:numel(Bout);
else
  % compute the projection per trial
  trllist = (1:numel(Bout))';
end  

Ae = cell(size(Bin));
for k = 1:size(trllist,1)
  indx = trllist(k,:); % this is either a scalar, or a vector
  [Uout,Sout,Vout] = svd(cat(2, Bout{indx}),'econ');
  [Uin, Sin, Vin]  = svd(cat(2, Bin{indx}), 'econ');
  Sout = diag(Sout);
  Sin  = diag(Sin);

  Nout = getN(Nout, Sout, 'outside');
  Nin  = getN(Nin,  Sin,  'inside');
  
  % compute unit-norm orthogonal time courses
  Qout = diag(1./Sout(1:Nout))*Uout(:,1:Nout)'*Bout(indx); % keep it in cell representation
  Qin  = diag(1./Sin(1:Nin)  )* Uin(:,1:Nin)' *Bin(indx);
  C    = Qin * cellfun(@transpose, Qout, 'UniformOutput', false);
  C    = sum(cat(3, C{:}), 3);

  % store the subspace information that is used in the next step
  subspace.Sin(k).U  = Uin;
  subspace.Sin(k).S  = Sin;
  subspace.Sin(k).n  = Nin;
  subspace.Sout(k).U = Uout;
  subspace.Sout(k).S = Sout;
  subspace.Sout(k).n = Nout;

  % covariance matrix of unit-norm 'components' -> how does this relate to
  % multivariate decomp? This is I guess equivalent mathematically
  [U,S] = svd(C);
  S     = diag(S);
  Nintersect = getN(Nintersect, S, 'intersection');
  
  Ae(indx) = U(:, 1:Nintersect)'*Qin;

  % keep the subspace information
  subspace.T(k).U = U;
  subspace.T(k).S = S;
  subspace.T(k).C = C; clear C;
  subspace.T(k).n = Nintersect;
end

function N = getN(N, S, name)

ttext = sprintf('enter the dimension for the %s field: ', name);
if isempty(N)
  N  = input(ttext);
elseif ischar(N) && isequal(N, 'interactive') && ~any(strcmp(name, {'outside' 'intersection'}))
  figure, plot(log10(S),'-o'); drawnow
  N = input(ttext);
elseif ischar(N) && isequal(N, 'interactive') && any(strcmp(name, {'outside' 'intersection'}))
  figure, plot(S, '-o'); drawnow
  N = input(ttext);
elseif ischar(N) && isequal(N, 'all')
  N = find(S./S(1)>1e5*eps, 1, 'last');
elseif isnumeric(N) && N<1
  N = find(S<=N, 1, 'last');
end
fprintf('Using %d dimensions for the %s field\n', N, name);
