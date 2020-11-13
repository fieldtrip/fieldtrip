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
%   cfg.sourcemodel      = structure, source model with precomputed leadfields, see FT_PREPARE_LEADFIELD
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check the input data
datain = ft_checkdata(datain, 'datatype', {'raw'}); % FIXME how about timelock and freq?

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});

% get the options
cfg.trials       = ft_getopt(cfg, 'trials',  'all', 1);
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.sourcemodel  = ft_getopt(cfg, 'sourcemodel');
cfg.dssp         = ft_getopt(cfg, 'dssp');         % sub-structure to hold the parameters
cfg.dssp.n_space = ft_getopt(cfg.dssp, 'n_space', 'all'); % number of spatial components to retain from the Gram matrix
cfg.dssp.n_in    = ft_getopt(cfg.dssp, 'n_in', 'all');    % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.dssp.n_out   = ft_getopt(cfg.dssp, 'n_out', 'all');   % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.dssp.n_intersect = ft_getopt(cfg.dssp, 'n_intersect', 0.9); % dimensionality of the intersection
cfg.output       = ft_getopt(cfg, 'output', 'original');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);

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

dat     = cat(2,datain.trial{:});
[dum, Ae, N, Nspace, Sout, Sin, Sspace, S] = dssp(dat, G, cfg.dssp.n_in, cfg.dssp.n_out, cfg.dssp.n_space, cfg.dssp.n_intersect);
datAe   = dat*Ae; % the projection is a right multiplication
% with a matrix (eye(size(Ae,1))-Ae*Ae'), since Ae*Ae' can become quite
% sizeable, it's computed slightly differently here.

% put some diagnostic information in the output cfg.
cfg.dssp.S_space        = Sspace;
cfg.dssp.n_space        = Nspace;
cfg.dssp.S_out          = Sout;
cfg.dssp.S_in           = Sin;
cfg.dssp.S_intersect    = S;
cfg.dssp.n_intersect    = N;

% compute the cleaned data and put in a cell-array
nsmp  = cellfun(@numel, datain.time);
csmp  = cumsum([0 nsmp]);
trial = cell(size(datain.trial));
switch cfg.output
  case 'original'
    for k = 1:numel(datain.trial)
      trial{k} = datain.trial{k} - datAe*Ae((csmp(k)+1):csmp(k+1),:)';
    end
  case 'complement'
    for k = 1:numel(datain.trial)
      trial{k} = datAe*Ae((csmp(k)+1):csmp(k+1),:)';
    end
  otherwise
    ft_error(sprintf('cfg.output = ''%s'' is not implemented',cfg.output));
end

% create the output argument
dataout       = keepfields(datain, {'label', 'time', 'fsample', 'trialinfo', 'sampleinfo', 'grad', 'elec', 'opto'}); % grad can be kept and does not need to be balanced, since the cleaned data is a mixture over time, not space.
dataout.trial = trial;

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrix
% kindly provided by Kensuke, and adjusted a bit by Jan-Mathijs
function [Bclean, Ae, Nee, Nspace, Sout, Sin, Sspace, S] = dssp(B, G, Nin, Nout, Nspace, Nee)

% Nc: number of sensors
% Nt: number of time points
% inputs
% B(Nc,Nt):  interference overlapped sensor data
% G(Nc,Nc): Gram matrix of voxel lead field
% Nout and Nin: dimensions of the two row spaces
% recom_Nspace: recommended value for the dimension of the pseudo-signal subspace
% outputs
% Bclean(Nc,Nt): cleaned sensor data
% Nee: dimension of the intersection
% Nspace: dimension of the pseudo-signal subspace
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%
% The code below is modified by Jan-Mathijs, no functional changes
% merely cosmetics

% eigen decomposition of the Gram matrix, matrix describing the spatial
% components
[U,S]   = eig(G);
Sspace  = abs(diag(S));

[Sspace, iorder] = sort(-Sspace);
Sspace           = -Sspace;
U(:,:)           = U(:,iorder);

if isempty(Nspace)
  ttext = 'enter the spatial dimension: ';
  Nspace    = input(ttext);
elseif ischar(Nspace) && isequal(Nspace, 'interactive')
  figure, plot(log10(Sspace),'-o');
  Nspace = input('enter spatial dimension of the ROI subspace: ');
elseif ischar(Nspace) && isequal(Nspace, 'all')
  Nspace = find(Sspace./Sspace(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions\n', Nspace);

% spatial subspace projector
Us   = U(:,1:Nspace);
USUS = Us*Us';

% Bin and Bout creations
Bin  =                  USUS  * B;
Bout = (eye(size(USUS))-USUS) * B;

% create the temporal subspace projector and apply it to the data
%[AeAe, Nee] = CSP01(Bin, Bout, Nout, Nin, Nee);
%Bclean      = B*(eye(size(AeAe))-AeAe);

[Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee);
Bclean    = B - (B*Ae)*Ae'; % avoid computation of Ae*Ae'


function [Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee)
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
%  ypost(1:Nc,1:Nt): denoised data
%  Nout: dimension of the interference subspace
%  Nin: dimension of the signal plus interference subspace
%  Nee: dimension of the intersection of the two subspaces
% outputs
% Ae = matrix from which the projector onto the intersection can
%      be obtained:
% AeAe: projector onto the intersection, which is equal to the
%       interference subspace.
% Nee: dimension of the intersection
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%

[dum,Sout,Vout] = svd(Bout,'econ');
[dum,Sin, Vin]  = svd(Bin, 'econ');
Sout = diag(Sout);
Sin  = diag(Sin);

if isempty(Nout)
  ttext = 'enter the spatial dimension for the outside field: ';
  Nout  = input(ttext);
elseif ischar(Nout) && isequal(Nout, 'interactive')
  figure, plot(Sout,'-o');
  Nout = input('enter dimension of the outside field: ');
elseif ischar(Nout) && isequal(Nout, 'all')
  Nout = find(Sout./Sout(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions for the outside field\n', Nout);

if isempty(Nin)
  ttext = 'enter the spatial dimension for the inside field: ';
  Nin  = input(ttext);
elseif ischar(Nin) && isequal(Nin, 'interactive')
  figure, plot(log10(Sin),'-o');
  Nin = input('enter dimension of the inside field: ');
elseif ischar(Nin) && isequal(Nin, 'all')
  Nin = find(Sin./Sin(1)>1e5*eps, 1, 'last');
end
fprintf('Using %d spatial dimensions for the inside field\n', Nin);

Qout = Vout(:,1:Nout);
Qin  = Vin(:, 1:Nin);

C     = Qin'*Qout;
[U,S] = svd(C);
S     = diag(S);
if (ischar(Nee) && strcmp(Nee, 'auto'))
  ft_error('automatic determination of intersection dimension is not supported');
elseif ischar(Nee) && strcmp(Nee, 'interactive')
  figure, plot(S,'-o');
  Nee  = input('enter dimension of the intersection: ');
elseif Nee<1
  % treat a numeric value < 1 as a threshold
  Nee = find(S>Nee,1,'last');
  if isempty(Nee), Nee = 1; end
end
fprintf('Using %d dimensions for the interaction\n', Nee);

Ae   = Qin*U;
Ae   = Ae(:,1:Nee);
%AeAe = Ae*Ae';
