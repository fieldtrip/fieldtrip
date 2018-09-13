function [dataout] = ft_denoise_dssp(cfg, datain)

% FT_DENOISE_DSSP 
%
% Use as 
%   dataout = ft_denoise_dssp(cfg, datain)
% where cfg is a configuration structure that contains
%   cfg.grid     = structure, source model with precomputed leadfields (see FT_PREPARE_LEADFIELD)
%   cfg.order    = integer, order of the source space
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

% get the options
cfg.grid = ft_getopt(cfg, 'grid');
cfg.dssp = ft_getopt(cfg, 'dssp'); %sub-structure to hold the parameters
cfg.dssp.n_space = ft_getopt(cfg.dssp, 'n_space'); % number of spatial components to retain from the Gram matrix
cfg.dssp.n_in    = ft_getopt(cfg.dssp, 'n_in');    % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.dssp.n_out   = ft_getopt(cfg.dssp, 'n_out');   % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.dssp.n_intersect = ft_getopt(cfg.dssp, 'n_intersect', 0.99); % dimensionality of the intersection

% check the input data
datain = ft_checkdata(datain, 'datatype', {'raw' 'timelock'});

% match the input data's channels with the labels in the leadfield
% TODO, and check whether the grid has leadfields

grid = cfg.grid;

% compute the Gram-matrix for the supplied forward model
lf = cat(2, grid.leadfield{:});
G  = lf*lf';

trial = cell(size(datain.trial));
for k = 1:numel(datain.trial)
  trial{k} = dssp(datain.trial{k}, G, cfg.dssp.n_in, cfg.dssp.n_out, cfg.dssp.n_space, cfg.dssp.n_intersect);
end

% do your stuff...
dataout = keepfields(datain, {'label','time','fsample','trialinfo','sampleinfo','grad', 'elec', 'opto'}); % grad can be kept and does not need to be balanced, since the cleaned data is a mixture over time, not space.
dataout.trial = trial;



ft_postamble debug              
ft_postamble trackconfig        
ft_postamble previous   datain  
ft_postamble provenance dataout 
ft_postamble history    dataout 
ft_postamble savevar    dataout 

function [yclean, Nee, mu] = dssp(yy, G, Nec, Net, mu, Nee)
%
% Nc: number of sensors
% Nt: number of time points
% inputs
% yy(Nc,Nt):  interference overlapped sensor data
% G(Nc,Nc): Gram matrix of voxel lead field
% Nec and Net: dimensions of the two row spaces
% recom_mu: recommended value for the dimension of the pseudo-signal subspace
% outputs
% yclean(Nc,Nt): cleaned sensor data
% Nee: dimension of the intersection 
% mu: dimension of the pseudo-signal subspace
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
eig_els = abs(diag(S));

[eig_el, iorder] = sort(-eig_els);
eig_el           = -eig_el;
U(:,:)           = U(:,iorder);

ev_max = length(eig_el);
figure(61)
semilogy(1:ev_max,eig_el,'-o');
axis([1, ev_max, 0, eig_el(1)*1.1]);
title('eigen value spectrum of the gram matrix')
xlabel('order of eigenvalues')
ylabel('relative magnitude')

if isempty(mu)
  ttext = 'enter the dimension: ';
  mu    = input(ttext);
end

% spatial subspace projector
Us = U(:,1:mu);  
USUS=Us*Us';

% Bin and Bout creations
yin  =                  USUS  * yy;
yout = (eye(size(USUS))-USUS) * yy;

% create the temporal subspace projector and apply it to the data
[AeAe, Nee] = CSP01(yin, yout, Nec, Net, Nee);
yclean      = yy*(eye(size(AeAe))-AeAe);


function [AeAe, Nee] = CSP01(ysig, yint, Net, Nec, Nee)
%
% interference rejection by removing the common temporal subspace of the two subspaces
% K. Sekihara,  March 28, 2012
% Golub and Van Loan, Matrix computations, The Johns Hopkins University Press, 1996
% 
%  Nc: number of channels
%  Nt: number of time points
% inputs
%  yint(1:Nc,1:Nt): interference data
%  ysig(1:Nc,1:Nt): signal plus interference data
%  ypost(1:Nc,1:Nt): denoised data
%  Nec: dimension of the interference subspace
%  Net: dimension of the signal plus interference subspace
%  Nee: dimension of the intersection of the two subspaces
% outputs
% AeAe: projector onto the intersection, which is equal to the interference subspace.
% Nee: dimension of the intersection
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%

[~,~,Vc] = svd(yint,'econ');
[~,~,Vt] = svd(ysig,'econ');

Qc = Vc(:,1:Nec);
Qt = Vt(:,1:Net);

C       = Qt'*Qc;
[U,S,~] = svd(C);
eivs    = diag(S);
if (ischar(Nee) && strcmp(Nee, 'auto'))
  keyboard
elseif ischar(Nee) && strcmp(Nee, 'interactive')
  figure, plot(eivs,'-o');
  Nee  = input('enter dimension of the intersection: ');
elseif Nee<1
  % treat a numeric value < 1 as a threshold
  Nee = find(eivs>Nee,1,'last');
end
A1   = Qt*U;
Ae   = A1(:,1:Nee);
AeAe = Ae*Ae';
