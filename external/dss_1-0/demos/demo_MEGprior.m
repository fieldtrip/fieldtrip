% This script demonstrates the use of the DSS MATLAB package in
% separation of MEG signals using prior information such as 
% reference signals (EOG, ECG) and visual inspection of the blind results
% 
%
%   See denss for description of the denoising source separation 
%   MATLAB package, see demo_MEGblind for the blind approach

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: demo_MEGprior.m,v 1.4 2005/12/07 14:04:04 jaakkos Exp $

load MEG_data.mat;                        % load the data

Xaux = X(123:end,:);                      % + 5 channels for EOG, etc.
Xaux(end,:) = -Xaux(end,:);               % make R-peaks positive
X = X(1:122,:);                           % exclude the reference channels from the data

% plot some channels and the references
% plot the results
figure
plots = 10;
for i = 1 : plots-3
  subplot(plots,1,i);
  plot(X(i,:));
end
subplot(plots,1,i+1); plot(Xaux(1,:));    % blinks
subplot(plots,1,i+2); plot(Xaux(2,:));    % saccades
subplot(plots,1,plots); plot(Xaux(5,:));  % ECG

state = dss_create_state(X)               % create the state structure 
% see dss_create_state for the definition of the different fields 
state.verbose = 3;                        % maximal information output
state.wdim = 50;                          % drop dimension to avoid overfitting
state = dss_preprocess(state.X,state);    % let's whiten for further use

% first the mascular artefacts around time indices 10000-13000
state = dss_set_denoising(state, @denoise_mask);   % mask-based denoising
% Now the following denoising function parameter has to be set: state.denf.params.mask

mask = zeros(1,length(state.X));      
mask(10000:13000)=1;                      % creating the on/off mask
state.denf.params.mask = mask;            % setting the mask for the denoising
state.algorithm = 'pca';                  % PCA approach (most efficient)
state.sdim = 5;                           % extract 5 components
state = denss(state);                     % run DSS

% add some labels (not used by DSS)
for i = 1 : state.sdim
  state.Slabels{i} = sprintf('muscular %d',i);
end
% then the ocular artefacts: let's use the EOG channels to generate the mask
for i = 1 : 2                             % blinks and saccades are Xaux(1:2,:)
  mask(i,:) = estimate_mask(Xaux(i,:));   % estimation of the mask
end

figure
for i = 1 : 2                             % plot the masks
  subplot(2,1,i);
  plot([Xaux(i,:) / max(Xaux(i,:));mask(i,:)]');
end

for i = 1 : 2
  state.denf.params.mask = mask(i,:);     % set the mask
  state.sdim = state.sdim + 1;            % extract one more component
  state.component = state.sdim;           % field has to be added because PCA DSS does not use it
  state.algorithm = 'defl';               % PCA algorithm cannot be used, use deflationary
  state = denss(state);                   % run DSS
end
state.Slabels{state.sdim-1}='blinks';
state.Slabels{state.sdim}='saccades';

% extraction of the cardiac: using ECG and ensemble averaging
trsig = (Xaux(5,:)>3*std(Xaux(5,:)));                        % simple peak detection
trsigdiff = trsig(2:end)-trsig(1:end-1);
triggers = find(trsigdiff==1)+1;                             % the triggers as indices
triggers = triggers(1:end-1);                                % last QRS close to the end of data
medianlength = median(triggers(2:end)-triggers(1:end-1));    % median length of the cycle
denf_params.tr = triggers(1:end-1);                          % set the denf parameters
denf_params.tr_begin = 0;
denf_params.tr_end = medianlength;
state = dss_set_denoising(state, @denoise_avg, denf_params); % set the denoising, parameters given too
state.sdim = state.sdim + 1;                                 % extract one more component

state = denss(state);                                        % run DSS
state.Slabels{state.sdim}='cardiac';

% plot the results
figure
for i = 1 : state.sdim
  subplot(state.sdim,1,i);
  plot(state.S(i,:));
  ylabel(state.Slabels{i});
end
