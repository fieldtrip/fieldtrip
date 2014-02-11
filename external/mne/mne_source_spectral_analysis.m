function [res] = mne_source_spectral_analysis(fname_rawdata, cfg)
%
% [res] = mne_source_spectral_analysis(fname_rawdata, fname_output, cfg)
%
% Estimate frequency spectra in the source space and optionally write out
% stc files which have frequencies along the "time" axis.
%
% fname_data     - Name of the data file
%
% MNE inversion
% cfg.inv        - Inverse operator structure or file name
% cfg.lambda2    - The regularization factor
% cfg.dSPM       - enable dSPM: 0 or 1
%
% Spectral estimation
% cfg.mode       - output quantity; 'amplitude', 'power', 'phase'
% cfg.window     - window type: 'hanning', 'hamming' or any other window
%                  function available on Matlab
% cfg.fft_length - FFT length in samples (half-overlapping windows used)
% cfg.foi        - Frequencies of interest as [f_min f_max]
%
% Output
% cfg.outfile    - The stem of the output STC file name holding the spectra
%
% (C)opyright Lauri Parkkonen, 2012
%
% $Log$
%

me='MNE:mne_spectral_analysis';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 2
   error(me,'Incorrect number of arguments'); 
end

% Set defaults
if ~isfield(cfg, 'lambda2')
    cfg.lambda2 = 1.0;
end
if ~isfield(cfg, 'dSPM')
    cfg.dSPM = 0;
end
if ~isfield(cfg, 'mode')
   cfg.mode = 'amplitude';
end
if ~isfield(cfg, 'window')
    cfg.window = 'hanning';
end
if ~isfield(cfg, 'fft_length')
    cfg.fft_length = 1024;
end

% Set up raw data input
rawfile = fiff_setup_read_raw(fname_rawdata);
if ~isfield(cfg, 'first_samp') || cfg.first_samp <= 0
    cfg.first_samp = rawfile.first_samp;
end
if ~isfield(cfg, 'last_samp') || cfg.last_samp == Inf
    cfg.last_samp = rawfile.last_samp;
end

% Make a selector for frequencies
faxis = ((1:cfg.fft_length/2) - 1) / cfg.fft_length * rawfile.info.sfreq;
if ~isfield(cfg, 'foi')
    cfg.foi = [faxis(1) faxis(end)];
end
foi_i = find(faxis >= cfg.foi(1) & faxis <= cfg.foi(2));
fprintf(1, 'Considering frequencies %g ... %g Hz\n', faxis(foi_i(1)), faxis(foi_i(end)));

% Compose the inverse operator
if isstruct(cfg.inv)
    inv = cfg.inv;
else
    inv = mne_read_inverse_operator(cfg.inv);
end
inv = mne_prepare_inverse_operator(inv, 1, cfg.lambda2, cfg.dSPM);
chsel = fiff_pick_channels(rawfile.info.ch_names, inv.noise_cov.names);
fprintf(1,'Using %d channels for the inverse solution\n', length(chsel));
fprintf(1,'Computing inverse...');
inverse_noleads = diag(sparse(inv.reginv))*inv.eigen_fields.data*inv.whitener*inv.proj;
if inv.eigen_leads_weighted
   %
   %     R^0.5 has been already factored in
   %
   fprintf(1,'(eigenleads already weighted)...');
   inverse = inv.eigen_leads.data*inverse_noleads;
else
   %
   %     R^0.5 has to factored in
   %
   fprintf(1,'(eigenleads need to be weighted)...');
   inverse = diag(sparse(sqrt(inv.source_cov.data)))*inv.eigen_leads.data*inverse_noleads;
end
fprintf(1, '[done]\n');

% Compute the average spectra
ave = [];
nspectra = 0;
windowfun =feval(cfg.window, cfg.fft_length);
for w = cfg.first_samp : cfg.fft_length/2 : cfg.last_samp
    data_sensor = fiff_read_raw_segment(rawfile, w, w + cfg.fft_length - 1, chsel);
    % Check for an incomplete last segment
    if size(data_sensor, 2) ~= cfg.fft_length
        fprintf(1, 'Skipping the remaining incomplete window\n');
        break;
    end
    nspectra = nspectra + 1;
    for ch = 1:size(data_sensor, 1)
        data_sensor(ch, :) = data_sensor(ch, :) .* windowfun';
    end
    data_fft = fft(data_sensor, cfg.fft_length, 2);
    data_fft = data_fft(:, foi_i);
    data_fft_source = inverse * double(data_fft);
    if isempty(ave)
        ave = zeros(size(data_fft_source));
    end
    ave = ave + abs(data_fft_source);
    fprintf('Window %d, %3.3g%% completed\n', nspectra, double(w - cfg.first_samp + cfg.fft_length/2) / double(cfg.last_samp - cfg.first_samp) * 100);
end
ave = ave / nspectra;

% Combine the three source orientations if needed
if inv.source_ori == FIFF.FIFFV_MNE_FREE_ORI
    % Do not do this for phase maps; just pick the normal component
    fprintf(1,'Combining the spectra of the current components\n');
    ave1 = zeros(size(ave,1)/3, size(ave,2));
    for k = 1:size(ave,2)
        ave1(:,k) = sqrt(mne_combine_xyz(ave(:,k)));
    end
    ave = ave1;
end

% Do dSPM weighting if requested
if cfg.dSPM
    fprintf(1,'Applying dSPM weighting\n');
    ave = inv.noisenorm * ave;
end

% Compose the result structure
res.inv     = inv;
res.rawfile = rawfile;
res.spectra = ave;
res.fmin    = faxis(foi_i(1));
res.fmax    = faxis(foi_i(end));
res.fstep   = faxis(2) - faxis(1);
res.faxis   = faxis;

% Write out the STC file
if isfield(cfg, 'outfile')
    fprintf('Writing out STC files...\n');
    mne_write_inverse_sol_stc(cfg.outfile, res.inv, res.spectra, 1e-3 * res.fmin, 1e-3 * res.fstep);
end

fprintf(1,'[done]\n');

return;
end