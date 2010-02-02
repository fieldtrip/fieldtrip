
%currently assumes that keep trials and keep tapers 
%data is in channels by samples
function [out] = specest_mtmfft(data, fsample, varargin) %maybe this should be called tapered fft as it can take any taper???

state         = keyval('state',         varargin);    if isempty(state),      stae       = [];  end
pad           = keyval('pad',           varargin);    if isempty(pad),        pad        = 'maxperlen'; end  %keep using maxperlen?
foilim        = keyval('foilim',        varargin);    if isempty(foilim),     foilim     = [0 fsample/2]; end
output        = keyval('output',        varargin);    if isempty(output),     output     = 'pow'; end
taper         = keyval('taper',         varargin);    if isempty(taper),      taper      = 'dpss'; end
tapsmofrq     = keyval('tapsmofrq',     varargin);    if isempty(tapsmofrq),  error('must specify number of tapers to use'); end

%only power spectrum for now
if isequal(output, 'pow')
    powflg = 1;
    fftflg = 0;
elseif isequal(output, 'fourier') 
    powflg = 0;
    fftflg = 1;
end

%we only have a single trial!!! Check the state variable to see if the
%current trial is the same length as the previous trial.

numdatbns = size(data,2);

% if cfg.pad is 'maxperlen', this is realized here:
if isequal(pad, 'maxperlen')
  pad = numdatbns./ fsample;
else
  % check that the specified padding is not too short
  if pad<(numdatbns/fsample)
    error('the padding that you specified is shorter than the longest trial in the data');
  end
end
numsmp = pad * fsample; % this used to be "cfg.pad .* data.fsample"
numsgn = size(data,1);
% doing the computation
boilim  = round(foilim ./ (fsample ./ numsmp)) + 1;
boi     = boilim(1):boilim(2);
numboi  = size(boi,2);
foi     = (boi-1) ./ pad;

if isequal(output, 'pow'),    powspctrm     = zeros(numsgn,numboi);   end
if isequal(output,'fourier'), fourierspctrm = complex(zeros(numsgn,numboi)); end

% trials are of equal length, compute the set of tapers only once . check
% the current number of data bins against the state variable!!!

if strcmp(taper, 'dpss')
    % create a sequence of DPSS (Slepian) tapers
    % ensure that the input arguments are double precision
    tap = double_dpss(numdatbns,numdatbns*(tapsmofrq./fsample))';
    elseif strcmp(taper, 'sine')
    tap = sine_taper(numdatbns, numdatbns*(tapsmofrq./fsample))';
    else
    % create a single taper according to the window specification as a
    % replacement for the DPSS (Slepian) sequence
    tap = window(cfg.taper, numdatbns)';
    tap = tap./norm(tap);
    % this function always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(2,:) = nan;
end

numtap = size(tap,1) - 1;

if (numtap < 1)
error(sprintf(...
  'datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
  numdatbns/fsample, tapsmofrq, fsample/numdatbns));
elseif (numtap < 2) && strcmp(taper, 'dpss')
fprintf('WARNING: using only one taper for specified smoothing\n');
end

pad = zeros(1,numsmp - numdatbns);
cumsumcnt = numdatbns; %???
numtap = size(tap,1) - 1;
cumtapcnt = numtap;


for taplop = 1:numtap

    autspctrmacttap = complex(zeros(numsgn,numboi));
    for sgnlop = 1:numsgn
        dum = fft([data(sgnlop,:) .* tap(taplop,:) , pad],[],2);
        autspctrmacttap(sgnlop,:) = dum(boi);
    end
    if taplop == 1
        fprintf('nfft: %d samples, taper length: %d samples, %d tapers\n',length(dum),size(tap,2),numtap);
    end
    if powflg
        powdum = 2 .* (autspctrmacttap .* conj(autspctrmacttap)) ./ numsmp; %cf Numercial Receipes 13.4.9
        powspctrm(:,:) = powdum;
    end
    if fftflg
        fourierdum = (autspctrmacttap) .* sqrt(2 ./ numsmp); %cf Numercial Receipes 13.4.9
        fourierspctrm(:,:) = fourierdum;
    end

end % taplop

if powflg, out = powspctrm; end
if fftflg, out = fourierspctrm; end
%add crsspctrm at later point

state.cumsumcnt = cumsumcnt;
stae.cumtapcnt  = cumtapcnt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin);
tap = dpss(double(a), double(b), varargin{:});





%     if csdflg
%       csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .* ...
%         conj(autspctrmacttap(cutdatindcmb(:,2),:))) ./ numsmp;
%       if keep == 1
%         crsspctrm(:,:) = crsspctrm(:,:) + csddum ./ numtap;
%       elseif keep == 2
%         crsspctrm(cnt,:,:) = crsspctrm(cnt,:,:) + permute(csddum,[3,1,2]) ./ numtap;
%       elseif keep == 4
%         crsspctrm(cnt,:,:) = csddum;
%       end
%     end

