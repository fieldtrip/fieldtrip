function [out,foi] = specest_mtmfft(dat, fsample, varargin) 

% SPECEST_MTMFFT computes a fast Fourier transform using many possible tapers
%
%
% Use as
%   [out,foi] = specest_mtmfft(dat,fsample...)       %%%% EITHER FSAMPLE OR TIME
%
%   dat      = matrix of chan*sample 
%   fsample  = number indicating sampling rate
%   out      = matrix of taper*chan*sample of fourier coefficients
%   foi      = vector of frequencies in out
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total number of samples after zero padding               %%% COULD ALSO BE GIVEN IN SECONDS....
%   freq       = vector, containing frequencies                                               
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%
%
%
%
%
%  TO DO:
%  - anchor phase to specific part of time-window (t=0 current suggestion), how to do it in the current format? need to do multiplication freq-domain (ala mtmconvol)?
%  - implement computation reduction by keeping tapers and such one way or another
%
%
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'dpss','pad','freq','tapsmofrq'});
taper     = keyval('dpss',        varargin); if isempty(taper),    taper   = 'dpss';     end
pad       = keyval('pad',         varargin); if isempty(pad),      pad     = 0;          end
freq      = keyval('freq',        varargin); if isempty(freq),     freq    = 'max';      end  
tapsmofrq = keyval('tapsmofrq',   varargin);


% zero padding
zeropad = zeros(1,pad-size(dat,2));

% set n's
[nchan,nsample] = size(dat);

% set boi and foi 
if isnumeric(freq) % if input is a vector
boi    = round(freq ./ (fsample ./ pad)) + 1;
numboi = size(boi,2);
foi    = (boi-1) ./ (pad / fsample); % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(freq,'max') % if input was 'max'
  boilim = round([0 fsample/2] ./ (fsample ./ pad)) + 1;
  boi    = boilim(1):boilim(2);
  numboi = size(boi,2);
  foi    = (boi-1) ./ (pad / fsample);
end

% create tapers
switch taper
  case 'dpss'
    % create a sequence of DPSS tapers, ensure that the input arguments are double precision
    tap = double_dpss(nsample,nsample*(tapsmofrq./fsample))';
    % remove the last taper
    tap = tap(1:(end-1), :);
    
    % give error/warning about number of tapers
    if isempty(tap)
      error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',nsample/fsample,tapsmofrq,fsample/fsample);
    elseif size(tap,1) == 1
      warning('using only one taper for specified smoothing')
    end
        
  case 'sine'
    tap = sine_taper(nsample, nsample*(tapsmofrq./fsample))';
    
  case 'alpha'
    error('not yet implemented');
    
  otherwise
    % create the taper and ensure that it is normalized
    tap = window(taper, nsample)';
    tap = tap ./ norm(tap);
end % switch taper
numtap = size(tap,1);


% compute fft per channel, keeping tapers automatically (per channel is about 40% faster than all channels at the same time)
out = [];
for itap = 1:numtap
  for ichan = 1:nchan
    dum = fft([dat(ichan,:) .* tap(itap,:) zeropad],[],2); % would be much better if fft could take boi as input (muuuuuch less computation)
    out(itap,ichan,:) = dum(boi);
  end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});








%%%%%%%%%% OLD CODE %%%%%%%%%%


% %maybe this should be called tapered fft as it can take any taper???
% 
% %currently assumes that keep trials and keep tapers
% %data is in channels by samples
% 
% state         = keyval('state',         varargin);    if isempty(state),      state       = [];  end
% pad           = keyval('pad',           varargin);    if isempty(pad),        pad        = 'maxperlen'; end  %keep using maxperlen?
% foilim        = keyval('foilim',        varargin);    if isempty(foilim),     foilim     = [0 fsample/2]; end
% output        = keyval('output',        varargin);    if isempty(output),     output     = 'pow'; end
% taper         = keyval('taper',         varargin);    if isempty(taper),      taper      = 'dpss'; end
% tapsmofrq     = keyval('tapsmofrq',     varargin);    if isempty(tapsmofrq),  error('must specify number of tapers to use'); end
% 
% [nchan, numdatbns] = size(dat);
% 
% % sheck whether the state can be reused
% options = varargin;
% if isfield(state.options) && ~isequal(state.options, options)
%   state = [];
% end
% 
% % set up the state
% if isempty(state)
% 
%   switch taper
%     case 'dpss'
%       % create a sequence of DPSS (Slepian) tapers
%       % ensure that the input arguments are double precision
%       tap = double_dpss(numdatbns,numdatbns*(tapsmofrq./fsample))';
%       % remove the last taper
%       tap = tap(1:(end-1), :);
% 
%     case 'sine'
%       tap = sine_taper(numdatbns, numdatbns*(tapsmofrq./fsample))';
% 
%     case 'alpha'
%       error('not yet implemented');
% 
%     otherwise
%       % create the taper and ensure that it is normalized
%       tap = window(taper, nsample);
%       tap = tap ./ norm(tap);
%   end % switch taper
% 
%   pad = zeros(nchan, padding-numdatbns);
% 
% else
%   tap = state.tap;
%   pad = state.pad;
% end % if previous state applies
% 
% 
% 
% numtap = size(tap,1);
% 
% if (numtap < 1)
%   error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz', numdatbns/fsample, tapsmofrq, fsample/numdatbns));
% elseif (numtap < 2) && strcmp(taper, 'dpss')
%   fprintf('WARNING: using only one taper for specified smoothing\n');
% end
% 
% numsmp = pad * fsample; % this used to be "cfg.pad .* data.fsample"
% numsgn = size(data,1);
% % doing the computation
% boilim  = round(foilim ./ (fsample ./ numsmp)) + 1;
% boi     = boilim(1):boilim(2);
% numboi  = size(boi,2);
% foi     = (boi-1) ./ pad;
% 
% % determine the time and frequency resolution
% dt = 1 ./ fsample;
% df = 1 ./ (nsample+npad)/fsample;
% 
% time = ...;
% freq = (1:(nsample+npad)) * df - df;
% 
% % trials are of equal length, compute the set of tapers only once . check
% % the current number of data bins against the state variable!!!
% 
% pad = zeros(1,numsmp - numdatbns);
% cumsumcnt = numdatbns; %???
% numtap = size(tap,1) - 1;
% cumtapcnt = numtap;
% 
% % pre-allocate memory that will contain the result
% spectrum = complex(zeros(numtap,numsgn,numboi));
% 
% for taplop = 1:numtap
% 
%   for sgnlop = 1:numsgn
%     dum = fft([data(sgnlop,:) .* tap(taplop,:), pad], [], 2);
%     spectrum(taplop,sgnlop,:) = dum(boi);
%   end
% 
%   if taplop == 1
%     fprintf('nfft: %d samples, taper length: %d samples, %d tapers\n',length(dum),size(tap,2),numtap);
%   end
% 
% end % taplop
% 
% % remember the state for the next call
% if isempty(state)
%   state.options = options;
%   state.tap     = tap;
%   state.pad     = pad;
%   state.cumsumcnt = cumsumcnt;
%   state.cumtapcnt = cumtapcnt;
% end






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

