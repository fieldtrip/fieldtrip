function [data_jittered] = ft_phasejitter(cfg, data)

% FT_PHASEJITTER Introduces phase jitter to some or all frequencies in
% raw timeseries data. This is useful for creating null data sets with the
% same spectral profile as the the original data. 

%
%
% The cfg.domain options determine which domains to shuffle along

% cfg.jittertype  = 'uniform' or 'vonmises'. Determines shape of noise
%                    distribution (default = 'uniform'). If using 'vonmises', 
%                    the following parameters also need ot be set. 
%                    Both paramaters can be single values or a range of values 
%                    for each frequency in cfg.foi or range in cfg.foilim.
%
% cfg.jittermean  = the mean of the phase jitter distribution in radians
%                    (default = 0). 
% cfg.jitterstd   = the standard deviation of the phase jitter distribution 
%                    radians (default = pi/4). NB Setting cfg.jitterstd = inf 
%                    is the same as setting cfg.jittertype = 'uniform'. 
% 
% frequency range of interst can be set using one of the following
% paramaters:
% cfg.foi          = a list of frequencies (default = 'all')
% cfg.foilim       = a 2-valued frequency range.


revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

% set defaults
cfg.jittertype            = ft_getopt(cfg, 'jittertype',    'uniform');
cfg.jittermean               = ft_getopt(cfg, 'jittermean',   0);
cfg.jitterstd           = ft_getopt(cfg, 'jitterstd',   pi/4);


% check freq paramaters
if isfield(cfg,'foi') && isfield(cfg,'foilim')
     warning('Both cfg.foi or cfg.foilim have been specified. Using cfg.foi.');
     cfg.foilim=[];
elseif isfield(cfg,'foi')
     cfg.foilim=[];
elseif isfield(cfg,'foilim')
    cfg.foi=[];
else
    cfg.foi='all';
end

% check von mises paramaters are ok
if ~strcmp(cfg.jittertype,'uniform') && ~strcmp(cfg.jittertype,'vonmises')
    error('Unknown jitter type specified. Valid options are ''uniform'' and ''vonmises''.');
end
if strcmp(cfg.jittertype,'vonmises')
    m=cfg.jittermean;
    s=cfg.jitterstd;
    if length(m)~=length(s)
        error('jittermean and jitterstd should be same size');
    end
elseif strcmp(cfg.jittertype,'uniform')
    m=0;
    s=inf;
end




% get frequency scale from data and match fois
fsample=data.fsample;

n_samples=length(data.time{1});
endtime=n_samples./fsample;

if strcmp(cfg.foi,'all');
    cfg.foilim=[0 fsample/2];
    cfg.foi=[];
end

if ~isempty(cfg.foi) % if input is foi
  freqboi   = round(cfg.foi ./ (fsample ./ n_samples)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif ~isempty(cfg.foilim) % if input is foilim
  freqboilim = round(cfg.foilim ./ (fsample ./ n_samples)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end

% get matching jitter paramaters

if length(m)==1 && length(s)==1
    m=m.*ones(size(freqboi));
    s=s.*ones(size(freqboi));
elseif length(m)==2 && length(s)==2 && ~isempty(cfg.foilim)
    m=linspace(m(1),m(2),length(freqoi));
    s=linspace(s(1),s(2),length(freqoi));
elseif length(m)~=length(freqoi) || length(s)~=length(freqoi)
    error('Size of jittermean and jitterstd should match that of foi.');
end

    
    
    
  
% now the main bit of code!

%create new output data structure
data_jittered=data;

% go through trials
n_trial=length(data.trial);
for i=1:n_trial
    
    % perform fft on data
    dat_temp=fft(data.trial{i},[],2);
    for f=freqboi

        if isinf(s(f))
            %generate uniformally distributed jitter
            rand_phase=rand(size(dat_temp(:,f))).*pi.*2;
        else
            %generate  von mises distriuted jitter
            rand_temp=randn(size(dat_temp(:,f)));
            rand_sin=rand_temp.*sin(s(f)) + sin(m(f));
            rand_cos=rand_temp.*cos(s(f)) + cos(m(f));
            
            rand_phase=atan2(rand_sin,rand_cos);
        end
        
        % add jitter to fft data
        dat_temp(:,f) = dat_temp(:,f).*exp(sqrt(-1).*rand_phase);
        
    end
        
    % invert the fft
    data_jittered.trial{i}=ifft(dat_temp,[],2,'symmetric');
end

