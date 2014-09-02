function [shuffled] = ft_shuffledata(cfg, data)

% FT_SHUFFLEDATA shuffles a dataset across across a a particular domain
% Use as
%   [shuffled] = ft_shuffledata(cfg, data)
%
% The cfg.domain options determine which domains to shuffle along

% cfg.domain  = 'time'       Shuffle data along time points (if time is present)
%               'freq_pow'   Shuffle amplitude data along frequency points  (if freq is present)
%               'freq_phase' Shuffle phase data along frequency points
%               'channel'    Shuffle data along channels (or sources or components)
%               A cell-array containing some or all the above to shuffle across multiple domains
%                    (default = {'time','freq_pow','freq_phase'})
%
% The following options allow selection of specific time/freq/channel
% windows to perfom suffling. 
%
% cfg.toi            = Time-points to shuffle (default = 'all')
% cfg.foi            = Frequency-points to shuffle (default = 'all')
% cfg.channel        = Channels (or sources or components) to shuffle (default = 'all')
% cfg.trial          = Trials to shuffle (default = 'all')
% cfg.keepunshuffled = 'yes' or 'no;. If any of the above options are specfied,
% this option will cut out any unshuffled parts of the data. (default =
% 'yes');
%
% The following options change the way randomisations are computed.
%
% cfg.fullrand       = 'yes' or 'no'. Perform full randomisatiion. 
%                       If 'yes', Sperate randomisations along the shuffling 
%                       domain are performed independently for each instance in the other domains. 
%                       If 'no', only one randomisation is performed along the
%                       shuffling domain and repeated for each instance in other domains.
%                       (default = 'no')
% cfg.randconn    = 'yes' or 'no'. (Only valid if data is contains
%                        connectivity data and 'channel' is a shuffling domain)
%                        If 'yes', randomisations on bivariate cconnectivity data are performed across
%                        channel pairs. If univariate data is present, this will not be affected.
%                       If 'no', randnomisations are performed on bivariate connectivity data acoss single
%                       channels. If univariate data is present, this be shuffled with the same radomisation. 
%                       shuffled with the (default = 'no')
%
% THIS IS DIFFICULT TO IMPLIMENT SO STILL NOT READY
% cfg.randlock      = A cell-array containing domains to lock to the same
%                      randomisation. If diasbled (cfg.randlock=[]), all domains 
%                     to be shuffled will be done to their own unique
%                     randaomisation. (default = []).
%                      e.g. to shuffle across time but to lock the randomisations 
%                       across channels, use {'time','channel'}. 
%                       To lock randomisations of time and trials to channels, use: {'time','channel','trial'}
%                       To lock randomisations of time to channels, and INDEPENDALTY lock trials
%                      to channels, use: {{'time','channel'},{'trial','channel'}}. 
    


revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

% set defaults
cfg.domain            = ft_getopt(cfg, 'domain',    {'time','freq_pow','freq_phase'});
cfg.toi               = ft_getopt(cfg, 'toi',    'all');
cfg.foi               = ft_getopt(cfg, 'foi',   'all');
cfg.channel           = ft_getopt(cfg, 'channel',   'all');
cfg.fullrand          = ft_getopt(cfg, 'fullrand',   'no');
cfg.randconn          = ft_getopt(cfg, 'randconn',   'no');
cfg.keepunshuffled    = ft_getopt(cfg, 'keepunshuffled',   'yes');

if strcmp(cfg.fullrand,'yes')
    fullrand=1;
else
    fullrand=0;
end

% ID which domains to shuffle
shufftime=0; 
shuffpow=0; 
shuffphase=0; 
shuffchan=0; 
for i=1:length(cfg.domain)
    switch cfg.domain{i}
        case 'time'
            shufftime=1;
        case 'freq_pow'
            shuffpow=1;
        case 'freq_phase'
            shuffphase=1;
        case 'channel'
            shuffchan=1;
        otherwise
            warning('Did not recognise shuffle domain ''%s''',cfg.domain{i});
    end
end

% identify conditions of data

if isfield(data, 'time')
    hastime=1;
else
    hastime=0;
end

if isfield(data, 'freq')
      hasfreq=1;
else
    hasfreq=0;
end  

if isfield(data, 'fourierspctrm') || isfield(data, 'crsspctrm')
      hascomplx=1;
else
    hascomplx=0;
end  

if shufftime && ~hastime
    warning('Time not found in data. Cannot shuffle across time');
    shufftime=0;
end

if shufffreq && ~hasfreq
    warning('Frequency not found in data. Cannot shuffle across frequency');
    shufffreq=0;
end

if isfield(data, 'dimord')
    israw=0;
    dimord = data.dimord;
    dimtok = tokenize(dimord, '_');
else
    israw=1;
    dimtok={};
end

if isfield(data, 'labelcmb')
    isconn=1;
else
isconn=0;
end

% identify channels
    cfg.channel = ft_channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    
    
% %% not implimenting this at the moment!
% % set up randlock as fully nested cell-array
% randlock=cfg.randlock;
% if ~iscell(randlock)
%     randlock={{randlock}};
% elseif ~iscell(randlock{1})
%     randlock={randlock};
% end
%    
% 
% % identify which domains to lock together
% randlock_mat=zeros(5,5,length(randlock));
% for i=1:length(randlock)
%     randlock_idx=match_str({'time','freq_pow','freq_phase','chan','trial'}, randlock{i});
%     randlock_mat(randlock_idx,randlock_idx)=1;
% end
% 
%    %% 
    
% now shuffle!

if israw 

    n_trial=length(data.trial);
    n_chan=length(data.label);


    if shufftime
        if fullrand
            
            % shuffle time points with full randomisation
            for i=1:n_trial
                n_time=length(data.time{i});
                for j=1:n_chan
                    [rubbish rand_idx]=sort(rand(1,n_time));
                    data.trial{i}(j,:)=data.trial{i}(j,rand_idx);
                end
            end
        else
            % shuffle time points with single randomisation
            for i=1:length(data.trial)
                n_time(i)=length(data.time{i});
            end
            if length(unique(n_time))>1
                error('Shuffling of time points with cfg.randconn=''no'' requires trials to be of equal length.');
            end
            
            [rubbish rand_idx]=sort(rand(1,n_time(1)));
            for i=1:length(data.trial)
                data.trial{i}=data.trial{i}(:,rand_idx);
            end
        end
    end
    
    
    if shuffchan
        if fullrand
            % shuffle channels with full randomisation
            for i=1:length(data.trial)
                n_time(i)=length(data.time{i});
                for j=1:n_time
                    [rubbish rand_idx]=sort(rand(1,n_chan));
                    data.trial{i}(:,j)=data.trial{i}(rand_idx,j);
                end
            end
        else
            % shuffle channels with single randomisation
            [rubbish rand_idx]=sort(rand(1,n_chan));
            for i=1:n_trial
                data.trial{i}=data.trial{i}(rand_idx,:);
            end
        end
        
    end
    
else
    
    if shufftime
        if fullrand
            tidx=match_str(dimtok,'time');
            if ~isempty(tidx)
                
    
end
    
    

                 
        
        
    


