function rt_onlineprocess(cfg)
% This function is for online BCI.
% it calculates the alpha lateralization index and based on the thresholding
% scheme, sends the outputs to the presentation machine as 1 for left, 2 for
% right and 3 for no move.

%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see CHANNELSELECTION (default = {'MLO' 'MRO'}')
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'last')
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string, default is shared memory
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%   cfg.ntraining     = number, the number of trials to be used in the training phase (default=inf) 
%
% To stop the realtime function, you have to press Ctrl-C
%  Copyright (C) 2009, Ali Bahramisharif, Marcel van Gerven, Robert Oostenveld



% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
if ~isfield(cfg, 'readevent'),      cfg.readevent = 'no';     end % capture events?
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end % jump to end of file at initialization
if ~isfield(cfg, 'nexample'),       cfg.nexample = inf;       end
if ~isfield(cfg, 'foi'),            cfg.foi = 10;             end % from alpha frequency
if ~isfield(cfg, 'ntraining'),      cfg.ntraining= inf;       end % the number of trials to be used in the training phase
if ~isfield(cfg, 'verbose'),        cfg.verbose=0;            end % verbose
if ~isfield(cfg, 'saveLI'),         cfg.saveLI=0;             end % save lateralization index on the local machine
if ~isfield(cfg, 'channel'),        cfg.channel = {'MLO' 'MRO'}';               end %channels to be used
if ~isfield(cfg, 'datafile'),       cfg.datafile='shm://';                      end %input stream
if ~isfield(cfg, 'ostream'),        cfg.ostream='tcp://presentation011:1976';   end %output stream
% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');
cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true, 'retry', true);

% define a subset of channels for reading
cfg.channel = channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
    error('no channels were selected');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

if strcmp(cfg.jumptoeof, 'yes')
    prevSample = hdr.nSamples * hdr.nTrials;
else
    prevSample  = 0;
end
count       = 0;
%accumulative command
ac_cmd=[];


cfg.count=0;

%%fixme! grad element is not produced in the online setting.
%load grad dimensions as it is not produced in online setting
load grad;

% to run basic functions only once
cfg.first=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while cfg.count<cfg.nexample
    % determine number of samples available in buffer
    hdr = read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);

    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

    if newsamples>=blocksize
        % determine the samples to process
        if strcmp(cfg.bufferdata, 'last')
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'first')
            begsample  = prevSample+1;
            endsample  = prevSample+blocksize ;
        else
            error('unsupported value for cfg.bufferdata');
        end

        % this allows overlapping data segments
        if overlap && (begsample>overlap)
            begsample = begsample - overlap;
            endsample = endsample - overlap;
        end

        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        %fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment from buffer
        dat = read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);


        % put the data in a fieldtrip-like raw structure
        data.trial{1} = dat;
        data.hdr      = hdr;
        data.time{1}  = (begsample+(0:(endsample-begsample)))/hdr.Fs;

        if cfg.first
            data.label    = hdr.label(chanindx);
            data.fsample  = hdr.Fs;
            data.grad     =data.hdr.grad;
            data.cfg      =cfg;
            % apply some preprocessing options
        end
        data.trial{1} = preproc_baselinecorrect(data.trial{1});
        
        % correct for timing
        tt=data.time{1};
        pt=tt;
        tt=tt-min(tt);
        tt=tt./max(tt);
        tt=tt.*cfg.blocksize;
        data.time{1}=tt-0.5*cfg.blocksize;
        
        if cfg.first
            % planar gradiant
            cfgpl = [];
            cfgpl.channel      = 'all';
            cfgpl.planarmethod = 'sincos';
            data.grad=grad;
            cfgpl.online=1;
        end
        plan = megplanar(cfgpl, data);
        if cfg.first
            cfgpl.onlineprocess=plan.cfg.onlineprocess;
        end
        
        if cfg.first
            % frequency analysis
            cfgf = [];
            cfgf.output       = 'pow';
            %the following can be changed for speeding up,it was 'all' for the demo
            cfgf.channel      = 'all';%{'MLO33_dH','MLO33_dV','MRO33_dH','MRO33_dV'};
            cfgf.keeptrials   = 'yes';
            cfgf.method       = 'mtmconvol';
            cfgf.foi          = cfg.foi;
            cfgf.toi          = 0;
            cfgf.t_ftimwin    = ones(1,length(cfgf.foi)) * (cfg.blocksize-0.001); % 0.5 s. timewindow
            cfgf.taper        = 'hanning';
            cfgf.online=1;
        end
        freq = freqanalysis(cfgf,plan);
        if cfg.first || ~isfield(cfgf,'onlineprocess') || ~isfield(cfgf.onlineprocess,'cfg')
            cfgf.onlineprocess=freq.cfg.onlineprocess;
        end
        
        comb=freq;
        clear freq
        if cfg.first
            %planar gradient
            comb.grad.type='ctf275_planar';
            planar    = planarchannelset(comb);
            sel_dH    = match_str(comb.label, planar(:,1));  % indices of the horizontal channels
            sel_dV    = match_str(comb.label, planar(:,2));  % indices of the vertical  channels
            [dum, sel_planar] = match_str(comb.label, planar(:,1));
            comb.label=planar(sel_planar,3);
        end
        comb.powspctrm=comb.powspctrm(:,sel_dH,:) + comb.powspctrm(:,sel_dV,:);

        if cfg.first
            %%fixme: just 'MLO' or 'MRO' is not working!
            %selecting occipital channel
            listL={'MLO33'};%{'MLO11'    'MLO13'    'MLO14'    'MLO21'    'MLO22'    'MLO23'   'MLO24' 'MLO31'   'MLO32'    'MLO33'   'MLO34'    'MLO41'    'MLO42'    'MLO43'  'MLO44'    'MLO51' 'MLO52'    'MLO53'};
            listR={'MRO33'};%{ 'MRO11' 'MRO12' 'MRO13'    'MRO14' 'MRO21'    'MRO22'    'MRO23'    'MRO24'    'MRO31'    'MRO32'    'MRO33' 'MRO34' 'MRO41'    'MRO42'    'MRO43'    'MRO44'    'MRO51'    'MRO52' 'MRO53'};
            sel_L=match_str(comb.label,channelselection(listL,comb.label));
            sel_R=match_str(comb.label,channelselection(listR,comb.label));
        end
        
        %the parameter of interest is calculated in the next line
        cmd1=log10(mean(mean(comb.powspctrm(:,sel_L,:),2),3)./mean(mean(comb.powspctrm(:,sel_R,:),2),3));

        cfg.count=cfg.count+1;
        if cfg.count<cfg.ntraining
            %accumulative command is saving for updating mean and std.
            ac_cmd=[ac_cmd,cmd1];
            if cfg.verbose
                disp(sprintf('trainig sample %d',cfg.count));
            end
            R_threshold=mean(ac_cmd)+1*std(ac_cmd);
            L_threshold=mean(ac_cmd)-1*std(ac_cmd);
        end
        if cmd1 > R_threshold
            cmd=1;
        elseif cmd1 < L_threshold
            cmd=2;
        else
            cmd=3;
        end
        if ~isempty(cfg.ostream)
            % send command
            cmdevent.type = 'uint';
            cmdevent.offset = [];
            cmdevent.duration = [];
            cmdevent.sample = abs(pt(1)*data.fsample);
            cmdevent.timestamp = pt(1);
            cmdevent.value = cmd;
            write_event(cfg.ostream,cmdevent);
        end
        if cfg.saveLI
            %for saving ac_cmd to check it.
            if (cfg.count/200)==floor(cfg.count/200)
                save('~/Desktop/ac_cmd','ac_cmd');
            end
        end
        cfg.first=0;
    end % if enough new samples
end % while true
