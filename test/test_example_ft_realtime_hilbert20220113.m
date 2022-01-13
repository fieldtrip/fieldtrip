function test_example_ft_realtime_hilbert

% MEM 4gb
% WALLTIME 00:10:00

%
%% Realtime neurofeedback application based on Hilbert phase estimation
%
% This page describes the neurofeedback application that was developed by Marco Rotonda for his MSc thesis **"Controllo volontario della sincronizzazione di fase intraemisferica nella banda gamma eeg mediante neurofeedback"** (Voluntary control of intrahemispheric phase sychronization in the EEG gamma band using neurofeedback).
%
% The thesis (in Italian) which describes all details and the results is available [here](/assets/pdf/example/rt_hilbert/rotonda_thesis_msc.pdf).
%
% The signal analysis details are available [here](/assets/pdf/example/gamma_analisys.pdf).
%
% Neurofeedback, as the word suggest, is feeding back to the subject it's neuronal activity analyzed by it's EEG activity in (almost) real time. In this training we give back to the subjest a visual stimulation (vertical bars) when it's gamma phase synchronization between F3-P3 and F4-P4 increased over the baseline (bars up) or decreased (bars down).
% Here is the feedback
%
%
% If you scroll the script you can see that this image is a particular point of view of a 3D visualization. I you change the POV you could see all the data from the last 3 minutes.
%
%% # Procedure
%
% Get the FieldTripBufferDemo from the workshop_bci2000 folder on the FieldTrip FTP server and . Start BCI2000 and modifiy the config. The script works originally with 40 channels (see the variable 'lab' in the realtime_hilbert script), the samplingrate used is 1000Hz. Use a blocksize of 1000 samples. Ensure BCI2000 runs for more than 2 minutes in the application tab. Set config and start BCI2000.
% Start MATLAB and in the shell type realtime_hilbert.
%
% Now it will start the function realtime_baseline that for 2 minutes will record the subject baseline and store the results for realtime_hilbert.
% Basically there are 12 channels: 4 for EEG (I used F3-P3-F4-P4), 4 for the EOG, 4 for the EMG (front and neck). In this script, those are selected from the forty originally recorded channels you can find in the variable 'lab'.
% The algorithm will care to take away both EOG and EMG artifacts.
%
%% # MATLAB code for ft_realtime_hilbert
%
function ft_realtime_hilbert()

% FT_REALTIME_HILBERT is a neurofeedback application based on Hilbert phase estimation.
%
% Use as
%   ft_realtime_hilbert()
% with the following configuration options that are coded inside the function
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'gui')
%   cfg.foilim     = [Flow Fhigh] (default = [0 120])
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'first')
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string ( es:  cfg.datafile='buffer://localhost:1972';)
%   cfg.headerfile    = string ( es:  cfg.headerfile='buffer://localhost:1972';)
%   cfg.filename      = string ( es:  cfg.filename = 'buffer://localhost:1972';)
%   cfg.eventfile     = string ( es:  cfg.eventfile = 'buffer://localhost:1972';)
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

% Take off the warning message to avoid problems with ATAN2 and hilbert.
warning off all;

% This is to save subject data
starttime= DATESTR(now, 30);
subjname = input('Insert the subject name.>>>>>', 's');
trialdata=strcat(starttime,subjname);

realtime_baseline(); %This it will launch the baseline script

load means;  % get means from file created by realtime_baseline

cfg = [];
cfg.datafile = 'buffer://localhost:1972';
cfg.headerfile = 'buffer://localhost:1972';
cfg.filename = 'buffer://localhost:1972';

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0.5;        end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = {'all', '-Fp1', '-Fp2',...
                                    '-F7','-Fz','-F8','-FT7', '-FC3', '-FCz',...
                                    '-FC4', '-FT8', '-T3', '-C3', '-Cz', '-C4','-T4','-TP7',...
                                    '-CP3', '-CPz', '-CP4','-TP8','-A1', '-T5', '-Pz',...
                                    '-T6','-A2', '-O1', '-Oz', '-O2'};
end                                                               % This will select F3 F4 P3 P4 for EEG
                                                                  % X1 X2 X3 X4 for EMG
                                                                  % X5 X6 X7 X8 for EOG
if ~isfield(cfg, 'foilim'),         cfg.foilim = [0.1 100];   end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata='last';   end

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading
lab=['X1 '; 'X2 '; 'Fp1'; 'Fp2'; 'X3 '; 'X4 '; 'F7 '; 'F3 '; 'Fz '; 'F4 ';...
    'F8 ';'FT7'; 'FC3'; 'FCz'; 'FC4'; 'FT8'; 'T3 '; 'C3 '; 'Cz '; 'C4 ';...
    'T4 ';'TP7'; 'CP3'; 'CPz'; 'CP4'; 'TP8'; 'A1 '; 'T5 '; 'P3 '; 'Pz ';...
    'P4 '; 'T6 ';'A2 '; 'O1 '; 'Oz '; 'O2 '; 'X5 '; 'X6 '; 'X7 '; 'X8 '];
label=cellstr(lab);
cfg.channel = ft_channelselection(cfg.channel, label);
chanindx    = match_str(label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

% determine the size of blocks to process
blocksize = cfg.blocksize * hdr.Fs;
overlap   = cfg.overlap * hdr.Fs;
prevSample  = 0;
count       = 0;

% Create arrays that contains the rhos found
vectordim = 360;         % vector dimention for the lasts 3 minutes of data
vl=zeros(vectordim,1);  % the visual vector of data for left synchrony
vr=zeros(vectordim,1);  % the visual vector of data for right synchrony
vl1=vl(1);              % the first element of the vector for the visual feedback
vr1=vr(1);              % the first element of the vector for the visual feedback

% This is the ratio above which accept the coherence
meanl = baselinel;           % It has to be calculated from the baseline
meanr = baseliner;           % It has to be calculated from the baseline

% This is the mean above which no EMG feedback has given
meanf = baselinef+(stf*2);             % It has to be calculated from the baseline
meann = baselinen+(stn*2);             % It has to be calculated from the baseline

% This is used to plot the feedback step
fullscreen = get(0,'ScreenSize');
fig1 = figure('NumberTitle','off', ...
    'MenuBar','none', ...
    'Units','pixels', ...
    'Position',[0 0 fullscreen(3) fullscreen(4)]);

% plot the feedback on the second monitor
% set(gcf,'position',[1025,1,1024,768]);

% fig2=figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

    % Create a matrix that contains the rhos found
    M = [vl vr]; %#ok`<NASGU>`

    % Create a matrix for the visual feedback
    K = [vl1 vr1];

    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'cache', true);

    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

    if newsamples>=blocksize

        % determine the samples to process
        if strcmp(cfg.bufferdata, 'last') && count==0
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'last')
          begsample  = prevSample+1;
          endsample  = prevSample+blocksize ;
        elseif strcmp(cfg.bufferdata, 'first')
          begsample  = prevSample+1;
          endsample  = prevSample+blocksize ;
        else
            error('unsupported value for cfg.bufferdata');
        end

        % this allows overlapping data segments
        if overlap && (begsample>overlap) %#ok`<BDLGI>`
          begsample = begsample - overlap;
          endsample = endsample - overlap;
        end

        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        % fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample,...
            'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to the hilbert phase sinchronisation from the data %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % put the data in a fieldtrip-like raw structure
        data.trial{1} = dat;
        data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;

        % correction of EOG based on algoritm fro
        % Author: German Gomez-Herrero
        %         german.gomezherrero@ieee.org
        %         http://www.cs.tut.fi/~gomezher/index.htm
        %         Institute of Signal Processing
        %         Tampere University of Technology, 2007
        % Reference
        % [1] P. He et al., Med. Biol. Comput. 42 (2004), 407-412
        % [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall

        data.trial{2}(9,:)=data.trial{1}(9,:)-data.trial{1}(10,:);
        data.trial{2}(11,:)=data.trial{1}(11,:)-data.trial{1}(12,:);
        opt.refdata=[data.trial{2}(9,:);data.trial{2}(11,:)];
        [data.trial{3}] = crls_regression(data.trial{1}(5:8,:), opt);

        % Build a FIR filter for EMG correction
        N    = 150;      % Order
        gammaband = [35 45];
        emgband = [60 80];
        emgfnband = [60 499];
        flag = 'scale';  % Sampling Flag
        Beta = 0.9;      % Window Parameter
        win = kaiser(N+1, Beta);

        % Correction between EMG and EEG based on Sheer D.E. "Biofeedback training
        % of 40-Hz eeg and behavior", pp. 325-362, on Behavior and
        % brain electrical activity (1975), Plenum Press. New York

        gammafilter  = fir1(N, gammaband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltgamma1 = filtfilt(gammafilter,1,data.trial{1}(5,:));
        datfiltgamma2 = filtfilt(gammafilter,1,data.trial{1}(6,:));
        datfiltgamma3 = filtfilt(gammafilter,1,data.trial{1}(7,:));
        datfiltgamma4 = filtfilt(gammafilter,1,data.trial{1}(8,:));
        datfiltgamma = [datfiltgamma1; datfiltgamma2; datfiltgamma3; datfiltgamma4];

        emgfilter = fir1(N, emgband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltemg1 = filtfilt(emgfilter,1,data.trial{1}(5,:));
        datfiltemg2 = filtfilt(emgfilter,1,data.trial{1}(6,:));
        datfiltemg3 = filtfilt(emgfilter,1,data.trial{1}(7,:));
        datfiltemg4 = filtfilt(emgfilter,1,data.trial{1}(8,:));
        datfiltemg = [datfiltemg1; datfiltemg2; datfiltemg3; datfiltemg4];

        datfiltemgsqr = datfiltemg.^2;
        datfiltgammasqr = datfiltgamma.^2;
        datfiltcrossqr = (datfiltemg.*datfiltgamma).^2;

        correction = datfiltgammasqr-(datfiltcrossqr./datfiltemgsqr);

        data.trial{3}(5:8,:) = datfiltgamma - correction;

        % Find the EMG on forehead and neck
        data.trial{2}(1,:)=data.trial{1}(1,:)-data.trial{1}(2,:);   % Frontal electrods
        data.trial{2}(3,:)=data.trial{1}(3,:)-data.trial{1}(4,:);   % Neck electrods
        emgfnfilter = fir1(N, emgfnband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltemgf = filter(emgfnfilter,1,data.trial{2}(1,:));
        datfiltemgn = filter(emgfnfilter,1,data.trial{2}(3,:));
        datfiltemgf = abs(datfiltemgf);
        datfiltemgn = abs(datfiltemgn);
        extrMaxValuef = datfiltemgf(find(diff(sign(diff(datfiltemgf)))==-2)+1);
        extrMaxValuen = datfiltemgn(find(diff(sign(diff(datfiltemgn)))==-2)+1);
        extrMaxIndexf =   find(diff(sign(diff(datfiltemgf)))==-2)+1;
        extrMaxIndexn =   find(diff(sign(diff(datfiltemgn)))==-2)+1;
        upf = extrMaxValuef;
        upn = extrMaxValuen;
        upf_t = data.time{1}(extrMaxIndexf);
        upn_t = data.time{1}(extrMaxIndexn);
        upf = interp1(upf_t,upf,data.time{1},'linear');
        upn = interp1(upn_t,upn,data.time{1},'linear');
        emgfmean = nanmean (upf'); %#ok`<UDIM>`
        emgnmean = nanmean (upn'); %#ok`<UDIM>`
        % plot(data.time{1},upf,'r')

%
        % Istantaneous (proto)phase difference found via Hilbert
        % Based on Pikovsky, A. R. (2001). Synchronization. A Universal
        % Concept In Nonlinear Sciences. Cambridge: Cambridge University
        % Press, pag. 368 A2.7

        % crate the data needed for phase coherence index
        chan1=data.trial{3}(5,:); % F3
        chan2=data.trial{3}(6,:); % F4
        chan3=data.trial{3}(7,:); % P3
        chan4=data.trial{3}(8,:); % P4
        chan1h = hilbert(chan1);
        chan2h = hilbert(chan2);
        chan3h = hilbert(chan3);
        chan4h = hilbert(chan4);
        chan1hi = imag(chan1h);
        chan2hi = imag(chan2h);
        chan3hi = imag(chan3h);
        chan4hi = imag(chan4h);

        % find the istantaneous left hemisphere (proto)phase difference
        phil = atan2(((chan1hi .* chan3)-(chan1 .* chan3hi)),...
                     ((chan1 .* chan3)+(chan1hi .* chan3hi)));

        % find the istantaneous right hemisphere (proto)phase difference
        phir = atan2(((chan2hi .* chan4)-(chan2 .* chan4hi)),...
                     ((chan2 .* chan4)+(chan2hi .* chan4hi)));

        % find the right hemisphere synchronization index
        sumsinr = sum(sin(phir))/blocksize;
        sumcosr = sum(cos(phir))/blocksize;
        rhor = sqrt(sumsinr.^2 + sumcosr.^2);

        % find the left hemisphere synchronization index
        sumsinl = sum(sin(phil))/blocksize;
        sumcosl = sum(cos(phil))/blocksize;
        rhol = sqrt(sumsinl.^2 + sumcosl.^2);

        % Give the visual feedback
        if (meanf>emgfmean) && (meann>emgnmean)
            clf;
            visual = bar3(K,0.3);
            view([-90 0]);
            grid off;
            shading interp;
            for i = 1:length(visual)
                zdata = get(visual(i),'Zdata');
                set(visual(i),'Cdata',zdata,'EdgeColor','none')
                colormap hot;
            end
            set(gca,'ZColor',[0.8 0.8 0.8],'Zlim',[-1 1],'YColor',[0.8 0.8 0.8],...
                'Ylim',[0 3],'XColor',[0.8 0.8 0.8],'Xlim',[0.85 1.15],...
                'Color',[0.8 0.8 0.8],'CLim', [-1 1]);

            % Create colorbar
            colorbar([0.5 0.148 0.02 0.73],'ZColor',[0.8 0.8 0.8],'YTick',[],...
                'YColor',[0.8 0.8 0.8],'XColor',[0.8 0.8 0.8]);

            % Create textboxes
            annotation('textbox',[0.496 0.89 0.03 0.04],'String',{'+'},...
                'HorizontalAlignment','center','FontSize',20,'FitBoxToText','off','EdgeColor','none');
            annotation('textbox',[0.50 0.11 0.02 0.04],'String',{'-'},...
                'HorizontalAlignment','center','FontSize',20,'FitBoxToText','off','EdgeColor','none');
            annotation('textbox',[0.25 0.49 0.05 0.07],'String',{'Emisfero','sinistro'},...
                'HorizontalAlignment','center','FontSize',14,'FitBoxToText','off','EdgeColor','none');
            annotation('textbox',[0.735 0.49 0.05 0.07],'String',{'Emisfero','destro'},...
                'HorizontalAlignment','center','FontSize',14,'FitBoxToText','off','EdgeColor','none');

            %create a copy of the plot for the experimenter
            % h1=gcf;
            % h2=figure;
            % objects=allchild(h1);
            % fig2=copyobj(get(h1,'children'),h2);

            % force MATLAB to update the figure
            drawnow ;

            % add the step to the array for the feedback and upgrade the M
            % for the final performance plot

            vl = vl([end 1:end-1]);
            vl(1)=rhol-meanl;
            vl1=vl(1);

            vr = vr([end 1:end-1]);
            vr(1)=rhor-meanr;
            vr1=vr(1);

        elseif  emgfmean>meanf
            clf;
            set(gca,'ZColor',[0.8 0.8 0.8],'Zlim',[0 1],'YColor',[0.8 0.8 0.8],...
                'Ylim',[0 3],'XColor',[0.8 0.8 0.8],'Xlim',[0.85 1.15],...
                'OuterPosition', [-0.0175 0.185 1 0.605],...
                'Color',[0.8 0.8 0.8],'CLim', [-50 50]);
            annotation(fig1,'textbox',[0.20 0.35 0.6292 0.1929],...
                'String',{'Rilassa i muscoli della fronte'},...
                'HorizontalAlignment','center','FontSize',20,'FitBoxToText','off','EdgeColor','none',...
                'Color',[1 0 0]);
            drawnow ;

        elseif  emgnmean>meann
            clf;
            set(gca,'ZColor',[0.8 0.8 0.8],'Zlim',[0 1],'YColor',[0.8 0.8 0.8],...
                'Ylim',[0 3],'XColor',[0.8 0.8 0.8],'Xlim',[0.85 1.15],...
                'OuterPosition', [-0.0175 0.185 1 0.605],...
                'Color',[0.8 0.8 0.8],'CLim', [-1 1]);
            annotation(fig1,'textbox',[0.20 0.35 0.6292 0.1929],...
                'String',{'Rilassa i muscoli del collo'},...
                'HorizontalAlignment','center','FontSize',20,'FitBoxToText','off','EdgeColor','none',...
                'Color',[1 0 0]);
            drawnow ;

        end % end feedbacks

    end % if enough new samples

    if count > vectordim
        cd 'C:\TEST\hilbert\risultati';
        save (trialdata);
        save hilbert.mat;
        close all hidden;
        cd 'C:\TEST';
        break;
    end

end % while true

%% # MATLAB code for ft_realtime_baseline
%
function ft_realtime_baseline()

% FT_REALTIME_BASELINE is a neurofeedback application based on Hilbert phase estimation.
%
% Use as
%   ft_realtime_baseline()
% with the following configuration options that are coded inside the function
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'gui')
%   cfg.foilim     = [Flow Fhigh] (default = [0 120])
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'first')
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string ( es:  cfg.datafile='buffer://localhost:1972';)
%   cfg.headerfile    = string ( es:  cfg.headerfile='buffer://localhost:1972';)
%   cfg.filename      = string ( es:  cfg.filename = 'buffer://localhost:1972';)
%   cfg.eventfile     = string ( es:  cfg.eventfile = 'buffer://localhost:1972';)
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

%Take off the warning message to avoid problems with ATAN2 and hilbert.
warning off all;

cfg = [];
cfg.datafile = 'buffer://localhost:1972';
cfg.headerfile = 'buffer://localhost:1972';
cfg.filename = 'buffer://localhost:1972';

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0.5;        end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = {'all', '-Fp1', '-Fp2',...
                                    '-F7','-Fz','-F8','-FT7', '-FC3', '-FCz',...
                                    '-FC4', '-FT8', '-T3', '-C3', '-Cz', '-C4','-T4','-TP7',...
                                    '-CP3', '-CPz', '-CP4','-TP8','-A1', '-T5', '-Pz',...
                                    '-T6','-A2', '-O1', '-Oz', '-O2'};
end                                                               % This will select F3 F4 P3 P4 for EEG
                                                                  % X1 X2 X3 X4 for EMG
                                                                  % X5 X6 X7 X8 for EOG
if ~isfield(cfg, 'foilim'),         cfg.foilim = [0.1 100];   end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata='last';   end

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');
cfg = checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading
lab=['X1 '; 'X2 '; 'Fp1'; 'Fp2'; 'X3 '; 'X4 '; 'F7 '; 'F3 '; 'Fz '; 'F4 ';...
    'F8 ';'FT7'; 'FC3'; 'FCz'; 'FC4'; 'FT8'; 'T3 '; 'C3 '; 'Cz '; 'C4 ';...
    'T4 ';'TP7'; 'CP3'; 'CPz'; 'CP4'; 'TP8'; 'A1 '; 'T5 '; 'P3 '; 'Pz ';...
    'P4 '; 'T6 ';'A2 '; 'O1 '; 'Oz '; 'O2 '; 'X5 '; 'X6 '; 'X7 '; 'X8 '];
label=cellstr(lab);
cfg.channel = ft_channelselection(cfg.channel, label);
chanindx    = match_str(label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end

% determine the size of blocks to process
blocksize = cfg.blocksize * hdr.Fs;
overlap   = cfg.overlap * hdr.Fs;
prevSample  = 0;
count       = 0;

% Create arrays that contains the rhos found
nummaxarray = 240; % the step is 500ms so 240 are 2 minutes
i=1;
vlmean=zeros(1,nummaxarray);
vrmean=zeros(1,nummaxarray);
frontmeans=zeros(1,nummaxarray);
neckmeans=zeros(1,nummaxarray);

% This is used to plot the screen for the subject
fullscreen = get(0,'ScreenSize');
fig1 = figure('NumberTitle','off', ...
    'MenuBar','none', ...
    'Units','pixels', ...
    'Position',[0 0 fullscreen(3) fullscreen(4)]);

% plot the feedback on the second monitor
% set(gcf,'position',[1025,1,1024,768]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true

    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'cache', true);

    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

    if newsamples>=blocksize

        % determine the samples to process
        if strcmp(cfg.bufferdata, 'last') && count==0
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'last')
          begsample  = prevSample+1;
          endsample  = prevSample+blocksize ;
        elseif strcmp(cfg.bufferdata, 'first')
          begsample  = prevSample+1;
          endsample  = prevSample+blocksize ;
        else
            error('unsupported value for cfg.bufferdata');
        end

        % this allows overlapping data segments
        if overlap && (begsample>overlap) %#ok`<BDLGI>`
          begsample = begsample - overlap;
          endsample = endsample - overlap;
        end

        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        % fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to the hilbert phase sinchronisation from the data %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % put the data in a fieldtrip-like raw structure
        data.trial{1} = dat;
        data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;

        % correction of EOG based on algoritm fro
        % Author: German Gomez-Herrero
        %         german.gomezherrero@ieee.org
        %         http://www.cs.tut.fi/~gomezher/index.htm
        %         Institute of Signal Processing
        %         Tampere University of Technology, 2007
        % Reference
        % [1] P. He et al., Med. Biol. Comput. 42 (2004), 407-412
        % [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall

        data.trial{2}(9,:)=data.trial{1}(9,:)-data.trial{1}(10,:);
        data.trial{2}(11,:)=data.trial{1}(11,:)-data.trial{1}(12,:);
        opt.refdata=[data.trial{2}(9,:);data.trial{2}(11,:)];
        [data.trial{3}] = crls_regression(data.trial{1}(5:8,:), opt);

        % Build a FIR filter for EMG correction
        N    = 150;      % Order
        gammaband = [35 45];
        emgband = [60 80];
        emgfnband = [60 499];
        flag = 'scale';  % Sampling Flag
        Beta = 0.9;      % Window Parameter
        win = kaiser(N+1, Beta);

        % Correction between EMG and EEG based on Sheer D.E. "Biofeedback training
        % of 40-Hz eeg and behavior", pp. 325-362, on Behavior and
        % brain electrical activity (1975), Plenum Press. New York

        gammafilter  = fir1(N, gammaband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltgamma1 = filtfilt(gammafilter,1,data.trial{1}(5,:));
        datfiltgamma2 = filtfilt(gammafilter,1,data.trial{1}(6,:));
        datfiltgamma3 = filtfilt(gammafilter,1,data.trial{1}(7,:));
        datfiltgamma4 = filtfilt(gammafilter,1,data.trial{1}(8,:));
        datfiltgamma = [datfiltgamma1; datfiltgamma2; datfiltgamma3; datfiltgamma4];

        emgfilter = fir1(N, emgband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltemg1 = filtfilt(emgfilter,1,data.trial{1}(5,:));
        datfiltemg2 = filtfilt(emgfilter,1,data.trial{1}(6,:));
        datfiltemg3 = filtfilt(emgfilter,1,data.trial{1}(7,:));
        datfiltemg4 = filtfilt(emgfilter,1,data.trial{1}(8,:));
        datfiltemg = [datfiltemg1; datfiltemg2; datfiltemg3; datfiltemg4];

        datfiltemgsqr = datfiltemg.^2;
        datfiltgammasqr = datfiltgamma.^2;
        datfiltcrossqr = (datfiltemg.*datfiltgamma).^2;

        correction = datfiltgammasqr-(datfiltcrossqr./datfiltemgsqr);

        data.trial{3}(5:8,:) = datfiltgamma - correction;

        % Find the EMG on forehead and neck
        data.trial{2}(1,:)=data.trial{1}(1,:)-data.trial{1}(2,:);   % Frontal electrods
        data.trial{2}(3,:)=data.trial{1}(3,:)-data.trial{1}(4,:);   % Neck electrods
        emgfnfilter = fir1(N, emgfnband/(hdr.Fs/2), 'bandpass', win, flag);
        datfiltemgf = filter(emgfnfilter,1,data.trial{2}(1,:));
        datfiltemgn = filter(emgfnfilter,1,data.trial{2}(3,:));
        datfiltemgf = abs(datfiltemgf);
        datfiltemgn = abs(datfiltemgn);
        extrMaxValuef = datfiltemgf(find(diff(sign(diff(datfiltemgf)))==-2)+1);
        extrMaxValuen = datfiltemgn(find(diff(sign(diff(datfiltemgn)))==-2)+1);
        extrMaxIndexf =   find(diff(sign(diff(datfiltemgf)))==-2)+1;
        extrMaxIndexn =   find(diff(sign(diff(datfiltemgn)))==-2)+1;
        upf = extrMaxValuef;
        upn = extrMaxValuen;
        upf_t = data.time{1}(extrMaxIndexf);
        upn_t = data.time{1}(extrMaxIndexn);
        upf = interp1(upf_t,upf,data.time{1},'linear');
        upn = interp1(upn_t,upn,data.time{1},'linear');
        emgfmean = nanmean (upf'); %#ok`<UDIM>`
        emgnmean = nanmean (upn'); %#ok`<UDIM>`
        % plot(data.time{1},upf,'r')

%
        % Istantaneous (proto)phase difference found via Hilbert
        % Based on Pikovsky, A. R. (2001). Synchronization. A Universal
        % Concept In Nonlinear Sciences. Cambridge: Cambridge University
        % Press, pag. 368 A2.7

        % crate the data needed for phase coherence index
        chan1=data.trial{3}(5,:); % F3
        chan2=data.trial{3}(6,:); % F4
        chan3=data.trial{3}(7,:); % P3
        chan4=data.trial{3}(8,:); % P4
        chan1h = hilbert(chan1);
        chan2h = hilbert(chan2);
        chan3h = hilbert(chan3);
        chan4h = hilbert(chan4);
        chan1hi = imag(chan1h);
        chan2hi = imag(chan2h);
        chan3hi = imag(chan3h);
        chan4hi = imag(chan4h);

        % find the istantaneous left hemisphere (proto)phase difference
        phil = atan2(((chan1hi .* chan3)-(chan1 .* chan3hi)),...
                     ((chan1 .* chan3)+(chan1hi .* chan3hi)));

        % find the istantaneous right hemisphere (proto)phase difference
        phir = atan2(((chan2hi .* chan4)-(chan2 .* chan4hi)),...
                     ((chan2 .* chan4)+(chan2hi .* chan4hi)));

        % find the right hemisphere synchronization index
        sumsinr = sum(sin(phir))/blocksize;
        sumcosr = sum(cos(phir))/blocksize;
        rhor = sqrt(sumsinr.^2 + sumcosr.^2);

        % find the left hemisphere synchronization index
        sumsinl = sum(sin(phil))/blocksize;
        sumcosl = sum(cos(phil))/blocksize;
        rhol = sqrt(sumsinl.^2 + sumcosl.^2);

        % update the array for the mean
        vlmean(i)=rhol;
        vrmean(i)=rhor;
        frontmeans(i)=emgfmean;
        neckmeans(i)=emgnmean;

        i=i+1;

    end % if enough new samples

    % screen for the subject
    clf;
    set(gca,'ZColor',[0.8 0.8 0.8],'Zlim',[-50 50],'YColor',[0.8 0.8 0.8],...
        'Ylim',[0 3],'XColor',[0.8 0.8 0.8],'Xlim',[0.85 1.15],...
        'DataAspectRatio',[0.2 0.1 3],'OuterPosition', [-0.0175 0.185 1 0.605],...
        'Color',[0.8 0.8 0.8],'CLim', [-50 50]);
    annotation(fig1,'textbox',[0.18 0.35 0.6292 0.1929],...
        'String',{'Rilassati mantenendo gli occhi aperti'},...
        'HorizontalAlignment','center','FontSize',20,'FitBoxToText','off','EdgeColor','none',...
        'Color',[1 0 0]);
    drawnow ;

    if count > (nummaxarray-1)
        close all hidden;
        break;
    end

end % while true

% print the mean baseline
baselinel = mean (vlmean);
fprintf ('mean baseline left= %d \n', baselinel);
baseliner = mean (vrmean);
fprintf ('mean baseline right= %d \n', baseliner);
baselinef = mean (frontmeans);
fprintf ('mean baseline forehead= %d \n', baselinef);
baselinen = mean (neckmeans);
fprintf ('mean baseline neck= %d \n', baselinen);
stf = std(frontmeans);
fprintf ('Std baseline front= %d \n', stf);
stn = std(neckmeans);
fprintf ('Std baseline neck= %d \n', stn);

save means.mat;   % save to the working folder - (edit jonaweber 2010)
