function ft_realtime_ouunpod(cfg)

% FT_REALTIME_OUUNPOD is an example realtime application for online power
% estimation and visualisation. It is designed for use with the OuUnPod, an
% OpenEEG based low cost EEG system with two channels, but in principle should
% work for any EEG or MEG system.
%
% Use as
%   ft_realtime_ouunpod(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.foilim     = [Flow Fhigh] (default = [0 80])
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.bufferdata = whether to start on the 'first or 'last' data that is available (default = 'last')
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C
%
% See also http://ouunpod.blogspot.com

% Copyright (C) 2008-2012, Robert Oostenveld
% Copyright (C) 2012, Stephen Whitmarsh

% setup MIDI, see http://en.wikipedia.org/wiki/General_MIDI
beatdrum = true;

if ~ispc
  eval('midiOut = @dummyfunction;');
end

m = midiOut; %'Microsoft GS Wavetable Synth = device number 2;
midiOut('O',2); % o for output; 2 for devicde nr 2
midiOut('.',1); % all off
midiOut('.',2); % all off

% midiOut('P',2,20); organ

midiOut('P',1,53);
midiOut('P',2,53);
midiOut('+',1,[64 65],[127 127]);         % command, channelnr, key, velocity
midiOut('+',2,[64 65 67],[127 127 127]);  % command, channelnr, key, velocity
midiOut(uint8([175+1,7,0]));              % change volume
midiOut(uint8([175+2,7,0]));

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat      = [];       end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat    = [];       end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat     = [];       end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize       = 0.05;     end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel         = 'all';    end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata      = 'last';   end % first or last
if ~isfield(cfg, 'dataset'),        cfg.dataset         = 'buffer:\\localhost:1972'; end;
if ~isfield(cfg, 'foilim'),         cfg.foilim          = [1 80];   end
if ~isfield(cfg, 'blockmem'),       cfg.blockmem        = 40;       end % how many blocks in the past will be taking into calc./disp
if ~isfield(cfg, 'scale'),          cfg.scale           = 1;        end

% schemerlamp = Lamp('com9');

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true, 'retry', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan~=2
  error('exactly two channels should be selected');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);

% this is used for scaling the figure
powmax(1:2) = 0;
prevSample  = 0;
count       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left_thresh_ampl = -100;
left_thresh_time = [cfg.blockmem*blocksize-blocksize*4 cfg.blockmem*blocksize];

% FIXME these are hardcoded, but might be incompatible with the cfg and data settings
right_freq          = [40 45];
right_offset        = 0.5;
right_mult          = 127/0.5;

close all;
figure

TFR = zeros(2,(cfg.foilim(2)-cfg.foilim(1)+1),100);

while true
  
  % determine number of samples available in buffer
  hdr = ft_read_header(cfg.headerfile, 'cache', true);
  
  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);
  
  if newsamples>=blocksize
    
    % determine the samples to process
    if strcmp(cfg.bufferdata, 'last')
      begsample  = hdr.nSamples*hdr.nTrials - blocksize*(cfg.blockmem) + 1;
      endsample  = hdr.nSamples*hdr.nTrials;
    elseif strcmp(cfg.bufferdata, 'first')
      begsample  = prevSample+1;
      endsample  = prevSample+blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end
    
    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
    
    if count > cfg.blockmem % memory
      
      % read the data segment from buffer
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
      dat = cfg.scale * dat;
      
      % construct a matching time axis
      time  = ((begsample:endsample)-1)/hdr.Fs;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % from here onward it is specific to the power estimation from the data
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % apply preprocessing options
      dat = ft_preproc_baselinecorrect(dat);
      dat = ft_preproc_highpassfilter(dat,hdr.Fs,5,1,'but','twopass');
      dat = ft_preproc_bandstopfilter(dat,hdr.Fs,[45  55],4,'but','twopass');
      dat = ft_preproc_bandstopfilter(dat,hdr.Fs,[95 115],4,'but','twopass');
            
      %             h = get(gca, 'children');
      %             hold on
      %
      %             if ~isempty(h)
      %                 % done on every iteration
      %                 delete(h);
      %             end
      %
      %             if isempty(h)
      %                 % done only once
      %                 powmax(1:2) = 0;
      %                 grid on
      %             end
      
      for i=1:nchan
        % pow(i) = psd(specest, dat(i,:), 'Fs', hdr.Fs);
        % [spec{i} ifreq{i}] = ft_specest_mtmfft(dat(i,:), time,'taper','dpss','tapsmofrq',4,'freqoi',[1:40]);
        [spec ifreq] = ft_specest_mtmfft(dat(i,:), time,'taper','dpss','tapsmofrq',2,'freqoi',[cfg.foilim(1) : cfg.foilim(2)]);
        % pow(i,:) = squeeze(abs(spec{i}));
        pow(i,:) = squeeze(mean(abs(spec)));
        % pow(i,:) = pow(i,:) ./ mean(pow(i,:));
        powmax(i) = max(max(pow(i,:)), powmax(i)*0.9); % this keeps a history
        for ii = 1 : size(TFR,3)-1
          TFR(i,:,ii) = TFR(i,:,ii+1);
        end
        TFR(i,:,end) = pow(i,:);
      end
      
      subplot(3,2,3);
      % plot(pow(1).Frequencies, pow(1).Data);
      % bar(pow(1).Frequencies,pow(1).Data);
      bar(1:length(ifreq),pow(1,:),0.5);
      axis([cfg.foilim(1) cfg.foilim(2) 0 powmax(1)*1.5]);
      axis([cfg.foilim(1) cfg.foilim(2) 0 10]);
      
      % str = sprintf('time = %d s\n', round(mean(time)));
      % title(str);
      xlabel('frequency (Hz)');
      ylabel('power');
      
      subplot(3,2,4);
      % plot(pow(2).Frequencies, pow(2).Data);
      % bar(pow(2).Frequencies,pow(2).Data);
      bar(1:length(ifreq),pow(2,:),0.5);
      ax = axis;
      % axis([cfg.foilim(1) cfg.foilim(2) 0 powmax(2)*1.5]);
      axis([cfg.foilim(1) cfg.foilim(2) 0 10]);
      line([right_freq(1) right_freq(1)],[ax(3) ax(4)]);
      line([right_freq(2) right_freq(2)],[ax(3) ax(4)]);
      line([right_freq(1) right_freq(2)],[right_offset right_offset]);
      
      % str = sprintf('time = %d s\n', round(mean(time)));
      % title(str);
      xlabel('frequency (Hz)');
      ylabel('power');
      
      subplot(3,2,1);
      plot(time,dat(1,:));
      axis('tight');
      axis([min(time) max(time) -300 300]);
      ax = axis;
      line([ax(1) ax(2)],[left_thresh_ampl left_thresh_ampl],'color','red');
      line([ax(1) + left_thresh_time(1)/hdr.Fs ax(1) + left_thresh_time(1)/hdr.Fs],[-300 300],'color','green');
      line([ax(1) + left_thresh_time(2)/hdr.Fs ax(1) + left_thresh_time(2)/hdr.Fs],[-300 300],'color','green');
      % axis('tight');
      ylabel('amplitude (uV)');
      xlabel(['time: ' int2str(cfg.blocksize*cfg.blockmem) 's']);
      grid;
      
      subplot(3,2,2);
      plot(time,dat(2,:));
      axis('tight');
      axis([min(time) max(time) -300 300]);
      % axis('tight');
      ylabel('amplitude (uV)');
      xlabel(['time: ' int2str(cfg.blocksize*cfg.blockmem) 's']);
      grid;
      
      subplot(3,2,5);
      h = surf(squeeze(TFR(1,:,:)));
      set(h, 'EdgeColor', 'none');
      shading interp
      ax = axis;
      axis([ax(1) ax(2) cfg.foilim(1) cfg.foilim(2) 0 ax(6)]);
      view(110,45);
      ylabel('Frequency (Hz)');
      xlabel('Time');
      zlabel('Power');
      box off;
      
      subplot(3,2,6);
      h = surf(squeeze(TFR(2,:,:)));
      set(h, 'EdgeColor', 'none');
      shading interp
      ax = axis;
      axis([ax(1) ax(2) cfg.foilim(1) cfg.foilim(2) 0 ax(6)]);
      view(110,45);
      ylabel('Frequency (Hz)');
      xlabel('Time');
      zlabel('Power');
      box off;
      
      if left_thresh_ampl < 0
        
        if min((dat(1,left_thresh_time(1):left_thresh_time(2)))) < left_thresh_ampl
          if beatdrum == true
            % midiOut('+',10,64,127);
            beatdrum = false;
          else
            % midiOut('+',10,31,127);
            beatdrum = true;
          end
        else
          midiOut('.',1);
        end;
        
      elseif max((dat(left_thresh_time(1):left_thresh_time(2)))) > left_thresh_ampl
        if beatdrum == true
          % midiOut('+',10,64,127);
          beatdrum = false;
        else
          % midiOut('+',10,31,127);
          beatdrum = true;
        end
      else
        % midiOut('.',1);
      end;
      
      % schemerlamp.setLevel(round(TFR(1,60,end) / mean(TFR(1,60,:)))*5);
      volume_right = round((mean(TFR(2,right_freq,end)) - right_offset) * right_mult);
      
      % midiOut(uint8([175+1,7, volume_left]));
      % midiOut(uint8([16*14+1-1,0, volume_left])); % ptich
      
      midiOut(uint8([175+2,7, volume_right]));
      midiOut(uint8([16*14+2-1,0, volume_right])); % ptich
      
      % schemerlamp.setLevel(9);

      % force an update of the figure
      drawnow
      
    end % if enough blocks
  end % if enough new samples
end % while true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this is a dropin for the midiOut function, which is not available
% on Mac OS X, and probably also not on Linux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = dummyfunction(varargin)
