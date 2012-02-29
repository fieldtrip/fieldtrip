function [cfg, artifact] = ft_artifact_manual(cfg);

% THIS FUNCTION IS DEPRECIATED, USE FT_REJECTVISUAL INSTEAD
%
% FT_ARTIFACT_MANUAL allows the user to detect artifacts manually using visual
% inspection.
%
% Use as:
%   [cfg, artifact] = ft_artifact_manual(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
%
% The configuration should also contain:
%   cfg.artfctdef.manual.channel        = cell-array with channels to be displayed.
%   (Be careful not to specify to much channels because the function then will be very slow.)
%   cfg.continuous                      = 'yes' or 'no' whether the file contains continuous data
%
% You can specify:
%   cfg.artfctdef.manual.pretrialtime   = time shown before trialstart (default 0)
%   cfg.artfctdef.manual.posttrialtime  = time shown after trialend    (default 0)
%   cfg.artfctdef.manual.fft            = 'no' (default) or 'yes' turns on FFT window
%   cfg.artfctdef.manual.padding        = 'no' (default) or FFT-padding in seconds
%
% The FFT will be executed on the complete time interval including pre and post
% trial times.
%
%   cfg.artfctdef.manual.timeaxrelative = 'yes' (default) or 'no'.
%
% Set to yes defines the time axes relative to trialstart. Set to no defines it
% relative to the beginning of the experiment.
%
%   cfg.artfctdef.manual.demean    = 'no' (default) or 'yes' apply baseline correction
%   cfg.artfctdef.manual.bpfilter  = 'no' (default) or 'yes' apply bandpass filter
%   cfg.artfctdef.manual.bpfreq    = [0.3 30] in Hz
%   cfg.artfctdef.manual.bpfiltord = 2
%
% See also FT_REJECTARTIFACT, FT_REJECTVISUAL

% Undocumented local options:
% cfg.artfctdef.manual.maxnumberofchannels = 20 (default)

% Copyright (C) 2004, Geerten Kramer, FCDC
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters if necessary.
if ~isfield(cfg, 'artfctdef'),                          cfg.artfctdef                            = [];       end
if ~isfield(cfg.artfctdef,'manual'),                    cfg.artfctdef.manual                     = [];       end
if ~isfield(cfg.artfctdef.manual,'zscale'),             cfg.artfctdef.manual.zscale              = 100;      end
if ~isfield(cfg.artfctdef.manual,'fft'),                cfg.artfctdef.manual.fft                 ='no';      end
if ~isfield(cfg.artfctdef.manual,'padding'),            cfg.artfctdef.manual.padding             ='no';      end
if ~isfield(cfg.artfctdef.manual,'pretrialtime'),       cfg.artfctdef.manual.pretrialtime        = 0;        end
if ~isfield(cfg.artfctdef.manual,'posttrialtime'),      cfg.artfctdef.manual.posttrialtime       = 0;        end
if ~isfield(cfg.artfctdef.manual,'timeaxrelative'),     cfg.artfctdef.manual.timeaxrelative      = 'yes';    end
if ~isfield(cfg.artfctdef.manual,'maxnumberofchannels'),cfg.artfctdef.manual.maxnumberofchannels = 20;       end
if ~isfield(cfg, 'headerformat'),                       cfg.headerformat                         = [];       end
if ~isfield(cfg, 'dataformat'),                         cfg.dataformat                           = [];       end

% for backward compatibility
if isfield(cfg.artfctdef.manual,'sgn')
  cfg.artfctdef.manual.channel = cfg.artfctdef.manual.sgn;
  cfg.artfctdef.manual         = rmfield(cfg.artfctdef.manual, 'sgn');
end
cfg.artfctdef = ft_checkconfig(cfg.artfctdef, 'renamed',    {'blc', 'demean'});
cfg.artfctdef = ft_checkconfig(cfg.artfctdef, 'renamed',    {'blcwindow' 'baselinewindow'});

if ~isfield(cfg.artfctdef.manual,'channel'),
  % set an unusual default because all crashes the program.
  cfg.artfctdef.manual.channel = [];
end

% read the header and do some preprocessing on the configuration
fprintf('Reading raw data...');
cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
cfg.artfctdef.manual.channel=ft_channelselection(cfg.artfctdef.manual.channel, hdr.label);
cfg.artfctdef.manual.trl=cfg.trl;
if(isempty(cfg.artfctdef.manual.channel))
  error(sprintf('\nNo channels selected for artifact_manual!\nSelect at least one channel in cfg.artfctdef.manual.channel'));
end;
channelindx=match_str(hdr.label, cfg.artfctdef.manual.channel);
ntrial=size(cfg.trl, 1);
nch=length(channelindx);
if(nch<1)
  error(sprintf('\nNo channels selected for artifact_manual!\nSelect at least one channel in cfg.artfctdef.manual.channel'));
elseif(nch>cfg.artfctdef.manual.maxnumberofchannels)
  error(sprintf('\nMore than %i channels selected in cfg.artfctdef.manual.channel',cfg.artfctdef.manual.maxnumberofchannels));
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
end

show=ft_read_data(cfg.datafile, 'header', hdr, 'begsample', 1, 'endsample', hdr.nTrials*hdr.nSamples, 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
show=show';

N=length(show);
S=floor(N./hdr.Fs); % In seconds
x=1:N;
x=x./hdr.Fs; % In seconds.

fprintf(' done.\n');
fprintf('Processing raw data for artifact_manual...');

% these elements are stored inside the figure so that the callback routines can modify them
dat.RejMarkList=zeros(length(cfg.trl),1);
dat.ResRejMarkList=zeros(length(cfg.trl),1);
dat.RejCount=0;
dat.ResRejCount=0;
dat.stop=0;
dat.trln=1;
dat.numtrl=length(cfg.trl);
dat.FFT=strcmp(cfg.artfctdef.manual.fft,'yes')|strcmp(cfg.artfctdef.manual.fft,'Yes')...
  |strcmp(cfg.artfctdef.manual.fft,'on')|strcmp(cfg.artfctdef.manual.fft,'On');
if(strcmp(cfg.artfctdef.manual.padding,'no')|strcmp(cfg.artfctdef.manual.padding,'No')...
    |strcmp(cfg.artfctdef.manual.padding,'Off')|strcmp(cfg.artfctdef.manual.padding,'off'))
  cfg.artfctdef.manual.padding=[];
end;
dat.IfRelOn=strcmp(cfg.artfctdef.manual.timeaxrelative,'yes')|strcmp(cfg.artfctdef.manual.timeaxrelative,'Yes')...
  |strcmp(cfg.artfctdef.manual.timeaxrelative,'on')|strcmp(cfg.artfctdef.manual.timeaxrelative,'On');

ival=cell(dat.numtrl,1);
dataX=cell(dat.numtrl,1);
dataY=cell(dat.numtrl,1);
dataFx=cell(dat.numtrl,1);
dataFy=cell(dat.numtrl,1);
dataXLim=cell(dat.numtrl,1);
dataYLim=cell(dat.numtrl,1);

% first we calculate everything and put it into varables.
for(i=1:dat.numtrl)
  begpadding = round(cfg.artfctdef.manual.pretrialtime.*hdr.Fs);
  endpadding = round(cfg.artfctdef.manual.posttrialtime.*hdr.Fs);
  ival{i}=(cfg.trl(i,1)-begpadding):(cfg.trl(i,2)+endpadding);
  dataY{i}=show(ival{i},:);
 
  % don't remove the padding
  [dataY{i}, dumlab, dumtime, dumcfg] = preproc(dataY{i}', cfg.artfctdef.manual.channel, offset2time(cfg.trl(i,3)-begpadding, hdr.Fs, size(dataY{i}',2)), cfg.artfctdef.manual, 0, 0);
  dataY{i} = dataY{i}';

  try, if(~rem(i,fix(dat.numtrl/10)))fprintf('.');end; end
end;
clear show; % Now, we don't need the complete dataset anymore...
for(i=1:dat.numtrl) % we go on with the calculations ....
  dataX{i}=x(ival{i})-dat.IfRelOn.*(x(ival{i}(1))+cfg.artfctdef.manual.pretrialtime);

  dataYLim{i}=[];
  for(j=1:nch)
    mini=min(dataY{i}(:,j));
    maxi=max(dataY{i}(:,j));

    dataYLim{i}=[dataYLim{i};[mini-.03.*(maxi-mini),maxi+.03.*(maxi-mini)]];
    dataXLim{i}=[min(dataX{i}),max(dataX{i})];
  end;

  if(dat.FFT) % if FFT is on, calculate the FFT of the trials too...
    tmp=[];
    for(k=1:nch) tmp=[tmp,dataY{i}(:,k)-mean(dataY{i}(:,k),1)];end;
    dataFy{i}=abs(fft(tmp,cfg.artfctdef.manual.padding,1));
    tmp=floor((max(ival{i})-min(ival{i}))/2);
    dataFx{i}=(1:tmp)./tmp.*hdr.Fs./2;
    dataFy{i}=dataFy{i}(1:length(dataFx{i}),:);
  end;
  if(~rem(i,fix(dat.numtrl/10)))fprintf('.');end;
end;

fprintf(' done.\n');

thisfig('fig_trace');% See the help of thisfig
set(fig_trace,'name','Traces');
if(dat.FFT);
  thisfig('fig_spec');
  set(fig_spec,'name','Spectra of traces');
end;
gt=[5 40 5 25 0 50 5 60 5 60 5 60 5 60 5 60 5 25 0 25];
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:1)),5,gt(2),18],'String','Done','Callback',@stop);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:3)),5,gt(4),18],'String','<','Callback',@prevtrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:5)),5,gt(6),18],'String','>','Callback',@nexttrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:7)),5,gt(8),18],'String','Reject >','Callback',@reject_next);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:9)),5,gt(10),18],'String','Reject','Callback',@reject);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:11)),5,gt(12),18],'String','Accept','Callback',@undo);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:13)),5,gt(14),18],'String','Accept All','Callback',@undoAll);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:15)),5,gt(16),18],'String','Redo All','Callback',@redoAll);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:17)),5,gt(18),18],'String','<<','Callback',@pptrial);
uicontrol(fig_trace,'units','pixels','position',[sum(gt(1:19)),5,gt(20),18],'String','>>','Callback',@nntrial);
set(fig_trace, 'KeyPressFcn',@keypress);
% set(get(fig_trace, 'children'),'KeyPressFcn',@keypress);
show_help;

thisfig('fig_trace'); % See the help of thisfig.

while(ishandle(fig_trace))

  if(dat.FFT) % Plots FFT when requested
    if(exist('HF')) % this is for remembering the zoom of the FFT-window
      if(ishandle(HF)) % while zapping trough the data.
        saveXLim=get(HF,'XLim');
        saveYLim=get(HF,'YLim');
      end;
    end;
    thisfig('fig_spec');
    plot(dataFx{dat.trln},dataFy{dat.trln});
    HF=get(fig_spec,'CurrentAxes');
    if(exist('saveXLim'))set(HF,'XLim',saveXLim);end;
    if(exist('saveYLim'))set(HF,'YLim',saveYLim);end;
    xlabel('Frequency [Hz]');
    ylabel('Power [A.U.]');
    thisfig('fig_trace');
  end;

  for(j=1:nch) % This loop fils the traces figure with the traces.
    H=subplot(nch,1,j);
    plot(dataX{dat.trln},dataY{dat.trln}(:,j),'b-');
    line([dataXLim{1}(1)+cfg.artfctdef.manual.pretrialtime dataXLim{1}(1)+cfg.artfctdef.manual.pretrialtime], dataYLim{dat.trln}(j,:), 'color', 'g');
    line([dataXLim{1}(2)-cfg.artfctdef.manual.posttrialtime dataXLim{1}(2)-cfg.artfctdef.manual.posttrialtime], dataYLim{dat.trln}(j,:), 'color', 'r');
    set(H,'YLim',dataYLim{dat.trln}(j,:) + [-eps +eps]);
    set(H,'XLim',dataXLim{dat.trln});
    xlabel('Time [s]','position',...
      [dataXLim{dat.trln}(2),dataYLim{dat.trln}(j,1)-(dataYLim{dat.trln}(j,2)-dataYLim{dat.trln}(j,1)).*.05]...
      ,'HorizontalAlignment','left');
    ylabel(sprintf('%s', cfg.artfctdef.manual.channel{j}));
  end;

  if(dat.RejMarkList(dat.trln)==0)
    tiet=sprintf('trial %i: ACCEPT',dat.trln);
  else
    tiet=sprintf('trial %i: REJECT',dat.trln);
  end;
  H=subplot(nch,1,1);
  title(tiet);

  guidata(fig_trace,dat);
  uiwait;
  if(ishandle(fig_trace)) dat=guidata(fig_trace);end;

  if(dat.stop)break;end;
end

if(~ishandle(fig_trace))
  error('Figure closed unexpectedly, no manual artifacts marked.');
else
  close(fig_trace)
end;

if(dat.FFT)
  if(ishandle(fig_spec))close(fig_spec);end;
end;

fprintf('Rejected %i trials.\n',dat.RejCount);
fprintf('Exported %i trials.\n',length(cfg.trl)-dat.RejCount);

artifact=cfg.trl(find((dat.RejMarkList)),[1,2]);

% remember the details that were used here
cfg.artfctdef.manual.artifact = artifact;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the SUBFUNCTIONS start that implement the gui callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=nexttrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln<dat.numtrl)
  dat.trln=dat.trln + 1;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;

function varargout=prevtrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln > 1),
  dat.trln=dat.trln - 1;
else
  dat.trln=dat.numtrl;
end;
guidata(h,dat);
uiresume;

function varargout=stop(h, eventdata, handles, varargin)
dat=guidata(h);
dat.stop=1;
guidata(h,dat);
uiresume;

function varargout=reject(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)==0)
  dat.RejCount=dat.RejCount+1;
  dat.RejMarkList(dat.trln)=1;
  if(~dat.ResRejMarkList(dat.trln))
    dat.ResRejCount=dat.ResRejCount+1;
    dat.ResRejMarkList(dat.trln)=1;
    fprintf('Rejected trial %i!\n',dat.trln);
  else
    fprintf('Rejected trial %i again!\n',dat.trln);
  end;
else
  fprintf('Trial %i is already marked for rejection.\n',dat.trln);
end;
guidata(h,dat);
uiresume;


function varargout=reject_next(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)==0)
  dat.RejCount=dat.RejCount+1;
  dat.RejMarkList(dat.trln)=1;
  if(~dat.ResRejMarkList(dat.trln))
    dat.ResRejCount=dat.ResRejCount+1;
    dat.ResRejMarkList(dat.trln)=1;
    fprintf('Rejected trial %i!\n',dat.trln);
  else
    fprintf('Rejected trial %i again!\n',dat.trln);
  end;
else
  fprintf('Trial %i is already marked for rejection.\n',dat.trln);
end;
if(dat.trln < dat.numtrl)
  dat.trln=dat.trln + 1;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;

function varargout=undo(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.RejMarkList(dat.trln)>0)
  dat.RejCount=dat.RejCount-1;
  dat.RejMarkList(dat.trln)=0;
  fprintf('Rejected trial %i recovered.\n',dat.trln);
end;
guidata(h,dat);
uiresume;

function varargout=undoAll(h, eventdata, handles, varargin)
dat=guidata(h);
dat.RejCount=0;
dat.RejMarkList=zeros(dat.numtrl,1);
fprintf('All rejected trials recovered.\n');
guidata(h,dat);
uiresume;

function varargout=redoAll(h, eventdata, handles, varargin)
dat=guidata(h);
dat.RejMarkList=dat.ResRejMarkList;
dat.RejCount=dat.ResRejCount;
fprintf('All recovered trials rejected again.\n');
guidata(h,dat);
uiresume;

function varargout=nntrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln<(dat.numtrl-10))
  dat.trln=dat.trln + 10;
else
  dat.trln=1;
end;
guidata(h,dat);
uiresume;


function varargout=pptrial(h, eventdata, handles, varargin)
dat=guidata(h);
if(dat.trln > 10),
  dat.trln=dat.trln - 10;
else
  dat.trln=dat.numtrl;
end;
guidata(h,dat);
uiresume;

function keypress(h, eventdata, handles, varargin)
dat=guidata(gcbf);
key=get(gcbf, 'CurrentCharacter');
if  key
  switch key
    case 28     % arrow left
      if(dat.trln > 1),
        dat.trln=dat.trln - 1;
      else
        dat.trln=dat.numtrl;
      end;
    case 29     % arrow right
      if(dat.trln<dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;
      %    case 30     % arrow up
      %      dat.zscale=data.zscale*1.5;
      %    case 31     % arrow down
      %      dat.zscale=data.zscale/1.5;
      %    case 'a'
      %      dat.reject(data.trlop)=0;
      %    case 'r'
      %      dat.reject(data.trlop)=1;
    case {'a', 'A'}
      if(dat.RejMarkList(dat.trln)==1)
        dat.RejCount=dat.RejCount-1;
        dat.RejMarkList(dat.trln)=0;
        fprintf('Accepted trial %i\n',dat.trln);
      end
      if(dat.trln < dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;
    case {'0', 'r', 'R'}
      if(dat.RejMarkList(dat.trln)==0)
        dat.RejCount=dat.RejCount+1;
        dat.RejMarkList(dat.trln)=1;
        if(~dat.ResRejMarkList(dat.trln))
          dat.ResRejCount=dat.ResRejCount+1;
          dat.ResRejMarkList(dat.trln)=1;
          fprintf('Rejected trial %i! \n',dat.trln);
        else
          fprintf('Rejected trial %i again! \n',dat.trln);
        end;
      else
        fprintf('Trial %i is already marked for rejection.\n',dat.trln);
      end;
      if(dat.trln < dat.numtrl)
        dat.trln=dat.trln + 1;
      else
        dat.trln=1;
      end;

    otherwise
      show_help;
  end
end
set(gcbf, 'CurrentCharacter', '-');
guidata(gcbf, dat);
uiresume(h);

function show_help;
fprintf('You can use the buttons in the window to navigate through the trials,   \n');
fprintf('reject trials or undo and redo your rejection selection. If you want to \n');
fprintf('use the keyboard, you first have to click in the figure every time you  \n');
fprintf('pressed a button with the mouse.                                        \n');
fprintf('                                                                        \n');
fprintf('Use the left and right arrow to walk through the trials                 \n');
fprintf('Press "a" to accept, "r" to reject.                                     \n');

function thisfig(NameOfHandle);
assignin('caller','This_Name_Will_Never_be_Used_gk',NameOfHandle);
if(~evalin('caller','exist(This_Name_Will_Never_be_Used_gk)')|~ishandle(evalin('caller',NameOfHandle)))
  H=figure;
  assignin('caller',NameOfHandle,H);
else
  figure(evalin('caller',NameOfHandle));
end;
evalin('caller','clear This_Name_Will_Never_be_Used_gk;');
