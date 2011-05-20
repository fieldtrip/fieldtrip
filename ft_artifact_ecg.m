function [cfg, artifact] = ft_artifact_ecg(cfg, data)

% FT_ARTIFACT_ECG performs a peak-detection on the ECG-channel. The
% heart activity can be seen in the MEG data as an MCG artifact and
% can be removed using independent component analysis.
%
% Use as
%   [cfg, artifact] = ft_artifact_ecg(cfg)
%   required configuration options:
%   cfg.dataset or both cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.artfctdef.ecg.channel = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.ecg.pretim  = 0.05; pre-artifact rejection-interval in seconds
%   cfg.artfctdef.ecg.psttim  = 0.3;  post-artifact rejection-interval in seconds
%   cfg.artfctdef.ecg.method  = 'zvalue'; peak-detection method
%   cfg.artfctdef.ecg.cutoff  = 3; peak-threshold
%   cfg.artfctdef.ecg.inspect = Nx1 list of channels which will be shown in a QRS-locked average
%   cfg.continuous            = 'yes' or 'no' whether the file contains continuous data
%
% The output artifact variable is an Nx2-matrix, containing the
% begin and end samples of the QRST-complexes in the ECG.
%
% See also FT_REJECTARTIFACT

% Undocumented local options:
% cfg.ft_datatype

% Copyright (c) 2005, Jan-Mathijs Schoffelen
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters for eog artifacts if necessary.
if ~isfield(cfg,'artfctdef'),            cfg.artfctdef               = [];            end
if ~isfield(cfg.artfctdef,'ecg'),        cfg.artfctdef.ecg           = [];            end
if ~isfield(cfg.artfctdef.ecg,'channel'),cfg.artfctdef.ecg.channel   = {'ECG'};       end
if ~isfield(cfg.artfctdef.ecg,'method'), cfg.artfctdef.ecg.method    = 'zvalue';      end
if ~isfield(cfg.artfctdef.ecg,'cutoff'), cfg.artfctdef.ecg.cutoff    = 3;             end
if ~isfield(cfg.artfctdef.ecg,'padding'),cfg.artfctdef.ecg.padding   = 0.5;           end
if ~isfield(cfg.artfctdef.ecg,'inspect'),cfg.artfctdef.ecg.inspect   = {'MLT' 'MRT'}; end
if ~isfield(cfg.artfctdef.ecg,'pretim'), cfg.artfctdef.ecg.pretim    = 0.05;          end
if ~isfield(cfg.artfctdef.ecg,'psttim'), cfg.artfctdef.ecg.psttim    = 0.3;           end
if ~isfield(cfg.artfctdef.ecg,'mindist'), cfg.artfctdef.ecg.mindist  = 0.5;           end
if ~isfield(cfg, 'headerformat'),         cfg.headerformat           = [];            end
if ~isfield(cfg, 'dataformat'),           cfg.dataformat             = [];            end

cfg.artfctdef = ft_checkconfig(cfg.artfctdef, 'renamed',    {'blc', 'demean'});
cfg.artfctdef = ft_checkconfig(cfg.artfctdef, 'renamed',    {'blcwindow' 'baselinewindow'});


% for backward compatibility
if isfield(cfg.artfctdef.ecg,'sgn')
  cfg.artfctdef.ecg.channel = cfg.artfctdef.ecg.sgn;
  cfg.artfctdef.ecg         = rmfield(cfg.artfctdef.ecg, 'sgn');
end

if ~strcmp(cfg.artfctdef.ecg.method, 'zvalue'),
  error('this method is not applicable');
end

if nargin == 1,
  cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  hdr = ft_read_header(cfg.headerfile,'headerformat', cfg.headerformat);
  trl = cfg.trl;
elseif nargin == 2,
  data = ft_checkdata(data, 'hastrialdef', 'yes');
  cfg  = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr  = ft_fetch_header(data);
  if isfield(data, 'sampleinfo'), 
    trl = data.sampleinfo;
    for k = 1:numel(data.trial)
      trl(k,3) = time2offset(data.time{k}, data.fsample);
    end
  else
    error('the input data does not contain a valid description of the sampleinfo');
  end  
end
artfctdef     = cfg.artfctdef.ecg;
padsmp        = round(artfctdef.padding*hdr.Fs);
ntrl          = size(trl,1);
artfctdef.trl = trl;
artfctdef.channel = ft_channelselection(artfctdef.channel, hdr.label);
artfctdef.demean  = 'yes';
sgnind        = match_str(hdr.label, artfctdef.channel);
numecgsgn     = length(sgnind);
fltpadding    = 0;

if numecgsgn<1
  error('no ECG channels selected');
elseif numecgsgn>1
  error('only one ECG channel can be selected');
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
end

% read in the ecg-channel and do demean and squaring
if nargin==2,
  tmpcfg = [];
  tmpcfg.channel = artfctdef.channel;
  ecgdata = ft_preprocessing(tmpcfg, data);
  ecg     = ecgdata.trial;
end

for j = 1:ntrl
  if nargin==1,
    ecg{j} = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(j,1), 'endsample', trl(j,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  ecg{j} = preproc(ecg{j}, artfctdef.channel, hdr.Fs, artfctdef, [], fltpadding, fltpadding);
  ecg{j} = ecg{j}.^2;
end

tmp   = cell2mat(ecg);
stmp  =  std(tmp, 0, 2);
mtmp  = mean(tmp, 2);
Nsmp  = max(trl(:,2));
trace = zeros(1,Nsmp);

% standardise the ecg
for j = 1:ntrl
  trace(trl(j,1):trl(j,2)) = (ecg{j}-mtmp)./stmp;
end

accept = 0;
while accept == 0,
  h = figure;
  plot(trace);zoom;
  hold on;
  plot([1 Nsmp], [artfctdef.cutoff artfctdef.cutoff], 'r:');
  hold off;
  xlabel('samples');
  ylabel('zscore');

  fprintf(['\ncurrent  ',artfctdef.method,' threshold = %1.3f'], artfctdef.cutoff);
  response = input('\nkeep the current value (y/n) ?\n','s');
  switch response
    case 'n'
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
    case 'y'
      oldcutoff = artfctdef.cutoff;
      accept = 1;
    otherwise
      warning('unrecognised response, assuming no');
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
  end;
  close
end

% detect peaks which are at least half a second apart and store
% the indices of the qrs-complexes in the artifact-configuration
mindist       = round(cfg.artfctdef.ecg.mindist.*hdr.Fs);
[pindx, pval] = peakdetect2(trace, artfctdef.cutoff, mindist);
%sel           = find(standardise(pval,2)<2);
%pindx         = pindx(sel);
%pval          = pval(sel);
artfctdef.qrs = pindx;

%---------------------------------------
% create trials for qrs-triggered average
trl = [];
trl(:,1) = pindx(:) - round(artfctdef.padding*(hdr.Fs))  ;
trl(:,2) = pindx(:) + round(artfctdef.padding*(hdr.Fs))-1;
trl(:,3) = -round(artfctdef.padding*(hdr.Fs));
trl(trl(:,1)<1,:) = [];
trl(trl(:,2)>hdr.nSamples.*hdr.nTrials,:) = [];
%------------

% ---------------------
% qrs-triggered average
% FIXME, at present this only works for continuous data: the assumption can
% be made that all trials are equally long.
sgn    = ft_channelselection(artfctdef.inspect, hdr.label);
megind = match_str(hdr.label, sgn);
sgnind = [megind(:); sgnind];
dat    = zeros(length(sgnind), trl(1,2)-trl(1,1)+1);
ntrl   = size(trl,1);

if ~isempty(sgnind)
  ntrlok = 0;
  for j = 1:ntrl
    fprintf('reading and preprocessing trial %d of %d\n', j, ntrl);
    if nargin==1,
      dum = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(j,1), 'endsample', trl(j,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
      dat = dat + ft_preproc_baselinecorrect(dum);
      ntrlok = ntrlok + 1;
    elseif nargin==2,
      dum = ft_fetch_data(data, 'header', hdr, 'begsample', trl(j,1), 'endsample', trl(j,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'docheck', 0);
      if any(~isfinite(dum(:))),
      else
        ntrlok = ntrlok + 1;
        dat    = dat + ft_preproc_baselinecorrect(dum);
      end
    end
  end
end

dat  = dat./ntrlok;
time = offset2time(trl(1,3), hdr.Fs, size(dat,2));
tmp  = dat(1:end-1,:);
mdat = max(abs(tmp(:)));

acceptpre = 0;
acceptpst = 0;
while acceptpre == 0 || acceptpst == 0,
  h = figure;
  subplot(2,1,1); plot(time, dat(end, :));
  abc = axis;
  axis([time(1) time(end) abc(3:4)]);
  subplot(2,1,2);
  axis([time(1) time(end) -1.1*mdat 1.1*mdat]);
  xpos   = -artfctdef.pretim;
  ypos   = -1.05*mdat;
  width  = artfctdef.pretim + artfctdef.psttim;
  height = 2.1*mdat;
  rectangle('Position', [xpos ypos width height], 'FaceColor', 'r');
  hold on; plot(time, dat(1:end-1, :), 'b');

  if acceptpre == 0,
    fprintf(['\ncurrent pre-peak interval = %1.3f'], artfctdef.pretim);
    response = input('\nkeep the current value (y/n) ?\n','s');
    switch response
      case 'n'
        oldpretim = artfctdef.pretim;
        artfctdef.pretim = input('\nenter new value \n');
      case 'y'
        oldpretim = artfctdef.pretim;
        acceptpre = 1;
      otherwise
        warning('unrecognised response, assuming no');
        oldpretim = artfctdef.pretim;
    end
  end
  if acceptpst == 0 && acceptpre == 1,
    fprintf(['\ncurrent post-peak interval = %1.3f'], artfctdef.psttim);
    response = input('\nkeep the current value (y/n) ?\n','s');
    switch response
      case 'n'
        oldpsttim = artfctdef.psttim;
        artfctdef.psttim = input('\nenter new value \n');
      case 'y'
        oldpsttim = artfctdef.psttim;
        acceptpst = 1;
      otherwise
        warning('unrecognised response, assuming no');
        oldpsttim = artfctdef.psttim;
    end
  end
  close
end

artifact(:,1) = trl(:,1) - trl(:,3) - round(artfctdef.pretim*hdr.Fs);
artifact(:,2) = trl(:,1) - trl(:,3) + round(artfctdef.psttim*hdr.Fs);

% remember the details that were used here
cfg.artfctdef.ecg          = artfctdef;
cfg.artfctdef.ecg.artifact = artifact;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();
