function [grad, elec] = mne2grad(hdr, dewar, coilaccuracy, coildeffile)

% MNE2GRAD converts a header from a fif file that was read using the MNE toolbox into
% a gradiometer structure that can be understood by the FieldTrip low-level forward
% and inverse routines.
%
% Use as
%   [grad, elec] = mne2grad(hdr, dewar, coilaccuracy)
% where
%   dewar        = boolean, whether to return it in dewar or head coordinates (default = false, i.e. head coordinates)
%   coilaccuracy = empty or a number (default = [])
%   coildeffile  = empty or a filename of a valid coil_def.dat file
%
% See also CTF2GRAD, BTI2GRAD

%%%%%%%%%%%%%%%%%% BEGIN REVISION HISTORY %%%%%%%%%%%%%%%%%%
%
% Robert Oostenveld 22/01/2016 reimplemented construction of integration
% points on basis of coil_def.dat file from MNE
%
% Laurence Hunt 03/12/2008 (with thanks to Joachim Gross's original script
% based on fiff_access). lhunt@fmrib.ox.ac.uk
%
% Teresa Cheung 09/24/2011 revision to Laurence Hunt's script. The coil
% geometry has been revised to better reflect the coil dimensions of the
% Vectorview (306) planar gradiometers and magnetometers. Coil geometry for
% the Neuromag-122 planar gradiometer is different and is not read with
% this script (coil_type should equal 2)
%
%%%%%%%%%%%%%%%%%% END REVISION HISTORY %%%%%%%%%%%%%%%%%%

% Copyrights (C) 2016, Robert Oostenveld
% Copyrights (C) 2011, Teresa Cheung
% Copyrights (C) 2008, Laurence Hunt
% Copyrights (C) 2005, Joachim Gross
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
% FieldTrip is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FieldTrip is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<2 || isempty(dewar)
  dewar = false;
end

if nargin<3 || isempty(coilaccuracy)
  % if empty it will use the original code
  % otherwise it will use the specified accuracy coil definition from the MNE coil_def.dat
  coilaccuracy = [];
end

if nargin<4 || isempty(coildeffile)
  % if coilaccuracy is not empty, it will use the default coil_def.dat
  % file, which is in external/mne
  coildeffile = [];
end

% orig = fiff_read_meas_info(filename);
% orig = fiff_setup_read_raw(filename);

if isfield(hdr, 'orig')
  orig = hdr.orig; % use the original MNE header, not the FieldTrip header
else
  orig = hdr; % assume that it is the original header
end

% start with empty gradiometer definition
grad = [];

% start with empty electrode definition
elec = [];

% begin by transforming all channel locations into the desired coordinate system, if possible
if ~dewar
  if ~isempty(orig.dev_head_t)
    ft_info('Transforming MEG channels according to the device to head transform from the data');
    orig.chs = fiff_transform_meg_chs(orig.chs,orig.dev_head_t);
    %orig.chs = fiff_transform_eeg_chs(orig.chs,orig.dev_head_t); % EEG channels are normally stored in head coordinates anyway, but what the heck
    if ~isempty(orig.ctf_head_t)
      orig.head_ctf_t.trans = inv(orig.ctf_head_t.trans);
      orig.head_ctf_t.from  = orig.ctf_head_t.to;
      orig.head_ctf_t.to    = orig.ctf_head_t.from;
      ft_info('Transforming MEG channels according to the head to ctf transform from the data');
      orig.chs = fiff_transform_meg_chs(orig.chs,orig.head_ctf_t);
    end
  else
    ft_warning('No device to head transform available in fif file');
    ft_warning('MEG channels will likely have coordinates in device frame, not head frame');
  end
else
  if ~isempty(orig.dev_head_t)
    % compute the transformation from head to device
    orig.head_dev_t.trans = inv(orig.dev_head_t.trans);
    orig.head_dev_t.from  = orig.dev_head_t.to;
    orig.head_dev_t.to    = orig.dev_head_t.from;
    orig.chs = fiff_transform_meg_chs(orig.chs,orig.head_dev_t); % MEG channels are normally stored in dewar coordinates anyway, but what the heck
    orig.chs = fiff_transform_eeg_chs(orig.chs,orig.head_dev_t);
  else
    ft_warning('No device to head transform available in fif file');
    ft_warning('EEG channels will likely have coordinates in head frame, not device frame');
  end
end

if ~isempty(coilaccuracy)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the coil definitions from the MNE coil_def.dat file
  % these allow for varying accuracy which is specified by
  % coilaccuracy = 0, 1 or 2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ft_hastoolbox('mne', 1);
  if isempty(coildeffile)
    [ftver, ftpath] = ft_version;
    coildeffile     = fullfile(ftpath, 'external', 'mne', 'coil_def.dat');
  end
  def = mne_load_coil_def(coildeffile);
  
  grad.chanpos = nan(length(orig.chs),3);
  grad.chanori = nan(length(orig.chs),3);
  grad.label   = cell(length(orig.chs),1);
  grad.tra      = nan(length(orig.chs),0);
  grad.chantype = repmat({'unknown'}, size(grad.label));
  grad.chanunit = repmat({'unknown'}, size(grad.label));
  
  k = 1; % counter for the coils
  for i=1:length(orig.chs)
    thisdef = def([def.id]==orig.chs(i).coil_type & [def.accuracy]==coilaccuracy);
    if isempty(thisdef)
      continue;
    end
    weight = thisdef.coildefs(:,1);
    pos = thisdef.coildefs(:,2:4);
    ori = thisdef.coildefs(:,5:7);
    H = orig.chs(i).coil_trans;
    pos = ft_warp_apply(H, pos);
    T = H;
    R = H;
    T(1:3,1:3) = 0; % remove the rotation, keep the translation
    R(1:3,4)   = 0; % remove the translation, keep the rotation
    ori = ft_warp_apply(R, ori); % only apply the rotation
    for j=1:thisdef.num_points
      grad.coilpos(k,:) = pos(j,:);
      grad.coilori(k,:) = ori(j,:);
      grad.tra(i,k)     = weight(j);
      k = k + 1;
    end
    grad.label{i} = orig.ch_names{i};
    grad.chanpos(i,:) = T(1:3,4);
    grad.chanori(i,:) = ft_warp_apply(R, [0 0 1]);
  end
  
  grad.label = grad.label(:);
  grad.unit  = 'm'; % the coil_def.dat file is in meter
  
  % all the chs.kinds and chs.coil_types are obtained from the MNE manual, p.210-211, 
  % and https://github.com/mne-tools/mne-matlab/blob/master/matlab/fiff_define_constants.m
  kind  = [orig.chs.kind];
  type  = [orig.chs.coil_type];
  logno = [orig.chs.logno]; 

  megplanars = [2 3011 3012 3013 3014 3015];
  megmags    = [3021 3022 3023 3024 3025 4001 8001 8002 8101 8201];
  meggrads   = [4002 5001 6001];
  megaxials  = 7001;
  refmags    = [4003 5002 6002];
  refgrads   = [4004 4005 5003 5004];

  grad.chantype(kind==1)   = {'meg'};
  grad.chantype(kind==301) = {'ref'};

  % overwrite with more specific values
  grad.chantype(kind==1 & ismember(type, megplanars)) = {'megplanar'}; % planar gradiometer
  grad.chantype(kind==1 & ismember(type, megmags))    = {'megmag'};    % magnetometer
  grad.chantype(kind==1 & ismember(type, meggrads))   = {'meggrad'};   % axial gradiometer
  grad.chantype(kind==1 & ismember(type, megaxials))  = {'megaxial'};  % axial gradiometer, historical exception, FIXME?
  grad.chantype(kind==301 & ismember(type, refmags))  = {'refmag'};    % reference magnetometer
  grad.chantype(kind==301 & ismember(type, refgrads)) = {'refgrad'};   % reference gradiometer
  
  grad.chantype(kind==2)   = {'eeg'}; %EEG channels
  grad.chantype(kind==201) = {'mcg'}; %MCG channels
  grad.chantype(kind==202) = {'eog'}; %EOG
  grad.chantype(kind==302) = {'emg'}; %EMG
  grad.chantype(kind==402) = {'ecg'}; %ECG
  grad.chantype(kind==502) = {'misc'}; %MISC
  grad.chantype(kind==602) = {'respiration'}; %Respiration  
  
  sel = find(kind==3); %Stim channels
    if any(ismember(logno(sel), [101 102])) % newer systems: 101 (and 102, if enabled) are digital; low numbers are 'pseudo-analog' (if enabled)
      grad.chantype(sel(logno(sel)==101)) = {'digital trigger'};
      grad.chantype(sel(logno(sel)==102)) = {'digital trigger'};
      grad.chantype(sel(logno(sel)<=32))  = {'analog trigger'};
      others = logno(sel)>32 & logno(sel)~=101 & logno(sel)~=102;
      grad.chantype(sel(others)) = {'other trigger'};
    elseif any(ismember(logno(sel), [14 15 16])) % older systems: STI 014/015/016 are digital; lower numbers 'pseudo-analog'(if enabled)
      grad.chantype(sel(logno(sel)==14)) = {'digital trigger'};
      grad.chantype(sel(logno(sel)==15)) = {'digital trigger'};
      grad.chantype(sel(logno(sel)==16)) = {'digital trigger'};
      grad.chantype(sel(logno(sel)<=13)) = {'analog trigger'};
      grad.chantype(sel(logno(sel)>16))  = {'other trigger'};
    else
      ft_warning('There does not seem to be a suitable trigger channel.');
      grad.chantype(sel) = {'other trigger'};
    end
  
  % FIFF.FIFF_UNIT_HZ  = 101;
  % FIFF.FIFF_UNIT_N   = 102;
  % FIFF.FIFF_UNIT_PA  = 103;
  % FIFF.FIFF_UNIT_J   = 104;
  % FIFF.FIFF_UNIT_W   = 105;
  % FIFF.FIFF_UNIT_C   = 106;
  % FIFF.FIFF_UNIT_V   = 107;
  % FIFF.FIFF_UNIT_F   = 108;
  % FIFF.FIFF_UNIT_OHM = 109;
  % FIFF.FIFF_UNIT_MHO = 110;
  % FIFF.FIFF_UNIT_WB  = 111;
  % FIFF.FIFF_UNIT_T   = 112;
  % FIFF.FIFF_UNIT_H   = 113;
  % FIFF.FIFF_UNIT_CEL = 114;
  % FIFF.FIFF_UNIT_LM  = 115;
  % FIFF.FIFF_UNIT_LX  = 116;
  
  for sel=find([orig.chs.unit]==107)'
    grad.chanunit(sel) = {'V'};
  end
  for sel=find([orig.chs.unit]==112)'
    grad.chanunit(sel) = {'T'};
  end
  for sel=find([orig.chs.unit]==201)'
    grad.chanunit(sel) = {'T/m'};
  end
  
  remove = cellfun(@isempty, grad.label);
  grad.label    = grad.label(~remove);
  grad.chantype = grad.chantype(~remove);
  grad.chanunit = grad.chanunit(~remove);
  grad.tra      = grad.tra(~remove,:);
  grad.chanpos  = grad.chanpos(~remove,:);
  grad.chanori  = grad.chanori(~remove,:);
  
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the original implementation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % how many Planar gradiometers?
  nPlaGrad = 0;
  for i = 1:orig.nchan
    nPlaGrad = nPlaGrad +(orig.chs(i).coil_type==3012 | orig.chs(i).coil_type==3013 | orig.chs(i).coil_type==3014 | orig.chs(i).coil_type==2);
  end
  
  % how many Magnetometers?
  nMag = 0;
  for i = 1:orig.nchan
    nMag = nMag + (orig.chs(i).coil_type==3022 | orig.chs(i).coil_type==3023 | orig.chs(i).coil_type==3024 | orig.chs(i).coil_type==8101);
  end
  
  % how many Axial gradiometers?
  nAxGrad = 0;
  for i = 1:orig.nchan
    nAxGrad = nAxGrad +(orig.chs(i).coil_type==7001); % babySQUID
  end
  
  % % how many IAS (internal active shielding) channels?
  % nIAS = 0;
  % for i = 1:orig.nchan;
  %   nIAS = nIAS + ~isempty(strmatch('IAS', orig.chs(i).ch_name));
  % end
  
  % how many magnetic field sensors in total?
  nSensors = nPlaGrad + nMag + nAxGrad;
  
  % how many coils in total?
  nCoils = nPlaGrad*2 + nMag + nAxGrad*2;
  
  % intialise grad structure
  grad.coilpos  = zeros(nCoils,3);
  grad.coilori  = zeros(nCoils,3);
  grad.tra      = zeros(nSensors,nCoils);
  grad.unit     = 'cm'; % see below for the conversion, the original fif units are in meter
  grad.label    = cell(nSensors,1);
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'neuromag';
  end
  
  % initialise elec structure, this can remain empty
  elec = [];
  
  % define coils
  kCoil = 1;
  kChan = 1;
  % cf. Joachim's original script - I've implemented it this way in case MEG
  % and EEG channels are not listed first in the .fif file; this shouldn't
  % ever be the case but acts as a safety net...
  
  for n = 1:orig.nchan
    if (orig.chs(n).coil_type==3022 || orig.chs(n).coil_type==3023 || orig.chs(n).coil_type==3024 || orig.chs(n).coil_type==8101) % magnetometer
      t = orig.chs(n).coil_trans;
      
      % TC 2011 09 24 I have changed the coil definition, the original was
      % grad.coilpos(kCoil,:) = t(1:3,4);
      
      grad.coilpos(kCoil,:) = t(1:3,4)+0.0003*t(1:3,3);
      grad.coilori(kCoil,:) = t(1:3,3);
      grad.tra(kChan,kCoil) = 1;
      kCoil = kCoil+1;
      grad.label{kChan} = deblank(orig.ch_names{n});
      grad.chantype{kChan,1}='megmag';
      kChan = kChan+1;
      
    elseif (orig.chs(n).coil_type==3012 || orig.chs(n).coil_type==3013 || orig.chs(n).coil_type==3014 || orig.chs(n).coil_type==2) % planar gradiometer
      t = orig.chs(n).coil_trans;
      
      % TC 2011 09 24 I have changed the coil definition, the original was
      % grad.coilpos(kCoil,:) = 0.0000*t(1:3,3)+t(1:3,4)-0.0084*t(1:3,1)); % for the 1st coil
      % grad.coilpos(kCoil,:) = 0.0000*t(1:3,3)+t(1:3,4)+0.0084*t(1:3,1)); % for the 2nd coil
      
      grad.coilpos(kCoil,:) = 0.0003*t(1:3,3)+t(1:3,4)-0.0084*t(1:3,1); % for the 1st coil
      grad.coilori(kCoil,:) = t(1:3,3);
      grad.tra(kChan,kCoil) = -1;
      kCoil = kCoil+1;
      
      grad.coilpos(kCoil,:) = 0.0003*t(1:3,3)+t(1:3,4)+0.0084*t(1:3,1); % for the 2nd coil
      grad.coilori(kCoil,:) = t(1:3,3);
      grad.tra(kChan,kCoil) = 1;
      kCoil = kCoil+1;
      
      grad.label{kChan} = deblank(orig.ch_names{n});
      grad.chantype{kChan,1}='megplanar';
      kChan = kChan+1;
      
    elseif (orig.chs(n).coil_type==7001)  % babySQUID axial gradiometer bottom coils
      t = orig.chs(n).coil_trans;
      
      grad.coilpos(kCoil,:)=t(1:3,4);  % for the 1st coil
      grad.coilori(kCoil,:)=t(1:3,3);
      grad.tra(kChan,kCoil)=1;
      kCoil=kCoil+1;
      
      grad.coilpos(kCoil,:)=t(1:3,4)+0.050*t(1:3,3);  % for the 2nd coil
      grad.coilori(kCoil,:)=t(1:3,3);
      grad.tra(kChan,kCoil)=-1;
      kCoil=kCoil+1;
      
      grad.label{kChan}=deblank(orig.ch_names{n});
      grad.chantype{kChan,1}='megaxial';
      kChan=kChan+1;
      
    else
      % do nothing - either an EEG channel or something else such as a stim channel
    end
  end
  
  % check we've got all the MEG channels:
  kChan = kChan-1;
  if kChan ~= (nPlaGrad + nMag +nAxGrad)
    ft_error('Number of MEG channels identified does not match number of channels in grad structure');
  end
  
  % determine the type of acquisition system
  if nAxGrad>0
    grad.type = 'babysquid74';
  elseif nPlaGrad>122 && nMag~=0
    grad.type = 'neuromag306';
  elseif nPlaGrad<=122 && nMag==0
    grad.type = 'neuromag122';
  else
    % do not specify type of acquisition system
  end
  
  % multiply by 100 to get cm
  grad.coilpos = 100*grad.coilpos;
end



% how many EEG channels?
nEEG = 0;
for i = 1:orig.nchan
  nEEG = nEEG +(orig.chs(i).kind==2);
end

% define EEG channels
if nEEG>0
  elec.elecpos  = zeros(nEEG,3);
  elec.unit     = 'cm';
  elec.label    = cell(nEEG,1);
  
  % Amendments to overcome problem with fiff_read_meas_info.m when >60 channels (thanks to Rik Henson)
  
  dig_eeg = find([orig.dig.kind]==3); % Find EEG digitisations
  
  if ~isempty(dig_eeg)
    dig_eeg(1) = []; % Remove reference
  end
  
  chn_eeg = find([orig.chs.kind]==2); % Find EEG channels
  
  kChan = 0;
  for n = chn_eeg
    kChan = kChan+1;
    if nEEG<=60
      elec.elecpos(kChan,1:3) = orig.chs(n).eeg_loc(1:3);
    else
      if kChan<=numel(dig_eeg)
        elec.elecpos(kChan,1:3) = orig.dig(dig_eeg(kChan)).r;
      else
        ft_warning('not all EEG channel positions have been digitized');
        elec.elecpos(kChan,1:3) = nan;
      end
    end
    elec.label{kChan} = deblank(orig.ch_names{n});
  end
  
  % check we've got all the EEG channels:
  if kChan~=(nEEG)
    ft_error('Number of EEG channels identified does not match number of channels in elec structure!!!!!');
  end
  
  % multiply by 100 to get cm
  elec.elecpos = 100*elec.elecpos;
end

% construct SSP projectors, these can be applied later to the data and grad
% structure using FT_DENOISE_SSP
grad.balance = [];
grad.balance.current = 'none';

for i = 1:length(hdr.projs)
  hdr.projs(i).active = 'true';
  montage = [];
  montage.tra = mne_make_projector(hdr.projs(i), hdr.ch_names, []);
  montage.labelorg = hdr.ch_names;
  montage.labelnew = hdr.ch_names;
  grad.balance.(fixname(hdr.projs(i).desc)) = montage;
end

% remove grad if completely empty
if size(grad.label,1) == 0
  grad = [];
end
