function [grad, elec] = mne2grad(hdr, dewar, coilaccuracy)

% MNE2GRAD converts a header from a fif file that was read using the MNE toolbox into
% a gradiometer structure that can be understood by the FieldTrip low-level forward
% and inverse routines.
%
% Use as
%   [grad, elec] = mne2grad(hdr, dewar, coilaccuracy)
% where
%   dewar        = boolean, whether to return it in dewar or head coordinates (default = false, i.e. head coordinates)
%   coilaccuracy = empty or a number (default = [])
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
    orig.chs = fiff_transform_meg_chs(orig.chs,orig.dev_head_t);
    orig.chs = fiff_transform_eeg_chs(orig.chs,orig.dev_head_t); % EEG channels are normally stored in head coordinates anyway, but what the heck
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
  [ftver, ftpath] = ft_version;
  def = mne_load_coil_def(fullfile(ftpath, 'external', 'mne', 'coil_def.dat'));
  
  k = 1;
  
  grad.chanpos = nan(length(orig.chs),3);
  grad.chanori = nan(length(orig.chs),3);
  grad.label = cell(length(orig.chs),1);
  grad.tra      = nan(length(orig.chs),0);
  grad.chantype = repmat({'unknown'}, size(grad.label));
  grad.chanunit = repmat({'unknown'}, size(grad.label));
  
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
  
  % all the chs.kinds and chs.coil_types are obtained from the MNE manual, p.210-211
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==2)' % planar gradiometers
    grad.chantype(sel) = {'megplanar'}; %Neuromag-122 planar gradiometer
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3012)' %planar gradiometers
    grad.chantype(sel) = {'megplanar'}; %Type T1 planar grad
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3013)' %planar gradiometers
    grad.chantype(sel) = {'megplanar'}; %Type T2 planar grad
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3014)' %planar gradiometers
    grad.chantype(sel) = {'megplanar'}; %Type T3 planar grad
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3022)' %magnetometers
    grad.chantype(sel) = {'megmag'};    %Type T1 magenetometer
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3023)' %magnetometers
    grad.chantype(sel) = {'megmag'};    %Type T2 magenetometer
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==3024)' %magnetometers
    grad.chantype(sel) = {'megmag'};    %Type T3 magenetometer
  end
  for sel=find([orig.chs.kind]==1 & [orig.chs.coil_type]==7001)' %axial gradiometer
    grad.chantype(sel) = {'megaxial'};
  end
  for sel=find([orig.chs.kind]==301)' %MEG reference channel, located far from head
    grad.chantype(sel) = {'ref'};
  end
  for sel=find([orig.chs.kind]==2)'   %EEG channels
    grad.chantype(sel) = {'eeg'};
  end
  for sel=find([orig.chs.kind]==201)' %MCG channels
    grad.chantype(sel) = {'mcg'};
  end
  for sel=find([orig.chs.kind]==3)' %Stim channels
    if any([orig.chs(sel).logno] == 101) %new systems: 101 (and 102, if enabled) are digital; low numbers are 'pseudo-analog' (if enabled)
      grad.chantype(sel([orig.chs(sel).logno] == 101)) = {'digital trigger'};
      grad.chantype(sel([orig.chs(sel).logno] == 102)) = {'digital trigger'};
      grad.chantype(sel([orig.chs(sel).logno] <= 32))  = {'analog trigger'};
      others = [orig.chs(sel).logno] > 32 & [orig.chs(sel).logno] ~= 101 & ...
        [orig.chs(sel).logno] ~= 102;
      grad.chantype(sel(others)) = {'other trigger'};
    elseif any([orig.chs(sel).logno] == 14) %older systems: STI 014/015/016 are digital; lower numbers 'pseudo-analog'(if enabled)
      grad.chantype(sel([orig.chs(sel).logno] == 14)) = {'digital trigger'};
      grad.chantype(sel([orig.chs(sel).logno] == 15)) = {'digital trigger'};
      grad.chantype(sel([orig.chs(sel).logno] == 16)) = {'digital trigger'};
      grad.chantype(sel([orig.chs(sel).logno] <= 13)) = {'analog trigger'};
      others = [orig.chs(sel).logno] > 16;
      grad.chantype(sel(others)) = {'other trigger'};
    else
      ft_warning('There does not seem to be a suitable trigger channel.');
      grad.chantype(sel) = {'other trigger'};
    end
  end
  for sel=find([orig.chs.kind]==202)' %EOG
    grad.chantype(sel) = {'eog'};
  end
  for sel=find([orig.chs.kind]==302)' %EMG
    grad.chantype(sel) = {'emg'};
  end
  for sel=find([orig.chs.kind]==402)' %ECG
    grad.chantype(sel) = {'ecg'};
  end
  for sel=find([orig.chs.kind]==502)' %MISC
    grad.chantype(sel) = {'misc'};
  end
  for sel=find([orig.chs.kind]==602)' %Resp
    grad.chantype(sel) = {'respiration'};
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
    nMag = nMag +(orig.chs(i).coil_type==3022 | orig.chs(i).coil_type==3023 | orig.chs(i).coil_type==3024);
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
    if (orig.chs(n).coil_type==3022 || orig.chs(n).coil_type==3023 || orig.chs(n).coil_type==3024) % magnetometer
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

% remove grad if completely empty
if size(grad.label,1) == 0
  grad = [];
end
