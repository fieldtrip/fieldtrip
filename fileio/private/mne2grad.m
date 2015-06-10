function [grad, elec] = mne2grad(hdr, dewar)

% MNE2GRAD creates gradiometer definition for FIFF dataset
%
% Use as
%   [grad, elec] = mne2grad(hdr, dewar)
%
% The optional second argument is a boolean that can be used to return the sensors
% in dewar coordinates (default is head coordinates)
%
% See also CTF2GRAD, BTI2GRAD

% Laurence Hunt 03/12/2008 (with thanks to Joachim Gross's original script 
% based on fiff_access). lhunt@fmrib.ox.ac.uk
%
% Teresa Cheung 09/24/2011 revision to Laurence Hunt's script. The coil
% geometry has been revised to better reflect the coil dimensions of the
% Vectorview (306) planar gradiometers and magnetometers. Coil geometry for
% the Neuromag-122 planar gradiometer is different and is not read with
% this script (coil_type should equal 2)

% Copyrights (C) 2011, Teresa Cheung
% Copyrights (C) 2008, Laurence Hunt
% Copyrights (C) 2005, Joachim Gross
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

grad = [];

% orig = fiff_read_meas_info(filename);
% orig = fiff_setup_read_raw(filename);

if isfield(hdr, 'orig')
  orig = hdr.orig; % use the original MNE header, not the FieldTrip header
else
  orig = hdr; % assume that it is the original header
end

% begin by transforming all channel locations into the desired coordinate system, if possible

if ~dewar
  if ~isempty(orig.dev_head_t)
    orig.chs = fiff_transform_meg_chs(orig.chs,orig.dev_head_t);
    orig.chs = fiff_transform_eeg_chs(orig.chs,orig.dev_head_t); % EEG channels are normally stored in head coordinates anyway, but what the heck
  else
    warning('No device to head transform available in fif file');
    warning('MEG channels will likely have coordinates in device frame, not head frame');
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
    warning('No device to head transform available in fif file');
    warning('EEG channels will likely have coordinates in head frame, not device frame');
  end
end

% how many Planar gradiometers?
nPlaGrad = 0;
for i = 1:orig.nchan;
  nPlaGrad = nPlaGrad +(orig.chs(i).coil_type==3012 | orig.chs(i).coil_type==3013 | orig.chs(i).coil_type==3014 | orig.chs(i).coil_type==2) ;
end

% how many Magnetometers?
nMag = 0;
for i = 1:orig.nchan;
  nMag = nMag +(orig.chs(i).coil_type==3022 | orig.chs(i).coil_type==3023 | orig.chs(i).coil_type==3024);
end

% how many Axial gradiometers?
nAxGrad = 0;
for i = 1:orig.nchan;
  nAxGrad = nAxGrad +(orig.chs(i).coil_type==7001); % babySQUID
end

% how many EEG channels?
nEEG = 0;
for i = 1:orig.nchan;
  nEEG = nEEG +(orig.chs(i).kind==2);
end

% % how many IAS (internal active shielding) channels?
% nIAS = 0;
% for i = 1:orig.nchan;
%   nIAS = nIAS + ~isempty(strmatch('IAS', orig.chs(i).ch_name));
% end

% how many sensors in total?
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

% multiply by 100 to get cm
grad.coilpos = 100*grad.coilpos;

% check we've got all the MEG channels:
kChan = kChan-1;
if kChan ~= (nPlaGrad + nMag +nAxGrad)
  error('Number of MEG channels identified does not match number of channels in grad structure');
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
    error('Number of EEG channels identified does not match number of channels in elec structure!!!!!');
  end
  
  % multiply by 100 to get cm
  elec.elecpos = 100*elec.elecpos;
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

% remove grad if empty
if size(grad.label,1) == 0
  grad = [];
end
