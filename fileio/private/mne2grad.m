function [grad,elec] = mne2grad(hdr)

% MNE2GRAD creates gradiometer definition for FIFF dataset
%
% Use as
%   [grad,elec] = mne2grad(hdr)
%
% See also CTF2GRAD, BTI2GRAD

% Laurence Hunt 03/12/2008 (with thanks to Joachim Gross's original script based on fiff_access). lhunt@fmrib.ox.ac.uk
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

grad = [];

% orig = fiff_read_meas_info(filename);
% orig = fiff_setup_read_raw(filename);

if isfield(hdr, 'orig')
  orig = hdr.orig;  % use the original MNE header, not the FieldTrip header
else
  orig = hdr;       % assume that it is the original header
end

%begin by transforming all channel locations into device frame,
%if possible
if isempty(orig.dev_head_t)==0
  orig.chs=fiff_transform_meg_chs(orig.chs,orig.dev_head_t);
  orig.chs=fiff_transform_eeg_chs(orig.chs,orig.dev_head_t); %EEG channels are normally stored in head coordinates anyway, but what the heck
else
  warning(['No device to head transform available in fif file']);
  warning(['MEG channels will likely have coordinates in device frame, not head frame']);
end

%how many Planar gradiometers?
nPlaGrad = 0;
for i = 1:orig.nchan;
  nPlaGrad = nPlaGrad +(orig.chs(i).coil_type==3012|orig.chs(i).coil_type==3013|orig.chs(i).coil_type==3014) ;
end

%how many Magnetometers?
nMag = 0;
for i = 1:orig.nchan;
  nMag = nMag +(orig.chs(i).coil_type==3022|orig.chs(i).coil_type==3023|orig.chs(i).coil_type==3024);
end

%how many EEG channels?
nEEG = 0;
for i = 1:orig.nchan;
  nEEG = nEEG +(orig.chs(i).kind==2);
end

%how many sensors in total?
nSensors = nPlaGrad + nMag;

%how many coils in total?
nCoils = nPlaGrad*2 + nMag;

%intialise grad structure
grad.pnt=zeros(nCoils,3);
grad.ori=zeros(nCoils,3);
grad.tra=zeros(nSensors,nCoils);
grad.unit='cm';
grad.label = cell(nSensors,1);

%initialise elec structure
% ??? %

%define coils
kCoil=1;
k=1; % cf. Joachim's original script - I've implemented it this way
% in case MEG and EEG channels are not listed first in the
% .fif file; this shouldn't ever be the case but acts as a
% safety net...
for n=1:orig.nchan
  if (orig.chs(n).coil_type==3022|orig.chs(n).coil_type==3023|orig.chs(n).coil_type==3024) %magnetometer
    t=orig.chs(n).coil_trans;
    grad.pnt(kCoil,:)=100*(t(1:3,4)); % multiply by 100 to get cm
    grad.ori(kCoil,:)=t(1:3,3);
    grad.tra(k,kCoil)=1;
    kCoil=kCoil+1;
    grad.label{k}=deblank(orig.ch_names{n});
    k=k+1;
  elseif (orig.chs(k).coil_type==3012|orig.chs(k).coil_type==3013|orig.chs(k).coil_type==3014) %planar gradiometer
    t=orig.chs(n).coil_trans;
    grad.pnt(kCoil,:)=100*(t(1:3,4)-0.008*t(1:3,1)); % multiply with 100 to get cm
    grad.ori(kCoil,:)=t(1:3,3);
    grad.tra(k,kCoil)= -1;
    kCoil=kCoil+1;
    grad.pnt(kCoil,:)=100*(t(1:3,4)+0.008*t(1:3,1));
    grad.ori(kCoil,:)=t(1:3,3);
    grad.tra(k,kCoil)=1;
    kCoil=kCoil+1;
    grad.label{k}=deblank(orig.ch_names{n});
    k=k+1;
  else
    %do nothing - either an EEG channel or something else (stim
    %channel etc.
  end
end

%check we've got all the MEG channels:
k=k-1;
if k ~= (nPlaGrad + nMag)
  error('Number of MEG channels identified does not match number of channels in grad structure!!!!!');
end

%%define EEG channels
elec = [];

if nEEG>0
    elec.pnt = zeros(nEEG,3);
    elec.unit = 'cm';
    elec.label = cell(nEEG,1);

    % Amendments to overcome problem with fiff_read_meas_info.m when >60
    % channels (thanks to Rik Henson)
    
    dig_eeg = find([orig.dig.kind]==3);  % Find EEG digitisations
    
    if ~isempty(dig_eeg)
        dig_eeg(1) = [];                     % Remove reference
    end
    
    chn_eeg = find([orig.chs.kind]==2);  % Find EEG channels

    k=0;
    for n = chn_eeg
        k=k+1;
        if nEEG<=60
            elec.pnt(k,1:3)=100*orig.chs(n).eeg_loc(1:3);  % multiply by 100 to get cm
        else
            elec.pnt(k,1:3)=100*orig.dig(dig_eeg(k)).r;    % multiply by 100 to get cm
        end
        elec.label{k}=deblank(orig.ch_names{n});
    end
    
    %check we've got all the EEG channels:
    if k ~= (nEEG)
        error('Number of EEG channels identified does not match number of channels in elec structure!!!!!');
    end
end

