function [grad,elec] = mne2grad(hdr)

% MNE2GRAD creates gradiometer definition for FIFF dataset
%
% Use as
%   [grad,elec] = mne2grad(hdr)
%
% See also CTF2GRAD, BTI2GRAD

% Laurence Hunt 03/12/2008 (with thanks to Joachim Gross's original script based on fiff_access). lhunt@fmrib.ox.ac.uk
% 
% $Log: mne2grad.m,v $
% Revision 1.8  2009/04/22 14:55:41  vlalit
% Major bug fix by Laurence who originally mixed up magnetometers and planar gradiometers
%  while building the grad.
%
% Revision 1.7  2009/04/20 17:19:01  vlalit
% Changed the MNE reader not to set hdr.elec when there are no EEG channels.
%
% Revision 1.6  2009/02/05 18:30:44  vlalit
% Updates by Laurence to recognize additional Neuromag sensor types
%
% Revision 1.5  2009/01/23 18:56:54  vlalit
% Another bug fix
%
% Revision 1.3  2009/01/23 16:29:40  roboos
% convert label to column
%
% Revision 1.2  2009/01/23 16:16:55  roboos
% changed the input to the function, removed unused fiff_setup_read_raw at begin of function, added cvs log
%

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
    grad.tra(k,kCoil)=1;
    kCoil=kCoil+1;
    grad.pnt(kCoil,:)=100*(t(1:3,4)+0.008*t(1:3,1));
    grad.ori(kCoil,:)=t(1:3,3);
    grad.tra(k,kCoil)=-1;
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

    k=1;
    for n=1:orig.nchan
        if (orig.chs(n).kind==2) %it's an EEG channel
            elec.pnt(k,1:3)=100*orig.chs(n).eeg_loc(1:3);  % multiply by 100 to get cm
            elec.label{k}=deblank(orig.ch_names{n});
            k=k+1;
        else
            %do nothing
        end
    end

    %check we've got all the EEG channels:
    k=k-1;
    if k ~= (nEEG)
        error('Number of EEG channels identified does not match number of channels in elec structure!!!!!');
    end
end

