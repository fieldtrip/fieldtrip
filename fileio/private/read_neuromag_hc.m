function [hc] = read_neuromag_hc(filename)

% READ_NEUROMAG_HC extracts the MEG headcoil marker positions from a neuromag
% fif file or from the FieldTrip buffer
%
% the definition of head coordinates is according to CTF standard:
% - Origin: Intersection of the line through LPA and RPA and a line orthogonal
%   to L passing through the nasion
% - X-axis from the origin towards the RPA point (exactly through)
% - Y-axis from the origin towards the nasion (exactly through)
% - Z-axis from the origin upwards orthogonal to the XY-plane
%
% hc = read_neuromag_hc(filename)
%
% returns a structure with the following fields
%   hc.dewar.nas    marker positions relative to dewar
%   hc.dewar.lpa
%   hc.dewar.rpa
%   hc.head.nas     marker positions relative to head (measured)
%   hc.head.lpa
%   hc.head.rpa
%   hc.standard.nas marker positions relative to head (expected)
%   hc.standard.lpa
%   hc.standard.rpa

% Copyright (C) 2013, Arjen Stolk
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% read neuromag fif file
hdr = ft_read_header(filename, 'checkmaxfilter', false);

% determine number of digitized points
nFid = size(hdr.orig.dig,2);

% extract the positions in head coordinates (default in header)
posN=1;
for i=1:nFid % loop over fiducials
  
  % 0 is unknown
  % 1 is device, i.e. dewar
  % 4 is fiducial system, i.e. head coordinates
  if hdr.orig.dig(i).coord_frame~=4
    warning(['Digitiser point (' num2str(i) ') not stored in head coordinates!']);
  end
  
  switch hdr.orig.dig(i).kind % constants defined in MNE - see p.215 of MNE manual
    case 1 % Cardinal point (Nasion, LPA or RPA)
      % get location of fiducial:
      hc.head.pos(posN,1:3) = hdr.orig.dig(i).r*100; % multiply by 100 to convert to cm
      switch hdr.orig.dig(i).ident
        case 1 % LPA
          hc.head.label{posN} = 'LPA';
        case 2 % nasion
          hc.head.label{posN} = 'Nasion';
        case 3 % RPA
          hc.head.label{posN} = 'RPA';
        otherwise
          error('Unidentified cardinal point in file');
      end
      posN = posN + 1;
      
    case 2 % HPI coil (up to 5)
      hc.head.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.head.label{posN} = strcat('hpi_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    case 3 % EEG electrode location (or ECG)
      hc.head.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.head.label{posN} = strcat('eeg_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    case 4 % Additional head point
      hc.head.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.head.label{posN} = strcat('extra_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    otherwise
      warning('Unidentified digitiser point in file!');
  end
  
end
hc.head.coordsys = 'neuromag';

% extract the positions in dewar/device coordinates
if ~isempty(hdr.orig.dev_head_t)
  % compute the transformation from head to device
  hdr.orig.head_dev_t.trans = inv(hdr.orig.dev_head_t.trans);
  hdr.orig.head_dev_t.from  = hdr.orig.dev_head_t.to;
  hdr.orig.head_dev_t.to    = hdr.orig.dev_head_t.from;
  
  % transform the digitized points from head to dewar coordinate system
  for k = 1:size(hdr.orig.dig,2)
    pos_dewar = hdr.orig.head_dev_t.trans*[hdr.orig.dig(k).r; 1];
    hdr.orig.dig(k).r = pos_dewar(1:3);
    hdr.orig.dig(k).coord_frame = 1;
  end
else
  warning('No device to head transform available in fif file');
  return
end

posN=1;
for i=1:nFid % loop over fiducials
  
  % 0 is unknown
  % 1 is device, i.e. dewar
  % 4 is fiducial system, i.e. head coordinates
  if hdr.orig.dig(i).coord_frame~=1
    warning(['Digitiser point (' num2str(i) ') not stored in head coordinates!']);
  end
  
  switch hdr.orig.dig(i).kind % constants defined in MNE - see p.215 of MNE manual
    case 1 % Cardinal point (Nasion, LPA or RPA)
      % get location of fiducial:
      hc.dewar.pos(posN,1:3) = hdr.orig.dig(i).r*100; % multiply by 100 to convert to cm
      switch hdr.orig.dig(i).ident
        case 1 % LPA
          hc.dewar.label{posN} = 'LPA';
        case 2 % nasion
          hc.dewar.label{posN} = 'Nasion';
        case 3 % RPA
          hc.dewar.label{posN} = 'RPA';
        otherwise
          error('Unidentified cardinal point in file');
      end
      posN = posN + 1;
      
    case 2 % HPI coil (up to 5)
      hc.dewar.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.dewar.label{posN} = strcat('hpi_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    case 3 % EEG electrode location (or ECG)
      hc.dewar.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.dewar.label{posN} = strcat('eeg_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    case 4 % Additional head point
      hc.dewar.pos(posN,1:3) = hdr.orig.dig(i).r*100;
      hc.dewar.label{posN} = strcat('extra_', num2str(hdr.orig.dig(i).ident));
      posN = posN + 1;
      
    otherwise
      warning('Unidentified digitiser point in file!');
  end
  
end
hc.dewar.coordsys = 'dewar';
