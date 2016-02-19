function grad = yokogawa2grad(hdr)

% YOKOGAWA2GRAD converts the position and weights of all coils that
% compromise a gradiometer system into a structure that can be used
% by FieldTrip. This implementation uses the old "yokogawa" toolbox.
%
% See also CTF2GRAD, BTI2GRAD, FIF2GRAD, MNE2GRAD, ITAB2GRAD,
% FT_READ_SENS, FT_READ_HEADER

% Copyright (C) 2005-2008, Robert Oostenveld
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

if ~ft_hastoolbox('yokogawa')
    error('cannot determine whether Yokogawa toolbox is present');
end

if isfield(hdr, 'label')
  label = hdr.label; % keep for later use
end

if isfield(hdr, 'orig')
  hdr = hdr.orig; % use the original header, not the FieldTrip header
end

% The "channel_info" contains
% 1  channel number, zero offset
% 2  channel type, type of gradiometer
% 3  position x (in m)
% 4  position y (in m)
% 5  position z (in m)
% 6  orientation of first coil (theta in deg)
% 7  orientation of first coil (phi in deg)
% 8  orientation from the 1st to 2nd coil for gradiometer (theta in deg)
% 9  orientation from the 1st to 2nd coil for gradiometer (phi in deg)
% 10 coil size (in m)
% 11 baseline (in m)

handles    = definehandles;
isgrad     = (hdr.channel_info(:,2)==handles.AxialGradioMeter | ...
              hdr.channel_info(:,2)==handles.PlannerGradioMeter | ...
              hdr.channel_info(:,2)==handles.MagnetoMeter | ...
              hdr.channel_info(:,2)==handles.RefferenceAxialGradioMeter);
% reference channels are excluded because the positions are not specified
%              hdr.channel_info(:,2)==handles.RefferencePlannerGradioMeter
%              hdr.channel_info(:,2)==handles.RefferenceMagnetoMeter
isgrad_handles = hdr.channel_info(isgrad,2);
ismag      = (isgrad_handles(:)==handles.MagnetoMeter | isgrad_handles(:)==handles.RefferenceMagnetoMeter);
grad.coilpos   = hdr.channel_info(isgrad,3:5)*100;    % cm
grad.unit  = 'cm';

% Get orientation of the 1st coil
ori_1st   = hdr.channel_info(find(isgrad),[6 7]);
% polar to x,y,z coordinates
ori_1st = ...
  [sin(ori_1st(:,1)/180*pi).*cos(ori_1st(:,2)/180*pi) ...
  sin(ori_1st(:,1)/180*pi).*sin(ori_1st(:,2)/180*pi) ...
  cos(ori_1st(:,1)/180*pi)];
grad.coilori = ori_1st;

% Get orientation from the 1st to 2nd coil for gradiometer
ori_1st_to_2nd   = hdr.channel_info(find(isgrad),[8 9]);
% polar to x,y,z coordinates
ori_1st_to_2nd = ...
  [sin(ori_1st_to_2nd(:,1)/180*pi).*cos(ori_1st_to_2nd(:,2)/180*pi) ...
  sin(ori_1st_to_2nd(:,1)/180*pi).*sin(ori_1st_to_2nd(:,2)/180*pi) ...
  cos(ori_1st_to_2nd(:,1)/180*pi)];
% Get baseline
baseline = hdr.channel_info(isgrad,size(hdr.channel_info,2));

% Define the location and orientation of 2nd coil
info = hdr.channel_info(isgrad,2); 
for i=1:sum(isgrad)
  if (info(i) == handles.AxialGradioMeter || info(i) == handles.RefferenceAxialGradioMeter )
    grad.coilpos(i+sum(isgrad),:) = [grad.coilpos(i,:)+ori_1st(i,:)*baseline(i)*100];
    grad.coilori(i+sum(isgrad),:) = -ori_1st(i,:);
  elseif (info(i) == handles.PlannerGradioMeter || info(i) == handles.RefferencePlannerGradioMeter)
    grad.coilpos(i+sum(isgrad),:) = [grad.coilpos(i,:)+ori_1st_to_2nd(i,:)*baseline(i)*100];
    grad.coilori(i+sum(isgrad),:) = -ori_1st(i,:);
  else
    grad.coilpos(i+sum(isgrad),:) = [0 0 0];
    grad.coilori(i+sum(isgrad),:) = [0 0 0];  
  end
end

% Define the pair of 1st and 2nd coils for each gradiometer
grad.tra = repmat(diag(ones(1,size(grad.coilpos,1)/2),0),1,2);

% for mangetometers change tra as there is no second coil
if any(ismag)
    sz_pnt = size(grad.coilpos,1)/2;
    % create logical variable
    not_2nd_coil = ([diag(zeros(sz_pnt),0)' ismag']~=0);
    grad.tra(ismag,not_2nd_coil) = 0;
end

% the gradiometer labels should be consistent with the channel labels in
% read_yokogawa_header, the predefined list of channel names in ft_senslabel
% and with ft_channelselection
% ONLY consistent with read_yokogawa_header as NO FIXED relation between
% channel index and type of channel exists for Yokogawa systems. Therefore
% all have individual label sequences: No useful support in ft_senslabel possible
if ~isempty(label)
    grad.label = label(isgrad);
else
    % this is only backup, if something goes wrong above.
    label = cell(size(isgrad));
    for i=1:length(label)
    label{i} = sprintf('AG%03d', i);
    end
    grad.label = label(isgrad);    
end
grad.unit = 'cm';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this defines some usefull constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = definehandles
handles.output = [];
handles.sqd_load_flag = false;
handles.mri_load_flag = false;
handles.NullChannel         = 0;
handles.MagnetoMeter        = 1;
handles.AxialGradioMeter    = 2;
handles.PlannerGradioMeter  = 3;
handles.RefferenceChannelMark = hex2dec('0100');
handles.RefferenceMagnetoMeter       = bitor( handles.RefferenceChannelMark, handles.MagnetoMeter );
handles.RefferenceAxialGradioMeter   = bitor( handles.RefferenceChannelMark, handles.AxialGradioMeter );
handles.RefferencePlannerGradioMeter = bitor( handles.RefferenceChannelMark, handles.PlannerGradioMeter );
handles.TriggerChannel      = -1;
handles.EegChannel          = -2;
handles.EcgChannel          = -3;
handles.EtcChannel          = -4;
handles.NonMegChannelNameLength = 32;
handles.DefaultMagnetometerSize       = (4.0/1000.0);       % ????4.0mm???????`
handles.DefaultAxialGradioMeterSize   = (15.5/1000.0);      % ???a15.5mm???~??
handles.DefaultPlannerGradioMeterSize = (12.0/1000.0);      % ????12.0mm???????`
handles.AcqTypeContinuousRaw = 1;
handles.AcqTypeEvokedAve     = 2;
handles.AcqTypeEvokedRaw     = 3;
handles.sqd = [];
handles.sqd.selected_start  = [];
handles.sqd.selected_end    = [];
handles.sqd.axialgradiometer_ch_no      = [];
handles.sqd.axialgradiometer_ch_info    = [];
handles.sqd.axialgradiometer_data       = [];
handles.sqd.plannergradiometer_ch_no    = [];
handles.sqd.plannergradiometer_ch_info  = [];
handles.sqd.plannergradiometer_data     = [];
handles.sqd.nullchannel_ch_no   = [];
handles.sqd.nullchannel_data    = [];
handles.sqd.selected_time       = [];
handles.sqd.sample_rate         = [];
handles.sqd.sample_count        = [];
handles.sqd.pretrigger_length   = [];
handles.sqd.matching_info   = [];
handles.sqd.source_info     = [];
handles.sqd.mri_info        = [];
handles.mri                 = [];

