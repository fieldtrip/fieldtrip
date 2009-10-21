function grad = yokogawa2grad(hdr)

% YOKOGAWA2GRAD converts the position and weights of all coils that
% compromise a gradiometer system into a structure that can be used
% by FieldTrip.
%
% See also READ_HEADER, CTF2GRAD, BTI2GRAD, FIF2GRAD

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: yokogawa2grad.m,v $
% Revision 1.1  2009/01/14 09:12:16  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2008/12/24 13:49:25  roboos
% added suggested changes by Kaoru Amano, see email 19 Dec 2008
%
% Revision 1.4  2008/12/18 11:53:20  roboos
% changed some comments and some slight cleanups, no functional change
%
% Revision 1.3  2008/12/18 08:31:23  roboos
% incorporated gradiometer, thanks to Kaoru Amano
%
% Revision 1.2  2008/05/15 13:20:36  roboos
% updated documentation
%
% Revision 1.1  2006/08/31 13:32:11  roboos
% moved from fieldtrip to fileio module
%
% Revision 1.1  2005/09/06 08:54:30  roboos
% new implementations for the Yokogawa 160 channel MEG system


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
% 7  orientation from the 1st to 2nd coil for gradiometer (theta in deg)
% 8  orientation of first coil (phi in deg)
% 9  orientation from the 1st to 2nd coil for gradiometer (phi in deg)
% 10 coil size (in m)
% 11 baseline (in m)

handles    = definehandles;
isgrad     = (hdr.channel_info(:,2)==handles.AxialGradioMeter | hdr.channel_info(:,2)==handles.PlannerGradioMeter);
grad.pnt   = hdr.channel_info(isgrad,3:5)*100;    % cm

% Get orientation of the 1st coil
ori_1st   = hdr.channel_info(find(isgrad),[6 8]);
% polar to x,y,z coordinates
ori_1st = ...
  [sin(ori_1st(:,1)/180*pi).*cos(ori_1st(:,2)/180*pi) ...
  sin(ori_1st(:,1)/180*pi).*sin(ori_1st(:,2)/180*pi) ...
  cos(ori_1st(:,1)/180*pi)];
grad.ori = ori_1st;

% Get orientation from the 1st to 2nd coil for gradiometer
ori_1st_to_2nd   = hdr.channel_info(find(isgrad),[7 9]);
% polar to x,y,z coordinates
ori_1st_to_2nd = ...
  [sin(ori_1st_to_2nd(:,1)/180*pi).*cos(ori_1st_to_2nd(:,2)/180*pi) ...
  sin(ori_1st_to_2nd(:,1)/180*pi).*sin(ori_1st_to_2nd(:,2)/180*pi) ...
  cos(ori_1st_to_2nd(:,1)/180*pi)];
% Get baseline
baseline = hdr.channel_info(isgrad,size(hdr.channel_info,2));

% Define the location and orientation of 2nd coil
for i=1:sum(isgrad)
  if hdr.channel_info(i,2) == handles.AxialGradioMeter
    grad.pnt(i+sum(isgrad),:) = [grad.pnt(i,:)+ori_1st(i,:)*baseline(i)*100];
    grad.ori(i+sum(isgrad),:) = -ori_1st(i,:);
  elseif hdr.channel_info(i,2) == handles.PlannerGradioMeter
    grad.pnt(i+sum(isgrad),:) = [grad.pnt(i,:)+ori_1st_to_2nd(i,:)*baseline(i)*100];
    grad.ori(i+sum(isgrad),:) = ori_1st(i,:);
  end
end

% Define the pair of 1st and 2nd coils for each gradiometer
grad.tra = repmat(diag(ones(1,size(grad.pnt,1)/2),0),1,2);

% Make the matrix sparse to speed up the multiplication in the forward
% computation with the coil-leadfield matrix to get the channel leadfield
grad.tra = sparse(grad.tra);

tmp = hdr.channel_info(isgrad,1);
for i=1:size(tmp,1)
  grad.label{i,1} = num2str(tmp(i)+1);
end
grad.unit='cm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this defines some usefull constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = definehandles;
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
handles.DefaultMagnetometerSize       = (4.0/1000.0);		% ????4.0mm???????`
handles.DefaultAxialGradioMeterSize   = (15.5/1000.0);		% ???a15.5mm???~??
handles.DefaultPlannerGradioMeterSize = (12.0/1000.0);		% ????12.0mm???????`
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

