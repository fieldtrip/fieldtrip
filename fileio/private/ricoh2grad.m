function grad = ricoh2grad(hdr)

%% RICOH2GRAD converts the position and weights of all coils that
%% compromise a gradiometer system into a structure that can be used
%% by FieldTrip. This implementation uses the "ricoh_meg_reader" toolbox.
%%
%% See also FT_READ_HEADER, CTF2GRAD, BTI2GRAD, FIF2GRAD, YOKOGAWA2GRAD_NEW
%

% The following line is only a safety measure:
% No function of the toolbox is actually called in this routine.
ft_hastoolbox('ricoh_meg_reader', 1);

if isfield(hdr, 'label')
  label = hdr.label; % keep for later use
end

if isfield(hdr, 'orig')
  hdr = hdr.orig; % use the original header, not the FieldTrip header
end

% The "channel_info.channel(i)" structure contains the details of coils.
% For the sepcifications, see Ricoh MEG Reader Toolbox 1.0 specifications.pdf.
% type, type of sensor
% data.x (in m, inner coil)
% data.y (in m, inner coil)
% data.z (in m, inner coil)
% data.size (in m, coil diameter)

% for gradiometers
% data.baseline (in m)

% for axial gradiometers and magnetometers
% data.zdir orientation of inner coil (theta in deg: angle to z-axis)
% data.xdir orientation of inner coil (phi in deg: angle to x-axis)

% shorten names
ch_info = hdr.channel_info.channel;
type = [ch_info.type];
handles    = definehandles;

% get all axial grads, planar grads, and magnetometers.
% All reference channels are excluded even if some of them have position information.
grad_ind = [1:hdr.channel_count];
isgrad   = (type==handles.AxialGradioMeter | type==handles.PlannerGradioMeter | type==handles.MagnetoMeter);

grad_ind = grad_ind(isgrad);
grad_nr = size(grad_ind,2);

grad = [];
grad.coilpos   = zeros(2*grad_nr,3);
grad.coilori   = zeros(2*grad_nr,3);

% define gradiometer and magnetometer
for i = 1:grad_nr
  ch_ind = grad_ind(i);
  grad.coilpos(i,1) =  ch_info(ch_ind).data.x*100; % cm
  grad.coilpos(i,2) =  ch_info(ch_ind).data.y*100; % cm
  grad.coilpos(i,3) =  ch_info(ch_ind).data.z*100; % cm
  grad.chanpos = grad.coilpos(1:grad_nr,:);

  if ch_info(ch_ind).type==handles.AxialGradioMeter || ch_info(ch_ind).type==handles.RefferenceAxialGradioMeter
      baseline = ch_info(ch_ind).data.baseline;

      ori_1st = [ch_info(ch_ind).data.zdir ch_info(ch_ind).data.xdir ];
      % polar to x,y,z coordinates
      ori_1st = ...
      [sin(ori_1st(:,1)/180*pi).*cos(ori_1st(:,2)/180*pi) ...
      sin(ori_1st(:,1)/180*pi).*sin(ori_1st(:,2)/180*pi) ...
      cos(ori_1st(:,1)/180*pi)];

      grad.coilori(i,:) = ori_1st;

      grad.coilpos(i+grad_nr,:) = [grad.coilpos(i,:)+ori_1st*baseline*100];
      grad.coilori(i+grad_nr,:) = -ori_1st;
  elseif ch_info(ch_ind).type==handles.PlannerGradioMeter || ch_info(ch_ind).type==handles.RefferencePlannerGradioMeter
      baseline = ch_info(ch_ind).data.baseline;

      ori_1st = [ch_info(ch_ind).data.zdir1 ch_info(ch_ind).data.xdir1 ];
      % polar to x,y,z coordinates
      ori_1st = ...
      [sin(ori_1st(:,1)/180*pi).*cos(ori_1st(:,2)/180*pi) ...
      sin(ori_1st(:,1)/180*pi).*sin(ori_1st(:,2)/180*pi) ...
      cos(ori_1st(:,1)/180*pi)];

      grad.coilori(i,:) = ori_1st;

      ori_1st_to_2nd   = [ch_info(ch_ind).data.zdir2 ch_info(ch_ind).data.xdir2 ];
      % polar to x,y,z coordinates
      ori_1st_to_2nd = ...
      [sin(ori_1st_to_2nd(:,1)/180*pi).*cos(ori_1st_to_2nd(:,2)/180*pi) ...
      sin(ori_1st_to_2nd(:,1)/180*pi).*sin(ori_1st_to_2nd(:,2)/180*pi) ...
      cos(ori_1st_to_2nd(:,1)/180*pi)];

      grad.coilpos(i+grad_nr,:) = [grad.coilpos(i,:)+ori_1st_to_2nd*baseline*100];
      grad.coilori(i+grad_nr,:) = -ori_1st;
  else % magnetometer
      ori_1st = [ch_info(ch_ind).data.zdir ch_info(ch_ind).data.xdir ];
      % polar to x,y,z coordinates
      ori_1st = ...
      [sin(ori_1st(:,1)/180*pi).*cos(ori_1st(:,2)/180*pi) ...
      sin(ori_1st(:,1)/180*pi).*sin(ori_1st(:,2)/180*pi) ...
      cos(ori_1st(:,1)/180*pi)];

      grad.coilori(i,:) = ori_1st;

      grad.coilpos(i+grad_nr,:) = [0 0 0];
      grad.coilori(i+grad_nr,:) = [0 0 0];
  end

  grad.chanori = grad.coilori(1:grad_nr,:);

end

% Define the pair of 1st and 2nd coils for each gradiometer
grad.tra = repmat(diag(ones(1,grad_nr),0),1,2);

% for mangetometers change tra as there is no second coil
for i = 1:grad_nr
  ch_ind = grad_ind(i);
  if ch_info(ch_ind).type==handles.MagnetoMeter
    grad.tra(i,grad_nr+i) = 0;
  end
end

% the gradiometer labels should be consistent with the channel labels in
% read_ricoh_header, the predefined list of channel names in ft_senslabel
% and with ft_channelselection:
% but it is ONLY consistent with read_ricoh_header as NO FIXED relation
% between channel index and type of channel exists for Ricoh systems.
% Therefore all have individual label sequences: Support in ft_senslabel is only partial.
if ~isempty(label)
    grad.label = label(grad_ind)';
else
    % this is only backup, if something goes wrong above.
    label = cell(grad_nr,1);
    for i=1:length(label)
    label{i,1} = sprintf('AG%03d', i);
    end
    grad.label = label;
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
handles.DefaultMagnetometerSize       = (4.0/1000.0);       % Square of 4.0mm in length
handles.DefaultAxialGradioMeterSize   = (15.5/1000.0);      % Circle of 15.5mm in diameter
handles.DefaultPlannerGradioMeterSize = (12.0/1000.0);      % Square of 12.0mm in length
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
