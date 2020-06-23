function varargout = xsens_mvnx(filename, hdr, begsample, endsample, chanindx)

% XSENS_MVNX reads motion tracking data from a file that was created by
% xsens MVN motion systems, see: https://www.xsens.com/motion-capture. The
% current function is designed to read in .mvnx files from release version
% 4.
%
% Use as
%   hdr = xsens_mvnx(filename);
%   dat = xsens_mvnx(filename, hdr, begsample, endsample, chanindx);
%   evt = xsens_mvnx(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV

% Copyright (C) 2020 Helena Cockx
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

persistent mvnx previous_fullname

% ensure that the external toolbox is available
ft_hastoolbox('xsens', 1); 

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);
if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(mvnx)
  % remember the full filename including path
  previous_fullname = fullname;
  % read the header and data
  mvnx = load_mvnx(fullname);
else
  % use the persistent variable to speed up subsequent read operations
end

if mvnx.version~=4
      ft_warning('this has only been tested with .mvnx files in the "version 4" format');
end

if needhdr
  %% parse the header  
  hdr = [];
  hdr.label = {}; 
  hdr.chanunit={};
  
  % loop over all labels for the segments
  if isfield(mvnx.subject.frames.frame,'orientation')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_orientation_Q0'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_orientation_Q1'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_orientation_Q2'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_orientation_Q3'];
          hdr.chanunit(end+1:end+4)=repmat({'quaternion'}, [1 4]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'position')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_position_X'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_position_Y'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_position_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'m'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'velocity')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_velocity_X'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_velocity_Y'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_velocity_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'m/s'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'acceleration')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_acceleration_X'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_acceleration_Y'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_acceleration_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'m/s^2'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'angularVelocity')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularVelocity_X'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularVelocity_Y'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularVelocity_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'rad/s'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'angularAcceleration')
      for seg=1:mvnx.subject.frames.segmentCount % loop over all segments
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularAcceleration_X'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularAcceleration_Y'];
          hdr.label{end+1}=['seg_' mvnx.subject.segments.segment(seg).label '_angularAcceleration_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'rad/s^2'}, [1 3]);
      end
  end
  
  % loop over the footContacts
  if isfield(mvnx.subject.frames.frame,'footContacts')
      for fc=1:numel(mvnx.subject.footContactDefinition.contactDefinition) % loop over foot contacts
          hdr.label{end+1}=['fc_' mvnx.subject.footContactDefinition.contactDefinition(fc).label '_footContacts'];
          hdr.chanunit(end+1)={'boolean'};
      end
  end
  
  % loop over the labels for the sensors
  if isfield(mvnx.subject.frames.frame,'sensorFreeAcceleration')
      for sen=1:mvnx.subject.frames.sensorCount % loop over all sensors
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorFreeAcceleration_X'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorFreeAcceleration_Y'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorFreeAcceleration_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'m/s^2'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'sensorMagneticField')
      for sen=1:mvnx.subject.frames.sensorCount % loop over all sensors
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorMagneticField_X'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorMagneticField_Y'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorMagneticField_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'a.u.'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'sensorOrientation')
      for sen=1:mvnx.subject.frames.sensorCount % loop over all sensors
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorOrientation_Q0'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorOrientation_Q1'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorOrientation_Q2'];
          hdr.label{end+1}=['sen_' mvnx.subject.sensors.sensor(sen).label '_sensorOrientation_Q3'];
          hdr.chanunit(end+1:end+4)=repmat({'quaternion'}, [1 4]);
      end
  end  

  % loop over the labels of the joints
  if isfield(mvnx.subject.frames.frame,'jointAngle')
      for jnt=1:mvnx.subject.frames.jointCount % loop over all joints
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngle_X'];
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngle_Y'];
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngle_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'deg'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'jointAngleXZY')
      for jnt=1:mvnx.subject.frames.jointCount % loop over all joints
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngleXZY_X'];
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngleXZY_Y'];
          hdr.label{end+1}=['jnt_' mvnx.subject.joints.joint(jnt).label '_jointAngleXZY_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'deg'}, [1 3]);
      end
  end
  
  % loop over the labels of the ergonomic joints
  if isfield(mvnx.subject.frames.frame,'jointAngleErgo')
      for jntx=1:numel(mvnx.subject.ergonomicJointAngles.ergonomicJointAngle) % loop over all ergonomic joints
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgo_X'];
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgo_Y'];
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgo_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'deg'}, [1 3]);
      end
  end
  if isfield(mvnx.subject.frames.frame,'jointAngleErgoXZY')
      for jntx=1:numel(mvnx.subject.ergonomicJointAngles.ergonomicJointAngle) % loop over all ergonomic joints
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgoXZY_X'];
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgoXZY_Y'];
          hdr.label{end+1}=['jntx_' mvnx.subject.ergonomicJointAngles.ergonomicJointAngle(jntx).label '_jointAngleErgoXZY_Z'];
          hdr.chanunit(end+1:end+3)=repmat({'deg'}, [1 3]);
      end
  end
  
  % loop over the labels of Center of mass
  if isfield(mvnx.subject.frames.frame,'centerOfMass')
      hdr.label{end+1}=['seg_COM_centerOfMass_X']; % add center of mass positions
      hdr.label{end+1}=['seg_COM_centerOfMass_Y'];
      hdr.label{end+1}=['seg_COM_centerOfMass_Z'];
      hdr.chanunit(end+1:end+3)=repmat({'m'}, [1 3]);
  end
  
  hdr.nChans      = numel(hdr.label);
  hdr.nSamples    = length(find(strcmp({mvnx.subject.frames.frame(:).type}, 'normal'))); %first three frames include 'identity', 'tpose' and 'tpose-isb' information.
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.Fs          = mvnx.subject.frameRate;
  hdr.chantype    = repmat({'motion'}, size(hdr.label));
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% parse the data  
  nchan   = hdr.nChans;
  nsample = hdr.nSamples;
  idx     = find(strcmp({mvnx.subject.frames.frame(:).type}, 'normal')); % indexes of all 'normal' frames
  dat = nan(nchan, nsample);
  i=1;
  
  if isfield(mvnx.subject.frames.frame,'orientation')
      nchan=length(find(contains(hdr.label, '_orientation')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).orientation], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'position')
      nchan=length(find(contains(hdr.label, '_position')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).position], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'velocity')
      nchan=length(find(contains(hdr.label, '_velocity')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).velocity], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'acceleration')
      nchan=length(find(contains(hdr.label, '_acceleration')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).acceleration], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'angularVelocity')
      nchan=length(find(contains(hdr.label, '_angularVelocity')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).angularVelocity], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'angularAcceleration')
      nchan=length(find(contains(hdr.label, '_angularAcceleration')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).angularAcceleration], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'footContacts')
      nchan=length(find(contains(hdr.label, '_footContacts')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).footContacts], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'sensorFreeAcceleration')
      nchan=length(find(contains(hdr.label, '_sensorFreeAcceleration')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).sensorFreeAcceleration], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'sensorMagneticField')
      nchan=length(find(contains(hdr.label, '_sensorMagneticField')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).sensorMagneticField], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'sensorOrientation')
      nchan=length(find(contains(hdr.label, '_sensorOrientation')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).sensorOrientation], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'jointAngle')
      nchan=length(find(contains(hdr.label, '_jointAngle_')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).jointAngle], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'jointAngleXZY')
      nchan=length(find(contains(hdr.label, '_jointAngleXZY')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).jointAngleXZY], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'jointAngleErgo')
      nchan=length(find(contains(hdr.label, '_jointAngleErgo_')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).jointAngleErgo], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'jointAngleErgoXZY')
      nchan=length(find(contains(hdr.label, '_jointAngleErgoXZY')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).jointAngleErgoXZY], [nchan, nsample]);
      i=i+nchan;
  end
  
  if isfield(mvnx.subject.frames.frame,'centerOfMass')
      nchan=length(find(contains(hdr.label, '_centerOfMass')));
      dat(i:i+nchan-1,:) = reshape([mvnx.subject.frames.frame(idx).centerOfMass], [nchan, nsample]);
      i=i+nchan;
  end
  
  % only select data that is asked for
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events
  ft_warning('reading of events is not yet implemented');
  
  evt = [];
  
  % return the events
  varargout = {evt};
  
end
