function [vol, sens, cfg] = prepare_headmodel(cfg, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that helps to prepare the electrodes/gradiometers and the volume
% this is used in sourceanalysis and dipolefitting
%
% This function will get the gradiometer/electrode definition from
%   cfg.channel
%   cfg.gradfile
%   cfg.elecfile
%   cfg.elec
%   cfg.grad
%   data.grad
%   data.elec
% and the volume conductor definition from
%   cfg.hdmfile
%   cfg.vol
%   data.vol
% Subsequently it will remove the gradiometers/electrodes that are not
% present in the data. Finally it with attach the gradiometers to a
% multi-sphere head model (if supplied) or attach the electrodes to
% a BEM head model.
%
% This function will return the electrodes/gradiometers in an order that is
% consistent with the order in cfg.channel, or in case that is empty in the
% order of the input electrode/gradiometer definition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2009, Robert Oostenveld
%
% $Log: prepare_headmodel.m,v $
% Revision 1.6  2009/03/11 11:26:21  roboos
% switched to using forwinv/prepare_vol_sens
%
% Revision 1.5  2009/01/21 12:47:46  sashae
% ensure vol is always a struct
%
% Revision 1.4  2009/01/19 12:09:00  roboos
% fixed typo in 4-sphere eeg code (thanks to John Iversen)
%
% Revision 1.3  2008/07/07 12:36:09  roboos
% fixed bug in skin-compartment determination for 4-sphere
%
% Revision 1.2  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.1  2008/04/10 07:57:37  roboos
% renamed from prepare_vol_sens into prepare_headmodel, based on rev
% 1.28. This is to avild a name clash with the lower level function
% that is part of the forwinv module (i.e. this is part of the
% rearrangement of high-level and low-level functionality between
% fieldtrip and teh seperate modules)
%
% Revision 1.28  2008/03/18 12:26:47  roboos
% use senstype function to add descriptive string to the output (sens.type=...)
%
% Revision 1.27  2008/03/05 10:50:06  roboos
% switched to read_sens for reading elec or grad structures from file
% moved convert_ama2vol to standalone function
%
% Revision 1.26  2007/04/19 17:15:56  roboos
% only initialize the nolte method when the gradiometer array is not empty (usefull for plotting)
%
% Revision 1.25  2006/07/20 15:04:23  roboos
% added instructive fprintf message
%
% Revision 1.24  2006/05/23 10:17:32  roboos
% changed some comments, no code changes
%
% Revision 1.23  2006/04/12 08:38:15  ingnie
% updated documentation
%
% Revision 1.22  2006/04/10 16:35:20  ingnie
% updated documentation
%
% Revision 1.21  2006/03/21 09:44:03  roboos
% implemented support for Guido Nolte's method using vol.type='nolte'
% for neuromag: store the surface normals in the vol.bnd.nrm (consistent with nolte)
%
% Revision 1.20  2005/12/14 10:42:26  roboos
% removed warning for synthetic gradiometers
% if topolabel is present in data, look at that instead of label (for ICA)
% always try to add vol.skin and vol.brain
%
% Revision 1.19  2005/11/16 09:14:06  roboos
% moved the conversion from dipoli/ama format to fieldtrip/vol format to subfunction
% changed the post-processing of the EEG BEM model: incorporate the tra and mat matrices into one matrix (for computational efficiency)
%
% Revision 1.18  2005/09/29 01:15:32  roboos
% changed construction of the localsphere per coil for CTF hdm file
%
% Revision 1.17  2005/09/29 00:47:05  roboos
% added support for mbfys_ama BEM volume conductor file
%
% Revision 1.16  2005/08/05 07:26:22  roboos
% switched to teh use of read_fcdc_elec for cfg.gradfile
%
% Revision 1.15  2005/07/18 10:13:23  roboos
% fixed reading gradiometer from cfg.gradfile for CTF using ctf2grad
% added reading gradiometer from fif file
%
% Revision 1.14  2005/06/28 19:50:42  roboos
% minor change to solve problem when compiling to c-code, functionality is the same
% load(cfg.filename) gives compilation error, so first copy the string with the filename from the structure in a plain variable and then load(...)
%
% Revision 1.13  2005/06/08 16:33:54  roboos
% changed reading of gradiometers from matlab file
% prune the coils for a higher-order gradiometer system (i.e. remove non-contributing coils) -> this temporary solves a problem with multisphere headmodels
%
% Revision 1.12  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.11  2005/03/03 10:52:37  roboos
% changed the handling of the channel selection. The input cfg now specifies
% the channels that should be kept in the sensor array. If no selection is
% specified in the input cfg (e.g. for dipolesimulation), the default is 'all'.
%
% Revision 1.10  2005/01/26 08:05:23  roboos
% added empty data if second input argument not given
%
% Revision 1.9  2005/01/17 14:52:50  roboos
% changed to use read_fcdc_elec
%
% Revision 1.8  2004/12/08 18:00:13  roboos
% implemented consistent method of selecting a subset of channels for
% forward and inverse computations using cfg.channel and updated the
% ducumentation
%
% Revision 1.7  2004/09/06 08:45:38  roboos
% moved reading of neuromag BEM bondary from find_inside_vol into prepare_vol_sens
%
% Revision 1.6  2004/09/03 09:15:28  roboos
% added channelselection to the volume for neuromag
%
% Revision 1.5  2004/09/03 06:39:09  roboos
% fixed bug in non-functional neuromag section, cfg->vol
%
% Revision 1.4  2004/09/01 17:30:09  roboos
% added some explanation to the chansel for neuromag
%
% Revision 1.3  2004/08/31 13:55:22  roboos
% added initialization of neuromag megmodel
%
% Revision 1.2  2004/08/06 08:54:53  roboos
% fixed bug for gradiometer info coming from matlab file
%
% Revision 1.1  2004/06/28 08:59:38  roboos
% moved files from fieldtrip to fieldtrip/private
%
% Revision 1.2  2004/05/08 21:06:25  roberto
% added support for reading gradiometer from matlab file
%
% Revision 1.1  2004/04/08 15:51:20  roberto
% initial submissiion into cvs
%

% set the defaults
if ~isfield(cfg, 'channel'), cfg.channel = 'all';   end
if ~isfield(cfg, 'order'),   cfg.order = 10;        end % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

if nargin<2
  data = [];
end

% get the volume conduction model
if isfield(cfg, 'hdmfile')
  fprintf('reading headmodel from file ''%s''\n', cfg.hdmfile);
  vol = read_vol(cfg.hdmfile);
elseif isfield(cfg, 'vol')
  fprintf('using headmodel specified in the configuration\n');
  vol = cfg.vol;
elseif isfield(data, 'vol')
  fprintf('using headmodel specified in the data\n');
  vol = data.vol;
else
  error('no headmodel specified');
end

% get the gradiometer or electrode definition
if isfield(cfg, 'gradfile')
  fprintf('reading gradiometers from file ''%s''\n', cfg.gradfile);
  sens = read_sens(cfg.gradfile);
elseif isfield(cfg, 'grad')
  fprintf('using gradiometers specified in the configuration\n');
  sens = cfg.grad;
elseif isfield(data, 'grad')
  fprintf('using gradiometers specified in the data\n');
  sens = data.grad;
elseif isfield(cfg, 'elecfile')
  fprintf('reading electrodes from file ''%s''\n', cfg.elecfile);
  sens = read_sens(cfg.elecfile);
elseif isfield(cfg, 'elec')
  fprintf('using electrodes specified in the configuration\n');
  sens = cfg.elec;
elseif isfield(data, 'elec')
  fprintf('using electrodes specified in the data\n');
  sens = data.elec;
else
  error('no electrodes or gradiometers specified');
end

if isfield(data, 'topolabel')
  % the data reflects a componentanalysis, where the topographic and the
  % timecourse labels are different
  cfg.channel = channelselection(cfg.channel, data.topolabel);
elseif isfield(data, 'label')
  % In the subsequent code, the matching channels in the sensor array and
  % in the configuration will be selected. To ensure that these channels
  % are also present in the data, update the configuration to match the data.
  cfg.channel = channelselection(cfg.channel, data.label);
else
  % update the selected channels based on the electrode/gradiometer definition
  cfg.channel = channelselection(cfg.channel, sens.label);
end

% ensure that these are a struct, which may be required in case configuration tracking is used
vol  = struct(vol);
sens = struct(sens);

% the prepare_vol_sens function from the forwinv module does most of the actual work
[vol, sens] = prepare_vol_sens(vol, sens, 'channel', cfg.channel, 'order', cfg.order);

% update the selected channels in the configuration
cfg.channel = sens.label;
