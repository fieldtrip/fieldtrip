function ft_trackusage(event, varargin)

% FT_TRACKUSAGE tracks the usage of specific FieldTrip components using a central
% tracking server. This involves sending a small snippet of information to the
% server. Tracking is only used to gather data on the usage of the FieldTrip
% toolbox, to get information on the number of users and on the frequency of use
% of specific toolbox functions. This allows the toolbox developers to improve the
% FIeldTrip toolbox source code, documentation and to provide better support.
%
% This function will NOT upload any information about the data, nor about the
% configuration that you are using in your analyses.
%
% This function will NOT upload any identifying details about you. Your username
% and computer name are "salted" and subsequently converted with the MD5
% cryptographic hashing function into a unique identifier. Not knowing the salt,
% it is impossible to decode these MD5 hashes and recover the original
% identifiers.
%
% It is possible to disable the tracking for all functions by specifying
% the following
%   global ft_defaults
%   ft_default.trackusage = 'no'
%
% See the following online documentation for more information
%   http://en.wikipedia.org/wiki/MD5
%   http://en.wikipedia.org/wiki/Salt_(cryptography)
%   http://www.fieldtriptoolbox.org/faq/tracking
%
% See also FT_DEFAULTS

% Copyright (C) 2015, Robert oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/donders/fieldtrip
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

global ft_default
persistent initialized

if isempty(initialized)
  initialized = false;
end

if nargin<1
  % there is nothing to track
  return
end

%% Since the functionality is still in beta testing, only enable the tracking for some developpers
knownuser = false;
knownuser = knownuser || (strcmp(getusername, 'roboos')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^mac011', 'once'))));
%knownuser = knownuser || (strcmp(getusername, 'jansch')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
knownuser = knownuser || (strcmp(getusername, 'jimher')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
%knownuser = knownuser || (strcmp(getusername, 'nielam')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
knownuser = knownuser || (strcmp(getusername, 'tzvpop')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
knownuser = knownuser || (strcmp(getusername, 'lucamb')  && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
knownuser = knownuser || (strcmp(getusername, 'elivzan') && (~isempty(regexp(gethostname, '^dccn', 'once')) || ~isempty(regexp(gethostname, '^fcdc', 'once'))));
if ~knownuser
  return
end

if ~strcmp(mfilename, 'ft_trackusage')
  % this function should not be used outside of the FieldTrip toolbox without updating the token (see below)
  return
end

if ~ft_platform_supports('rng')
  % this function should not (yet) be used on Octave
  return
end

%% The first part pertains to keeping the tracking settings consistent over multiple MATLAB sessions

% This functionality overlaps in part with what normally would be done using
% ft_defaults, but is replicated here to make the tracking independent from the path
% settings in ft_defaults.

% locate the file that contains the persistent FieldTrip preferences
fieldtripprefs  = fullfile(prefdir, 'fieldtripprefs.mat');

if ~isfield(ft_default, 'trackusage')
  % read options from the preferences file
  if exist(fieldtripprefs, 'file')
    prefs      = load(fieldtripprefs); % the file contains multiple fields
    ft_default = mergeconfig(ft_default, prefs);
  end
end

if ~isfield(ft_default, 'trackusage')
  % the default is to allow tracking
  % create a salt for one-way encryption of identifying information
  rng('shuffle');
  trackusage = dec2hex(intmax('uint32')*rand(1));  % create a secret salt, this is never shared
  warning('enabling online tracking of FieldTrip usage, see http://www.fieldtriptoolbox.org/faq/tracking');
  if exist(fieldtripprefs, 'file')
    % update the existing preferences file
    save(fieldtripprefs, 'trackusage', '-append');
  end
  % keep it in the global variable
  ft_default.trackusage = trackusage;
  clear trackusage
end

if ~exist(fieldtripprefs, 'file')
  % save it to a new preferences file
  trackusage = ft_default.trackusage;
  save(fieldtripprefs, 'trackusage');
  clear trackusage
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The second part pertains to the actual tracking

if isequal(ft_default.trackusage, false) || isequal(ft_default.trackusage, 'no') || isequal(ft_default.trackusage, 'off')
  return
end

% this are the default properties to track
properties.token       = '1187d9a6959c39d0e733d6273d1658a5'; % this is specific for the FieldTrip project
properties.user        = ft_hash(sprintf('%s%s', ft_default.trackusage, getusername)); % hash it with a secret salt
properties.host        = ft_hash(sprintf('%s%s', ft_default.trackusage, gethostname)); % hash it with a secret salt
properties.matlab      = version('-release');
properties.fieldtrip   = ft_version;
properties.computer    = lower(computer);
properties.distinct_id = properties.user; % this links the event to the profile

% add the custom properties, these come in key-value pairs
for i=1:2:numel(varargin)
  properties.(varargin{i}) = varargin{i+1};
end

% construct the HTTP request for Mixpanel, see https://mixpanel.com/help/reference/http
event_json   = sprintf('{"event": "%s", "properties": {%s}}', event, ft_struct2json(properties));
event_base64 = base64encode(event_json);
event_http   = sprintf('http://api.mixpanel.com/track/?data=%s', event_base64);


[output, status] = ft_urlread(event_http);
if ~status
  disp(output);
  warning('could not send tracker information for "%s"', event);
end

if ~initialized
  % this only gets send once
  user_json   = sprintf('{"$token": "%s", "$distinct_id": "%s", "$ip": "%s", "$set": {} }',  properties.token, properties.user, getaddress());
  user_base64 = base64encode(user_json);
  user_http   = sprintf('http://api.mixpanel.com/engage/?data=%s', user_base64);

  [output, status] = ft_urlread(user_http);
  if ~status
    disp(output);
    warning('could not send tracker information for "%s"', event);
  end

  initialized = true;
end % if initialized


