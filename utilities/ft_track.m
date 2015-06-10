function ft_track(event)

% FT_TRACK tracks the usage of specific FieldTrip components using a central
% tracking server. This involves sending a small snippet of information to the
% server. Tracking is only used to gather data on the usage of the FieldTrip
% toolbox, to get information on the number of users and on the frequency of use
% of specific toolbox functions. This allows the toolbox developers to improve the
% toolbox and to provide better support.
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
% It is possible to enable or disable the tracking with the option
%   ft_default.track.allow = true or false (default = true)
%
% See also
%   http://en.wikipedia.org/wiki/MD5
%   http://en.wikipedia.org/wiki/Salt_(cryptography)
%   http://www.fieldtriptoolbox.org/faq/tracking

global ft_default

%% Since the functionality is still in beta testing, disable the tracking for most users
if ~strcmp(getusername, 'roboos')
  return
end

if ~(isempty(regexp(gethostname, '^mac011', 'once')) || isempty(regexp(gethostname, '^dccn', 'once')))
  return
end

if ~strcmp(mfilename, 'ft_track')
  return
end

%% The first part pertains to keeping the tracking settings consistent over multiple MATLAB sessions

% this functionality overlaps in part with ft_defaults but is replicated here to
% make the tracking independent from the path settings in ft_defaults

% locate the file that contains the persistent FieldTrip preferences
fieldtripprefs  = fullfile(prefdir, 'fieldtripprefs.mat');

if ~isfield(ft_default, 'track') || ~isfield(ft_default.track, 'salt')
  % try to read the tracking options from the preferences file
  if exist(fieldtripprefs, 'file')
    prefs      = load(fieldtripprefs); % the file contains multiple fields
    ft_default = mergeconfig(ft_default, prefs);
  end
end

if ~isfield(ft_default, 'track') || ~isfield(ft_default.track, 'salt')
  % the tracking options still don't exist, create them on the fly
  rng('shuffle');
  track.salt  = dec2hex(intmax('uint32')*rand(1));  % create a secret salt, this is never shared
  track.allow = true;                               % default is to allow tracking
  warning('enabling online tracking of FieldTrip usage, see http://www.fieldtriptoolbox.org/faq/tracking');
  if exist(fieldtripprefs, 'file')
    % update the existing preferences file
    save(fieldtripprefs, 'track', '-append');
  end
  % keep it in the global variable
  ft_default.track = track;
  clear track
end

if ~exist(fieldtripprefs, 'file')
  % save it to a new preferences file
  track = ft_default.track;
  save(fieldtripprefs, 'track');
  clear track
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The second part pertains to the actual tracking

if ~ft_default.track.allow
  return
end

% this contains the CalcMD5 function
ft_hastoolbox('fileexchange', 1);

properties.token     = '1187d9a6959c39d0e733d6273d1658a5'; % this is specific for the FieldTrip project
properties.user      = CalcMD5(sprintf('%s%s', ft_default.track.salt, getusername)); % hash it with a secret salt
properties.host      = CalcMD5(sprintf('%s%s', ft_default.track.salt, gethostname)); % hash it with a secret salt
properties.matlab    = version;
properties.fieldtrip = ft_version;
properties.computer  = lower(computer);

% construct the HTTP request for Mixpanel, see https://mixpanel.com/help/reference/http
event_json   = sprintf('{"event": "%s", "properties": {%s}}', event, struct2json(properties));
event_base64 = base64encode(event_json);
event_http   = sprintf('http://api.mixpanel.com/track/?data=%s', event_base64);

[output, status] = urlread(event_http, 'TimeOut', 15);

if ~status
  error('could not send tracker event for "%s"', event);
else
  fprintf('tracked the usage of "%s"\n', event);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function j = struct2json(s)
fn = fieldnames(s);
fv = cell(size(fn));
for i=1:numel(fn)
  val = s.(fn{i});
  switch class(val)
    case 'char'
      fv{i} = val;
    case 'double'
      fv{i} = num2str(val);
    otherwise
      error('class %s is not supported\n', type(val));
  end
end
f = cat(1, fn', fv');
j = sprintf('"%s": "%s", ', f{:});
j = j(1:end-2); % remove the last comma and space
