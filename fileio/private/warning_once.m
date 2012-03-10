function [ws warned] = warning_once(varargin)
%
% Use as one of the following
%   warning_once(string)
%   warning_once(string, timeout)
%   warning_once(id, string, timeout)
% where timeout should be inf if you don't want to see the warning ever
% again. The default timeout value is 60 seconds.
%
% It can be used instead of the MATLAB built-in function WARNING, thus as
%   s = warning_once(...)
% or as
%   warning_once(s)
% where s is a structure with fields 'identifier' and 'state', storing the
% state information. In other words, warning_once accepts as an input the
% same structure it returns as an output. This returns or restores the
% states of warnings to their previous values.
%
% It can also be used as
%    [s w] = warning_once(...)
% where w is a boolean that indicates whether a warning as been thrown or not.
%
% Please note that you can NOT use it like this
%   warning_once('the value is %d', 10)
% instead you should do
%   warning_once(sprintf('the value is %d', 10))

persistent stopwatch previous

if nargin < 1
  error('You need to specify at least a warning message');
end

warned = false;
if isstruct(varargin{1})
  warning(varargin{1});
  return;
end

if nargin==3
  msgid = varargin{1};
  msgstr = varargin{2};
  timeout = varargin{3};
elseif nargin==2
  msgstr= ''; % this becomes irrelevant
  msgid = varargin{1}; % this becomes the real msgstr
  timeout = varargin{2};
elseif nargin==1
  msgstr= ''; % this becomes irrelevant
  msgid = varargin{1}; % this becomes the real msgstr
  timeout = 60; % default timeout in seconds
end

if isempty(timeout)
  error('Timeout ill-specified');
end

if isempty(stopwatch)
  stopwatch = tic;
end
if isempty(previous)
  previous = struct;
end

now = toc(stopwatch); % measure time since first function call
fname = fixname([msgid '_' msgstr]); % make a nice string that is allowed as structure fieldname, copy the subfunction from  ft_hastoolbox
fname = decomma(fname);

if length(fname) > 63 % MATLAB max name
  fname = fname(1:63);
end

if isfield(previous, fname) && now>previous.(fname).timeout
  % it has timed out, give the warning again
  ws = warning(msgid, msgstr);
  previous.(fname).timeout = now+timeout;
  previous.(fname).ws = ws;
  warned = true;
elseif ~isfield(previous, fname)
  % the warning has not been issued before
  ws = warning(msgid, msgstr);
  previous.(fname).timeout = now+timeout;
  previous.(fname).ws = ws;
  warned = true;
else
  % the warning has been issued before, but has not timed out yet
  ws = previous.(fname).ws;
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function name = decomma(name)
name(name==',')=[];
end % function
