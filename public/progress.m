function progress(varargin)

% PROGRESS shows a graphical or non-graphical progress indication similar
% to the standard Matlab WAITBAR function, but with the extra option of
% printing it in the command window as a plain text string or as a rotating
% dial. Alternatively, you can also specify it not to give feedback on the
% progress.
%
% Prior to the for-loop, you should call either
%   progress('init', 'none',    'Please wait...')
%   progress('init', 'gui',     'Please wait...')
%   progress('init', 'etf',     'Please wait...')      % estimated time to finish
%   progress('init', 'dial',    'Please wait...')      % rotating dial
%   progress('init', 'textbar', 'Please wait...')      % ascii progress bar
%   progress('init', 'text',    'Please wait...')
%   progress('init', 'textcr',  'Please wait...')      % force cariage return
%   progress('init', 'textnl',  'Please wait...')      % force newline
%
% In each iteration of the for-loop, you should call either
%   progress(x)                                       % only show percentage
%   progress(x, 'Processing event %d from %d', i, N)  % show string, x=i/N
%
% After finishing the for-loop, you should call
%   progress('close')
%
% Here is an example for the use of a progress indicator
%    progress('init', 'etf',     'Please wait...');
%    for i=1:42
%      progress(i/42, 'Processing event %d from %d', i, 42);
%      pause(0.1);
%    end
%    progress('close')

% Copyright (C) 2004-2008, Robert Oostenveld
%
% $Log: progress.m,v $
% Revision 1.2  2008/11/13 10:59:19  roboos
% added estimated time to finish (etf) and example
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.8  2007/05/08 20:55:03  roboos
% replaced all single & by double &&, hopefully resulting in a small speedup
%
% Revision 1.7  2006/12/11 10:52:32  roboos
% removed "feature accel off", since it was not solving the problem.
% The acceleration feature has to be disabled BEFORE this function
% is called the first time.
%
% Revision 1.6  2006/10/26 10:44:00  roboos
% fixed bug in double variable use
%
% Revision 1.5  2006/10/26 09:27:06  roboos
% added "feature accel off" to prevent matlab 7.3 from crashing
%
% Revision 1.4  2006/06/12 11:07:39  roboos
% modified the 1% update criterium, fixed a bug in the length of the textbar
%
% Revision 1.3  2005/03/31 12:22:38  roboos
% added ascii-art "textbar" option using the full screen width
%
% Revision 1.2  2004/10/22 07:23:39  roboos
% improved help, comments and code layout
%
% Revision 1.1  2004/10/21 17:37:17  roboos
% new implementation, to be used to replace subfunction in  beamformer scan
% and in other fieldtrip functions
%

persistent p        % the previous value of the progress
persistent c        % counter for the number of updates that is done
persistent t0       % initial time, required for ETF
persistent p0       % initial percentage, required for ETF
persistent t        % type of feedback, string with none, gui, text, textcr, textnl
persistent h        % the handle of the dialog (in case of type=gui)
persistent a        % the angle in degrees, for dial or textbar
persistent s        % the string containing the title

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>1 && ischar(varargin{1}) && strcmp(varargin{1}, 'init')
  a = 0;
  p = 0;
  h = 0;
  c = 0;
  % determine the type of feedback
  t = varargin{2};
  % determine the title of the dialog
  if nargin>2
    s = varargin{3};
  else
    s = '';
  end
  switch t
  case 'gui'
    % initialise the waitbar dialog
    if ~isempty(s)
      h = waitbar(0, s);
    else
      h = waitbar(0, 'Please wait');
    end
  case {'text', 'textnl', 'textcr'}
    if ~isempty(s)
      % print the title to the screen and go to the next line
      fprintf('%s\n', s)
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin==1 && ischar(varargin{1}) && strcmp(varargin{1}, 'close')
  switch t
  case 'gui'
    % close the waitbar dialog
    close(h);
  case {'textcr', 'dial', 'textbar'}
    % finish by going to the next line
    fprintf('\n');
  end
  % reset these to the defaults
  a  = 0;
  h  = 0;
  p  = 0;
  t  = 'none';
  s  = '';
  t0 = [];
  p0 = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  if strcmp(t, 'dial')
    % display should always be updated for the dial
    % continue;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'gui')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'textbar')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'etf')
    % display should not be updated it the difference is less than one percent
    return;
  end

  % count the number of updates, for debugging
  c = c+1;

  % remember the current value for the next function call
  p = varargin{1};

  switch t
  case 'gui'
    % update the the length of the bar in the waitbar dialog
    waitbar(varargin{1}, h);

  case 'etf'
    % compute the estimated time that the computation still needs to finish
    if isempty(t0) || isempty(p0)
      t0 = clock;
      p0 = p;
    end
    elapsed = etime(clock, t0);
    if nargin>1 && ~isempty(varargin{2})
      % include the specified string
      fprintf(varargin{2:end});
      fprintf(' - estimated time to finish is %d seconds\n', round(elapsed*(1-p)/(p-p0)));
    else
      % only print the estimated time to finish
      fprintf(' - estimated time to finish is %d seconds\n', round(elapsed*(1-p)/(p-p0)));
    end

  case 'dial'
    dial = '|/-\|/-\';
    if ~isempty(s)
      % print the title and draw a new hand of the rotating dial
      fprintf('\r%s %s', s, dial(1+a/45));
    else
      % draw a new hand of the rotating dial
      fprintf('\r%s', dial(1+a/45));
    end
    % increment the angle with 45 degrees
    a = a + 45;
    if a==360
      % reset the angle to 0 degrees
      a = 0;
    end

  case 'textbar'
    dial = '|/-\|/-\';
    % construct the line looking like [------/          ]
    len  = 75 - length(s) - 3;
    len1 = round(p*len);            % number of '-' characters before the dial
    len2 = len - len1;              % number of ' ' characters after the dial
    line = [s, ' [' repmat('-',1,len1), dial(1+a/45), repmat(' ',1,len2) ,']'];
    fprintf('\r%s', line);
    % increment the angle with 45 degrees
    a = a + 45;
    if a==360
      % reset the angle to 0 degrees
      a = 0;
    end

  case 'text'
    if nargin>1
      % print the string as it is
      fprintf(varargin{2:end});
    else
      fprintf('%6.2f %%\n', 100*varargin{1});
    end

  case 'textnl'
    if nargin>1
      % ensure that the string ends with a newline
      if length(varargin{2})>1 && all(varargin{2}((end-1):end) == '\r')
        varargin{2}((end-1):end) = '\n';
      elseif length(varargin{2})>1 && ~all(varargin{2}((end-1):end) == '\n')
        varargin{2}((end+1):(end+2)) = '\n';
      elseif length(varargin{2})<2
        varargin{2}((end+1):(end+2)) = '\n';
      end
      fprintf(varargin{2:end});
    else
      fprintf('%6.2f %%\n', 100*varargin{1});
    end

  case 'textcr'
    if nargin>1
      % ensure that the string ends with a cariage return
      if length(varargin{2})>1 && all(varargin{2}((end-1):end) == '\n')
        varargin{2}((end-1):end) = '\r';
      elseif length(varargin{2})>1 && ~all(varargin{2}((end-1):end) == '\r')
        varargin{2}((end+1):(end+2)) = '\r';
      elseif length(varargin{2})<2
        varargin{2}((end+1):(end+2)) = '\r';
      end
      fprintf(varargin{2:end});
    else
      fprintf('%6.2f %%\r', 100*varargin{1});
    end

  end % case gui, dial, text, textnl, textcr
end % updating the displayed value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some test code follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% progress('init', 'gui');  for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')
% progress('init', 'dial'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')
% progress('init', 'none'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')


