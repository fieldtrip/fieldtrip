function ft_progress(varargin)

% FT_PROGRESS shows a graphical or non-graphical progress indication similar to the
% standard WAITBAR function, but with the extra option of printing it in the command
% window as a plain text string or as a rotating dial. Alternatively, you can also
% specify it not to give feedback on the progress.
%
% Prior to the for-loop, you should call either
%   ft_progress('init', 'none',    'Please wait...')
%   ft_progress('init', 'text',    'Please wait...')
%   ft_progress('init', 'textbar', 'Please wait...')      % ascii progress bar
%   ft_progress('init', 'dial',    'Please wait...')      % rotating dial
%   ft_progress('init', 'etf',     'Please wait...')      % estimated time to finish
%   ft_progress('init', 'gui',     'Please wait...')
%
% In each iteration of the for-loop, you should call either
% ft_progress(x)                                          % only show percentage
% ft_progress(x, 'Processing event %d from %d', i, N)     % show string, x=i/N
%
% After finishing the for-loop, you should call
%   ft_progress('close')
%
% Here is an example for the use of a progress indicator
%   ft_progress('init', 'etf', 'Please wait...');
%   for i=1:100
%     ft_progress(i/100, 'Processing event %d from %d', i, 100);
%     pause(0.03);
%   end
%   ft_progress('close')
%
% See also WAITBAR

% Copyright (C) 2004-2022, Robert Oostenveld
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

persistent p         % the previous value of the progress
persistent c         % counter for the number of updates that is done
persistent h         % the handle of the dialog (in case of type=gui)
persistent t0        % initial time, required for ETF
persistent p0        % initial percentage, required for ETF
persistent type      % type of feedback, string with none, gui, text, textcr, textnl
persistent title     % the string containing the title
persistent angle     % the angle in degrees, for dial or textbar
persistent strlen    % the length of the previously printed string, used to remove it by \b
persistent stopwatch % the time of previous invocation, used to restrict number of updates
persistent closing
persistent lastargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>1 && ischar(varargin{1}) && strcmp(varargin{1}, 'init')
  % reset these to the defaults
  p = 0;
  c = 0;
  h = 0;
  t0 = [];
  p0 = [];
  type = 'none'; % set below
  title = ''; % set below
  angle = 0;
  strlen = 0;
  stopwatch = tic();
  lastargin = [];
  closing = false;

  % determine the type of feedback
  type = varargin{2};
  if strcmp(type, 'textcr') || strcmp(type, 'textnl')
    type = 'text';
  end
  % determine the title of the dialog
  if nargin>2
    title = varargin{3};
  else
    title = '';
  end
  switch type
    case 'gui'
      % initialise the waitbar dialog
      if ~isempty(title)
        h = waitbar(0, title);
      else
        h = waitbar(0, 'Please wait');
      end
    case {'text', 'textnl', 'textcr'}
      if ~isempty(title)
        % print the title to the screen and go to the next line
        fprintf('%s\n', title)
      end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin==1 && ischar(varargin{1}) && strcmp(varargin{1}, 'close')

  if ~isempty(lastargin)
    % the last input argument is used when ft_progress('close') is called but the
    % previous invocation was not processed yet due to the restriction in the number of
    % updates to once every 100 ms
    closing = true;
    ft_progress(lastargin{:});
  end

  switch type
    case 'gui'
      % close the waitbar dialog
      close(h);
    case {'text', 'etf', 'dial', 'textbar'}
      % finish by going to the next line
      fprintf('\n');
  end

  % reset these to the defaults
  p = 0;
  c = 0;
  h = 0;
  t0 = [];
  p0 = [];
  type = 'none';
  title = '';
  angle = 0;
  strlen = 0;
  stopwatch = [];
  lastargin = [];
  closing = false;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

  % make sure we don't update more than once every 100 ms, otherwise certain
  % conditions result in a significant performance hit
  if ~isempty(stopwatch) && toc(stopwatch) < 0.1 && ~closing
    lastargin = varargin;
    return;
  end
  stopwatch = tic();
  lastargin = [];

  if strcmp(type, 'dial')
    % display should always be updated for the dial
  elseif (varargin{1}-p)<0.01 && strcmp(type, 'gui')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(type, 'textbar')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(type, 'etf')
    % display should not be updated it the difference is less than one percent
    return;
  end

  % count the number of updates, for debugging
  c = c+1;

  % remember the current value for the next function call
  p = varargin{1};

  switch type
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

        varargin{2} = [repmat(sprintf('\b'), [1 strlen]) varargin{2}];

        part1 = sprintf(varargin{2:end});
        part2 = sprintf(' - estimated time to finish is %d seconds', round(elapsed*(1-p)/(p-p0)));

        fprintf([part1 part2]);

        % record actual string length that was printed (subtracting all the \b's)
        strlen = length(part1) + length(part2) - strlen;
      else
        % only print the estimated time to finish
        part1 = sprintf([repmat(sprintf('\b'), [1 strlen]) ' - estimated time to finish is %d seconds'], round(elapsed*(1-p)/(p-p0)));
        fprintf(part1);
        strlen = length(part1) - strlen + 1;
      end

    case 'dial'
      dial = '|/-\|/-\';
      if ~isempty(title)
        % print the title and draw a new hand of the rotating dial
        part1 = sprintf([repmat(sprintf('\b'), [1 strlen]) '%s %s'], title, dial(1+angle/45));
        strlen = length(title) + 2;
      else
        % draw a new hand of the rotating dial
        part1 = sprintf([repmat(sprintf('\b'), [1 strlen]) '%s'], dial(1+angle/45));
        strlen = 1;
      end
      if part1(end) == '\'
        % a single backslash is not valid, add another one
        part1(end+1) = '\';
      end
      fprintf(part1);
      % increment the angle with 45 degrees
      angle = angle + 45;
      if angle==360
        % reset the angle to 0 degrees
        angle = 0;
      end

    case 'textbar'
      dial = '|/-\|/-\';
      % construct the line looking like this [------/ ]
      len = 75 - length(title) - 3;
      len1 = round(p*len); % number of '-' characters before the dial
      len2 = len - len1; % number of ' ' characters after the dial
      line = [title, ' [' repmat('-', 1, len1), dial(1+angle/45), repmat(' ', 1, len2), ']'];
      if c~=1
        backline = repmat('\b', [1 length(line)]);
      else
        backline = '';
      end
      % don't use carriage return, it sometimes leads to a new line
      fprintf([backline '%s'], line);
      % increment the angle with 45 degrees
      angle = angle + 45;
      if angle==360
        % reset the angle to 0 degrees
        angle = 0;
      end

    case 'text'
      if nargin>1
        % ensure the string does not end with a newline or carriage return, since
        % either would break compatibility with a -nodesktop matlab environment
        if length(varargin{2})>1 && (all(varargin{2}((end-1):end) == '\r')...
            || all(varargin{2}((end-1):end) == '\n'))
          varargin{2} = varargin{2}(1:end-2);
        end

        varargin{2} = [repmat(sprintf('\b'), [1 strlen]) varargin{2}];
        if usejava('desktop')
          % a newline is appropriate when using the desktop environment
          varargin{2} = [varargin{2} '\n'];
        end

        part1 = sprintf(varargin{2:end});
        fprintf(part1);
        strlen = length(part1) - strlen;

      else
        part1 = sprintf([repmat(sprintf('\b'), [1 strlen]) '%6.2f %%'], 100*varargin{1});
        fprintf(part1);
        strlen = length(part1) - strlen;
      end

  end % case gui, etf, dial, textbar, text
end % updating the displayed value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some test code follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft_progress('init', 'gui'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'textbar'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'dial'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'none'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')

% ft_progress('init', 'text'); for i=1:100, ft_progress(i/100, '%d/%d' , i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'text'); for i=1:100, ft_progress(i/100, '%d/%d\n', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'text'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')

% ft_progress('init', 'textnl'); for i=1:100, ft_progress(i/100, '%d/%d' , i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'textnl'); for i=1:100, ft_progress(i/100, '%d/%d\n', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'textnl'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')

% ft_progress('init', 'textcr'); for i=1:100, ft_progress(i/100, '%d/%d' , i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'textcr'); for i=1:100, ft_progress(i/100, '%d/%d\n', i, 100); pause(0.03); end; ft_progress('close')
% ft_progress('init', 'textcr'); for i=1:100, ft_progress(i/100, '%d/%d\r', i, 100); pause(0.03); end; ft_progress('close')
