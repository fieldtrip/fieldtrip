function pausejava(tPauseSec)

% PAUSEJAVA uses the Java Virtual Machine to pause for a specified amount of time.
% If the JVM is not running, it defaults to the builtin MATLAB pause.
%
% Use as
%   pause(seconds)
%
% The builtin MATLAB pause function has a known memory leak in R2011b and R2012a.
% Whenever pause is called, the graphics event queue (EDT) is flushed, thereby
% updating all Matlab figure windows.
%
% http://undocumentedmatlab.com/blog/pause-for-the-better/
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1997

if usejava('jvm')
  th = java.lang.Thread.currentThread();  % Get current Thread
  th.sleep(1000*tPauseSec)                % Pause thread, conversion to milliseconds
else
  pause(tPauseSec)
end
