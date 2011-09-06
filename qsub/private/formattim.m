function str = formattim(tim)

% FORMATTIM returns the time as hh:mm:ss
%
% Use as
%   str = formattim(tim)

hour   = 60*60;
minute = 60;

h = floor(tim/hour  ); tim = tim - h*hour;
m = floor(tim/minute); tim = tim - m*minute;  % note small m

str = sprintf('%02d:%02d:%02d',h,m,tim);

