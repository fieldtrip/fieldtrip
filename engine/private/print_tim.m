%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for pretty-printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = print_tim(tim)
% partition the time in seconds into years, months, etc.
year   = 60*60*24*7*365.25;
month  = 60*60*24*365.25/12;
week   = 60*60*24*7;
day    = 60*60*24;
hour   = 60*60;
minute = 60;
org = tim; % remember the original time in seconds
Y = floor(tim/year  ); tim = tim - Y*year;
M = floor(tim/month ); tim = tim - M*month;   % note capital M
w = floor(tim/week  ); tim = tim - w*week;
d = floor(tim/day   ); tim = tim - d*day;
h = floor(tim/hour  ); tim = tim - h*hour;
m = floor(tim/minute); tim = tim - m*minute;  % note small m
s = tim;
if Y>=1
  str = sprintf('%.1f years', org/year);
elseif M>=1
  str = sprintf('%.1f months', org/month);
elseif w>=1
  str = sprintf('%.1f weeks', org/week);
elseif d>=1
  str = sprintf('%.1f days', org/day);
elseif h>=1
  str = sprintf('%.1f hours', org/hour);
elseif m>=1
  str = sprintf('%.1f minutes', org/minute);
else
  % note that timreq and timavail are implemented as integers
  str = sprintf('%.0f seconds', org);
end

