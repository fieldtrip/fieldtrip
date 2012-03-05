%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for pretty-printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = print_mem(val)
if val<1024
  str = sprintf('%d bytes', val);
elseif val<1024^2
  str = sprintf('%.1f KB', val/1024);
elseif val<1024^3
  str = sprintf('%.1f MB', val/1024^2);
else
  str = sprintf('%.1f GB', val/1024^3);
end

