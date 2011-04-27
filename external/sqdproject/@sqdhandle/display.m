function display(t)
% DISPLAY Display object properties


newline = sprintf('\n');
for i = 1:length(t)
    disp('Sqdhandle object properties: ');
    disp('---------------------------');
    get(t(i));
end;