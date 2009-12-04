function om_write_cond(condfile,c)
%   OM_WRITE_COND   
%       [] = OM_WRITE_COND(CONDFILE,C)
% 
%   Write conductivity file for OpenMEEG
%   
%   Created by Alexandre Gramfort on 2009-03-30.
%   Copyright (c)  Alexandre Gramfort. All rights reserved.

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information


if length(c) == 3
    cfid = fopen(condfile, 'w');
    fprintf(cfid,'# Properties Description 1.0 (Conductivities)\n');
    fprintf(cfid,'                                             \n');
    fprintf(cfid,'Air         0.0                              \n');
    fprintf(cfid,'Skin        %f                               \n',c(3));
    fprintf(cfid,'Brain       %f                               \n',c(1));
    fprintf(cfid,'Skull       %f                               \n',c(2));
    fclose(cfid);
elseif length(c) == 2
    cfid = fopen(condfile, 'w');
    fprintf(cfid,'# Properties Description 1.0 (Conductivities)\n');
    fprintf(cfid,'                                             \n');
    fprintf(cfid,'Air         0.0                              \n');
    fprintf(cfid,'Skin        %f                               \n',c(2));
    fprintf(cfid,'Skull       %f                               \n',c(1));
    fclose(cfid);
elseif length(c) == 1
    cfid = fopen(condfile, 'w');
    fprintf(cfid,'# Properties Description 1.0 (Conductivities)\n');
    fprintf(cfid,'                                             \n');
    fprintf(cfid,'Air         0.0                              \n');
    fprintf(cfid,'Head        %f                               \n',c(1));
    fclose(cfid);
else
    error('OpenMEEG only supports geometry with 2 or 3 interfaces')
end

end %  function