function version = om_getgccversion
% checks for gcc compiler version (works if superior to gcc3)
% Copyright (C) 2010, Alexandre Gramfort, INRIA

% $Id$
% $LastChangedBy: alegra $
% $LastChangedDate: 2010-04-19 11:03:39 +0200 (Mon, 19 Apr 2010) $
% $Revision$

tmpdir = pwd;
cd /tmp
[junk,tname] = fileparts(tempname);
txtfile  = [tname '.txt'];
dos(['gcc -v >& ' txtfile]);
efid = fopen(txtfile);

tmp = ''; cnt = 1;
vec = [];
while ~isnumeric(tmp)
    tmp = fgetl(efid);
    vec{cnt} = tmp;
    cnt = cnt+1;
end
fclose(efid);
delete(txtfile);
cd(tmpdir);
tmp = deblank(vec{cnt-2});
num = findstr('gcc version ',tmp);
version = str2num(tmp(num+11:num+12));
if isempty(version) && findstr('clang',vec{2}) > 0
    version = 5;
end
