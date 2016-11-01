function p=getvarfrom(ws,name)
%
% p=getvarfrom(ws,name)
%
% get variable value by name from specified work-space
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    ws: name of the work-space, for example, 'base'
%    name: name string of the variable
%
% output:
%    p: the value of the specified variable, if the variable does not
%       exist, return empty array
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

wsname=ws;
if(~iscell(ws))
   wsname=cell(1);
   wsname{1}=ws;
end

p=[];
for i=1:length(wsname)
    isdefined=evalin(wsname{i},['exist(''' name ''')']);
    if(isdefined==1)
        p=evalin(wsname{i},name);
        break;
    end
end
