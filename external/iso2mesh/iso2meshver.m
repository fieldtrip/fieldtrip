function [major,minor,patchnum,extra]=iso2meshver
%
% [major,minor,patchnum,extra]=iso2meshver
%      or
% v=iso2meshver
%
% get the version number of iso2mesh toolbox
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% output:
%    if you ask for a single output:
%      v: a string denotes the current version number; the string is 
%       typically in the following format: "major.minor.patch-extra"
%       where major/minor/patch are typically integers, and extra can
%       be an arbitrary string and is optional
%    if you ask for 4 outputs:
%     [major,minor,patchnum,extra] are each field of the version string
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

major=1;
minor=5;
patchnum=0;
extra='$Rev:: 500 $';
extra=regexprep(extra,'[\s$:]', '');

iso2meshvstr=sprintf('%d.%d.%d',major,minor,patchnum);
if(~isempty(extra))
   iso2meshvstr=[iso2meshvstr '-' extra];
end

if(nargout==0)
   fprintf(1,'iso2mesh toolbox version: %s\n',iso2meshvstr);
   clear major;
elseif(nargout==1)
   major=iso2meshvstr;
elseif(nargout~=4)
   error('you need to return either 1 or 4 output variables');
end
