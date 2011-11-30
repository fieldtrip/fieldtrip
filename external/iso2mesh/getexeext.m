function exesuff=getexeext()
%
% exesuff=getexeext()
%
% get meshing external tool extension names for the current platform
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% output:
%     exesuff: file extension for iso2mesh tool binaries
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff='.exe';
if(isunix) 
	exesuff=['.',mexext];
end
if(isoctavemesh)
   if(isempty(strfind(computer,'msdos')) & isempty(strfind(computer,'mingw')))
	   if(isempty(regexp(computer,'86_64')))
	      exesuff='.mexglx';
	   else
              exesuff='.mexa64';
	   end
   else
      exesuff='.exe';
   end
end
