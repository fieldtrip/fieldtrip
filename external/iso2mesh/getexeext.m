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
   if(~ispc)
      if(~ismac)
	   if(isempty(regexp(computer,'86_64')))
	      exesuff='.mexglx';
	   else
              exesuff='.mexa64';
	   end
      else
           if(isempty(regexp(computer,'86_64')))
              exesuff='.mexmaci';
           else
              exesuff='.mexmaci64';
           end
      end
   else
      exesuff='.exe';
   end
end
