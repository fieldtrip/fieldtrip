function exesuff=fallbackexeext(exesuffix, exename)
%
% exesuff=fallbackexeext(exesuffix, exename)
%
% get the fallback external tool extension names for the current platform
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%     exesuffix: the output executable suffix from getexeext
%     exename: the executable name
%
% output:
%     exesuff: file extension for iso2mesh tool binaries
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff=exesuffix;
if(strcmp(exesuff,'.mexa64') & exist([mcpath(exename) exesuff],'file')==0) % fall back to i386 linux
        exesuff='.mexglx';
	return;
end
if(strcmp(exesuff,'.mexmaci64') & exist([mcpath(exename) exesuff],'file')==0) % fall back to i386 mac
        exesuff='.mexmaci';
end
if(strcmp(exesuff,'.mexmaci') & exist([mcpath(exename) exesuff],'file')==0) % fall back to ppc mac
        exesuff='.mexmac';
end

