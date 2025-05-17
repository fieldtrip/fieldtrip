function exesuff = fallbackexeext(exesuffix, exename)
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

exesuff = exesuffix;
if (strcmp(exesuff, '.mexa64') && exist([mcpath(exename) exesuff], 'file') == 0) % fall back to i386 linux
    exesuff = '.mexglx';
end
if (strcmp(exesuff, '.mexmaci64') && exist([mcpath(exename) exesuff], 'file') == 0) % fall back to i386 mac
    exesuff = '.mexmaci';
end
if (strcmp(exesuff, '.mexmaci') && exist([mcpath(exename) exesuff], 'file') == 0) % fall back to ppc mac
    exesuff = '.mexmac';
end
if (exist([mcpath(exename) exesuff], 'file') == 0) % fall back to OS native package
    exesuff = '';
end

if (exist([mcpath(exename) exesuff], 'file') == 0)
    if (system(['which ' exename]) == 0)
        return
    end
    error(['The following executable:\n' ...
           '\t%s%s\n' ...
           'is missing. Please download it from ' ...
           'https://github.com/fangq/iso2mesh/tree/master/bin/ ' ...
           'and save it to the above path, then rerun the script.\n' ...
          ], mcpath(exename), getexeext);
end
