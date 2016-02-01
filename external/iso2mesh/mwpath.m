function tempname=mwpath(fname)
%
% tempname=meshtemppath(fname)
%
% get full temp-file name by prepend working-directory and current session name
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%    fname: input, a file name string
%
% output:
%    tempname: output, full file name located in the working directory
%
%    if global variable ISO2MESH_TEMP is set in 'base', it will use it
%    as the working directory; otherwise, will use matlab function tempdir
%    to return a working directory.
%
%    if global variable ISO2MESH_SESSION is set in 'base', it will be
%    prepended for each file name, otherwise, use supplied file name.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

p=getvarfrom({'caller','base'},'ISO2MESH_TEMP');
session=getvarfrom({'caller','base'},'ISO2MESH_SESSION');

username=getenv('USER'); % for Linux/Unix/Mac OS

if(isempty(username))
   username=getenv('UserName'); % for windows
end

if(~isempty(username))
   username=['iso2mesh-' username];
end

tempname=[];
if(isempty(p))
      if(isoctavemesh & tempdir=='\')
		tempname=['.'  filesep session fname];
	else
		tdir=tempdir;
		if(tdir(end)~=filesep)
			tdir=[tdir filesep];
		end
		if(~isempty(username))
                    tdir=[tdir username filesep];
                    if(exist(tdir)==0) mkdir(tdir); end
        end
        if(nargin==0)
            tempname=tdir;
        else
            tempname=[tdir session fname];
        end
	end
else
	tempname=[p filesep session fname];
end
