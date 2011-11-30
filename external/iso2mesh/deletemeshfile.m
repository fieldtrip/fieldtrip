function flag=deletemeshfile(fname)
%
% flag=deletemeshfile(fname)
%
% delete a given work mesh file under the working directory
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input: 
%     fname: specified file name (without path)
%
% output:
%     flag: not used
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

try
    if(exist(fname)) 
	delete(fname); 
    end
catch
    error(['You do not have permission to delete temporary files. If you are working in a multi-user ',...
         'environment, such as Unix/Linux and there are other users using iso2mesh, ',...
         'you may need to define ISO2MESH_SESSION=''yourstring'' to make your output ',...
         'files different from others; if you do not have permission to ',mwpath(''),...
         ' as the temporary directory, you have to define ISO2MESH_TEMP=''/path/you/have/write/permission'' ',...
         'in matlab/octave base workspace.']);
end

