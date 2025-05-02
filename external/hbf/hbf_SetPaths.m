function hbf_SetPaths
% HBF_SETPATHS initializes the hbf library, setting needed paths
%
% v200331 Matti Stenroos
p=mfilename('fullpath');
[basedirectory,name,ext]=fileparts(p);
dirlist={'hbf_calc','hbf_mesh'};
fprintf('hbf: adding paths...\n');
for I=1:numel(dirlist)
    thisdir=fullfile(basedirectory,dirlist{I});
    addpath(thisdir)
    fprintf('\t%s\n',thisdir);
end
    