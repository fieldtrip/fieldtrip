function status = om_check_vol(vol)
%   OM_CHECK_VOL   Check meshes of volume conductor for BEM modeling
%       [STATUS] = OM_CHECK_VOL(VOL)
%
%   returns 1 if there is a problem with geometry
%   else returns 0
%
% Copyright (C) 2010-2017, OpenMEEG developers

openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries

% the first compartment should be the skin, the last the source
% flip the order of the compartments if necessary
if vol.skin==length(vol.bnd) && vol.source==1
    vol.bnd    = fliplr(vol.bnd(:)');
    vol.skin   = 1;
    vol.source = length(vol.bnd);
end

assert(vol.skin == 1)
assert(vol.source == length(vol.bnd))

% Flip faces for openmeeg convention
for ii=1:length(vol.bnd)
    vol.bnd(ii).tri = fliplr(vol.bnd(ii).tri);
end

try
    % store the current path and change folder to the temporary one
    tmpfolder = cd;

    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for ii=1:length(vol.bnd)
        [junk,tname] = fileparts(tempname);
        bndfile{ii} = [tname '.tri'];
        om_save_tri(bndfile{ii}, vol.bnd(ii).pos, vol.bnd(ii).tri);
    end

    % these will hold the shell script and the inverted system matrix
    [junk,tname] = fileparts(tempname);
    if ~ispc
      exefile = [tname '.sh'];
    else
      exefile = [tname '.bat'];
    end

    [junk,tname] = fileparts(tempname);
    geomfile  = [tname '.geom'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);

    % Exe file
    status = system([prefix 'om_check_geom -g ' geomfile]);
    cleaner(vol,bndfile,geomfile)
    cd(tmpfolder)
catch
    cleaner(vol,bndfile,geomfile)
    cd(tmpfolder)
    rethrow(lasterror)
end

function cleaner(vol,bndfile,geomfile)

% delete the temporary files
for i=1:length(vol.bnd)
    delete(bndfile{i})
end

delete(geomfile);
return
