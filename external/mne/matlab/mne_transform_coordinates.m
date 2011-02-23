function [trans_pos] = mne_transform_coordinates(filename,pos,from,to)
%
% [trans_pos] = mne_transform_coordinates(filename,pos,from,to)
%
% Transform locations between various MRI-related coordinate frames
%
% filename   - Name of a fif file containing the coordinate transformations
%              This file can be conveniently created with mne_collect_transforms
% pos        - N x 3 array of locations to transform (in meters)
% from       - Coordinate frame of the above locations
%              Allowed choices are: FIFFV_COORD_MRI (surface RAS coordinates)
%              and FIFFV_COORD_HEAD (MEG head coordinates)
% to         - Coordinate frame of the result
%              Allowed choices are: FIFFV_COORD_MRI, FIFFV_COORD_HEAD,
%              FIFFV_MNE_COORD_MNI_TAL (MNI Talairach), and
%              FIFFV_MNE_COORD_FS_TAL (FreeSurfer Talairach)
%
%              All of the above constants are define in fiff_define_constants
%
% trans_pos  - The transformed locations
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   No part of this program may be photocopied, reproduced,
%   or translated to another program language without
%   prior written consent of the author.
%   Revision 1.1  2008/11/16 21:31:24  msh
%   Added mne_transform_coordinates and new coordinate frame definitions
%
%

me = 'mne:mne_transform_coordinates';
global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin ~= 4
    error(me,'Some input arguments are missing');
end
%
%   Read the fif file containing all necessary transformations
%
[ fid, tree, dir ] = fiff_open(filename);

T0      = [];
T1      = [];
T2      = [];
T3plus  = [];
T3minus = [];
for k = 1:length(dir)
    if dir(k).kind == FIFF.FIFF_COORD_TRANS
        tag = fiff_read_tag(fid,dir(k).pos);
        trans = tag.data;
        if trans.from == FIFF.FIFFV_COORD_MRI && trans.to == FIFF.FIFFV_COORD_HEAD
            T0 = fiff_invert_transform(trans);
        elseif trans.from == FIFF.FIFFV_COORD_MRI && trans.to == FIFF.FIFFV_MNE_COORD_RAS
            T1 = trans;
        elseif trans.from == FIFF.FIFFV_MNE_COORD_RAS && trans.to == FIFF.FIFFV_MNE_COORD_MNI_TAL
            T2 = trans;
        elseif trans.from == FIFF.FIFFV_MNE_COORD_MNI_TAL
            if trans.to == FIFF.FIFFV_MNE_COORD_FS_TAL_GTZ
                T3plus = trans;
            elseif trans.to == FIFF.FIFFV_MNE_COORD_FS_TAL_LTZ
                T3minus = trans;
            end
        end
    end
end
fclose(fid);
%
%   Check we have everything we need
%
if (from == FIFF.FIFFV_COORD_HEAD && isempty(T0)) || isempty(T1) || isempty(T2) || ...
        (to == FIFF.FIFFV_MNE_COORD_FS_TAL && (isempty(T3minus) || isempty(T3minus)))
    error(me,'All required coordinate transforms not found');
end
%
%   Go ahead and transform the data
%
if ~isempty(pos)
    if size(pos,2) ~= 3
        error(me,'Coordinates must be given in a N x 3 array');
    end
    if to == from
        trans_pos = pos;
    else
        pos = [ pos ones(size(pos,1),1) ]';
        if from == FIFF.FIFFV_COORD_HEAD
            pos = T0.trans*pos;
        elseif from ~= FIFF.FIFFV_COORD_MRI
            error(me,'Input data must be in MEG head or surface RAS coordinates');
        end
        if to == FIFF.FIFFV_COORD_HEAD
            pos = inv(T0.trans)*pos;
        elseif to ~= FIFF.FIFFV_COORD_MRI
            pos = (T2.trans*T1.trans)*pos;
            if to ~= FIFF.FIFFV_MNE_COORD_MNI_TAL
                if to == FIFF.FIFFV_MNE_COORD_FS_TAL
                    for k = 1:size(pos,2)
                        if pos(3,k) > 0
                            pos(:,k) = T3plus.trans*pos(:,k);
                        else
                            pos(:,k) = T3minus.trans*pos(:,k);
                        end
                    end
                else
                    error(me,'Illegal choice for the output coordinates');
                end
            end
        end
        trans_pos = pos(1:3,:)';
    end
end
end

