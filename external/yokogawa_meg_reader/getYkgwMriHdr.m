% Get header of mri
% [ Yokogawa MEG Reader toolbox for MATLAB ]
% 
% brief:
%   This function retrieves header information of specified mri file (*.mri).
%
% usage:
%   mri_header = getYkgwMriHdr(filepath)
%
% arguments:
%   filepath                  : file path of *.mri file
%
% return values:
%   mri_header                : structure of mri header information
%     .data_style             : double, data style (0 : DICOM, others : Polhemus)
%     .model                  : structure of conductor model
%       .done                 :   bool, is conductor model defined ? ( true -> defined )
%       .type                 : double, type of conductor model as follows:
%                                UNKNOWN_MODEL   = -1;
%                                NO_MODEL        =  0;
%                                SPHERICAL_MODEL =  1;
%                                LAYERED_MODEL   =  2;
%     [when model_type is SPHERICAL_MODEL]
%       .cx                   : double, x coordinate [meter] of spherical center position on MRI coordinate
%       .cy                   : double, y coordinate [meter] of spherical center position on MRI coordinate
%       .cz                   : double, z coordinate [meter] of spherical center position on MRI coordinate
%       .radius               : double,       radius [meter] of spherical conductor on MRI coordinate
%     [when model_type is LAYERED_MODEL]
%       .ax                   : double, coefficient 'ax' of planar equation 'ax * x + ay * y + az * z = c'
%       .ay                   : double, coefficient 'ay' of planar equation 'ax * x + ay * y + az * z = c'
%       .az                   : double, coefficient 'az' of planar equation 'ax * x + ay * y + az * z = c'
%       .c                    : double, coefficient 'c'  of planar equation 'ax * x + ay * y + az * z = c'
%
%     .hpi                    : The structure of point data about picked HPI.
%       .done                 : boolean, Is pick-up of a HPI point done ? (true : done)
%       .mri_pos              : double, HPI position [x, y, z] on MRI coordinate [meter]
%       .label                : string, HPI label as follows:
%                                'LPA'   : Left   PreAuricular
%                                'RPA'   : Right  PreAuricular
%                                'CPF'   : Center PreFrontal
%                                'LPF'   : Left   PreFrontal
%                                'RPF'   : Right  PreFrontal
%
%     .image_parameter        : structure of image parameters
%       .intensity            : 1 x 2 row vector(double), minimum and maximum of image values
%       .initial_color        : 1 x 2 row vector(double), minimum and maximum of initial brightness
%       .color                : 1 x 2 row vector(double), minimum and maximum of current brightness
%
%     .normalize              : structure of HEAD coordinate system ( LPA(x-), RPA(x+), nasion(y+) )
%       .done                 :   bool, is HEAD coordinate system defined ? ( true -> defined )
%       .mri2normalize        : 4 x 4 matrix(double), matrix transforming MRI coordinate to HEAD coordinate [meter]
%                               usage: [xhead, yhead, zhead, 1]' = mri_header.normalize.mri2normalize * [xmri, ymri, zmri, 1]'
%       .point                : structure array of HEAD fiducial points
%         .done               :   bool, is pick-up of a HEAD fiducial point done ? ( true -> done )
%         .name               : string, name of HEAD fiducial points
%         .x                  : double, x coordinate [meter] of a HEAD fiducial point on MRI coordinate
%         .y                  : double, y coordinate [meter] of a HEAD fiducial point on MRI coordinate
%         .z                  : double, z coordinate [meter] of a HEAD fiducial point on MRI coordinate
%
%     .besa_fiducial          : structure of BESA fiducial information
%       .point                : structure array of BESA fiducial points
%         .done               :   bool, is pick-up of a BESA fiducial point done ? ( true -> done )
%         .x                  : double, x coordinate [meter] of a BESA fiducial point on MRI coordinate
%         .y                  : double, y coordinate [meter] of a BESA fiducial point on MRI coordinate
%         .z                  : double, z coordinate [meter] of a BESA fiducial point on MRI coordinate
%
% rivision history
%   2 : 2011.04.27 : add HPI(marker) information
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
