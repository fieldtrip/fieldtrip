% Get header of the system information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information about coregistration in the specified file.
%
% usage:
%   coregist = getYkgwHdrCoregist(filepath)
%
% arguments:
%   filepath       : file path
%
% return values:
%   coregist       : structure of coregistration
%     .done        :   bool, Is coregistration done ? (true : done)
%     .mri_type    : string, type of MRI as follows:
%                      NoMriFile       = 0;
%                      NormalMriFile   = 1;
%                      VirtualMriFile  = 2;
%     .mri_file    : string, file path of MRI file (*.mri)
%     .hpi_file    : string, file path of HPI file (*.mrk)
%     .meg2mri     : 4 x 4 matrix(double), matrix transforming MEG device coordinate to MRI coordinate [meter]
%                    usage: [xmri, ymri, zmri, 1]' = coregist.meg2mri * [xmeg, ymeg, zmeg, 1]'
%     .mri2meg     : 4 x 4 matrix(double), matrix transforming MRI coordinate to MEG device coordinate [meter]
%                    usage: [xmeg, ymeg, zmeg, 1]' = coregist.meg2mri * [xmri, ymri, zmri, 1]'
%     .hpi         : structure of HPI
%       .meg_pos   : 1 x 3 matrix(double), HPI position on MEG device coordinate [meter]
%       .mri_pos   : 1 x 3 matrix(double), HPI position on MRI coordinate [meter]
%       .label     : string, HPI label as follows:
%                    'LPA'  Left   PreAuricular
%                    'RPA'  Right  PreAuricular
%                    'CPF'  Center PreFrontal
%                    'LPF'  Left   PreFrontal
%                    'RPF'  Right  PreFrontal
%     .model       : structure of conductor model
%       .type      : double, type of conductor model as follows:
%                                UNKNOWN_MODEL   = -1;
%                                NO_MODEL        =  0;
%                                SPHERICAL_MODEL =  1;
%                                LAYERED_MODEL   =  2;
%     [when model.type is SPHERICAL_MODEL]
%       .cx        : double, x coordinate [meter] of spherical center position on MRI coordinate
%       .cy        : double, y coordinate [meter] of spherical center position on MRI coordinate
%       .cz        : double, z coordinate [meter] of spherical center position on MRI coordinate
%       .radius    : double,       radius [meter] of spherical conductor on MRI coordinate
%     [when model.type is LAYERED_MODEL]
%       .ax        : double, coefficient 'ax' of planar equation 'ax * x + ay * y + az * z = c'
%       .ay        : double, coefficient 'ay' of planar equation 'ax * x + ay * y + az * z = c'
%       .az        : double, coefficient 'az' of planar equation 'ax * x + ay * y + az * z = c'
%       .c         : double, coefficient 'c'  of planar equation 'ax * x + ay * y + az * z = c'
%
% rivision history
%   1 : 2011.02.14 : add HPI(marker) information
%                    1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
