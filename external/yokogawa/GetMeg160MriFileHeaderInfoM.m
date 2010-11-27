%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    [data_style, model, marker, image_parameter, normalize, besa_fiducial_point] = GetMeg160MriFileHeaderInfoM( fid )
% 
%  arguments:
%    fid             : MRI file ID
% 
%  return values:
%    data_style : string ('DICOM' or 'NOT DICOM')
%    model   : When spherical conductor model is defined, following fields return. If not, empty will return.
%                .cx, .cy, .cz: coordinate of sphere's center [meter] MRIcoordinate
%                .r: radius of sphere [meter] MRIcoordinate
%    marker  : When pick-up of mri marker is done, following fields return. If not, empty will return.
%                .count   : count of picked up marker
%                .pos(8)  : Each structure has following fields.
%                            .done       : Flag for whether picked up or not.
%                            .x, .y, .z  : Coordinate of marker. If done is false, these are empty. [meter]
%    image_parameter : This is mri imaging parameter which has following fields.
%                        .intensity      : 1x2 matrix of min/max intensity of image. ([min max])
%                        .initial_color  : 1x2 matrix of min/max initial color. ([min max])
%                        .color          : 1x2 matrix of min/max color. ([min max])
%    normalize       : When head coordinate is defined, following fields return. If not, empty will return.
%                        .mri_to_normalize   : 4x4 matrix of affine transformation for mri to head coordinate.
%                        .normalize_data(3)  : Normalize parameter includes following fields.
%                                                .name       : Normalize parameter name
%                                                .done       : Flag for whether normalize done or not.
%                                                .x, .y, .z  : Coordinate of normalize parameter. If done is false, these are empty. [meter]
%    besa_fiducial_point(5)  : These structures are besa fiducial points information which has following fields.
%                                .done       : Flag for whether picked up or not.
%                                .x, .y, z   : Coordinate of fiducial point. If done is false, these are empty.
%  
%  confirmation of revision:
%   GetMeg160MriFileHeaderInfoM( Inf ) will show and return revision of this function.
%
