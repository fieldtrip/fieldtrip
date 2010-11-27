%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    source_info = GetMeg160SourceInfoM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    source_info        : (N) structures of source information. (N: number of sample point when ECD is estimated)
%                        .sample_no          : sample no when source is estimated
%                        .gof                : GOF
%                        .correlation        : correlation
%                        .dipole_count       : dipole count at the sample_no
%                        .current_dipole     : (dipole count) structures of estimated dipole infomation
%                                .x, .y, .z      : dipole position [m] in MEG coordinate
%                                .zdir, xdir     : dipole direction [deg]
%                                .intensity      : dipole intensity [Am]
%  
%  confirmation of revision:
%   GetMeg160SourceInfoM( Inf ) will show and return revision of this function.
%
