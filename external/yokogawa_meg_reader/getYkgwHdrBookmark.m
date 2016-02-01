% Get header of bookmark
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information about bookmark in the specified file.
%
% usage:
%   bookmark = getYkgwHdrBookmark(filepath)
%
% arguments:
%   filepath            : file path
%
% return values:
%   bookmark            : structure array of bookmark
%     .sample_no        : double, sample number of bookmark
%     .label            : double, label of bookmark
%     .comment          : double, comment of bookmark
%
% rivision history
%   2 : 2011.03.03 : Structure fields (reference_no, type) which had not been used were removed.
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
