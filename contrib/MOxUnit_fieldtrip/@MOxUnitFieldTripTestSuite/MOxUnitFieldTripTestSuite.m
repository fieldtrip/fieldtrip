function obj=MOxUnitFieldTripTestSuite(varargin)
% Initialize empty FieldTrip test suite
%
% Input:
%   name            Optional name of the test suite
% Output:
%   obj             MOxUnitFieldTripTestSuite instance with no tests.
%
% Notes:
%   MOxUnitFieldTripTestSuite is a subclass of MOxUnitTestSuite,
%   which is is a subclass of MoxUnitTestNode.
%
% See also: MoxUnitTestNode, MOxUnitTestSuite
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    class_name='MOxUnitFieldTripTestSuite';
    obj=class(struct(),class_name,MOxUnitTestSuite(varargin{:}));

