function obj=addFromFile(obj,fn)
% Add FieldTrip test from file
%
% obj=addFromFile(obj,fn)
%
% Inputs:
%   obj             MoxUnitFieldTripTestSuite instance.
%   fn              name of file that contains a FieldTrip test function
%                   (typically in FieldTrip's 'test' directory)
%
% Output:
%   obj             MoxUnitFieldTripTestSuite instance with the test
%                   added
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    test_case=MOxUnitFieldTripTestCase(fn);
    obj=addTest(obj,test_case);