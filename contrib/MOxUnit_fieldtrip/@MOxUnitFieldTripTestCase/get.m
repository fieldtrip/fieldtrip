function value=get(obj, key)
% Get test property
%
% Inputs:
%   obj                 MOxUnitFieldTripTestCase instance
%   key                 test property key, one of:
%                       'mem'       Maximum memory the test uses (bytes)
%                       'walltime'  Maximum time the test takes (seconds)
%                       'dccnpath'  Set to true if this function requires
%                                   dccnpath
%
% Output:
%   value               The value associated with key, or empty if the
%                       property is not set for this test case
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    properties=obj.properties;

    if ~isfield(properties,key)
        value=[];
        return;
    end

    value=properties.(key);

