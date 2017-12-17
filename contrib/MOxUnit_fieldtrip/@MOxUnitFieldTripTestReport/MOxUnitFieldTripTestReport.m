function obj=MOxUnitFieldTripTestReport(verbosity,varargin)
% Initialize empty MOxUnitFieldTripTestReport instance
%
% obj=MOxUnitFieldTripTestReport([verbosity]
%
% Inputs:
%   verbosity       Integer indicating how verbose the output is when
%                   tests are run using this instance (default: 2).
%                   Use 1 for less output, or 0 for no output.
%
% Returns:
%   obj             Empty MOxUnitTestReport instance, with no test errors,
%                   failures, skips, or successes stored.
%
% Notes:
%   - MOxUnitFieldTripTestReport is a subclass of MOxUnitTestReport. It
%     adds a function submit_to_dashboard that sends test reports to the
%     FieldTrip dashboard.
%


    if nargin<1
        verbosity=2; % be verbose by default
    end

    parent=MOxUnitTestReport(verbosity,varargin{:});

    s=struct();
    obj=class(s,'MOxUnitFieldTripTestReport',parent);