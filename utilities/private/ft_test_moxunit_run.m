function passed = ft_test_moxunit_run(unused,varargin)

% FT_TEST_MOXUNIT_RUN

% Copyright (C) 2017, Nikolaas N. Oosterhof
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

    % ensure path is set for MOxUnit fieldtrip test functions, but
    % set to original state after running this function


    orig_path=path();
    path_resetter=onCleanup(@()path(orig_path));
    ensure_moxunit_fieldtrip_path_is_set();
    check_dependencies();

    % By default, running FieldTrip excludes files if their name
    % starts with 'failed'. Here this default behaviour is mimicked.
    override_default_arg={'exclude_if_prefix_equals_failed',true};
    arg=cat(2,override_default_arg,varargin);

    passed=moxunit_fieldtrip_runtests(arg{:});



function check_dependencies()
% throw an error if MOxUnit is not available

    if isempty(which('moxunit_runtests'))
        error(['MOxUnit is required; see '...
                    'https://github.com/moxunit/moxunit']);
    end

function fieldtrip_root_dir=get_fieldtrip_root_dir()
    fieldtrip_root_dir=fileparts(which('ft_defaults'));
    if isempty(fieldtrip_root_dir)
        error('Fieldtrip path is not set');
    end



function ensure_moxunit_fieldtrip_path_is_set()
    if isempty(which('MOxUnitFieldTripTestSuite'))
        moxunit_fieldtrip_dir=fullfile(get_fieldtrip_root_dir(),...
                                       'contrib','MOxUnit_fieldtrip');
        addpath(moxunit_fieldtrip_dir);
    end
