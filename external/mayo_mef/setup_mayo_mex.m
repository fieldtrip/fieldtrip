function setup_mayo_mex(options)
    % SETUP_MAYO_MEX make mex binary necessary for reading MEF dataset
    %
    % Syntax:
    %   setup_mayo_mex()
    %
    % Input(s):
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also make_mex_mef, test_mayo_mef.

    % Copyright 2020 Richard J. Cui. Created: Fri 05/15/2020 10:33:00.474 AM
    % $ Revision: 0.2 $  $ Date: Thu 01/26/2023 10:35:01.233 PM $
    %
    % Mayo Foundation for Medical Education and Research
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905, USA
    %
    % Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
    end % positional

    arguments
        options.ForceBuildMex (1, 1) logical = false
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    force_build_mex = options.ForceBuildMex;

    % get current directory
    % ---------------------
    cur_dir = pwd;

    % check mex binary
    % ----------------
    % directory of setup_mayo_mex.m assumed in mayo_mef
    mayo_mef = fileparts(mfilename('fullpath'));

    if force_build_mex
        cd([mayo_mef, filesep, 'mex_mef'])
        make_mex_mef
    else % check mex files in mayo_mef
        valid_mex = check_mex_files(mayo_mef);

        if valid_mex == false
            cd([mayo_mef, filesep, 'mex_mef'])
            make_mex_mef
        end % if

    end % if

    % return to original directory
    % ----------------------------
    cd(cur_dir)

end

% ==========================================================================
% subroutines
% ==========================================================================
function valid_mex = check_mex_files(mex_path)
    % SETUP_MAYO_MEX.CHECK_MEX_FILES check if mex files are valid

    arguments
	mex_path (1, 1) string % full path to mex files
    end % positional

end % function

% [EOF]
