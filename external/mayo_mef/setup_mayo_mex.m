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

    % Copyright 2020-2023 Richard J. Cui. Created: Fri 05/15/2020 10:33:00.474 AM
    % $ Revision: 0.3 $  $ Date: Sun 10/08/2023 12:42:06.711 PM $
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
        options.DHNRootPath (1, :) char = '' % full path to DHN root directory
        options.ForceBuildMex (1, 1) logical = false
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    dhn_root = options.DHNRootPath; % if empty, use default root directory
    force_build_mex = options.ForceBuildMex;

    % * set default DHN root directory
    if isempty(dhn_root)

        switch computer
            case 'MACI64' % Mac
                user = getenv('USER');
                dhn_root = fullfile(filesep, 'Users', user, 'DHN');
            case 'GLNXA64' % Linux
                user = getenv('USER');
                dhn_root = fullfile(filesep, 'home', user, 'DHN');
            case 'PCWIN64' % Windows
                driver = getenv('HOMEDRIVE');
                user = getenv('HOMEPATH');
                dhn_root = fullfile(driver, user, 'DHN');
            otherwise
                ft_error('MAYO_MEF:setup_mayo_mex', ...
                    'Unknown computer type %s. MED is not supported.\n', ...
                    computer)
        end % switch

    end % if

    if isfolder(dhn_root)
        % add DHN root directory to MATLAB path
        addpath(genpath(dhn_root))
        med_mex_path = fullfile(dhn_root, 'read_MED', 'Resources');
    else
        ft_warning('MAYO_MEF:setup_mayo_mex', ...
            'DHN root directory %s does not exist. please install read_MED package (http://darkhorseneuro.com) or manually set DHN root directory\n', dhn_root)
        med_mex_path = string(1, 0);
    end % if

    % install mayo_mef package
    % ------------------------
    ft_hastoolbox('mayo_mef', 1);

    % get current directory
    % ---------------------
    cur_dir = pwd;

    % check MEF mex binary
    % --------------------
    fprintf('***************************\n')
    fprintf('* Checking MEF mex binary *\n')
    fprintf('***************************\n')

    % directory of setup_mayo_mex.m assumed in mayo_mef
    mayo_mef = fileparts(mfilename('fullpath')); % store mef mex here
    mef_mex_path = fullfile(mayo_mef, 'mex_mef');

    if force_build_mex
        cd(mef_mex_path)
        make_mex_mef
    else % check mex files in mayo_mef
        valid_mex = check_mex_files(mayo_mef, MexType = "MEF");

        if valid_mex == false
            cd(mef_mex_path)
            make_mex_mef
        end % if

    end % if

    % check MED mex binary
    % --------------------
    if ~isempty(med_mex_path)
        fprintf('\n')
        fprintf('***************************\n')
        fprintf('* Checking MED mex binary *\n')
        fprintf('***************************\n')

        % check mex files in read_MED
        valid_mex = check_mex_files(med_mex_path, MexType = 'MED');

        if valid_mex == false
            ft_warning('MAYO_MEF:setup_mayo_mex', ...
            'None or not all MED mex files can be found. Please install read_MED package (http://darkhorseneuro.com) or manually set DHN root directory\n')
        end % if

    end % if

    % return to original directory
    % ----------------------------
    cd(cur_dir)

end

% ==========================================================================
% subroutines
% ==========================================================================
function valid_mex = check_mex_files(mex_path, options)
    % SETUP_MAYO_MEX.CHECK_MEX_FILES check if mex files are valid

    arguments
        mex_path (1, 1) string % full path to mex files
    end % positional

    arguments
        options.MexType (1, 1) string {mustBeMember(options.MexType, ["MEF", "MED"])} ...
            = "MEF" % MEF or MED
    end % optional

    mex_type = options.MexType;

    switch mex_type
        case 'MEF'
            valid_mex_1 = isfile(fullfile(mex_path, "read_mef_header_2p1" + "."+mexext));
            valid_mex_2 = isfile(fullfile(mex_path, "decompress_mef_2p1" + "."+mexext));
            valid_mex_3 = isfile(fullfile(mex_path, "read_mef_session_metadata" + "."+mexext));
            valid_mex_4 = isfile(fullfile(mex_path, "read_mef_ts_data" + "."+mexext));

            valid_mex = valid_mex_1 && valid_mex_2 && valid_mex_3 && valid_mex_4;

            if valid_mex
                fprintf('MEF mex files are found.\n')
            else
                ft_warning('None or not all MEF mex files can be found.\n')
            end % if

        case 'MED'
            valid_mex_1 = isfile(fullfile(mex_path, "load_session" + "."+mexext));
            valid_mex_2 = isfile(fullfile(mex_path, "matrix_MED_exec" + "." + mexext));
            valid_mex_3 = isfile(fullfile(mex_path, "read_MED_exec" + "." + mexext));

            valid_mex = valid_mex_1 && valid_mex_2 && valid_mex_3;

            if valid_mex
                fprintf('MED mex files are found.\n')
            else
                ft_warning('None or not all MED mex files can be found.\n')
            end % if

    end % switch

end % function

% [EOF]
