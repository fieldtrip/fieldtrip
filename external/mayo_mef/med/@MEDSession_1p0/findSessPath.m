function [sesspath, channames] = findSessPath(this, filename)
    % MEDSESSION_1P0.FINDSESSPATH find the path of the session and channel names
    %
    % Syntax:
    %  [sesspath, channames] = findSessPath(this, filename)
    %
    % Input(s):
    %   this 		- [obj] MEDsession_1p0 object
    %   filename 	- [char] filename of the session
    %
    % Output(s):
    %   sesspath        - [char] session path
    %   channames       - [string] channel names included
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 02/14/2023  9:43:31.226 PM
    % $Revision: 0.2 $  $Date: Sun 02/26/2023 11:09:45.566 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDSession_1p0
        filename (1, :) char = ''
    end % positional

    % =========================================================================
    % main
    % =========================================================================
    if isempty(filename)
        sesspath = '';
        channames = "";
        return
    end % if

    % find session path and channel name
    % ----------------------------------
    [~, ~, f_ext] = fileparts(filename);

    switch f_ext
        case '.medd' % session folder

            if isfolder(filename)
                sesspath = filename; % session path
                listing = dir(fullfile(filename, '*.ticd')); % channel names

                if isempty(listing)
                    channames = ""; % no channel
                else
                    [~, channames] = fileparts(string({listing.name}));
                end % if

            else
                error('MEDsession_1p0:findSessPath:InvalidInput', ...
                'Session folder does not exist. Can only process MED 1.0 files.');
            end % if

        case '.ticd' % channel folder

            if isfolder(filename)
                chanpath = fileparts(filename); % channel path
                [sesspath, channames] = fileparts(chanpath); % session path and channel names
            else
                error('MEDsession_1p0:findSessPath:InvalidInput', ...
                'Channel folder does not exist. Can only process MED 1.0 files.');
            end % if

        case '.recd' % record folder

            if isfolder(filename)
                recpath = fileparts(filename); % record path
                sesspath = fileparts(recpath); % session path
                listing = dir(fullfile(filename, '*.ticd')); % channel names

                if isempty(listing)
                    channames = ""; % no channel
                else
                    channames = string({listing.name});
                end % if

            else
                error('MEDsession_1p0:findSessPath:InvalidInput', ...
                'Record folder does not exist. Can only process MED 1.0 files.');
            end % if

        otherwise
            error('MEDsession_1p0:findSessPath:InvalidInput', ...
            'Invalid file type. Can only process MED 1.0 files.')

    end % switch

end % function findSessPath

% [EOF]
