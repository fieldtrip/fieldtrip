function dhn_out = read_dhn_med10(filename, password, sortchannel, hdr, begsample, endsample, chanindx)
    % READ_DHN_MED10 read header, event, and waveform data formated in Dark Horse Neuron MED 1.0
    %
    % Syntax:
    %   hdr = read_dhn_med10(filename)
    %   hdr = read_dhn_med10(filename, password)
    %   hdr = read_dhn_med10(filename, password, sortchannel)
    %   evt = read_dhn_med10(filename, password, sortchannel, hdr)
    %   dat = read_dhn_med10(filename, password, sortchannel, hdr, begsample, endsample, chanindx)
    %
    % Input(s):
    %   filename        - [char] name of the file or folder of the dataset
    %   password        - [struct] (opt) password structure of MED 1.0 data (see
    %                     MEDSession_1p0)
    %   sortchannel     - [char] (opt) sort channel order either alphabetically
    %                     'alphabet' or numerically 'number' (default = 'alphabet')
    %   hdr             - [struct] (opt) header structure of the dataset (see FT_READ_HEADER; default = struct([]))
    %   begsample       - [num] (opt) first sample to read (default = [])
    %   endsample       - [num] (opt) last smaple to read (default = [])
    %   chanindx        - [num] (opt) list of channel indices to read (default = [])
    %
    % Output(s):
    %   hdr             - [struct] header structure of the dataset (see FT_READ_HEADER)
    %   evt             - [struct] event structure of the dataset (see FT_READ_EVENT)
    %   dat             - [num] data read in
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also FT_FILETYPE, FT_READ_HEADER, FT_READ_EVENT, FT_READ_DATA.

    % Copyright 2023 Richard J. Cui. Created: Sat 02/11/2023  5:47:28.254 PM
    % $Revision: 0.3 $  $Date: Sun 10/08/2023 01:56:00.468 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        filename (1, :) char
        password (1, 1) struct = struct('Level1Password', '', 'Level2Password', '', ...
            'AccessLevel', 1); % example_data password =='L1_password' or 'L2_password'
        sortchannel (1, 1) logical = false
        hdr (1, 1) struct = struct()
        begsample (1, 1) double = nan
        endsample (1, 1) double = nan
        chanindx (1, :) double = nan
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    % check the consistency of SortChannel
    % ------------------------------------
    if ~isempty(fieldnames(hdr)) && hdr.SortChannel ~= sortchannel
        warning('off', 'backtrace')
        warning('dhn_med10:invalidSortChannel', ...
            'SortChannel provided -%s- is not consistent with -%s- in header. use that in header', ...
            sortchannel, hdr.SortChannel)
        warning('on', 'backtrace')

        sortchannel = hdr.SortChannel;
    end % if

    % setup the instance of MEDFieldTrip_1p0
    % -------------------------------------
    med_ft = MEDFieldTrip_1p0(filename, password, sortchannel); % dealing MED 1.0 data for FieldTrip
    channames = med_ft.SelectedChannel;

    % get the desired information
    % ---------------------------
    switch nargin
        case {1, 2, 3} % get header
            dhn_out = med_ft.getHeader(channames);
            dhn_out.SortChannel = sortchannel;
        case 4 % get event
            dhn_out = med_ft.getEvent(channames);
        case 7 % get data
            dhn_out = med_ft.getData(channames, hdr.sampleunit, begsample, endsample, ...
                chanindx);
        otherwise
            % error
            error('FieldTrip:dhn_med10:invalidInput', ...
            'invalid number of inputs of the function')
    end % switch

end % function read_dhn_med10

% [EOF]
