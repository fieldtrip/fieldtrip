function hdr = getHeader(this, channames)
    % MEDFIELDTRIP_1P0.GETHEADER get header information of MED 1.0 session for fieldtrip
    %
    % Syntax:
    %   hdr = getHeader(this)
    %   hdr = getHeader(this, channames)
    %
    % Input(s):
    %   this            - [obj] MEDFieldTrip_1p0 object
    %   channames       - [string] channel names
    %
    % Output(s):
    %   hdr             - [struct] structure of header information (from
    %                     session metadata information)
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Mon 10/02/2023 12:43:05.396 AM
    % $Revision: 0.1 $  $Date: Wed 10/04/2023 01:14:37.925 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDFieldTrip_1p0
        channames (1, :) string = []
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    % check channel names if provided
    % -------------------------------
    if ~isempty(channames)
        sess_chan = this.ChannelName; % all channels in the session

        if all(ismember(channames, sess_chan)) == false
            error('MEFFieldTrip_3p0:getHeader:invalidChannel', ...
            'invalid channel names')
        end % if

    end % if

    % construct header structure
    % --------------------------
    metadata = this.MetaData; % session metadata

    % sampling frequency
    hdr.Fs = this.SamplingFrequency;

    % number of channels
    hdr.nChans = length(this.ChannelName);

    % number of samples per trial
    hdr.nSamples = metadata.absolute_end_sample_number;

    % number of trials
    hdr.nTrials = 1;

    % label of channel
    l = convertStringsToChars(this.ChannelName);
    hdr.label = l(:);

    % channel type
    num_ch = hdr.nChans;
    ct = cell(num_ch, 1);
    ct(:) = {'unknown'};
    hdr.chantype = ct;

    % channel unit
    unit_des = metadata.amplitude_units_description;
    cu = cell(num_ch, 1);
    cu(:) = {unit_des};
    hdr.chanunit = cu;

    % number of pre-trigger samples in each trial
    hdr.nSamplesPre = this.Samples;

    % save the metadata in the header
    hdr.metadata = metadata;

end % function getHeader

% [EOF]
