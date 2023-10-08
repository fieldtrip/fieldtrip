function dat = getData(this, channames, sampleunit, begsample, endsample, chanindx)
    % MEDFILEDTRIP_1P0.GETDATA read data from MED 1.0 dataset for FieldTrip
    %
    % Syntax:
    %   dat = getData(this)
    %   dat = getData(__, channames, sampleunit, begsample, endsample, chanindx)
    %
    % Input(s):
    %   this            - [obj] MEDFieldTrip_1p0 object
    %   channames       - [string] (opt) channel names (default = all channels)
    %   sampleunit      - [char] (opt) unit of begsample and endsample (default
    %                     = index)
    %   begsample       - [num] (opt) time point of begin sample (default = 1st one)
    %   endsample       - [num] (opt) time point of end sample (default = last one)
    %   chanindx        - [num] (opt) index of selection from channames
    %
    % Output(s):
    %   dat             - [num] data array channel x samples
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 10/03/2023 12:25:01.918 AM
    % $Revision: 0.2 $  $Date: Sun 10/08/2023 03:04:54.209 PM $
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
        channames (1, :) string = this.ChannelName
        sampleunit (1, :) char = 'index'
        begsample (1, 1) double = 1
        endsample (1, 1) double = this.Samples
        chanindx (1, :) double = []
    end % positional

    if isempty(chanindx)
        chanindx = 1:numel(channames);
    end % if

    % ======================================================================
    % main
    % ======================================================================
    sel_chan = channames(chanindx);
    dat = this.importSession([begsample, endsample], sampleunit, this.SessionPath, ...
        'SelectedChannel', sel_chan);

end % function getData

% [EOF]
