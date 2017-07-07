function failed_old_besa2fieldtrip

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_besa2fieldtrip

% Test the besa2fieldtrip function, which reads an object from a BESA file
% and returns it as FieldTrip object (e.g. timelocked, freq, source).
%
% The besa2fieldtrip function both should work on original besa files,
% as well as on MATLAB files containing a particular structure (for the
% direct besa->matlab interface).

basedir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/besa/file');

filename = {
  'besa_avr/Rolandic-Segment.avr'
  'besa_avr/Rolandic-Segment.elp'
  'besa_avr/Rolandic-Spike-Child-export.avr'
  'besa_avr/Rolandic-Spike-Child-export.elp'
  'besa_image/LAURA.dat'
  'besa_image/MinimumNorm.dat'
  'besa_image/MSBF.dat'
  'besa_image/Rolandic_CLARA.dat'
  'besa_image/Rolandic_LAURA.dat'
  'besa_mul/Rolandic-Spike-Child-export.elp'
  'besa_mul/Rolandic-Spike-Child-export.mul'
  'besa_mul/Rolandic_Segment.elp'
  'besa_mul/Rolandic_Segment.mul'
  'besa_sb/EpochedData.atf'
  'besa_sb/EpochedData.dat'
  'besa_sb/EpochedData.elp'
  'besa_sb/EpochedData.generic'
  'besa_sb/Rolandic-Segment.dat'
  'besa_sb/Rolandic-Segment.elp'
  'besa_sb/Rolandic-Segment.generic'
  'besa_sb/Rolandic-Spike-Child-export.dat'
  'besa_sb/Rolandic-Spike-Child-export.elp'
  'besa_sb/Rolandic-Spike-Child-export.evt'
  'besa_sb/Rolandic-Spike-Child-export.generic'
  'besa_sourcewaveforms/Spikes-Child1.swf'
  'besa_sourcewaveforms/SWF_Rolandic.swf'
  'besa_tfc/AC_Osc20_CO2_F9_rfr.tfc'
  'besa_tfc/AC_Osc20_ERA.tfc'
  'besa_tfc/AC_Osc20_UsrMtgRC0.tfc'
  'besa_tfc/AC_Osc20_VirtualRef.tfc'
};

success = false(size(filename));
for i=1:length(filename)
  fname = fullfile(basedir, filename{i});
  try
    ft_struct = besa2fieldtrip(fname);
    % FIXME a call to checkdata should be inserted here
    success(i) = true;
  catch
    success(i) = false;
  end
end

if ~all(success)
  error('problem reading file %s\n',  filename{~success});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = dccnpath('/home/common/matlab/fieldtrip/testdata/original/besa/struct');

filename = {
  'besa_channels/EpochedData.mat'
  'besa_channels/Rolandic-Spike-Child-export.mat'
  'besa_channels/Rolandic_Segment.mat'
  'besa_image/LAURA.mat'
  'besa_image/MinimumNorm.mat'
  'besa_image/MSBF.mat'
  'besa_image/Rolandic_LAURA.mat'
  'besa_image/Rolandic_LAURAandCLARA.mat'
  'besa_sourcewaveforms/Spikes-Child1.mat'
  'besa_sourcewaveforms/SWF_Rolandic.mat'
  'besa_tfc/3Items_UsrMtgRC0.mat'
  'besa_tfc/3Items_VirtRef.mat'
  'besa_tfc/AC_Osc20_CO2_F9_rfr.mat'
  'besa_tfc/AC_OSC20_ERA.mat'
};

success = false(size(filename));
for i=1:length(filename)
  fname = fullfile(basedir, filename{i});
  try
    % loat the mat file
    content = load(fname);
    % determine the content of the mat file
    fn = fieldnames(content);
    % get the besa struct and convert to FieldTrip struct
    besa_struct = getfield(content, fn{i});
    ft_struct = besa2fieldtrip(besa_struct);
    % FIXME a call to checkdata should be inserted here
    success(i) = true;
  catch
    success(i) = false;
  end
end

if ~all(success)
  error('problem reading file %s\n',  filename{~success});
end

