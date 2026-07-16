function asc = read_eyelink_asc(filename, asc, begsample, endsample, chanindx)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file. The output events are
% represented as matlab tables (after Aug 2022)
%
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010-2015, Robert Oostenveld
% Copyright (C) 2022-2026, Jan-Mathijs Schoffelen
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

needhdr = (nargin==1); % always
needevt = (nargin==1); % always
needdat = (nargin==5);

if needhdr || needevt
  % skip this
  % read the whole file at once
  fid = fopen_or_error(filename, 'rt');
  aline = fread(fid, inf, 'char=>char');
  fclose(fid);

  aline = splitlines(string(aline(:)'));
  sel   = startsWith(aline, string(0:9));

  datline = aline(sel);
  aline   = aline(~sel);

  % according to the eyelink documentation, the sample lines have different
  % flavours, depending on how the data was collected (monocular or
  % binocular), and possibly the conversion settings from edf2asc
  % Monocular,                              <time> <xp>  <yp>  <ps>
  % Monocular, with velocity                <time> <xp>  <yp>  <ps>  <xv>  <yv>
  % Monocular, with resolution              <time> <xp>  <yp>  <ps>  <xr>  <yr>
  % Monocular, with velocity and resolution <time> <xp>  <yp>  <ps>  <xv>  <yv>  <xr>  <yr>
  % Binocular,                              <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr>
  % Binocular, with velocity                <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr>
  % Binocular, with resolution              <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xr> <yr>
  % Binocular, with velocity and resolution <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr><xr> <yr>
  %
  % Then, there may be additional columns, if recorded with cornea reflection
  % mode:
  %
  % MONOCULAR Corneal Reflection (CR) Samples
  % "..." if no warning for sample
  % first character is "I" if sample was interpolated
  % second character is "C" if CR missing
  % third character is "R" if CR recovery in progress
  %
  % BINOCULAR Corneal Reflection (CR) Samples?
  % "....." if no warning for sample
  % first character is "I" if sample was interpolated
  % second character is "C" if LEFT CR missing
  % third character is "R" if LEFT CR recovery in progress
  % fourth character is "C" if RIGHT CR missing
  % fifth character is "R" if RIGHT CR recovery in progress
  %
  % Then, there may be additional columns, if data collection was done in
  % remote mode:
  %
  % Data files recorded using the Remote Mode have extra columns to encode the
  % target distance, position, and eye/target status information. The first three
  % columns are:
  % <target x>: X position of the target in camera coordinate (a value from 0 to 10000).
  % Returns "MISSING_DATA" (-32768) if target is missing.
  % <target y>: Y position of the target in camera coordinate (a value from 0 to 10000).
  % Returns "MISSING_DATA" (-32768) if target is missing.
  % <target distance>: Distance between the target and camera (in millimeters).
  % Returns "MISSING_DATA" (-32768) if target is missing.
  % The next thirteen fields represent warning messages for that sample relating to
  % the target and eye image processing."............." if no warning for target and eye image
  % first character is "M" if target is missing
  % second character is "A" if extreme target angle occurs
  % third character is "N" if target is near eye so that the target window and eye window overlap
  % fourth character is "C" if target is too close
  % fifth character is "F" if target is too far
  % sixth character is "T" if target is near top edge of the camera image
  % seventh character is "B" if target is near bottom edge of the camera image
  % eighth character is "L" if target is near left edge of the camera image
  % ninth character is "R" if target is near right edge of the camera image
  % tenth character is "T" if eye is near top edge of the camera image
  % eleventh character is "B" if eye is near bottom edge of the camera image
  % twelfth character is "L" if eye is near left edge of the camera image
  % thirteenth character is "R" if eye is near right edge of the camera image
  % For a binocular recording, there will be seventeen target/eye status columns,
  % with the last eight columns reporting the warning messages for the left and
  % right eyes separately.

  % So, anecdotally, parsing the file is a bit tedious, and not very robust.
  % Reading and converting per line is super slow, so here I chose to pipe
  % the datalines through a checker that removes all '...' '.....'
  % '.............' and '.................' (or containing warning) stuff
  % so that we end up with more or less well behaved data.


  % check if the formatting for all datlines is similar
  ntab   = count(datline, sprintf('\t'));

  % check whether all lines have the same number of columns
  assert(all(ntab==ntab(1)));
  tstamps = double(extractBefore(datline, sprintf('\t')));
  udelta  = unique(diff(tstamps));
  if udelta(1)==0 && udelta(2)==1
    % anecdotally, the edf2asc has converted a 2kHz file, but rounded the
    % time stamps...
    ft_warning('it looks like a 2 kHz file with rounded off time stamps, this may lead to a 0.5 ms inaccuracy in the timing of the events');
    Fs1 = 2000;

    ft_info('updating timestamps');
    sel = [false;diff(tstamps)==0];
    tstamps(sel) = tstamps(sel)+0.5;
  elseif udelta(1)>=1
    Fs1 = 1000/udelta(1); % this is a guess
  end

  % create a mapping from the time stamps to samples in the data
  tstamps_samples = (Fs1./1000).*(tstamps-tstamps(1)) + 1;
  samples         = nan(1, max(tstamps_samples));
  samples(tstamps_samples) = 1:numel(tstamps_samples);

  % create a mapping from samples in the data to discontinuous trials
  trialidx = cumsum([1;double(diff(tstamps_samples)>1)]);

  asc = [];
  asc.filename = filename;
  asc.datline  = datline; % don't parse it for now
  asc.tstamps  = tstamps;
  asc.samples  = samples(:);
  asc.trialidx = trialidx;

  % make tables for the events sfix, ssacc, sblink etc
  selsfix   = startsWith(aline, 'SFIX');
  selssacc  = startsWith(aline, 'SSACC');
  selsblink = startsWith(aline, 'SBLINK');
  selection = {selsfix selssacc selsblink};
  fn        = {'sfix' 'ssacc' 'sblink'};
  for i  = find(cellfun(@sum, selection))
    s        = extractAfter(aline(selection{i}), " ");
    s        = split(s, " ");
    stime    = double(s(:,end)); % in ms, possibly rounded off (if Fs=2000)
    sample   = timestamp2samples(stime, tstamps, samples, Fs1);
    eye      = s(:,1);
    asc.(fn{i}) = table(eye, stime, sample);
  end

  selefix   = startsWith(aline, 'EFIX');
  if sum(selefix)
    s   = extractAfter(aline(selefix), " ");
    eye = extractBefore(s, " ");
    s   = double(split(extractAfter(s, " "), sprintf('\t')));
    if size(s,2)==6
      varnames = {'stime', 'etime', 'dur', 'axp', 'ayp', 'aps'};
    elseif size(s,2)==8
      varnames = {'stime', 'etime', 'dur', 'axp', 'ayp', 'aps', 'xr', 'yr'};
    else
      ft_error('unknown number of columns in efix');
    end
    s_smp    = timestamp2samples(s(:,1:3), tstamps, samples, Fs1);
    asc.efix = cat(2, table(eye), array2table(s, 'VariableNames', varnames), ...
      array2table(s_smp, 'VariableNames', {'sample' 'sample_end' 'sample_duration'}));
  end

  seleblink  = startsWith(aline, 'EBLINK');
  if sum(seleblink)
    s   = extractAfter(aline(seleblink), " ");
    eye = extractBefore(s, " ");
    s   = double(split(extractAfter(s, " "), sprintf('\t')));
    s_smp    = timestamp2samples(s(:,1:3), tstamps, samples, Fs1);
    asc.eblink = cat(2, table(eye), array2table(s, 'VariableNames', {'stime', 'etime', 'dur'}), ...
      array2table(s_smp, 'VariableNames', {'sample' 'sample_end' 'sample_duration'}));
  end

  selesacc   = startsWith(aline, 'ESACC');
  if sum(selesacc)
    s   = extractAfter(aline(selesacc), " ");
    eye = extractBefore(s, " ");
    s   = double(split(extractAfter(s, " "), sprintf('\t')));
    if size(s,2)==9
      varnames = {'stime', 'etime', 'dur', 'sxp', 'syp', 'exp', 'eyp', 'ampl', 'pv'};
    elseif size(s,2)==11
      varnames = {'stime', 'etime', 'dur', 'sxp', 'syp', 'exp', 'eyp', 'ampl', 'pv', 'xr', 'yr'};
    else
      ft_error('unknown number of columns in esacc');
    end
    s_smp    = timestamp2samples(s(:,1:3), tstamps, samples, Fs1);
    asc.esacc = cat(2, table(eye), array2table(s, 'VariableNames', varnames), ...
      array2table(s_smp, 'VariableNames', {'sample' 'sample_end' 'sample_duration'}));
  end

  selinput   = startsWith(aline, 'INPUT');
  if sum(selinput)
    s      = extractAfter(aline(selinput), sprintf('\t'));
    s      = double(split(s, sprintf('\t')));
    sample = timestamp2samples(s(:,1), tstamps, samples, Fs1);
    asc.input = array2table([s sample],'VariableNames', {'stime' 'value' 'sample'});
  end

  selhdr     = startsWith(aline, '**');
  if sum(selhdr)
    asc.header = char(aline(selhdr));
  end

  selmsg     = startsWith(aline, 'MSG');
  if sum(selmsg)
    s = extractAfter(aline(selmsg), sprintf('\t'));
    stime   = double(extractBefore(s, " "));
    sample  = timestamp2samples(stime, tstamps, samples, Fs1);
    message = extractAfter(s, " ");
    asc.msg = cat(2, array2table(stime), table(message), array2table(sample));
  end
end

if needhdr
  hdr.nChans              = ntab(1)+1;
  hdr.nSamples            = numel(tstamps);
  hdr.nSamplesPre         = 0;
  hdr.nTrials             = 1;
  hdr.FirstTimeStamp      = tstamps(1);

  sel = startsWith(asc.msg.message,"RECCFG");
  msg = split(asc.msg.message(sel), " ");
  if sum(sel)==1
    % if there's only a single msg, then the above yields a column
    msg = msg(:)';
  end
  Fs  = double(msg(:,3));
  assert(all(Fs==Fs(1)) && Fs(1)==Fs1);

  hdr.TimeStampPerSample  = 1;
  hdr.Fs                  = Fs(1); 
  
  % give this warning only once
  ft_warning('creating fake channel names');
  hdr.label{1} = 'timestamps';
  for i=2:hdr.nChans
    hdr.label{i} = sprintf('chan%02d', i-1);
  end
  asc.hdr = hdr;
end

if needevt
  % convert all metadata in to events
  %type, timestamp, value, duration, duration, offset
  
  % events with 'eye' as value, but no duration
  fn = {'sfix' 'ssacc' 'sblink'};
  event = cell(0,1);
  for i = 1:numel(fn)
    if isfield(asc, fn{i})
      t = asc.(fn{i});
      t = renamevars(t, {'eye', 'stime'}, {'value', 'timestamp'});
      t.value = convertStringsToChars(t.value);
      t = t(:, contains(t.Properties.VariableNames, {'value' 'timestamp' 'sample'}));
      n = size(t,1);
      type   = repmat(fn(i), n, 1);
      offset = repmat({[]},  n, 1);
      duration = repmat({[]}, n, 1);
      t = cat(2, t, table(type), table(offset), table(duration));
      event{end+1} = table2struct(t);
    end
  end

  % events with 'eye' as value, but with duration
  fn = {'efix' 'eblink' 'esacc'};
  for i = 1:numel(fn)
    if isfield(asc, fn{i})
      t = asc.(fn{i});
      t = renamevars(t, {'eye', 'stime' 'sample_duration'}, {'value', 'timestamp' 'duration'});
      t.value = convertStringsToChars(t.value);
      t = t(:, ismember(t.Properties.VariableNames, {'value' 'timestamp' 'duration' 'sample'}));
      n = size(t,1);
      type   = repmat(fn(i), n, 1);
      offset = repmat({[]},  n, 1);
      t = cat(2, t, table(type), table(offset));
      event{end+1} = table2struct(t);
    end
  end
  
  % other type of events: input
  if isfield(asc, 'input')
    t = asc.input;
    t = renamevars(t, {'stime'}, {'timestamp'});
    t = t(:, contains(t.Properties.VariableNames, {'value' 'timestamp' 'sample'}));
    n = size(t,1);
    type   = repmat({'input'}, n, 1);
    offset = repmat({[]},  n, 1);
    duration = repmat({[]}, n, 1);
    t = cat(2, t, table(type), table(offset), table(duration));
    event{end+1} = table2struct(t);
  end

  % other type of events: msg
  if isfield(asc, 'msg')
    t = asc.msg;
    t = renamevars(t, {'message' 'stime'}, {'value' 'timestamp'});
    t = t(:, contains(t.Properties.VariableNames, {'value' 'timestamp' 'sample'}));
    sel = isfinite(t.sample) & ~startsWith(t.value, '!');
    t = t(sel,:);
    t.value = convertStringsToChars(t.value);

    n = size(t,1);
    type   = repmat({'msg'}, n, 1);
    offset = repmat({[]},  n, 1);
    duration = repmat({[]}, n, 1);
    t = cat(2, t, table(type), table(offset), table(duration));
    event{end+1} = table2struct(t);
  end

  event = cat(1, event{:});
  [srt, ix] = sort([event.sample]);
  event = event(ix);

  asc.event = event;
end

if needdat
  datline = asc.datline(begsample:endsample);
  ntab    = count(datline, sprintf('\t'));

  % check whether all lines have the same number of columns
  assert(all(ntab==ntab(1)));

  % identify the chunks of consecutive '...' and comments
  chunks = 0:100000:numel(datline);
  asc.dat = nan(ntab(1)+1, numel(datline));
  idy = true(ntab(1)+1,1);
  for i = 1:numel(chunks)
    fprintf('parsing the numeric data in chunk %d/%d\n', i, numel(chunks));
    idx = [(chunks(i)+1) chunks(i)+100000];
    idx(idx>numel(datline)) = numel(datline);
    d   = split(datline(idx(1):idx(2)), sprintf('\t'));
    asc.dat(:,idx(1):idx(2)) = double(d(:,idy))';
  end
  asc.dat = asc.dat(chanindx,:);
end


function out = timestamp2samples(in, tstamps, samples, fs)

out = nan(size(in));
for i = 1:min(size(in,2), 2)
  indx = (fs/1000).*(in(:, i) - tstamps(1)) + 1;
  if ~all(indx==round(indx))
    ft_warning('rounding off time stamps to the nearest sample indices');
    indx = round(indx);
  end
  out(indx>0 & indx<=numel(samples), i) = samples(indx(indx>0 & indx<=numel(samples)));
end
% assume 3d column always to be duration
if size(in,2)==3
  out(:, 3)   = (fs/1000).*in(:, 3);
end
