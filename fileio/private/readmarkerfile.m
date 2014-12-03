function [marker] = readmarkerfile(folder)

% Read the MarkerFile.mrk file in a CTF dataset.
%
% Use as
%   marker = readmarkerfile(folder)
%
% Creates a marker structure which contains number_markers,
% number_samples, marker_names, and trial_times.

% Contributed by Tom Holroyd by email on 21 May 2005

name = fullfile(folder, 'MarkerFile.mrk');
if ~exist(name, 'file')
  error('%s not found', name);
end

f = fopen(name, 'rt');
markfile = {};
while true
  l = fgetl(f);
  if ~ischar(l)
    break
  end
  markfile{end + 1} = l;
end
fclose(f);

% Only one of these.
i = strmatch('NUMBER OF MARKERS:', markfile, 'exact') + 1;
nmarkers = str2num(markfile{i});

% Get all the marker names.
i = strmatch('NAME:', markfile, 'exact') + 1;
names = markfile(i);

% Get the number of samples for each.
i = strmatch('NUMBER OF SAMPLES:', markfile, 'exact') + 1;
nsamples = str2num(char(markfile(i)));

for i = 1:length(nsamples)
  if nsamples(i) == 0
    warning('marker %s in %s has zero samples', names{i}, folder);
  end
end

% Get the samples.  Each is trial and time in seconds.
j = strmatch('LIST OF SAMPLES:', markfile, 'exact') + 2;
for i = 1:nmarkers
  marks{i} = str2num(char(markfile(j(i):j(i) + nsamples(i))));
  
  % Convert from index origin 0 to 1
  if nsamples(i) ~= 0
    marks{i}(:, 1) = marks{i}(:, 1) + 1;
  end
end

marker = struct('number_markers', {nmarkers}, 'number_samples', {nsamples}, ...
  'marker_names', {names}, 'trial_times', {marks});

