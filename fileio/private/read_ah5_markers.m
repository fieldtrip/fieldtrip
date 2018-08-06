function [ event ] = read_ah5_markers(hdr, filename )

if isempty(hdr)
    orig = read_ahdf5_hdr(filename);
    hdr.orig = orig;
    hdr.Fs = orig.channels(1).samplingRate;
    hdr.nChans = numel(orig.channels);
    hdr.nSamples = orig.numberOfSamples;
    hdr.nTrials = orig.numberOfBlocks;
    hdr.nSamplesPre = 0;
    hdr.label = orig.label;
end
    
marks = h5read(filename, '/Blocks/1/markers');
labels = marks.label';
positions = marks.position';
durations = marks.duration';
values = marks.value';

for i=1:size(labels, 1)
    event(i).type = strcat(labels(i, 1:256));
    event(i).value = values(i);
    event(i).sample = positions(i) * hdr.Fs;
    event(i).duration = durations(i) * hdr.Fs;
end
end

