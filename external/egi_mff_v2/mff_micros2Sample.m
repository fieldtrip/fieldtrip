%% mff_micros2Sample.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012, 2013 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%
%  Converts from microseconds to samples, given the sampling rate. 
%%
function [sampleNum, remainder] = mff_micros2Sample(microsecs, sampRate)
microsecs = double(microsecs);
sampDuration = 1000000/sampRate;
sampleNum = microsecs/sampDuration;
remainder = uint64(rem(microsecs, sampDuration));
sampleNum = fix(sampleNum);
