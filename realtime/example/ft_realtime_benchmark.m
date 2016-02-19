function [TW, TPW, TR, TPR] = ft_realtime_benchmark(target)

% FT_REALTIME_BENCHMARK times the reading and writing of data
%
% Use as
%   ft_realtime_benchmark(target)
% where target is the location of the buffer or the file, for
% example 'buffer://localhost:1972'
%
% Please note that any data in the target will be overwritten.

% Copyright (C) 2010, Stefan Klanke
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

N = 100;

hdr        = [];
hdr.Fs     = 1000; % this is not really used here
hdr.nChans = 64; 

fprintf(1,'Writing header and flushing data in %s...\n', target);
ft_write_data(target, single([]), 'header', hdr, 'append', false);

blocklen = [8, 16, 32, 64, 128, 256, 512, 1024, 2048];
blocksize = hdr.nChans * blocklen;

TW  = zeros(N, length(blocklen));	% Times for writing
TPW = zeros(1, length(blocklen));	% Throughput writing
TR  = zeros(N, length(blocklen));	% Times for reading
TPR = zeros(1, length(blocklen));	% Throughput reading


for k=1:length(blocklen);
  fprintf(1,'\nDetermining WRITE throughput at %ix%i samples per block...\n', hdr.nChans, blocklen(k));
  % generate random data for writing
  X = single(randn(hdr.nChans, blocklen(k)));
  for n=1:N;
    tic;
    ft_write_data(target, X, 'header', hdr, 'append', true);
    TW(n,k) = toc;
  end
  mt = mean(TW(:,k));
  st = std(TW(:,k));
  ss = blocksize(k) * N / sum(TW(:,k));
  fprintf(1,'Time per block %f +/- %f  => %d samples/sec\n', mt, st, ss);
  TPW(k) = ss;
  
  hdr = ft_read_header(target);
  
  fprintf(1,'\nDetermining READ throughput at %ix%i samples per block...\n', hdr.nChans, blocklen(k));
  begS = 1;
  for n=1:N;
    endS = begS+blocklen(k)-1;
    tic;
    X = ft_read_data(target, 'header', hdr, 'begsample', begS, 'endsample', endS);
    TR(n,k) = toc;
    begS = endS+1;
  end
  mt = mean(TR(:,k));
  st = std(TR(:,k));
  ss = blocksize(k) * N / sum(TR(:,k));
  fprintf(1,'Time per block %f +/- %f  => %d samples/sec\n', mt, st, ss);
  
  TPR(k) = ss;
end

loglog(blocksize, TPW,'r+-', blocksize, TPR, 'b+-');
legend('Writing','Reading','Location','NorthWest');
xlabel('block size');
ylabel('samples per sec');

  
  
  
  
