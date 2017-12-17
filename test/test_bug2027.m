function test_bug2027

% MEM 1500mb
% WALLTIME 00:20:00

% TEST read_4d_hdr
% TEST ft_read_header
% TEST ft_channelselection

% Bug: ft_read_header returns only 152 MEG channels in the hdr.label
% This goes wrong at a low level, i.e. in read_4d_hdr

% reproduce
datasets = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2027/colorado/e,rfhp1.0Hz,COH')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2027/glasgow/e,rfDC')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2027/marseille/e,rfhp1.0Hz,COH')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2027/stlouis/e,rfhp1.0Hz,COH')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2027/konstanz/c,rfhp0.1Hz')
};
  
nummeg = [248 248 248 248 148];          
for k = 1:numel(datasets)
  hdr = ft_read_header(dccnpath(datasets{k}));
  assert(numel(ft_channelselection('MEG',hdr.label))==nummeg(k));
end
