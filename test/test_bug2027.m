% TEST test_bug2027
% TEST read_4d_hdr
% TEST ft_read_header

% Bug: ft_read_header returns only 152 MEG channels in the hdr.label
% This goes wrong at a low level, i.e. in read_4d_hdr

% reproduce
dataset = '/home/common/matlab/fieldtrip/data/test/bug2027/e,rfhp1.0Hz,COH';
hdr     = ft_read_header(dataset);

assert(numel(ft_channelselection('MEG',hdr.label))==248);

