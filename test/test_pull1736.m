function test_pullXXX

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY openpose_keypoints

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/motion/openpose/S_Npred_02_crop_000000010521_keypoints.json');

assert(ft_filetype(filename, 'openpose_keypoints'));

%%

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

%%

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.viewmode = 'vertical';
cfg.colorgroups = 'labelchar8'; % this is the "people" number, i.e. 1 or 2
cfg.preproc.demean = 'yes';
cfg.blocksize = 30; % seconds
cfg.channel = find(contains(data.label, 'q'));
ft_databrowser(cfg, data);

%%

close all
figure

for i=1:25
  chanx = (i-1)*3 + 1;
  chany = (i-1)*3 + 2;
  chanq = (i-1)*3 + 3;
  
  datx = data.trial{1}(chanx,:);
  daty = data.trial{1}(chany,:);
  datq = data.trial{1}(chanq,:);
  
  datx(datq<0.8) = nan;
  daty(datq<0.8) = nan;
  
  subplot(5,5,i)
  plot(datx, daty, 'b')
  axis equal
  axis([1 704 1 576])
end

for i=1:25
  chanx = (i-1)*3 + 1 + 75;
  chany = (i-1)*3 + 2 + 75;
  chanq = (i-1)*3 + 3 + 75;
  
  datx = data.trial{1}(chanx,:);
  daty = data.trial{1}(chany,:);
  datq = data.trial{1}(chanq,:);
  
  datx(datq<0.8) = nan;
  daty(datq<0.8) = nan;
  
  subplot(5,5,i)
  hold on
  plot(datx, daty, 'r')
  axis equal
  axis([1 704 1 576])
end
