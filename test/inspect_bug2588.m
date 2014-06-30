function inspect_bug2588

% TEST inspect_bug2588
% TEST ft_databrowser

data = [];
for i=1:10
  data.label{i} = num2str(i);
end

data.trial{1} = zeros(10,300);
data.time{1}  = (1:300)/300;
for i=1:10
  sel = ((10*(i-1)+1):(10*i)) + i*10;
  data.trial{1}(i,sel) = 1;
end

cfg = [];
ft_databrowser(cfg, data);
