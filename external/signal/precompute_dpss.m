function w = precompute_dpss

% assume 1 second and 1000 Hz sampling rate

s = 1:0.5:50;
n = 1000;
w = cell(size(s));
for i=1:length(s)
  w{i} = dpss(n, n*s(i)/1000);
end


