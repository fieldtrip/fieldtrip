function tests = test_ft_preproc_lowpassfilter

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_lowpassfilter

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  fn = localfunctions;
  for i=1:numel(fn)
    feval(fn{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testOptions(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan   = 8;
nsample = 1000;

dat = randn(nchan, nsample) + 1;
Fs  = 1000;
Flp = 35;

instabilityfix  = [];
df              = []; 
wintype         = []; 
dev             = []; 
plotfiltresp    = [];
usefftfilt      = [];

filttypes = {'brickwall', 'firws' 'fir' 'firls' 'but'};
filtdirs  = {'onepass' 'onepass-reverse' 'twopass' 'twopass-reverse' 'twopass-average' 'onepass-zerophase' 'onepass-reverse-zerophase' 'onepass-minphase'};
filtorders = {[] 1 2 4 8};

result = {};
opts   = {};
for i1 = 1:numel(filttypes)
  if isequal(filttypes{i1}, 'brickwall')
    opts{end+1}   =[filttypes{i1}]; 
    result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, [], filttypes{i1}, [], instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
    continue;
  elseif isequal(filttypes{i1}, 'but')
    sel2 = 1:5;
    sel3 = 2:4; % order 8 may become unstable, and thus requires an instabilityfix
  elseif isequal(filttypes{i1}, 'firws')
    sel2 = find(~contains(filtdirs, 'twopass')); % this leads to a filter that is not stable
    sel3 = 1; % the low number orders do not make sense at all for the finite impulse filters.
  else
    sel2 = 1:numel(filtdirs);
    sel3 = 1;
  end
  for i2 = sel2
    for i3 = sel3
      tmp = filtorders{i3};
      if isempty(tmp)
        tmp = 'defaultorder';
      elseif isnumeric(tmp)
        tmp = sprintf('order %d', tmp);
      end
      opts{end+1}   = sprintf('%s_%s_%s', filttypes{i1}, filtdirs{i2}, tmp); 
      result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, filtorders{i3}, filttypes{i1}, filtdirs{i2}, instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
    end
  end
end

% result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, 8, 'but'      , 'onepass'                  , 'split', [], [], [], [], []); % this one is instable
% result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, 8, 'but'      , 'onepass-reverse'          , 'split', [], [], [], [], []);
% result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, 8, 'but'      , 'twopass'                  , 'split', [], [], [], [], []);
% result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, 8, 'but'      , 'twopass-reverse'          , 'split', [], [], [], [], []);
% result{end+1} = ft_preproc_lowpassfilter(dat, Fs, Flp, 8, 'but'      , 'twopass-average'          , 'split', [], [], [], [], []);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

for i=1:size(result{1},1)
  for k=1:numel(result)
    for m=(k+1):numel(result)
      b(k,m) = result{k}(i,:)/result{m}(i,:);
      b(m,k) = result{m}(i,:)/result{k}(i,:);
    end
    B(:,:,i) = b;
  end
end
n = numel(result);
figure;imagesc((mean(B,3)+mean(B,3)')/2);
caxis([0.82 1.01]);
set(gca, 'xtick', 1:n, 'ytick', 1:n, 'xticklabel', opts', 'yticklabel', opts', 'ticklabelinterpreter', 'none');
set(gcf, 'position', [230 47 993 750]);
colorbar
axis equal;axis tight
