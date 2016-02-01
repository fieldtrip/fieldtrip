function test_bug946

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_appenddata

% make some dummy data
data1 = [];
data1.trial = repmat({rand(10,100)},1,10);
data1.time = repmat({(1:100)/100},1,10);
data1.label = {'CH01','CH02','CH03','CH04','CH05','CH06','CH07','CH08','CH09','CH10'};
data1.fsample = 100;
data1.sampleinfo(:,1) = 1:100:901;
data1.sampleinfo(:,2) = 100:100:1000;
data1.trialinfo = (1:10)';

data2 = data1;
data2.sampleinfo = data2.sampleinfo + 1000;
data2.trialinfo = (11:20)';

% perform concatenation
cfg = [];
concat = ft_appenddata(cfg, data1, data2);

% check whether trialinfo and sampleinfo are correct
testSampleInfo = concat.sampleinfo==cat(1,data1.sampleinfo,data2.sampleinfo);
testTrialInfo = concat.trialinfo==cat(1,data1.trialinfo,data2.trialinfo);

wrong = [];
if ~all(testSampleInfo(:))
  wrong = [wrong 'sampleinfo,'];
end
if ~all(testTrialInfo(:))
  wrong = [wrong 'trialinfo,'];
end
if ~isempty(wrong)
  error([wrong(1:end-1) ' not concatenated properly']);
end

end
