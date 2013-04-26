function test_bug2104

% TEST test_bug2104 
% TEST ft_redefinetrial ft_datatype_sens ft_datatype_raw
% TEST ft_checkdata channelposition

% % ORIGINAL BUGREPORT:
% Hi all, 
% 
% I would like to do a source analysis with MEG data (ctf) which was preprocessed
% with an older fieldtrip version. For this, I need to redefine the trials of the
% preprocessed data, which I did with the command from the fieldtrip tutorial:
% 
% cfg = [];                                           
% cfg.toilim = [-.5 0];                       
% dataPre = ft_redefinetrial(cfg, data_clean);
% 
% 
% 
% However, I get following error message: 
% 
% 
% Attempted to access numcoils(1); index out of bounds because
% numel(numcoils)=0.

% Note: this is probably caused by mismatch between old and new sensor definition


% load data
load test_bug2104.mat

cfg = [];                                           
cfg.toilim = [-.5 0];                       
data = ft_redefinetrial(cfg, data);




