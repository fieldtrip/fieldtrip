function test_bug3144

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_checkdata ft_datatype_sens
% DATA private

%%

load(dccnpath('/project/3031000.02/test/bug3144a.mat'))

ft_checkdata(data);

ft_datatype_sens(data.grad);

%% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3144#c14

load(dccnpath('/project/3031000.02/test/bug3144b.mat'))

cfg=[];
timelock_cgra_A       = ft_combineplanar(cfg, timelock_gra_A);
