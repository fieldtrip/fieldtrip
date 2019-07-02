function test_issue1159

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_senstype undobalancing

%%
% see https://github.com/fieldtrip/fieldtrip/issues/1159

load issue1159.mat

%%
% this works

ft_senstype(grad)

%%
% this fails

ft_senstype(rmfield(grad, 'type'))