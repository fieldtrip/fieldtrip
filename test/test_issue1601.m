function test_issue1601

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY bids_sidecar ft_read_header ft_read_sens

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1601'));

%%

hdr  = ft_read_header('sub-RESP0521/ses-1/ieeg/sub-RESP0521_ses-1_task-Sleep_run-060307_ieeg.vhdr');
dat  = ft_read_data  ('sub-RESP0521/ses-1/ieeg/sub-RESP0521_ses-1_task-Sleep_run-060307_ieeg.vhdr');
evt  = ft_read_event ('sub-RESP0521/ses-1/ieeg/sub-RESP0521_ses-1_task-Sleep_run-060307_ieeg.vhdr');
elec = ft_read_sens  ('sub-RESP0521/ses-1/ieeg/sub-RESP0521_ses-1_task-Sleep_run-060307_ieeg.vhdr');
