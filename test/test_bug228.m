% this script tests whether the change in name for baselinecorrection is
% done appropriately

% strategy: rename blc into demean in preproc, and blcwindow into
% baselinewindow.
% also change the options in ft_checkconfig
% 
% then grep all occurrences of blc and blcwindow. use the renamed option of
% ft_checkconfig in all those functions to ensure backward compatibility
% this is only needed for those functions which directly call
% ft_preproc_baselinecorrect. functions calling preproc should be fine for
% backward compatibility. in that case, only change the default naming

% don't forget to change the wiki FIXME

tname = [tempname,'.txt'];
system(['grep blc *.m > ',tname]);
delete(tname);