function test_example_checkconfig

% MEM 4gb
% WALLTIME 00:10:00

%
%% How to use ft_checkconfig
%
% The function **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** checks the input configuration (cfg) for the main FieldTrip functions. This is similar to what **[ft_checkdata](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkdata.m)** does for the input data. The **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** function is automatically called when you use a FieldTrip function. You may not even notice this, unless it gives you feedback about your cfg (i.e., warning or error messages).
%
%
%% # Introduction to ft_checkconfig
%
% This function checks the input cfg of the main FieldTrip functions in three ways:
%
% a warning when renamed or deprecated options are used, and it makes sure
% no forbidden options are used. If necessary and possible, this function
% will adjust the cfg to the input requirements. If the input cfg does NOT
% correspond to the requirements, this function gives an elaborate warning
% message.
%
% functions, by putting them into substructures or converting them into the
% required format.
%
% relevant and used fields. The size of fields in the output cfg is also
% controlled: fields exceeding a certain maximum size are emptied.
%
%% # How to control the behavior of ft_checkconfig
%
% When you use a FieldTrip function, this automatically calls **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** to check the cfg you supplied. If necessary, **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** will give you feedback. How can you control this feedback? As explained in the help documentation
%
% The behavior of checkconfig can be controlled by the following cfg options,
% which can be set as global FieldTrip defaults (see FT_DEFAULTS
%   cfg.checkconfig = 'pedantic', 'loose' or 'silent' (control the feedback behavior of checkconfig)
%   cfg.trackconfig = 'cleanup', 'report' or 'off'
%   cfg.checksize   = number in bytes, can be inf (set max size allowed for output cfg fields)

% When you use a FieldTrip function, this automatically calls the function `ft_defaults`, which takes care of path setting, plus it sets defaults to be used throughout FieldTrip. It does this by creating a global variable called `ft_default`, i.e. a variable that is available to all functions, but not directly visible to the user. You can make it visible by typing `global ft_default`. The variable `ft_default` has the following fields and default settings that pertain to **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)**.
%
ft_default.checkconfig = 'loose';
ft_default.trackconfig = 'off';
ft_default.checksize   = 1e5;

% These settings control the behavior of **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)**. If you want to change these settings, either do it via the global variable (this way they will apply to all FieldTrip functions and automatically be added to the cfg), or do it directly via the cfg when you call a specific function. What are the available options?
%
%
% This setting determines the type of feedback **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** gives about the input cfg. An important function of **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** is to check whether the input cfg contains all the required options, and no renamed, unused, deprecated or forbidden options. If possible, **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** will adjust the cfg to the input requirements and give feedback to the user. You can control the type of feedback given. This can either be **'silent'**, which means no feedback is given at all, or **'loose'** which means warnings are given for all inconsistencies, or **'pedantic'**, which means errors are given for each inconsistency in your input cfg. Note that a missing required field in the cfg will always lead to an error, because FieldTrip simply will not run without it.
%
% To give an example using **[ft_freqdescriptives](https://github.com/fieldtrip/fieldtrip/blob/release/ft_freqdescriptives.m)**, if your cfg contains the field 'jacknife' (which was used in a previous version of FieldTrip but has since been renamed to 'jackknife'
%

% JM added this to make the function functional:
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg/freq_mtmfft_trl_ctf151.mat'));

cfg = [];
cfg.jacknife = 'no';
test = ft_freqdescriptives(cfg, freq)

% You will get the following feedback
%
% Warning: use cfg.jackknife instead of cfg.jacknife

% In this case, you don't have to do anything: **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** will rename the field for you, and **[ft_freqdescriptives](https://github.com/fieldtrip/fieldtrip/blob/release/ft_freqdescriptives.m)** can do its job. Of course the idea is, that you will use this feedback to improve your scripts!
%
%
%
%
%
% This determines the maximum size allowed for output cfg fields (i.e. the data.cfg). Some fields in the output cfg can be very large, e.g., the cfg.grid field when you do sourceanalysis. To avoid that several MBs or even GBs of your disk space are taken up by the data.cfg, you can set a maximum and **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** will empty all the cfg fields that are too big, before adding the cfg to the data. Crucial fields such as the cfg.trl and cfg.event will never be removed (currently, the following fields are ignored: 'checksize', 'trl', 'trlold', 'event', 'artifact', 'artfctdef', 'previous'). The default is set to *100000* bytes, but you can change this to anything you want. If you do not want any fields to be removed, set checksize to *inf*.
%
%% # How to use trackconfig to figure out which options were and were not used
%
% Ok, let's see this in action.
% In this example the tutorial dataset is used, but you can of course use this on any dataset.
% If you are using the tutorial dataset, first get the trial definitio
%
cfg = [];
cfg.dataset              = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.trialdef.eventtype   = 'backpanel trigger';
cfg.trialdef.prestim     = 1;
cfg.trialdef.poststim    = 2;
cfg.trialdef.eventvalue  = 3;
cfg = ft_definetrial(cfg);

% Now, run **[ft_preprocessing](https://github.com/fieldtrip/fieldtrip/blob/release/ft_preprocessing.m)** with trackconfig set to cleanup. This will give you both a report (on screen) and a cleaned data.cfg:
%
cfg.trackconfig = 'cleanup';
cfg.unused = 1; % to test whether it works!
data = ft_preprocessing(cfg);

% % The report will look like this:
% %
% The following config fields were specified by YOU and were USED
%   cfg.dataset
%   cfg.trackconfig
%   cfg.checkconfig
%   cfg.datafile
%   cfg.headerfile
% 
% The following config fields were specified by YOU and were NOT USED
%   cfg.trialdef
%   cfg.trialfun
%   cfg.unused
% 
% The following config fields were set to DEFAULTS and were USED
%   cfg.channel
%   cfg.removemcg
%   cfg.feedback
%   cfg.precision
%   cfg.padding
%   cfg.headerformat
%   cfg.dataformat
%   cfg.dftfilter
%   cfg.lpfilter
%   cfg.hpfilter
%   cfg.bpfilter
%   cfg.bsfilter
%   cfg.medianfilter
%   cfg.reref
%   cfg.refchannel
%   cfg.implicitref
%   cfg.continuous
%   cfg.polyremoval
%   cfg.detrend
%   cfg.blc
%   cfg.hilbert
%   cfg.derivative
%   cfg.rectify
%   cfg.boxcar
%   cfg.absdiff
%   cfg.conv
%   cfg.montage
% 
% The following config fields were set to DEFAULTS and were NOT USED
%   cfg.removeeog
%   cfg.polyorder
%   cfg.blcwindow
%   cfg.lpfiltord
%   cfg.hpfiltord
%   cfg.bpfiltord
%   cfg.bsfiltord
%   cfg.lpfilttype
%   cfg.hpfilttype
%   cfg.bpfilttype
%   cfg.bsfilttype
%   cfg.lpfiltdir
%   cfg.hpfiltdir
%   cfg.bpfiltdir
%   cfg.bsfiltdir
%   cfg.medianfiltord
%   cfg.dftfreq
%   cfg.dftinvert

% Thus, it specifies which options were set by you and whether they were used. As you can see, the cfg.unused indeed ends up as "set by you, not used". Furthermore, the report shows the options that were added by the **[ft_preprocessing](https://github.com/fieldtrip/fieldtrip/blob/release/ft_preprocessing.m)** function, and whether they were actually used for your analysis. There is quite a list of unused defaults. When we now look at the output cfg we see that it is nicely cleaned. Compare this with a call to **[ft_preprocessing](https://github.com/fieldtrip/fieldtrip/blob/release/ft_preprocessing.m)** without using trackconfig.
%
cfg.trackconfig = 'off';
data2 = ft_preprocessing(cfg);

data.cfg  % cleaned
data2.cfg % not cleaned

% This example showed how trackconfig can be used when doing preprocessing. The same approach can be applied for all the main FieldTrip functions. (Currently, trackconfig has been implemented in about 20 FieldTrip functions, this will be expanded.)
%
%% # How to use checksize to save disk space
%
% As explained above, you can use cfg.checksize to set a limit to the size of fields in the output cfg, **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** then empties fields that are too big. Crucial fields such as the trl will never be removed. This all pertains to the output cfg, i.e. the cfg that comes out of the function (depending on the FieldTrip function you are using, this is either cfg or data.cfg). This does _not_ change any data.cfg.previous fields.
%
% However, you may have lots of analysed data on disk, with data.cfgs that might be taking up quite some of your disk space. Especially after doing beamforming (sourceanalysis) the output cfg can be large, since the grid is always kept in the data.cfg. If you would like to free some disk space (and are sure you can do without these fields), the following trick can be applie
%
%%% script to downsize cfgs of stored data
%%% this can free up significant amounts of disk space

% JM commented this out because this is based on dummy data
% downsizefiles={
%   '/mydir/dataset1'
%   '/mydir/dataset2'
%   '/mydir/etc.'};
% 
% for k=1:length(downsizefiles)
%   load(downsizefiles{k})
% 
%   data.cfg.checksize=100000;
%   data.cfg=ft_checkconfig(data.cfg, 'checksize', 'yes');
% 
%   save(fullfile(tempdir, downsizefiles{k}), 'data')
% end

% This way **[ft_checkconfig](https://github.com/fieldtrip/fieldtrip/blob/release/ft_checkconfig.m)** will run recursively through the entire data.cfg, including all the previous fields, and empty the fields that are larger than the specified maximum. This can be used on all FieldTrip data.
