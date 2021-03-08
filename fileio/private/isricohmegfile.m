function val = isricohmegfile(filename)

%% The extentions, .con, .ave, and .mrk are common between Ricoh and Yokogawa systems.
%% isricohmegfile idetifies whether the file is generated from Ricoh system or not.
%% This function uses a function in YOKOGAWA_MEG_READER toolbox, getYkgwHdrSystem.

if ft_hastoolbox('yokogawa_meg_reader', 3)
	sys_info = getYkgwHdrSystem(filename);
	ver = sys_info.version;
	if ver > 2   % The file is Ricoh one.
		val = true;
		ft_hastoolbox('ricoh_meg_reader', 3);
	else
		val = false;
	end
else
	msg = 'yokogawa_meg_reader toolbox is not installed.';
	error(msg)
end
