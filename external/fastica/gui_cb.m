function gui_cb(action)
%
% This file is used by FASTICAG

% This file holds the callbacks to the main window

% @(#)$Id: gui_cb.m,v 1.5 2003/09/10 10:33:41 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the main figure
global hf_FastICA_MAIN;

% Handles for needed controls in main figure;
global ht_FastICA_mixedStatus;
global ht_FastICA_dim;
global ht_FastICA_numOfSamp;
global ht_FastICA_newDim;
global ht_FastICA_whiteStatus;
global ht_FastICA_icaStatus;
global hpm_FastICA_approach;
global he_FastICA_numOfIC;
global hpm_FastICA_g;
global hpm_FastICA_stabilization;

% Main values are stored here
global g_FastICA_mixedsig;
global g_FastICA_mixedmean;
global g_FastICA_pca_D;
global g_FastICA_pca_E;
global g_FastICA_white_sig;
global g_FastICA_white_wm;
global g_FastICA_white_dwm;
global g_FastICA_ica_sig;
global g_FastICA_ica_A;
global g_FastICA_ica_W;
global g_FastICA_initGuess;
global g_FastICA_approach;
global g_FastICA_numOfIC;
global g_FastICA_g;
global g_FastICA_finetune;
global g_FastICA_a1;
global g_FastICA_a2;
global g_FastICA_myy;
global g_FastICA_stabilization;
global g_FastICA_epsilon;
global g_FastICA_maxNumIte;
global g_FastICA_maxFinetune;
global g_FastICA_sampleSize;
global g_FastICA_initState;
global g_FastICA_displayMo;
global g_FastICA_displayIn;
global g_FastICA_verbose;

% String values are here
global c_FastICA_appr_strV;
global c_FastICA_g1_strD;
global c_FastICA_g1_strV;
global c_FastICA_g2_strD;
global c_FastICA_g2_strV;
global c_FastICA_finetune_strD;
global c_FastICA_finetune_strV;
global c_FastICA_stabili_strV;
global c_FastICA_iSta_strV;
global c_FastICA_dMod_strV;
global c_FastICA_verb_strV;

% What is the load type of load dialog
global g_FastICA_loadType;

% Global variable for stopping the ICA calculations
global g_FastICA_interrupt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What ever we do, it will take some time... not much, but some :-)
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'InitAll'

 % If the data is already loaded, then get the information from data
 % and show to the user (also set g_FastICA_numOfIC)
 if ~isempty(g_FastICA_mixedsig)
   set(ht_FastICA_mixedStatus, 'String', '');
   [dim, numofsamp] = size(g_FastICA_mixedsig);
   set(ht_FastICA_dim, 'String', int2str(dim));
   set(ht_FastICA_numOfSamp, 'String', int2str(numofsamp));
   set(ht_FastICA_newDim, 'String', int2str(dim));
   set(he_FastICA_numOfIC, 'String', int2str(dim));
   g_FastICA_numOfIC = dim;
   g_FastICA_mixedmean = [];
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'LoadData'

 handle = findobj('Tag','f_FastICALoad');  % Check if the window is already
 if isempty(handle)                        % open. If not then open it.
   pos = get(hf_FastICA_MAIN, 'Position');
   % Based on the feedback obtained from some users, it seems
   % that at least in some systems, pos can sometimes be empty. A
   % similar check is done a few lines below.
   if ~isempty (pos),
     gui_l(pos(1), pos(2));
   else
     gui_l (0, 0);
   end
 else
   if strcmp(g_FastICA_loadType, 'data')   % Check if it was the same load
     figure(handle);                       % window. If it wasn't then
   else                                    % close the other window first
     close(handle);                        % and then open the load window
     fprintf('''Load initial guess'' -dialog closed!\n');
     pos = get(hf_FastICA_MAIN, 'Position');
     if ~isempty (pos),
       gui_l(pos(1), pos(2));
     else
       gui_l (0, 0);
     end
   end
 end
 
 % gui_cb NewData; - is called from the load function if not canceled...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'NewData'

 % New data is loaded or the old data changed. We need to find out
 % somethings about the new data... and do some other stuff also...
 [dim, numofsamp] = size(g_FastICA_mixedsig);
 set(ht_FastICA_dim, 'String', dim);
 set(ht_FastICA_newDim, 'String', dim);
 set(ht_FastICA_numOfSamp, 'String', numofsamp);
 set(he_FastICA_numOfIC, 'String', int2str(dim));
 
 g_FastICA_numOfIC = dim;    % Default for numOfIC = the new dimension
			     % PCA needs to be calculated again.
 g_FastICA_pca_E = [];       % We use this to check if PCA is calculated
 gui_cb NullWhite;           % Whitening needs to be done again also
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'NullWhite'

 % Whitening needs to done again next time it's needed
 g_FastICA_white_sig = [];   % We use this to check if whitening is calculated
 gui_cb NullICA;             % The FPICA must be calculated again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'NullICA'

 % If IC's are needed they have to bee calculated again
 g_FastICA_ica_sig = [];     % We use this to check if FPICA is calculated
 set(ht_FastICA_icaStatus,'String','Not yet done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Transpose'

 if isempty(g_FastICA_mixedmean)
   g_FastICA_mixedsig = g_FastICA_mixedsig';
 else
   g_FastICA_mixedsig = (g_FastICA_mixedsig + ...
			 g_FastICA_mixedmean * ...
			 ones(1,size(g_FastICA_mixedsig, 2)))';
   g_FastICA_mixedmean = [];
 end
 gui_cb NewData;             % Data has been changed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'DoPCA'

 if ~isempty(g_FastICA_mixedsig)
   % We'll remove mean of the data here also just in case...
   if isempty(g_FastICA_mixedmean)
     [g_FastICA_mixedsig, g_FastICA_mixedmean] = remmean(g_FastICA_mixedsig);
   end
   
   % Do PCA interactively: ask the user for eigenvalues
   [g_FastICA_pca_E, g_FastICA_pca_D] = pcamat(g_FastICA_mixedsig, ...
					       0, 0, 'gui', ...
					       deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));
   
   newdim = size(g_FastICA_pca_D, 1);
   set(ht_FastICA_newDim, 'String', int2str(newdim));
   set(he_FastICA_numOfIC, 'String', int2str(newdim));
   g_FastICA_numOfIC = newdim;
   gui_cb NullWhite;           % Whitening needs to be done again also
			       % but we'll do it when it's needed.
 else
   fprintf('Data not loaded yet!\n\n');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'OrigDim'

 gui_cb NewData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ShowMixed'

 if ~isempty(g_FastICA_mixedsig)
   handle = findobj('Tag','f_FastICA_mix');  % Check if the window is already
   if isempty(handle)                        % open. If not then open it.
     figure('Tag', 'f_FastICA_mix', ...
	    'Name', 'FastICA: Plot data', ...
	    'NumberTitle', 'off');
   else
     figure(handle);
     clf;		% clear the figure for next plots
   end
   if isempty(g_FastICA_mixedmean)
     icaplot('dispsig',g_FastICA_mixedsig, 0, 0, 0, 'Mixed signals');
   else
     icaplot('dispsig',g_FastICA_mixedsig + g_FastICA_mixedmean * ...
	     ones(1, size(g_FastICA_mixedsig, 2)), 0, 0, 0, 'Mixed signals');
   end
 else
   fprintf('Data not loaded yet!\n\n');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ShowWhite'

 if ~isempty(g_FastICA_mixedsig)
   if isempty(g_FastICA_white_sig)     % if whitening is not done, we need to
     gui_cb Whiten;                    % do it before we can display the
   end                                 % whitened signals
   
   handle = findobj('Tag','f_FastICA_white');  % Check if the window is already
   if isempty(handle)                          % open. If not then open it.
     figure('Tag', 'f_FastICA_white', ...
	    'Name', 'FastICA: Plot whitened', ...
	    'NumberTitle', 'off');
   else
     figure(handle);
     clf;		% clear the figure for next plots
   end
   icaplot('dispsig',g_FastICA_white_sig,0,0,0,'Whitened signals');
 else
   fprintf('Data not loaded yet!\n\n');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Whiten'

 set(ht_FastICA_whiteStatus,'String','Computing...');
 
 % If PCA is not calculated, we'll have to calculate it now,
 % we'll do it without guestions - we don't reduce the dimension
 % here - but PCAMAT might reduce the dimension automatically.
 if isempty(g_FastICA_pca_E)
   % We'll remove mean of the data here also just in case...
   if isempty(g_FastICA_mixedmean)
     [g_FastICA_mixedsig, g_FastICA_mixedmean] = remmean(g_FastICA_mixedsig);
   end
   
   [g_FastICA_pca_E, g_FastICA_pca_D] = pcamat(g_FastICA_mixedsig, 1, ...
					       size(g_FastICA_mixedsig, ...
						    1), 'off', ...
					       deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));
   
   % Check if the dimension was reduced automatically
   newdim = size(g_FastICA_pca_D, 1);
   set(ht_FastICA_newDim, 'String', int2str(newdim));
   % Check if the numOfIC now has illegal value entered 
   % We do that by telling the program that there is new value 
   % entered for NumOfIC.
   gui_cb ChangeNumOfIC;
 end
 
 % And now we can calculate whitening...
 [g_FastICA_white_sig, g_FastICA_white_wm, g_FastICA_white_dwm] = ...
								  whitenv(g_FastICA_mixedsig, g_FastICA_pca_E, g_FastICA_pca_D, deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));
 
 set (ht_FastICA_whiteStatus,'String','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ChangeApproach'

 % Get the old value for g
 eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);
 old_g = deblank(g_str(g_FastICA_g,:));
 
 % Get and set the new value for approach
 g_FastICA_approach = get(hpm_FastICA_approach, 'Value');
 
 % The possible values for g depend on the value of approach...
 eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strD;']);
 set(hpm_FastICA_g, 'String', g_str);
 
 % Match the old g value from the new g values so that if the 
 % old_g can be found from the new values (anywhere), then set new g
 % to that value, and if it's not found then set the new value to 1.
 match = 0;
 eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);
 for i=1:size(g_str,1)
   if strcmp(old_g, deblank(g_str(i,:)))
     match = i;
   end
 end
 if match == 0
   match = 1;   % the old g is not availabe anymore, set g = 1.
 end
 g_FastICA_g = match;
 set(hpm_FastICA_g, 'Value', match);
 
 gui_cb NullICA;      % The options are changed so we must calculate ICA again
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ChangeNumOfIC'

 % Get the new value... and store it later on after some checks
 numofic = str2num(get(he_FastICA_numOfIC, 'String'));
 
 % The number of IC can't be less than 1 or more than the reduced dimension.
 numoficmax = str2num(get(ht_FastICA_newDim, 'String'));
 if numofic < 1
   set(he_FastICA_numOfIC, 'String', '1');
   g_FastICA_numOfIC = 1;
 elseif numofic > numoficmax
   set(he_FastICA_numOfIC, 'String', int2str (numoficmax));
   g_FastICA_numOfIC = numoficmax;
 else
   g_FastICA_numOfIC = numofic;
 end
 
 gui_cb NullICA;      % The options are changed so we must calculate ICA again
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ChangeG'

 % Get the new value for g.
 g_FastICA_g = get(hpm_FastICA_g, 'Value');
 
 gui_cb NullICA;      % The options are changed so we must calculate ICA again
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ChangeStab'

 % Get the new value for g.
 g_FastICA_stabilization = get(hpm_FastICA_stabilization, 'Value');
 
 gui_cb NullICA;      % The options are changed so we must calculate ICA again
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'AdvOpt'

 handle = findobj('Tag','f_FastICAAdvOpt');
 if isempty(handle)                        % Check to see if the window is
   pos = get(hf_FastICA_MAIN, 'Position');     % already open...
   if ~isempty (pos),
     gui_adv(pos(1), pos(2));
   else
     gui_adv(0, 0);
   end
 else
   figure(handle)
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ShowICASig'

 if ~isempty(g_FastICA_mixedsig)
   % If the IC's are not already calculated, we'll do it now
   if isempty(g_FastICA_ica_sig)
     gui_cb DoFPICA;
   end
   
   % The signals may have been already displaued by the FPICA function
   % BUT the FPICA may also have shown either the basis of the filters
   % so the signals still need to be shown - besides the mean was added
   % in only later after FPICA
   
   % Also notice that in this version if there was something wrong in FPICA
   % Then the results are []. - We don't try to plot them!
   if ~isempty(g_FastICA_ica_sig')
     handle = findobj('Tag','f_FastICA_ica');  % Check if the window is already
     if isempty(handle)                        % open. If not then open it.
       figure('Tag', 'f_FastICA_ica', ...
	      'Name', 'FastICA: Plot ICs', ...
	      'NumberTitle', 'off');
     else
       figure(handle);
       clf;		% clear the figure for next plots
     end
     
     icaplot('dispsig', g_FastICA_ica_sig, 0, ...
	     0, 0, 'Independent components');
   end
 else
   fprintf('Data not loaded yet!\n\n');
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'DoFPICA'

 gui_cb DisableButtons;
 g_FastICA_interrupt = 0;
 
 if ~isempty(g_FastICA_mixedsig)
   if isempty(g_FastICA_white_sig)     % We need the whitened signal here
     gui_cb Whiten;
   end
   
   set(ht_FastICA_icaStatus,'String','Computing...');
   
   % The possible values for g depend on approach
   eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);
   
   % We'll contruct a command string which we'll later on evaluate...
   % This is where the Fixed point algorithm is used.
   command_str = ['[g_FastICA_ica_A,g_FastICA_ica_W]=' ...
		  'fpica(g_FastICA_white_sig,' ... 
		  'g_FastICA_white_wm,' ... 
		  'g_FastICA_white_dwm,' ...
		  '''' deblank(c_FastICA_appr_strV(g_FastICA_approach,:)) ...
		  ''',' ...
		  'g_FastICA_numOfIC,' ...
		  '''' deblank(g_str(g_FastICA_g,:)) ''',' ...
		  '''' ...
		  deblank(c_FastICA_finetune_strV(g_FastICA_finetune,:)) ...
		  ''',' ...
		  'g_FastICA_a1,' ...
		  'g_FastICA_a2,' ...
		  'g_FastICA_myy,' ...
		  '''' ...
		  deblank(c_FastICA_stabili_strV(g_FastICA_stabilization,:)) ...
		  ''',' ...
		  'g_FastICA_epsilon,' ...
		  'g_FastICA_maxNumIte,' ...
		  'g_FastICA_maxFinetune,' ...
		  '''' deblank(c_FastICA_iSta_strV(g_FastICA_initState,:)) ...
		  ''',' ...
		  'g_FastICA_initGuess,' ...
		  'g_FastICA_sampleSize,' ...
		  '''' deblank(c_FastICA_dMod_strV(g_FastICA_displayMo,:)) ...
		  ''',' ...
		  'g_FastICA_displayIn,' ...
		  '''' deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)) ...
		  ''');'];
   

   % If the user wants to plot while computing...
   % let's at least plot it to the right figure then
   if ~strcmp(deblank(c_FastICA_dMod_strV(g_FastICA_displayMo,:)),'off')
     handle = findobj('Tag','f_FastICA_ica');  % Check if the window is already
     if isempty(handle)                        % open. If not then open it.
       figure('Tag', 'f_FastICA_ica', ...
	      'Name', 'FastICA: Plot ICs', ...
	      'NumberTitle', 'off');
     else
       figure(handle);
       clf;		% clear the figure for next plots
     end
   end
   
   % ... and so let's do it...
   eval(command_str);
   
   % Also notice that in this version if there was something wrong in FPICA
   % Then the results are [].
   if ~isempty(g_FastICA_ica_W)
     % Add the mean back in.
     g_FastICA_ica_sig = g_FastICA_ica_W * g_FastICA_mixedsig ...
			 + (g_FastICA_ica_W * g_FastICA_mixedmean) ...
			 * ones(1,size(g_FastICA_mixedsig, 2));
     set (ht_FastICA_icaStatus,'String','Done');
   else
     gui_cb NullICA;  % set icasig=[] and do what ever needs to be done then
   end
   
   if ~(g_FastICA_interrupt)
     gui_cb EnableButtons;
   end
 else
   fprintf('Data not loaded yet!\n\n');
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Interrupt'
 g_FastICA_interrupt = 1;
 set(ht_FastICA_icaStatus,'String','Interrupted');
 gui_cb EnableButtons;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'DisableButtons'
 set(findobj('Tag','b_Transpose'),'Enable','off');
 set(findobj('Tag','b_ShowMixed'),'Enable','off');
 set(findobj('Tag','b_DoPCA'),'Enable','off');
 set(findobj('Tag','b_OrigDim'),'Enable','off');
 set(findobj('Tag','b_ShowWhite'),'Enable','off');
 set(findobj('Tag','b_advOpt'),'Enable','off');
 set(findobj('Tag','b_ShowICASig'),'Enable','off');
 set(findobj('Tag','b_LoadData'),'Enable','off');
 set(findobj('Tag','b_DoFPICA'),'Enable','off');
 set(findobj('Tag','b_SaveData'),'Enable','off');
 set(findobj('Tag','b_Quit'),'Enable','off');
 set(findobj('Tag','b_Interrupt'),'Visible','on');
 drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'EnableButtons'
 set(findobj('Tag','b_Transpose'),'Enable','on');
 set(findobj('Tag','b_ShowMixed'),'Enable','on');
 set(findobj('Tag','b_DoPCA'),'Enable','on');
 set(findobj('Tag','b_OrigDim'),'Enable','on');
 set(findobj('Tag','b_ShowWhite'),'Enable','on');
 set(findobj('Tag','b_advOpt'),'Enable','on');
 set(findobj('Tag','b_ShowICASig'),'Enable','on');
 set(findobj('Tag','b_LoadData'),'Enable','on');
 set(findobj('Tag','b_DoFPICA'),'Enable','on');
 set(findobj('Tag','b_SaveData'),'Enable','on');
 set(findobj('Tag','b_Quit'),'Enable','on');
 set(findobj('Tag','b_Interrupt'),'Visible','off');
 drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'SaveData'

 handle = findobj('Tag','f_FastICASave');  % Check if the window is already
 if isempty(handle)                        % open. If not then open it.
   pos = get(hf_FastICA_MAIN, 'Position');
   if ~isempty (pos),
     gui_s(pos(1), pos(2));
   else
     gui_s(0, 0);
   end
 else
   figure(handle);                       % window. If it wasn't then
 end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Quit'

 % We'll close the other dialogs if they are open.
 Tags = ['f_FastICALoad  '
	 'f_FastICAAdvOpt'
	 'f_FastICASave  '
	 'f_FastICA_mix  '
	 'f_FastICA_white'
	 'f_FastICA_ica  '];
 for i=1:size(Tags,1)
   handle = findobj('Tag', deblank(Tags(i,:)));
   if ~isempty(handle)
     close(handle);
   end
 end

 % Close this window
 close(hf_FastICA_MAIN);
 
 % Clear the used global variables.
 gui_cg;

 % Use return to avoid reaching the watchoff statement at the end
 % (There used to be a 'break' statement here, but it resulted in
 % errors in more recent version of Matlab -- jarmo)
 return;
 % ... and we're done.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'About'

 gui_help('gui_cb_about');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Help'

 gui_help('gui_cb_help');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end    % switch

watchoff (watchonInFigure);
