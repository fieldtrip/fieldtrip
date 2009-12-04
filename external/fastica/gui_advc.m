function gui_advc (action)
%
% This file is needed by FASTICAG

% This file holds the callbacks for advanced options -dialog

% @(#)$Id: gui_advc.m,v 1.3 2003/09/08 11:28:58 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global hf_FastICA_AdvOpt;

% Handles to some of the controls in window
global hpm_FastICA_finetune;
global he_FastICA_a1;
global he_FastICA_a2;
global he_FastICA_myy;
global he_FastICA_epsilon;
global he_FastICA_maxIterations;
global he_FastICA_maxFinetune;
global he_FastICA_sampleSize;
global hpm_FastICA_initState;
global hb_FastICA_initGuess;
global ht_FastICA_initGuess;
global hpm_FastICA_displayMode;
global he_FastICA_displayInterval;
global hpm_FastICA_verbose;

% Needed handles to the main window
global hf_FastICA_MAIN;
global ht_FastICA_icaStatus;
global ht_FastICA_numOfSamp;
global hpm_FastICA_stabilization;

% Some of the main variables needed
global g_FastICA_initGuess;
global g_FastICA_numOfIC;
global g_FastICA_finetune;
global g_FastICA_a1;
global g_FastICA_a2;
global g_FastICA_myy;
global g_FastICA_epsilon;
global g_FastICA_maxNumIte;
global g_FastICA_maxFinetune;
global g_FastICA_initState;
global g_FastICA_sampleSize;
global g_FastICA_displayMo;
global g_FastICA_displayIn;
global g_FastICA_verbose;

global c_FastICA_iSta_strV;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Checka1'

 e_a1_val = str2num(get(he_FastICA_a1, 'String'));
 if e_a1_val <= 0
   set(he_FastICA_a1, 'String', '0.1');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Checka2'
  
 e_a2_val = str2num(get(he_FastICA_a2, 'String'));
 if e_a2_val <= 0
   set(he_FastICA_a2, 'String', '0.1');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'CheckMyy'

 e_myy_val = str2num(get(he_FastICA_myy, 'String'));
 if e_myy_val <= 0
   set(he_FastICA_myy, 'String', '0.1');
 elseif e_myy_val > 1
   set(he_FastICA_myy, 'String', '1');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'CheckSampleSize'
  
 e_sampleSize_val = str2num(get(he_FastICA_sampleSize, 'String'));
 if e_sampleSize_val > 1
   set(he_FastICA_sampleSize, 'String', '1');
 else 
   numOfSamp = str2num(get(ht_FastICA_numOfSamp, 'String'));
   if numOfSamp < 1000
     set(he_FastICA_sampleSize, 'String', '1');
     fprintf('Can''t reduce sample size. Already less than 1000 samples.\n');
   elseif (e_sampleSize_val * numOfSamp) < 1000
     fprintf('Can''t reduce sample size to less than 1000 samples.\n');
     set(he_FastICA_sampleSize, 'String', sprintf('%0.3f',1000/numOfSamp));
   end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'loadGuess'
 
 handle = findobj('Tag','f_FastICALoad');  % Check if the window is already
 if isempty(handle)                        % open. If not then open it.
   pos = get(hf_FastICA_MAIN, 'Position');
   gui_l(pos(1), pos(2));
 else
   if strcmp(g_FastICA_loadType, 'guess')  % Check if it was the same load
     figure(handle);                       % window. If it wasn't then
   else                                    % close the other window first
     close(handle);                        % and then open the load window
     fprintf('''Load data'' -dialog closed!\n');
     pos = get(hf_FastICA_MAIN, 'Position');
     gui_l(pos(1), pos(2));
   end
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'OK'
  
 gui_advc Apply;

 close(hf_FastICA_AdvOpt);

 % Use return to avoid reaching the watchoff statement at the end
 % (There used to be a 'break' statement here, but it resulted in
 % errors in more recent version of Matlab -- jarmo)
 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Apply'

 newValues = 0;

 val = g_FastICA_finetune;
 g_FastICA_finetune = get(hpm_FastICA_finetune, 'Value');
 if (g_FastICA_finetune ~= val)
   newValues = 1;
 end
 
 val = g_FastICA_a1;
 g_FastICA_a1 = str2num(get(he_FastICA_a1, 'String'));
 if (g_FastICA_a1 ~= val)
   newValues = 1;
 end

 val = g_FastICA_a2;
 g_FastICA_a2 = str2num(get(he_FastICA_a2, 'String'));
 if (g_FastICA_a2 ~= val)
   newValues = 1;
 end

 val = g_FastICA_myy;
 g_FastICA_myy = str2num(get(he_FastICA_myy, 'String'));
 if (g_FastICA_myy ~= val)
   newValues = 1;
 end

 % If myy < 1 then will use stabilazed code, and we don't care 
 % about the parameter stabilization :-)
 if (g_FastICA_myy == 1)
   set(hpm_FastICA_stabilization, 'Enable', 'on');
 else
   set(hpm_FastICA_stabilization, 'Enable', 'off');
 end

 val = g_FastICA_epsilon;
 g_FastICA_epsilon = str2num(get(he_FastICA_epsilon, 'String'));
 if (g_FastICA_epsilon ~= val)
   newValues = 1;
 end

 val = g_FastICA_maxNumIte;
 g_FastICA_maxNumIte = str2num(get(he_FastICA_maxIterations, 'String'));
 if (g_FastICA_maxNumIte ~= val)
   newValues = 1;
 end
 
 val = g_FastICA_maxFinetune;
 g_FastICA_maxFinetune = str2num(get(he_FastICA_maxFinetune, 'String'));
 if (g_FastICA_maxFinetune ~= val)
   newValues = 1;
 end
 
 val = g_FastICA_sampleSize;
 g_FastICA_sampleSize = str2num(get(he_FastICA_sampleSize, 'String'));
 if (g_FastICA_sampleSize ~= val)
   newValues = 1;
 end

 val = g_FastICA_initState;
 g_FastICA_initState = get(hpm_FastICA_initState, 'Value');
 if (g_FastICA_initState ~= val)
   newValues = 1;
 end

 val = g_FastICA_initGuess;
 g_FastICA_initGuess = get(hb_FastICA_initGuess, 'UserData');
 if min(size(val) == size(g_FastICA_initGuess)) == 0
   newValues = 1;
 else
   if (g_FastICA_initGuess ~= val)
     newValues = 1;
   end
 end
 
 g_FastICA_displayMo = get(hpm_FastICA_displayMode, 'Value');
 g_FastICA_displayIn = str2num(get(he_FastICA_displayInterval, 'String'));
 g_FastICA_verbose = get(hpm_FastICA_verbose, 'Value');
 
 if newValues == 1
   gui_cb NullICA;
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Cancel'
  
 close(hf_FastICA_AdvOpt);

 % Use return to avoid reaching the watchoff statement at the end
 % (There used to be a 'break' statement here, but it resulted in
 % errors in more recent version of Matlab -- jarmo)
 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Default'
  
 % set default values to controls
 set(hpm_FastICA_finetune, 'Value',5);
 set(he_FastICA_a1, 'String','1');
 set(he_FastICA_a2, 'String','1');
 set(he_FastICA_myy, 'String','1');
 set(he_FastICA_epsilon, 'String','0.0001');
 set(he_FastICA_maxIterations, 'String','1000');
 set(he_FastICA_sampleSize, 'String','1');
 set(hpm_FastICA_initState, 'Value',1);
 set(hpm_FastICA_displayMode, 'Value',1);
 set(he_FastICA_displayInterval, 'String','1');
 set(hpm_FastICA_verbose, 'Value',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Help'

 gui_help('gui_advc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % switch

watchoff (watchonInFigure);