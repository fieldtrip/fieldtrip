function dss_gui_help()
% DSS_GUI_HELP
% This file includes all the text behind the help button of the DSS GUI.
% This file is used by other DSS files. Calling this function will
% display all the help texts in Matlab's help browser.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% $Id: dss_gui_help.m,v 1.8 2005/04/20 10:19:24 kosti Exp $

% Help topics.
% If cell array notation is used, rows don't have to be equally long.

title1 = 'DSS GUI: General information';
help_text1 = {...
'This is the help file for the DSS MATLAB package GUI.                  '
'Also see README_GUI.txt and README.txt for more information about      '
'the DSS and the DSS GUI.                                               '
'                                                                       '
'DSS GUI can be started by simply typing:                               '
'>> dss_gui                                                             '
'                                                                       '
'States, signals, and parameter strings can be given as a parameter     '
'in any order. No more than two parameters can be given.                '
'Currently there is only one parameter string which can be given.       '
'String "pro" will make all panels visible at start up.                 '
'                                                                       '
'Using the GUI is quite straightforward. There are three'
'panels: "Load data", "Preprocessing" and "Choose parameters".'
'So first you load some data. Then you choose some parameters.'
'DSS start with the "Start DSS" button. After that results can be saved'
'and examined under the "Results" button.'
'See later sections for more information on each panel.'
''
'Things to note:                                                        '
'Reset buttons reset selections to system defaults if signals           '
'are loaded, and to state defaults if a state is loaded.                '
'System reset button always resets to system defaults.                  '};

title2 = 'Panel 1: Load data';
help_text2 = {...
'Use "Load data" -panel to load a set of input signals or               '
'a state structure. Data can be loaded from workspace or from           '
'a file. Data can be plotted and transposed in this panel.              '
'                                                                       '
'When a state is loaded, the selections in the GUI are changed to       '
'reflect the parameters in the state. If a signal set is loaded to      '
'overwrite a previous state, settings are changed to defaults.          '
'Signals overwriting signals do not change settings.                    '};
    
title3 = 'Panel 2: Preprocessing';
help_text3 = {...
'Use "Preprocessing" -panel to choose a preprocessing function          '
'for the loaded data. Input data dimensions can also be reduced         '
'here. Preprocessed data can be plotted from this panel.                '
'If preprocessed data is plotted, preprocessing is not calculated       '
'again when DSS is started.                                             '
'                                                                       '
'Currently there is only one preprocessing function available with      '
'the package. It is Default whiten and it spheres the data...           '
'Custom preprocessing functions can be made.                            '};

title4 = 'Panel 3: Parameters';
help_text4 = {...
'Use "Choose parameters" -panel to choose a suitable set of'
'parameters. The chosen algorithmic approach affects other'
'selections in this panel (available denoising functions) and'
'in advanced options panel (available alpha/beta/gamma functions).'
'Chosen denoising function also affects available beta functions.'
''
'Changing the denoising function will disable any previously'
'selected beta function. Changing the approach will disable any'
'previously selected alpha, beta and gamma functions.'
''
'Descriptions of different parameters: ...'};

title5 = 'Advanced options';
help_text5 = {...
'Advanced options window is opened from panel 3. It contains some less'
'frequently used parameters.'
''
'Verbosity affects the amount of command line output the DSS generates.'
'If interruptability is selected, calculations can be interrupted from'
'a separate window that appears during the runs.'
'Alpha/Beta/Gamma...'};

    
% Show entries inside the cell notation in the Matlab help browser with
% the entry of the second parameter in view.
helpwin({title1 help_text1;title2 help_text2;title3 help_text3;title4 help_text4;...
    title5 help_text5;title6 help_text6}, title1);