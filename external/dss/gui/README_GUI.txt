This is the readme file for the GUI of the DSS Matlab Package.
If you are looking for information about the actual DSS
functionality, refer to README.txt in the root directory.

$Id: README_GUI.txt,v 1.8 2005/04/20 10:19:24 kosti Exp $


------
About:
------

  Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
  Distributed by Laboratory of Computer and Information Science,
  Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


--------------
Using the GUI:
--------------

Start the DSS GUI by simply typing:
>> dss_gui
.

States, signals, and parameter strings can be given as a parameter
in any order. No more than two parameters can be given.
Currently there is only one parameter string which can be given.
String 'pro' will make all panels visible at start up.

Examples:

% Only state
>> dss_gui(state)

% Only signals
>> dss_gui(X)

% Only pro mode
>> dss_gui('pro')

% State and pro-mode
>> dss_gui('pro', state)
or
>> dss_gui(state, 'pro')

% Signals and pro-mode
>> dss_gui('pro', X)
or
>> dss_gui(X, 'pro')


Using the GUI is quite straightforward. There are three
panels: 'Load data', 'Preprocessing' and 'Choose parameters'.

Use 'Load data' -panel to load a set of input signals or
a state structure. Data can be loaded from workspace or from
a file. Data can be plotted and transposed in this panel.

Use 'Preprocessing' -panel to choose a preprocessing function
for the loaded data. Input data dimensions can also be reduced
here. Preprocessed data can be plotted from this panel.
If preprocessed data is plotted, preprocessing is not calculated
again when DSS is started.

Use 'Choose parameters' -panel to choose a suitable set of
parameters. The chosen algorithmic approach affects other
selections in this panel (available denoising functions) and
in advanced options panel (available alpha/beta/gamma functions).
Changing the denoising function will disable any previously
selected beta function. Changing the approach will disable any
previously selected alpha, beta and gamma functions.

Things to note:
Reset buttons reset selections to system defaults if signals
are loaded, and to state defaults if a state is loaded.
System reset button always resets to system defaults. 

Access the actual help files from the GUI's help button for
more information.


------
Files:
------

The GUI is made with Matlab's GUIDE editor, so there are
two of most files. One of the files is an m-file containing
the callback functions and other functionality and the other
one is a fig-file containing the visual layout and the graphical
objects.

dss_gui.m               Functionality of the main window. The GUI starts 
                        with this one.
dss_gui.fig             -

dss_gui_browse.m        Functionality of the file/workspace browsing.
dss_gui_browse.fig      -

dss_gui_advOptions.m	Functionality of the advanced options window.
dss_gui_advOptions.fig 	-

dss_gui_funcParams.m	Functionality of the function parameters window.
dss_gui_funcParams.fig	-

dss_gui_insertFunc.m	Functionality of the insert function window.
dss_gui_insertFunc.fig	-

dss_gui_modalDlg.m      Functionality for a general purpose Yes/No
                        question window.
dss_gui_modalDlg.fig	-

dss_gui_save.m          Functionality of the save results window.
dss_gui_save.fig        -

dss_gui_about.m         Functionality of the about window.
dss_gui_about.fig       -

dss_gui_help.m          File containing the help texts.


---------------
Known problems:
---------------

Result reporting is not completed. (Very simple
at the moment.)

Documentation is not completed.

If a state with whitened data is given to DSS, wdim and sdim
in parameters are ignored. -> Not possible to reduce dimensions.