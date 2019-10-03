Lag Extraction EEGLAB plugin
----------------------------

Copyright (c) 2008-2011 Alexandre Gramfort

This code implements the method detailed in this paper::

    Graph-based variability estimation in single-trial event-related neural responses.
    Gramfort A, Keriven R, Clerc M.
    IEEE Trans Biomed Eng. 2010 May;57(5):1051-61. Epub 2010 Feb 5.

URL : http://www.ncbi.nlm.nih.gov/pubmed/20142163

Notes
-----

This plugin provides a set of functions used to extract lags of
evoked response on single trial M/EEG data

compile mex files using :

>> lagextraction_setup

then run the demo with synthetic and real data :

>> lagextraction_demo.m

This work uses the max flow algorithm provided by Vladimir Kolmogorov :

http://www.adastral.ucl.ac.uk/~vladkolm/software.html

A small documentation is available in 'doc/lag_extraction_demo.html'

License
-------

Copyright (c) 2008-2011 Alexandre Gramfort

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
