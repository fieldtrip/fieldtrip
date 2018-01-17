This directory contains some code that assists in accessing the
data of the CTF system in real time. You should compile this code
on your CTF MEG acquisition computer and start it prior to starting
Acq.

This directory contains three versions of "ctf2ft" for historical
reasons. It is recommended that you use the latest version, which
copies the data from CTF shared memory on the acquisition computer
to the FieldTrip buffer. You can use MATLAB on another computer to
read the header, data and events.

The "ctf2ft_v1" is the very first implementation and does not copy
the data to the FieldTrip buffer; it only maintains the CTF shared
memory. It keeps as many data records as possible in shared memory,
but alsways ensures that there are enough empty data records for
Acq to write new data towards. It is to be used with MATLAB running
on the acquisition computer in conjunction with the "ctf_shm"
fileformat option in the low-level reading functions.

For more information please visit http://fieldtriptoolbox.org/development/realtime

This software is free but copyrighted software, distributed
under the terms of the GNU General Public Licence as published by
the Free Software Foundation (either version 2, or at your option
any later version). See the file COPYING for more details.

Copyright (C) 2005-2008, F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, The Netherlands
Copyright (C) 2008-2018, Donders Institute for Brain, Cognition and Behaviour, Nijmegen, The Netherlands

