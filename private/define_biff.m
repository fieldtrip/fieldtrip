% DEFINE_BIFF script to define the known BIFF chunk types 
%
% This function should not be called separately, it is intended to
% be used from withing read_biff and write_biff.
%
% See also READ_BIFF, WRITE_BIFF

% Copyright (C) 2000, Robert Oostenveld
% 
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

CAT_.desc = {'CAT ', 'cat', 'unknown'};
LIST.desc = {'LIST', 'list',    'unknown'};
BIFF.desc = {'BIFF', 'biff',    'group'};

BIFF.SEMG.desc = {'SEMG', 'semg',   'group'};

BIFF.SEMG.VERS.desc = {'VERS', 'version',       'group'};
BIFF.SEMG.DINF.desc = {'DINF', 'data_info',     'group'};
BIFF.SEMG.EXPI.desc = {'EXPI', 'experiment',        'group'};
BIFF.SEMG.MEAS.desc = {'MEAS', 'measurement',       'group'};
BIFF.SEMG.SGRD.desc = {'SGRD', 'semggrid',      'group'};
BIFF.SEMG.PATI.desc = {'PATI', 'patient',       'group'};

BIFF.SEMG.VERS.ID__.desc = {'ID  ', 'id',       'string'};
BIFF.SEMG.VERS.LIC_.desc = {'LIC ', 'license',      'string'};

BIFF.SEMG.DATI.CHAN.desc = {'CHAN', 'unknown'       'unknown'};
BIFF.SEMG.DATI.MODE.desc = {'MODE', 'unknown'       'unknown'};
BIFF.SEMG.DATI.TYPE.desc = {'TYPE', 'unknown'       'unknown'};
BIFF.SEMG.DATI.FS__.desc = {'FS  ', 'unknown'       'unknown'};
BIFF.SEMG.DATI.LABL.desc = {'LABL', 'unknown'       'unknown'};
BIFF.SEMG.DATI.UNIT.desc = {'UNIT', 'unknown'       'unknown'};
BIFF.SEMG.DATI.FILT.desc = {'FILT', 'unknown'       'unknown'};

BIFF.SEMG.EXPI.NAME.desc = {'NAME', 'name',     'string'};
BIFF.SEMG.EXPI.EXAM.desc = {'EXAM', 'examination',  'string'};
BIFF.SEMG.EXPI.EXRC.desc = {'EXRC', 'exercise',     'string'};
BIFF.SEMG.EXPI.RPHS.desc = {'RPHS', 'ref_physician',    'string'};
BIFF.SEMG.EXPI.DATE.desc = {'DATE', 'date',     'unknown'};
BIFF.SEMG.EXPI.TIME.desc = {'TIME', 'time',     'unknown'};

BIFF.SEMG.MEAS.TIME.desc = {'TIME', 'time',     'string'};
BIFF.SEMG.MEAS.ORFN.desc = {'ORFN', 'filename',     'string'};
BIFF.SEMG.MEAS.MUSC.desc = {'MUSC', 'muscle',       'string'};
BIFF.SEMG.MEAS.SIDE.desc = {'SIDE', 'side',     'unknown'};
BIFF.SEMG.MEAS.FNUM.desc = {'FNUM', 'filenumber',   'unknown'};
BIFF.SEMG.MEAS.TEMP.desc = {'TEMP', 'temperature',  'unknown'};
BIFF.SEMG.MEAS.ANOT.desc = {'ANOT', 'annotation',   'string'};
BIFF.SEMG.MEAS.FORC.desc = {'FORC', 'force',        'string'};

BIFF.SEMG.SGRD.NPD_.desc = {'NPD ', 'npd',      'unknown'};
BIFF.SEMG.SGRD.NML_.desc = {'NML ', 'nml',      'unknown'};
BIFF.SEMG.SGRD.IED_.desc = {'IED ', 'ied',      'unknown'};
BIFF.SEMG.SGRD.FIBD.desc = {'FIBD', 'fibre_dir',    'unknown'};
BIFF.SEMG.SGRD.POS_.desc = {'POS ', 'tle_pos',      'unknown'};
BIFF.SEMG.SGRD.CABL.desc = {'CABL', 'cable_pos',    'unknown'};
BIFF.SEMG.SGRD.GRID.desc = {'GRID', 'grid',     'int16vec'};

BIFF.SEMG.PATI.UPID.desc = {'UPID', 'unique_id',    'unknown'};
BIFF.SEMG.PATI.GEND.desc = {'GEND', 'gendre',       'unknown'};
BIFF.SEMG.PATI.BDAT.desc = {'BDAT', 'birthdate',    'unknown'};
BIFF.SEMG.PATI.HPID.desc = {'HPID', 'hospital_id',  'unknown'};
BIFF.SEMG.PATI.ANAM.desc = {'ANAM', 'anamnesis',    'unknown'};

