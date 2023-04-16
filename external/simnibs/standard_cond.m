function S=standard_cond
% USAGE: S=standard_cond
% 
%        returns a structure with standard conductivity values
%
% A. Thielscher, A. Antunes 2013
% updated 2018; A. Thielscher, G. Saturnino


%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2018 Axel Thielscher, Andre Antunes, Guilherme Saturnino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.




% WM
S(1)=sim_struct('COND');
S(1).name = 'WM';
S(1).value = 0.126;
S(1).descrip = 'brain white matter (from Wagner 2004)';

% GM
S(2)=sim_struct('COND');
S(2).name = 'GM';
S(2).value = 0.275;
S(2).descrip = 'brain gray matter (from Wagner 2004)';

% CSF
S(3)=sim_struct('COND');
S(3).name = 'CSF';
S(3).value = 1.654;
S(3).descrip = 'cerebrospinal fluid (from Wagner 2004)';

% Bone
S(4)=sim_struct('COND');
S(4).name = 'Bone';
S(4).value = 0.010;
S(4).descrip = 'average bone (from Wagner 2004)';

% Scalp
S(5)=sim_struct('COND');
S(5).name = 'Scalp';
S(5).value = 0.465;
S(5).descrip = 'average scalp (from Wagner 2004)';

% Eye balls (vitreous humour)
S(6)=sim_struct('COND');
S(6).name = 'Eye_balls';
S(6).value = 0.5;
S(6).descrip = 'vitreous humour (from Opitz, ..., Thielscher, NI, 2015)';

% Compact bone
S(7)=sim_struct('COND');
S(7).name = 'Compact_bone';
S(7).value = 0.008;
S(7).descrip = 'compact bone (from Opitz, ..., Thielscher, NI, 2015)';

% Spongy bone
S(8)=sim_struct('COND');
S(8).name = 'Spongy_bone';
S(8).value = 0.025;
S(8).descrip = 'spongy bone (from Opitz, ..., Thielscher, NI, 2015)';

% Blood
S(9)=sim_struct('COND');
S(9).name = 'Blood';
S(9).value = 0.6;
S(9).descrip = 'Blood (from Gabriel et al., 2009)';

% Muscle
S(10)=sim_struct('COND');
S(10).name = 'Muscle';
S(10).value = 0.16;
S(10).descrip = 'Muscle (from Gabriel et al., 2009)';

% Rubber
S(100)=sim_struct('COND');
S(100).name = 'Electrode_rubber';
S(100).value = 29.4;
S(100).descrip = 'for tDCS rubber electrodes (for neuroConn rubber: Wacker Elastosil R 570/60 RUSS)';


% Saline
S(500)=sim_struct('COND');
S(500).name = 'Saline';
S(500).value = 1.0;
S(500).descrip = 'for tDCS sponge electrodes';


