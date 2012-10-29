function factor = scalingfactor(old, new)

% SCALINGFACTOR determines the scaling factor from old to new units.
%
% Use as
%   factor = scalingfactor(old, new)
% where old and new are strings with the units.
%
% For example
%   scalingfactor('m', 'cm')          % returns 100
%   scalingfactor('V', 'uV')          % returns 1000
%   scalingfactor('T/cm', 'fT/m')     % returns 10^15 divided by 10^-2, which is 10^17
%
% The following fundamental units are supported
%   metre       m   length  l (a lowercase L), x, r L
%   kilogram    kg  mass    m   M
%   second      s   time    t   T
%   ampere      A   electric current    I (an uppercase i)  I
%   kelvin      K   thermodynamic temperature   T   #
%   mole        mol amount of substance n   N
%   candela     cd  luminous intensity  Iv (an uppercase i with lowercase non-italicized v subscript)   J
%
% The following derived units are supported
%   hertz       Hz  frequency   1/s T-1
%   radian      rad angle   m/m dimensionless
%   steradian   sr  solid angle m2/m2   dimensionless
%   newton      N   force, weight   kg#m/s2 M#L#T-2
%   pascal      Pa  pressure, stress    N/m2    M#L-1#T-2
%   joule       J   energy, work, heat  N#m = C#V = W#s M#L2#T-2
%   coulomb     C   electric charge or quantity of electricity  s#A T#I
%   volt        V   voltage, electrical potential difference, electromotive force   W/A = J/C   M#L2#T-3#I-1
%   farad       F   electric capacitance    C/V M-1#L-2#T4#I2
%   siemens     S   electrical conductance  1/# = A/V   M-1#L-2#T3#I2
%   weber       Wb  magnetic flux   J/A M#L2#T-2#I-1
%   tesla       T   magnetic field strength V#s/m2 = Wb/m2 = N/(A#m)    M#T-2#I-1
%   henry       H   inductance  V#s/A = Wb/A    M#L2#T-2#I-2
%   lumen       lm  luminous flux   cd#sr   J
%   lux         lx  illuminance lm/m2   L-2#J
%   becquerel   Bq  radioactivity (decays per unit time)    1/s T-1
%   gray        Gy  absorbed dose (of ionizing radiation)   J/kg    L2#T-2
%   sievert     Sv  equivalent dose (of ionizing radiation) J/kg    L2#T-2
%   katal       kat catalytic activity  mol/s   T-1#N
%
% The following derived units are not supported due to potential confusion
% between their ascii character representation
%   ohm             #   electric resistance, impedance, reactance   V/A M#L2#T-3#I-2
%   watt            W   power, radiant flux J/s = V#A   M#L2#T-3
%   degree Celsius	°C	temperature relative to 273.15 K	K	?
%
% See also http://en.wikipedia.org/wiki/International_System_of_Units

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% convert the old and new to SI units

m   = 1;
g   = 0.001; % rather than representing kg, here g gets represented, the k is added further down
s   = 1;
A   = 1;
K   = 1;
mol = 1;
cd  = 1;

Hz  = 1;
rad = 1;
sr  = 1;
N   = 1;
Pa  = 1;
J   = 1;
C   = 1;
V   = 1;
F   = 1;
S   = 1;
Wb  = 1;
T   = 1;
H   = 1;
lm  = 1;
lx  = 1;
Bq  = 1;
Gy  = 1;
Sv  = 1;
kat = 1;

unit = {'m' 'g' 's' 'A' 'K' 'mol' 'cd' 'Hz' 'rad' 'sr' 'N' 'Pa' 'J' 'C' 'V' 'F' 'S' 'Wb' 'T' 'H' 'lm' 'lx' 'Bq' 'Gy' 'Sv' 'kat'};

% deal with the possible prefixes
for i=1:length(unit)
  eval(sprintf('d%s = 1e-1  * %s;', unit{i}, unit{i}));
  eval(sprintf('c%s = 1e-2  * %s;', unit{i}, unit{i}));
  eval(sprintf('m%s = 1e-3  * %s;', unit{i}, unit{i}));
  eval(sprintf('u%s = 1e-6  * %s;', unit{i}, unit{i})); % note that u is used for micro
  eval(sprintf('n%s = 1e-9  * %s;', unit{i}, unit{i}));
  eval(sprintf('p%s = 1e-12 * %s;', unit{i}, unit{i}));
  eval(sprintf('f%s = 1e-15 * %s;', unit{i}, unit{i}));
  eval(sprintf('a%s = 1e-18 * %s;', unit{i}, unit{i}));
  eval(sprintf('z%s = 1e-21 * %s;', unit{i}, unit{i}));
  eval(sprintf('y%s = 1e-24 * %s;', unit{i}, unit{i}));
  
  eval(sprintf('da%s = 1e1  * %s;', unit{i}, unit{i}));
  eval(sprintf('h%s  = 1e2  * %s;', unit{i}, unit{i}));
  eval(sprintf('k%s  = 1e3  * %s;', unit{i}, unit{i}));
  eval(sprintf('M%s  = 1e6  * %s;', unit{i}, unit{i}));
  eval(sprintf('G%s  = 1e9  * %s;', unit{i}, unit{i}));
  eval(sprintf('T%s  = 1e12 * %s;', unit{i}, unit{i}));
  eval(sprintf('P%s  = 1e15 * %s;', unit{i}, unit{i}));
  eval(sprintf('E%s  = 1e18 * %s;', unit{i}, unit{i}));
  eval(sprintf('Z%s  = 1e21 * %s;', unit{i}, unit{i}));
  eval(sprintf('Y%s  = 1e24 * %s;', unit{i}, unit{i}));
end

eval(sprintf('old2si = %s;', old));
eval(sprintf('new2si = %s;', new));

factor = old2si/new2si;
