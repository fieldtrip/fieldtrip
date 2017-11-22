function factor = ft_scalingfactor(old, new)

% FT_SCALINGFACTOR determines the scaling factor from old to new units, i.e. it
% returns a number with which the data in the old units needs to be multiplied
% to get it expressed in the new units. 
%
% Use as
%   factor = ft_scalingfactor(old, new)
% where old and new are strings that specify the units.
%
% For example
%   ft_scalingfactor('m', 'cm')          % returns 100
%   ft_scalingfactor('V', 'uV')          % returns 1000
%   ft_scalingfactor('T/cm', 'fT/m')     % returns 10^15 divided by 10^-2, which is 10^17
%   ft_scalingfactor('cm^2', 'mm^2')     % returns 100
%   ft_scalingfactor('1/ms', 'Hz')       % returns 1000
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
% The following alternative units are supported
%   inch        inch  length
%   feet        feet  length
%   gauss       gauss magnetic field strength
%
% The following derived units are not supported due to potential confusion
% between their ascii character representation
%   ohm             #   electric resistance, impedance, reactance   V/A M#L2#T-3#I-2
%   watt            W   power, radiant flux J/s = V#A   M#L2#T-3
%   degree Celsius	?C	temperature relative to 273.15 K	K	?
%
% See also http://en.wikipedia.org/wiki/International_System_of_Units

% Copyright (C) 2012-2015, Robert Oostenveld
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

persistent previous_old previous_new previous_factor

if isempty(previous_old)
  previous_old = {};
  previous_new = {};
  previous_factor = [];
end

if ~isequal(class(old), class(new))
  ft_error('the input arguments should be of the same class');
end

if iscell(old)
  factor = cellfun(@ft_scalingfactor, old(:), new(:));
  return
end

cachehit = strcmp(old, previous_old) & strcmp(new, previous_new);
if any(cachehit)
  factor = previous_factor(cachehit);
  return
end

if isequal(old, new)
  % this applies regardless of the units
  factor = 1;
  return
end

unit = {'m' 'g' 's' 'A' 'K' 'mol' 'cd' 'Hz' 'rad' 'sr' 'N' 'Pa' 'J' 'C' 'V' 'F' 'S' 'Wb' 'T' 'H' 'lm' 'lx' 'Bq' 'Gy' 'Sv' 'kat' 'inch' 'feet' 'gauss'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following section pertains to checking that units are compatible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unknown = -1;

% each of the fundamental units is represented by a prime number
m   = 2;
kg  = 3; g = 3; % besides representing kg, also represent g to facilitate the processing further down
s   = 5;
A   = 7;
K   = 11;
mol = 13;
cd  = 17;

% each of the derives units is represented by a product and/or ratio of prime numbers
Hz  = 1/s;
rad = nan; % this is dimensionless, cannot be converted
sr  = nan; % this is dimensionless, cannot be converted
N   = kg*m/s^2;
Pa  = N/m^2;
J   = N*m;
C   = s*A;
V   = J/C;
F   = C/V;
S   = A/V;
Wb  = J/A;
T   = Wb/m^2;
H   = V*s/A;
lm  = cd*sr;
lx  = lm/m^2;
Bq  = 1/s;
Gy  = J/kg;
Sv  = J/kg;
kat = mol/s;

% these are the same units as more fundamental units, and only differ with respect to multiplying them with a constant
inch = m;     % it is actually 0.0254 m
feet = inch;  % it is actually 12 inch
gauss = T;    % it is actually 10^-4 T

% deal with all possible prefixes
for i=1:length(unit)
  eval(sprintf('d%s = %s;', unit{i}, unit{i}));
  eval(sprintf('c%s = %s;', unit{i}, unit{i}));
  eval(sprintf('m%s = %s;', unit{i}, unit{i}));
  eval(sprintf('u%s = %s;', unit{i}, unit{i})); % note that u is used for micro
  eval(sprintf('n%s = %s;', unit{i}, unit{i}));
  eval(sprintf('p%s = %s;', unit{i}, unit{i}));
  eval(sprintf('f%s = %s;', unit{i}, unit{i}));
  eval(sprintf('a%s = %s;', unit{i}, unit{i}));
  eval(sprintf('z%s = %s;', unit{i}, unit{i}));
  eval(sprintf('y%s = %s;', unit{i}, unit{i}));

  eval(sprintf('da%s = %s;', unit{i}, unit{i}));
  eval(sprintf('h%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('k%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('M%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('G%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('T%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('P%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('E%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('Z%s  = %s;', unit{i}, unit{i}));
  eval(sprintf('Y%s  = %s;', unit{i}, unit{i}));
end

eval(sprintf('oldunit = %s;', old));
eval(sprintf('newunit = %s;', new));
if ~isequal(oldunit, newunit)
  ft_error('cannot convert %s to %s', old, new);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following pertains to determining the scaling factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fundamental units
m   = 1;
kg  = 1; g = 0.001; % besides representing kg, also represent g to facilitate the processing further down
s   = 1;
A   = 1;
K   = 1;
mol = 1;
cd  = 1;

% derived units
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

% these are the same units as more fundamental units, and only differ with respect to multiplying them with a constant
inch  = 0.0254;
feet  = 0.0254*12;
gauss = 1e-4;

% deal with all possible prefixes
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

% remember the input args and the result, this will speed up the next call if the input is the same
previous_old{end+1}    = old;
previous_new{end+1}    = new;
previous_factor(end+1) = factor;
