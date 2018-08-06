function [elc, lab] = elec1020_locate(pnt, dhk, nas, ini, lpa, rpa, feedback)

% ELEC1020_LOCATE determines 10-20 (20%, 10% and 5%) electrode positions
% on a scalp surface that is described by its surface triangulation

% Copyright (C) 2003, Robert Oostenveld
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

if nargin<7
  feedback = false;
end


% determine the approximate location of the vertex
ori = (lpa+rpa+nas+ini)/4;      % center of head
ver =  cross(rpa-lpa, nas-ini); % orientation
ver = ver /sqrt(norm(ver));     % make correct length
ver = ori + 0.7*ver;            % location from center of head

if feedback
  figure
  ft_plot_mesh(struct('pos', pnt, 'tri', dhk), 'edgecolor', 'none', 'facecolor', 'skin')
  lighting gouraud
  material dull
  lightangle(0, 90);
  alpha 0.9
  ft_plot_mesh(nas, 'vertexsize', 30)
  ft_plot_mesh(lpa, 'vertexsize', 30)
  ft_plot_mesh(ini, 'vertexsize', 30)
  ft_plot_mesh(rpa, 'vertexsize', 30)
  ft_plot_mesh(ver, 'vertexsize', 30)
  grid on
  hold on
  view([1 1 0.5])
end


% point near LPA that is at 50% of left lower contour
[cnt1, cnt2] = elec1020_follow(pnt, dhk, nas, lpa, ini, feedback);
mle = elec1020_fraction(cnt1, cnt2, 0.5);

% point near RPA that is at 50% of right lower contour
[cnt1, cnt2] = elec1020_follow(pnt, dhk, nas, rpa, ini, feedback);
mre = elec1020_fraction(cnt1, cnt2, 0.5);

% determine two points that approximate the vertex
[cnt1, cnt2] = elec1020_follow(pnt, dhk, nas, ver, ini, feedback);
ver1 = elec1020_fraction(cnt1, cnt2, 0.5);
[cnt1, cnt2] = elec1020_follow(pnt, dhk, mle, ver, mre, feedback);
ver2 = elec1020_fraction(cnt1, cnt2, 0.5);

% refined estimate is the average of these two
ver = (ver1+ver2)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start contouring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ant-post contour through vertex
fprintf('constructing vertical ant-post contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, nas, ver, ini, feedback);
Nz   = elec1020_fraction(cnt1, cnt2,  0/20);
NFpz = elec1020_fraction(cnt1, cnt2,  1/20);
Fpz  = elec1020_fraction(cnt1, cnt2,  2/20);
AFpz = elec1020_fraction(cnt1, cnt2,  3/20);
AFz  = elec1020_fraction(cnt1, cnt2,  4/20);
AFFz = elec1020_fraction(cnt1, cnt2,  5/20);
Fz   = elec1020_fraction(cnt1, cnt2,  6/20);
FFCz = elec1020_fraction(cnt1, cnt2,  7/20);
FCz  = elec1020_fraction(cnt1, cnt2,  8/20);
FCCz = elec1020_fraction(cnt1, cnt2,  9/20);
Cz   = elec1020_fraction(cnt1, cnt2, 10/20);
CCPz = elec1020_fraction(cnt1, cnt2, 11/20);
CPz  = elec1020_fraction(cnt1, cnt2, 12/20);
CPPz = elec1020_fraction(cnt1, cnt2, 13/20);
Pz   = elec1020_fraction(cnt1, cnt2, 14/20);
PPOz = elec1020_fraction(cnt1, cnt2, 15/20);
POz  = elec1020_fraction(cnt1, cnt2, 16/20);
POOz = elec1020_fraction(cnt1, cnt2, 17/20);
Oz   = elec1020_fraction(cnt1, cnt2, 18/20);
OIz  = elec1020_fraction(cnt1, cnt2, 19/20);
Iz   = elec1020_fraction(cnt1, cnt2, 20/20);

% left-right through vertex
fprintf('constructing C contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, mle, ver, mre, feedback);
T9   = elec1020_fraction(cnt1, cnt2,  0/20);
T9h  = elec1020_fraction(cnt1, cnt2,  1/20);
T7   = elec1020_fraction(cnt1, cnt2,  2/20);
T7h  = elec1020_fraction(cnt1, cnt2,  3/20);
C5   = elec1020_fraction(cnt1, cnt2,  4/20);
C5h  = elec1020_fraction(cnt1, cnt2,  5/20);
C3   = elec1020_fraction(cnt1, cnt2,  6/20);
C3h  = elec1020_fraction(cnt1, cnt2,  7/20);
C1   = elec1020_fraction(cnt1, cnt2,  8/20);
C1h  = elec1020_fraction(cnt1, cnt2,  9/20);
Cz   = elec1020_fraction(cnt1, cnt2, 10/20);
C2h  = elec1020_fraction(cnt1, cnt2, 11/20);
C2   = elec1020_fraction(cnt1, cnt2, 12/20);
C4h  = elec1020_fraction(cnt1, cnt2, 13/20);
C4   = elec1020_fraction(cnt1, cnt2, 14/20);
C6h  = elec1020_fraction(cnt1, cnt2, 15/20);
C6   = elec1020_fraction(cnt1, cnt2, 16/20);
T8h  = elec1020_fraction(cnt1, cnt2, 17/20);
T8   = elec1020_fraction(cnt1, cnt2, 18/20);
T10h = elec1020_fraction(cnt1, cnt2, 19/20);
T10  = elec1020_fraction(cnt1, cnt2, 20/20);

% horizontal ant-post through T7
fprintf('constructing horizontal left contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, Fpz, T7, Oz, feedback);
Fp1h = elec1020_fraction(cnt1, cnt2,  1/20);
Fp1  = elec1020_fraction(cnt1, cnt2,  2/20);
AFp7 = elec1020_fraction(cnt1, cnt2,  3/20);
AF7  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF7 = elec1020_fraction(cnt1, cnt2,  5/20);
F7   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT7 = elec1020_fraction(cnt1, cnt2,  7/20);
FT7  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT7 = elec1020_fraction(cnt1, cnt2,  9/20);
T7   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP7 = elec1020_fraction(cnt1, cnt2, 11/20);
TP7  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP7 = elec1020_fraction(cnt1, cnt2, 13/20);
P7   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO7 = elec1020_fraction(cnt1, cnt2, 15/20);
PO7  = elec1020_fraction(cnt1, cnt2, 16/20);
POO7 = elec1020_fraction(cnt1, cnt2, 17/20);
O1   = elec1020_fraction(cnt1, cnt2, 18/20);
O1h  = elec1020_fraction(cnt1, cnt2, 19/20);

% horizontal ant-post through T8
fprintf('constructing horizontal right contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, Fpz, T8, Oz, feedback);
Fp2h = elec1020_fraction(cnt1, cnt2,  1/20);
Fp2  = elec1020_fraction(cnt1, cnt2,  2/20);
AFp8 = elec1020_fraction(cnt1, cnt2,  3/20);
AF8  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF8 = elec1020_fraction(cnt1, cnt2,  5/20);
F8   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT8 = elec1020_fraction(cnt1, cnt2,  7/20);
FT8  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT8 = elec1020_fraction(cnt1, cnt2,  9/20);
T8   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP8 = elec1020_fraction(cnt1, cnt2, 11/20);
TP8  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP8 = elec1020_fraction(cnt1, cnt2, 13/20);
P8   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO8 = elec1020_fraction(cnt1, cnt2, 15/20);
PO8  = elec1020_fraction(cnt1, cnt2, 16/20);
POO8 = elec1020_fraction(cnt1, cnt2, 17/20);
O2   = elec1020_fraction(cnt1, cnt2, 18/20);
O2h  = elec1020_fraction(cnt1, cnt2, 19/20);

fprintf('constructing AFp contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, AFp7, AFpz, AFp8, feedback);
AFp7h  = elec1020_fraction(cnt1, cnt2,  1/16);
AFp5   = elec1020_fraction(cnt1, cnt2,  2/16);
AFp5h  = elec1020_fraction(cnt1, cnt2,  3/16);
AFp3   = elec1020_fraction(cnt1, cnt2,  4/16);
AFp3h  = elec1020_fraction(cnt1, cnt2,  5/16);
AFp1   = elec1020_fraction(cnt1, cnt2,  6/16);
AFp1h  = elec1020_fraction(cnt1, cnt2,  7/16);
AFp2h  = elec1020_fraction(cnt1, cnt2,  9/16);
AFp2   = elec1020_fraction(cnt1, cnt2, 10/16);
AFp4h  = elec1020_fraction(cnt1, cnt2, 11/16);
AFp4   = elec1020_fraction(cnt1, cnt2, 12/16);
AFp6h  = elec1020_fraction(cnt1, cnt2, 13/16);
AFp6   = elec1020_fraction(cnt1, cnt2, 14/16);
AFp8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing AF contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, AF7, AFz, AF8, feedback);
AF7h  = elec1020_fraction(cnt1, cnt2,  1/16);
AF5   = elec1020_fraction(cnt1, cnt2,  2/16);
AF5h  = elec1020_fraction(cnt1, cnt2,  3/16);
AF3   = elec1020_fraction(cnt1, cnt2,  4/16);
AF3h  = elec1020_fraction(cnt1, cnt2,  5/16);
AF1   = elec1020_fraction(cnt1, cnt2,  6/16);
AF1h  = elec1020_fraction(cnt1, cnt2,  7/16);
AF2h  = elec1020_fraction(cnt1, cnt2,  9/16);
AF2   = elec1020_fraction(cnt1, cnt2, 10/16);
AF4h  = elec1020_fraction(cnt1, cnt2, 11/16);
AF4   = elec1020_fraction(cnt1, cnt2, 12/16);
AF6h  = elec1020_fraction(cnt1, cnt2, 13/16);
AF6   = elec1020_fraction(cnt1, cnt2, 14/16);
AF8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing AFF contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, AFF7, AFFz, AFF8, feedback);
AFF7h  = elec1020_fraction(cnt1, cnt2,  1/16);
AFF5   = elec1020_fraction(cnt1, cnt2,  2/16);
AFF5h  = elec1020_fraction(cnt1, cnt2,  3/16);
AFF3   = elec1020_fraction(cnt1, cnt2,  4/16);
AFF3h  = elec1020_fraction(cnt1, cnt2,  5/16);
AFF1   = elec1020_fraction(cnt1, cnt2,  6/16);
AFF1h  = elec1020_fraction(cnt1, cnt2,  7/16);
AFF2h  = elec1020_fraction(cnt1, cnt2,  9/16);
AFF2   = elec1020_fraction(cnt1, cnt2, 10/16);
AFF4h  = elec1020_fraction(cnt1, cnt2, 11/16);
AFF4   = elec1020_fraction(cnt1, cnt2, 12/16);
AFF6h  = elec1020_fraction(cnt1, cnt2, 13/16);
AFF6   = elec1020_fraction(cnt1, cnt2, 14/16);
AFF8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing F contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, F7, Fz, F8, feedback);
F7h  = elec1020_fraction(cnt1, cnt2,  1/16);
F5   = elec1020_fraction(cnt1, cnt2,  2/16);
F5h  = elec1020_fraction(cnt1, cnt2,  3/16);
F3   = elec1020_fraction(cnt1, cnt2,  4/16);
F3h  = elec1020_fraction(cnt1, cnt2,  5/16);
F1   = elec1020_fraction(cnt1, cnt2,  6/16);
F1h  = elec1020_fraction(cnt1, cnt2,  7/16);
F2h  = elec1020_fraction(cnt1, cnt2,  9/16);
F2   = elec1020_fraction(cnt1, cnt2, 10/16);
F4h  = elec1020_fraction(cnt1, cnt2, 11/16);
F4   = elec1020_fraction(cnt1, cnt2, 12/16);
F6h  = elec1020_fraction(cnt1, cnt2, 13/16);
F6   = elec1020_fraction(cnt1, cnt2, 14/16);
F8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing FFC contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, FFT7, FFCz, FFT8, feedback);
FFT7h  = elec1020_fraction(cnt1, cnt2,  1/16);
FFC5   = elec1020_fraction(cnt1, cnt2,  2/16);
FFC5h  = elec1020_fraction(cnt1, cnt2,  3/16);
FFC3   = elec1020_fraction(cnt1, cnt2,  4/16);
FFC3h  = elec1020_fraction(cnt1, cnt2,  5/16);
FFC1   = elec1020_fraction(cnt1, cnt2,  6/16);
FFC1h  = elec1020_fraction(cnt1, cnt2,  7/16);
FFC2h  = elec1020_fraction(cnt1, cnt2,  9/16);
FFC2   = elec1020_fraction(cnt1, cnt2, 10/16);
FFC4h  = elec1020_fraction(cnt1, cnt2, 11/16);
FFC4   = elec1020_fraction(cnt1, cnt2, 12/16);
FFC6h  = elec1020_fraction(cnt1, cnt2, 13/16);
FFC6   = elec1020_fraction(cnt1, cnt2, 14/16);
FFT8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing FC contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, FT7, FCz, FT8, feedback);
FT7h  = elec1020_fraction(cnt1, cnt2,  1/16);
FC5   = elec1020_fraction(cnt1, cnt2,  2/16);
FC5h  = elec1020_fraction(cnt1, cnt2,  3/16);
FC3   = elec1020_fraction(cnt1, cnt2,  4/16);
FC3h  = elec1020_fraction(cnt1, cnt2,  5/16);
FC1   = elec1020_fraction(cnt1, cnt2,  6/16);
FC1h  = elec1020_fraction(cnt1, cnt2,  7/16);
FC2h  = elec1020_fraction(cnt1, cnt2,  9/16);
FC2   = elec1020_fraction(cnt1, cnt2, 10/16);
FC4h  = elec1020_fraction(cnt1, cnt2, 11/16);
FC4   = elec1020_fraction(cnt1, cnt2, 12/16);
FC6h  = elec1020_fraction(cnt1, cnt2, 13/16);
FC6   = elec1020_fraction(cnt1, cnt2, 14/16);
FT8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing FCC contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, FTT7, FCCz, FTT8, feedback);
FTT7h  = elec1020_fraction(cnt1, cnt2,  1/16);
FCC5   = elec1020_fraction(cnt1, cnt2,  2/16);
FCC5h  = elec1020_fraction(cnt1, cnt2,  3/16);
FCC3   = elec1020_fraction(cnt1, cnt2,  4/16);
FCC3h  = elec1020_fraction(cnt1, cnt2,  5/16);
FCC1   = elec1020_fraction(cnt1, cnt2,  6/16);
FCC1h  = elec1020_fraction(cnt1, cnt2,  7/16);
FCC2h  = elec1020_fraction(cnt1, cnt2,  9/16);
FCC2   = elec1020_fraction(cnt1, cnt2, 10/16);
FCC4h  = elec1020_fraction(cnt1, cnt2, 11/16);
FCC4   = elec1020_fraction(cnt1, cnt2, 12/16);
FCC6h  = elec1020_fraction(cnt1, cnt2, 13/16);
FCC6   = elec1020_fraction(cnt1, cnt2, 14/16);
FTT8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing CCP contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, TTP7, CCPz, TTP8, feedback);
TTP7h  = elec1020_fraction(cnt1, cnt2,  1/16);
CCP5   = elec1020_fraction(cnt1, cnt2,  2/16);
CCP5h  = elec1020_fraction(cnt1, cnt2,  3/16);
CCP3   = elec1020_fraction(cnt1, cnt2,  4/16);
CCP3h  = elec1020_fraction(cnt1, cnt2,  5/16);
CCP1   = elec1020_fraction(cnt1, cnt2,  6/16);
CCP1h  = elec1020_fraction(cnt1, cnt2,  7/16);
CCP2h  = elec1020_fraction(cnt1, cnt2,  9/16);
CCP2   = elec1020_fraction(cnt1, cnt2, 10/16);
CCP4h  = elec1020_fraction(cnt1, cnt2, 11/16);
CCP4   = elec1020_fraction(cnt1, cnt2, 12/16);
CCP6h  = elec1020_fraction(cnt1, cnt2, 13/16);
CCP6   = elec1020_fraction(cnt1, cnt2, 14/16);
TTP8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing CP contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, TP7, CPz, TP8, feedback);
TP7h  = elec1020_fraction(cnt1, cnt2,  1/16);
CP5   = elec1020_fraction(cnt1, cnt2,  2/16);
CP5h  = elec1020_fraction(cnt1, cnt2,  3/16);
CP3   = elec1020_fraction(cnt1, cnt2,  4/16);
CP3h  = elec1020_fraction(cnt1, cnt2,  5/16);
CP1   = elec1020_fraction(cnt1, cnt2,  6/16);
CP1h  = elec1020_fraction(cnt1, cnt2,  7/16);
CP2h  = elec1020_fraction(cnt1, cnt2,  9/16);
CP2   = elec1020_fraction(cnt1, cnt2, 10/16);
CP4h  = elec1020_fraction(cnt1, cnt2, 11/16);
CP4   = elec1020_fraction(cnt1, cnt2, 12/16);
CP6h  = elec1020_fraction(cnt1, cnt2, 13/16);
CP6   = elec1020_fraction(cnt1, cnt2, 14/16);
TP8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing CPP contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, TPP7, CPPz, TPP8, feedback);
TPP7h  = elec1020_fraction(cnt1, cnt2,  1/16);
CPP5   = elec1020_fraction(cnt1, cnt2,  2/16);
CPP5h  = elec1020_fraction(cnt1, cnt2,  3/16);
CPP3   = elec1020_fraction(cnt1, cnt2,  4/16);
CPP3h  = elec1020_fraction(cnt1, cnt2,  5/16);
CPP1   = elec1020_fraction(cnt1, cnt2,  6/16);
CPP1h  = elec1020_fraction(cnt1, cnt2,  7/16);
CPP2h  = elec1020_fraction(cnt1, cnt2,  9/16);
CPP2   = elec1020_fraction(cnt1, cnt2, 10/16);
CPP4h  = elec1020_fraction(cnt1, cnt2, 11/16);
CPP4   = elec1020_fraction(cnt1, cnt2, 12/16);
CPP6h  = elec1020_fraction(cnt1, cnt2, 13/16);
CPP6   = elec1020_fraction(cnt1, cnt2, 14/16);
TPP8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing P contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, P7, Pz, P8, feedback);
P7h  = elec1020_fraction(cnt1, cnt2,  1/16);
P5   = elec1020_fraction(cnt1, cnt2,  2/16);
P5h  = elec1020_fraction(cnt1, cnt2,  3/16);
P3   = elec1020_fraction(cnt1, cnt2,  4/16);
P3h  = elec1020_fraction(cnt1, cnt2,  5/16);
P1   = elec1020_fraction(cnt1, cnt2,  6/16);
P1h  = elec1020_fraction(cnt1, cnt2,  7/16);
P2h  = elec1020_fraction(cnt1, cnt2,  9/16);
P2   = elec1020_fraction(cnt1, cnt2, 10/16);
P4h  = elec1020_fraction(cnt1, cnt2, 11/16);
P4   = elec1020_fraction(cnt1, cnt2, 12/16);
P6h  = elec1020_fraction(cnt1, cnt2, 13/16);
P6   = elec1020_fraction(cnt1, cnt2, 14/16);
P8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing PPO contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, PPO7, PPOz, PPO8, feedback);
PPO7h  = elec1020_fraction(cnt1, cnt2,  1/16);
PPO5   = elec1020_fraction(cnt1, cnt2,  2/16);
PPO5h  = elec1020_fraction(cnt1, cnt2,  3/16);
PPO3   = elec1020_fraction(cnt1, cnt2,  4/16);
PPO3h  = elec1020_fraction(cnt1, cnt2,  5/16);
PPO1   = elec1020_fraction(cnt1, cnt2,  6/16);
PPO1h  = elec1020_fraction(cnt1, cnt2,  7/16);
PPO2h  = elec1020_fraction(cnt1, cnt2,  9/16);
PPO2   = elec1020_fraction(cnt1, cnt2, 10/16);
PPO4h  = elec1020_fraction(cnt1, cnt2, 11/16);
PPO4   = elec1020_fraction(cnt1, cnt2, 12/16);
PPO6h  = elec1020_fraction(cnt1, cnt2, 13/16);
PPO6   = elec1020_fraction(cnt1, cnt2, 14/16);
PPO8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing PO contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, PO7, POz, PO8, feedback);
PO7h  = elec1020_fraction(cnt1, cnt2,  1/16);
PO5   = elec1020_fraction(cnt1, cnt2,  2/16);
PO5h  = elec1020_fraction(cnt1, cnt2,  3/16);
PO3   = elec1020_fraction(cnt1, cnt2,  4/16);
PO3h  = elec1020_fraction(cnt1, cnt2,  5/16);
PO1   = elec1020_fraction(cnt1, cnt2,  6/16);
PO1h  = elec1020_fraction(cnt1, cnt2,  7/16);
PO2h  = elec1020_fraction(cnt1, cnt2,  9/16);
PO2   = elec1020_fraction(cnt1, cnt2, 10/16);
PO4h  = elec1020_fraction(cnt1, cnt2, 11/16);
PO4   = elec1020_fraction(cnt1, cnt2, 12/16);
PO6h  = elec1020_fraction(cnt1, cnt2, 13/16);
PO6   = elec1020_fraction(cnt1, cnt2, 14/16);
PO8h  = elec1020_fraction(cnt1, cnt2, 15/16);

fprintf('constructing POO contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, POO7, POOz, POO8, feedback);
POO7h  = elec1020_fraction(cnt1, cnt2,  1/16);
POO5   = elec1020_fraction(cnt1, cnt2,  2/16);
POO5h  = elec1020_fraction(cnt1, cnt2,  3/16);
POO3   = elec1020_fraction(cnt1, cnt2,  4/16);
POO3h  = elec1020_fraction(cnt1, cnt2,  5/16);
POO1   = elec1020_fraction(cnt1, cnt2,  6/16);
POO1h  = elec1020_fraction(cnt1, cnt2,  7/16);
POO2h  = elec1020_fraction(cnt1, cnt2,  9/16);
POO2   = elec1020_fraction(cnt1, cnt2, 10/16);
POO4h  = elec1020_fraction(cnt1, cnt2, 11/16);
POO4   = elec1020_fraction(cnt1, cnt2, 12/16);
POO6h  = elec1020_fraction(cnt1, cnt2, 13/16);
POO6   = elec1020_fraction(cnt1, cnt2, 14/16);
POO8h  = elec1020_fraction(cnt1, cnt2, 15/16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start contouring the low electrode locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low horizontal ant-post through T9
fprintf('constructing low horizontal left contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, Nz, T9, Iz, feedback);
AFp9 = elec1020_fraction(cnt1, cnt2,  3/20);
AF9  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF9 = elec1020_fraction(cnt1, cnt2,  5/20);
F9   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT9 = elec1020_fraction(cnt1, cnt2,  7/20);
FT9  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT9 = elec1020_fraction(cnt1, cnt2,  9/20);
T9   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP9 = elec1020_fraction(cnt1, cnt2, 11/20);
TP9  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP9 = elec1020_fraction(cnt1, cnt2, 13/20);
P9   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO9 = elec1020_fraction(cnt1, cnt2, 15/20);
PO9  = elec1020_fraction(cnt1, cnt2, 16/20);
POO9 = elec1020_fraction(cnt1, cnt2, 17/20);
I1   = elec1020_fraction(cnt1, cnt2, 18/20);
I1h  = elec1020_fraction(cnt1, cnt2, 19/20);

[cnt1, cnt2] = elec1020_follow(pnt, dhk, NFpz, T9h, OIz, feedback);
AFp9h = elec1020_fraction(cnt1, cnt2,  3/20);
AF9h  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF9h = elec1020_fraction(cnt1, cnt2,  5/20);
F9h   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT9h = elec1020_fraction(cnt1, cnt2,  7/20);
FT9h  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT9h = elec1020_fraction(cnt1, cnt2,  9/20);
T9h   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP9h = elec1020_fraction(cnt1, cnt2, 11/20);
TP9h  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP9h = elec1020_fraction(cnt1, cnt2, 13/20);
P9h   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO9h = elec1020_fraction(cnt1, cnt2, 15/20);
PO9h  = elec1020_fraction(cnt1, cnt2, 16/20);
POO9h = elec1020_fraction(cnt1, cnt2, 17/20);
OI1   = elec1020_fraction(cnt1, cnt2, 18/20);
OI1h  = elec1020_fraction(cnt1, cnt2, 19/20);

% low horizontal ant-post through T10
fprintf('constructing low horizontal right contour\n');
[cnt1, cnt2] = elec1020_follow(pnt, dhk, Nz, T10, Iz, feedback);
AFp10 = elec1020_fraction(cnt1, cnt2,  3/20);
AF10  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF10 = elec1020_fraction(cnt1, cnt2,  5/20);
F10   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT10 = elec1020_fraction(cnt1, cnt2,  7/20);
FT10  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT10 = elec1020_fraction(cnt1, cnt2,  9/20);
T10   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP10 = elec1020_fraction(cnt1, cnt2, 11/20);
TP10  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP10 = elec1020_fraction(cnt1, cnt2, 13/20);
P10   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO10 = elec1020_fraction(cnt1, cnt2, 15/20);
PO10  = elec1020_fraction(cnt1, cnt2, 16/20);
POO10 = elec1020_fraction(cnt1, cnt2, 17/20);
I2    = elec1020_fraction(cnt1, cnt2, 18/20);
I2h   = elec1020_fraction(cnt1, cnt2, 19/20);

[cnt1, cnt2] = elec1020_follow(pnt, dhk, NFpz, T10h, OIz, feedback);
AFp10h = elec1020_fraction(cnt1, cnt2,  3/20);
AF10h  = elec1020_fraction(cnt1, cnt2,  4/20);
AFF10h = elec1020_fraction(cnt1, cnt2,  5/20);
F10h   = elec1020_fraction(cnt1, cnt2,  6/20);
FFT10h = elec1020_fraction(cnt1, cnt2,  7/20);
FT10h  = elec1020_fraction(cnt1, cnt2,  8/20);
FTT10h = elec1020_fraction(cnt1, cnt2,  9/20);
T10h   = elec1020_fraction(cnt1, cnt2, 10/20);
TTP10h = elec1020_fraction(cnt1, cnt2, 11/20);
TP10h  = elec1020_fraction(cnt1, cnt2, 12/20);
TPP10h = elec1020_fraction(cnt1, cnt2, 13/20);
P10h   = elec1020_fraction(cnt1, cnt2, 14/20);
PPO10h = elec1020_fraction(cnt1, cnt2, 15/20);
PO10h  = elec1020_fraction(cnt1, cnt2, 16/20);
POO10h = elec1020_fraction(cnt1, cnt2, 17/20);
OI2    = elec1020_fraction(cnt1, cnt2, 18/20);
OI2h   = elec1020_fraction(cnt1, cnt2, 19/20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect all the computed electrode locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lab = {
  'LPA'
  'RPA'
  'NAS'
  'INI'
  'Nz'
  'Fp1'
  'Fpz'
  'Fp2'
  'AF9'
  'AF7'
  'AF5'
  'AF3'
  'AF1'
  'AFz'
  'AF2'
  'AF4'
  'AF6'
  'AF8'
  'AF10'
  'F9'
  'F7'
  'F5'
  'F3'
  'F1'
  'Fz'
  'F2'
  'F4'
  'F6'
  'F8'
  'F10'
  'FT9'
  'FT7'
  'FC5'
  'FC3'
  'FC1'
  'FCz'
  'FC2'
  'FC4'
  'FC6'
  'FT8'
  'FT10'
  'T9'
  'T7'
  'C5'
  'C3'
  'C1'
  'Cz'
  'C2'
  'C4'
  'C6'
  'T8'
  'T10'
  'TP9'
  'TP7'
  'CP5'
  'CP3'
  'CP1'
  'CPz'
  'CP2'
  'CP4'
  'CP6'
  'TP8'
  'TP10'
  'P9'
  'P7'
  'P5'
  'P3'
  'P1'
  'Pz'
  'P2'
  'P4'
  'P6'
  'P8'
  'P10'
  'PO9'
  'PO7'
  'PO5'
  'PO3'
  'PO1'
  'POz'
  'PO2'
  'PO4'
  'PO6'
  'PO8'
  'PO10'
  'O1'
  'Oz'
  'O2'
  'I1'
  'Iz'
  'I2'
  'AFp9h'
  'AFp7h'
  'AFp5h'
  'AFp3h'
  'AFp1h'
  'AFp2h'
  'AFp4h'
  'AFp6h'
  'AFp8h'
  'AFp10h'
  'AFF9h'
  'AFF7h'
  'AFF5h'
  'AFF3h'
  'AFF1h'
  'AFF2h'
  'AFF4h'
  'AFF6h'
  'AFF8h'
  'AFF10h'
  'FFT9h'
  'FFT7h'
  'FFC5h'
  'FFC3h'
  'FFC1h'
  'FFC2h'
  'FFC4h'
  'FFC6h'
  'FFT8h'
  'FFT10h'
  'FTT9h'
  'FTT7h'
  'FCC5h'
  'FCC3h'
  'FCC1h'
  'FCC2h'
  'FCC4h'
  'FCC6h'
  'FTT8h'
  'FTT10h'
  'TTP9h'
  'TTP7h'
  'CCP5h'
  'CCP3h'
  'CCP1h'
  'CCP2h'
  'CCP4h'
  'CCP6h'
  'TTP8h'
  'TTP10h'
  'TPP9h'
  'TPP7h'
  'CPP5h'
  'CPP3h'
  'CPP1h'
  'CPP2h'
  'CPP4h'
  'CPP6h'
  'TPP8h'
  'TPP10h'
  'PPO9h'
  'PPO7h'
  'PPO5h'
  'PPO3h'
  'PPO1h'
  'PPO2h'
  'PPO4h'
  'PPO6h'
  'PPO8h'
  'PPO10h'
  'POO9h'
  'POO7h'
  'POO5h'
  'POO3h'
  'POO1h'
  'POO2h'
  'POO4h'
  'POO6h'
  'POO8h'
  'POO10h'
  'OI1h'
  'OI2h'
  'Fp1h'
  'Fp2h'
  'AF9h'
  'AF7h'
  'AF5h'
  'AF3h'
  'AF1h'
  'AF2h'
  'AF4h'
  'AF6h'
  'AF8h'
  'AF10h'
  'F9h'
  'F7h'
  'F5h'
  'F3h'
  'F1h'
  'F2h'
  'F4h'
  'F6h'
  'F8h'
  'F10h'
  'FT9h'
  'FT7h'
  'FC5h'
  'FC3h'
  'FC1h'
  'FC2h'
  'FC4h'
  'FC6h'
  'FT8h'
  'FT10h'
  'T9h'
  'T7h'
  'C5h'
  'C3h'
  'C1h'
  'C2h'
  'C4h'
  'C6h'
  'T8h'
  'T10h'
  'TP9h'
  'TP7h'
  'CP5h'
  'CP3h'
  'CP1h'
  'CP2h'
  'CP4h'
  'CP6h'
  'TP8h'
  'TP10h'
  'P9h'
  'P7h'
  'P5h'
  'P3h'
  'P1h'
  'P2h'
  'P4h'
  'P6h'
  'P8h'
  'P10h'
  'PO9h'
  'PO7h'
  'PO5h'
  'PO3h'
  'PO1h'
  'PO2h'
  'PO4h'
  'PO6h'
  'PO8h'
  'PO10h'
  'O1h'
  'O2h'
  'I1h'
  'I2h'
  'AFp9'
  'AFp7'
  'AFp5'
  'AFp3'
  'AFp1'
  'AFpz'
  'AFp2'
  'AFp4'
  'AFp6'
  'AFp8'
  'AFp10'
  'AFF9'
  'AFF7'
  'AFF5'
  'AFF3'
  'AFF1'
  'AFFz'
  'AFF2'
  'AFF4'
  'AFF6'
  'AFF8'
  'AFF10'
  'FFT9'
  'FFT7'
  'FFC5'
  'FFC3'
  'FFC1'
  'FFCz'
  'FFC2'
  'FFC4'
  'FFC6'
  'FFT8'
  'FFT10'
  'FTT9'
  'FTT7'
  'FCC5'
  'FCC3'
  'FCC1'
  'FCCz'
  'FCC2'
  'FCC4'
  'FCC6'
  'FTT8'
  'FTT10'
  'TTP9'
  'TTP7'
  'CCP5'
  'CCP3'
  'CCP1'
  'CCPz'
  'CCP2'
  'CCP4'
  'CCP6'
  'TTP8'
  'TTP10'
  'TPP9'
  'TPP7'
  'CPP5'
  'CPP3'
  'CPP1'
  'CPPz'
  'CPP2'
  'CPP4'
  'CPP6'
  'TPP8'
  'TPP10'
  'PPO9'
  'PPO7'
  'PPO5'
  'PPO3'
  'PPO1'
  'PPOz'
  'PPO2'
  'PPO4'
  'PPO6'
  'PPO8'
  'PPO10'
  'POO9'
  'POO7'
  'POO5'
  'POO3'
  'POO1'
  'POOz'
  'POO2'
  'POO4'
  'POO6'
  'POO8'
  'POO10'
  'OI1'
  'OIz'
  'OI2'
  'T3' % this is now called T7
  'T4' % this is now called T8
  'T5' % this is now called P7
  'T6' % this is now called P8
  'M1' % left mastoid
  'M2' % right mastoid
  'A1' % left ear lobe
  'A2' % right ear lobe
  };

% allow for upper case electrode position labels
NAS = nas;
INI = ini;
LPA = lpa;
RPA = rpa;

% it would be possible to assign locations to these old locations
% but there are no locations determined for M1/2 and A1/2
if false
  T3 = T7;
  T4 = T8;
  T5 = P7;
  T6 = P8;
end

% assign the known local electrode positions to the output matrix
nlab = numel(lab);
elc = ones(nlab,3) * nan;
for i=1:nlab
  if exist(lab{i}, 'var')
    eval(sprintf('elc(%d,:) = %s;', i, lab{i}));
  else
    fprintf('not placing electrode %s\n', lab{i});
  end
end

% remove unknown electrode positions
sel = ~isnan(elc(:,1));
elc = elc(sel, :);
lab = lab(sel);

if feedback
  elec = [];
  elec.elecpos = elc;
  elec.label = lab;
  ft_plot_sens(elec)
end
