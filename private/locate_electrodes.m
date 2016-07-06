function [elc, lab] = locate_electrodes(pnt, dhk, nas, ini, lpa, rpa)

% LOCATE_ELECTORDES determines 10-20 (10% and 5%) electrode positions on a head
% that is described by its surface triangulation
%
% [elc, lab] = locate_electrodes(pnt, dhk, nas, ini, lpa, rpa)

% Copyright (C) 2003, Robert Oostenveld

% determine the direction where the vertex is
ver = (lpa+rpa+nas+ini)/4 + cross(rpa-lpa, nas-ini);

% point near LPA that is at 50% of left lower contour
[cnt1, cnt2] = contour_follow(pnt, dhk, nas, lpa, ini);
mle = contour_fraction(cnt1, cnt2, 0.5);

% point near RPA that is at 50% of right lower contour
[cnt1, cnt2] = contour_follow(pnt, dhk, nas, rpa, ini);
mre = contour_fraction(cnt1, cnt2, 0.5);

% determine two points that approximate the vertex 
[cnt1, cnt2] = contour_follow(pnt, dhk, nas, ver, ini);
ver1 = contour_fraction(cnt1, cnt2, 0.5);
[cnt1, cnt2] = contour_follow(pnt, dhk, mle, ver, mre);
ver2 = contour_fraction(cnt1, cnt2, 0.5);

% refined estimate is the average of these two
ver = (ver1+ver2)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start contouring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ant-post contour through vertex
fprintf('vertical ant-post contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, nas, ver, ini);
Nz   = contour_fraction(cnt1, cnt2,  0/20);
NFpz = contour_fraction(cnt1, cnt2,  1/20);
Fpz  = contour_fraction(cnt1, cnt2,  2/20);
AFpz = contour_fraction(cnt1, cnt2,  3/20);
AFz  = contour_fraction(cnt1, cnt2,  4/20);
AFFz = contour_fraction(cnt1, cnt2,  5/20);
Fz   = contour_fraction(cnt1, cnt2,  6/20);
FFCz = contour_fraction(cnt1, cnt2,  7/20);
FCz  = contour_fraction(cnt1, cnt2,  8/20);
FCCz = contour_fraction(cnt1, cnt2,  9/20);
Cz   = contour_fraction(cnt1, cnt2, 10/20);
CCPz = contour_fraction(cnt1, cnt2, 11/20);
CPz  = contour_fraction(cnt1, cnt2, 12/20);
CPPz = contour_fraction(cnt1, cnt2, 13/20);
Pz   = contour_fraction(cnt1, cnt2, 14/20);
PPOz = contour_fraction(cnt1, cnt2, 15/20);
POz  = contour_fraction(cnt1, cnt2, 16/20);
POOz = contour_fraction(cnt1, cnt2, 17/20);
Oz   = contour_fraction(cnt1, cnt2, 18/20);
OIz  = contour_fraction(cnt1, cnt2, 19/20);
Iz   = contour_fraction(cnt1, cnt2, 20/20);

% left-right through vertex
fprintf('C contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, mle, ver, mre);
T9   = contour_fraction(cnt1, cnt2,  0/20);
T9h  = contour_fraction(cnt1, cnt2,  1/20);
T7   = contour_fraction(cnt1, cnt2,  2/20);
T7h  = contour_fraction(cnt1, cnt2,  3/20);
C5   = contour_fraction(cnt1, cnt2,  4/20);
C5h  = contour_fraction(cnt1, cnt2,  5/20);
C3   = contour_fraction(cnt1, cnt2,  6/20);
C3h  = contour_fraction(cnt1, cnt2,  7/20);
C1   = contour_fraction(cnt1, cnt2,  8/20);
C1h  = contour_fraction(cnt1, cnt2,  9/20);
Cz   = contour_fraction(cnt1, cnt2, 10/20);
C2h  = contour_fraction(cnt1, cnt2, 11/20);
C2   = contour_fraction(cnt1, cnt2, 12/20);
C4h  = contour_fraction(cnt1, cnt2, 13/20);
C4   = contour_fraction(cnt1, cnt2, 14/20);
C6h  = contour_fraction(cnt1, cnt2, 15/20);
C6   = contour_fraction(cnt1, cnt2, 16/20);
T8h  = contour_fraction(cnt1, cnt2, 17/20);
T8   = contour_fraction(cnt1, cnt2, 18/20);
T10h = contour_fraction(cnt1, cnt2, 19/20);
T10  = contour_fraction(cnt1, cnt2, 20/20);

% horizontal ant-post through T7
fprintf('horizontal left contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, Fpz, T7, Oz);
Fp1h = contour_fraction(cnt1, cnt2,  1/20);
Fp1  = contour_fraction(cnt1, cnt2,  2/20);
AFp7 = contour_fraction(cnt1, cnt2,  3/20);
AF7  = contour_fraction(cnt1, cnt2,  4/20);
AFF7 = contour_fraction(cnt1, cnt2,  5/20);
F7   = contour_fraction(cnt1, cnt2,  6/20);
FFT7 = contour_fraction(cnt1, cnt2,  7/20);
FT7  = contour_fraction(cnt1, cnt2,  8/20);
FTT7 = contour_fraction(cnt1, cnt2,  9/20);
T7   = contour_fraction(cnt1, cnt2, 10/20);
TTP7 = contour_fraction(cnt1, cnt2, 11/20);
TP7  = contour_fraction(cnt1, cnt2, 12/20);
TPP7 = contour_fraction(cnt1, cnt2, 13/20);
P7   = contour_fraction(cnt1, cnt2, 14/20);
PPO7 = contour_fraction(cnt1, cnt2, 15/20);
PO7  = contour_fraction(cnt1, cnt2, 16/20);
POO7 = contour_fraction(cnt1, cnt2, 17/20);
O1   = contour_fraction(cnt1, cnt2, 18/20);
O1h  = contour_fraction(cnt1, cnt2, 19/20);

% horizontal ant-post through T8
fprintf('horizontal right contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, Fpz, T8, Oz);
Fp2h = contour_fraction(cnt1, cnt2,  1/20);
Fp2  = contour_fraction(cnt1, cnt2,  2/20);
AFp8 = contour_fraction(cnt1, cnt2,  3/20);
AF8  = contour_fraction(cnt1, cnt2,  4/20);
AFF8 = contour_fraction(cnt1, cnt2,  5/20);
F8   = contour_fraction(cnt1, cnt2,  6/20);
FFT8 = contour_fraction(cnt1, cnt2,  7/20);
FT8  = contour_fraction(cnt1, cnt2,  8/20);
FTT8 = contour_fraction(cnt1, cnt2,  9/20);
T8   = contour_fraction(cnt1, cnt2, 10/20);
TTP8 = contour_fraction(cnt1, cnt2, 11/20);
TP8  = contour_fraction(cnt1, cnt2, 12/20);
TPP8 = contour_fraction(cnt1, cnt2, 13/20);
P8   = contour_fraction(cnt1, cnt2, 14/20);
PPO8 = contour_fraction(cnt1, cnt2, 15/20);
PO8  = contour_fraction(cnt1, cnt2, 16/20);
POO8 = contour_fraction(cnt1, cnt2, 17/20);
O2   = contour_fraction(cnt1, cnt2, 18/20);
O2h  = contour_fraction(cnt1, cnt2, 19/20);

fprintf('AFp contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, AFp7, AFpz, AFp8);
AFp7h  = contour_fraction(cnt1, cnt2,  1/16);
AFp5   = contour_fraction(cnt1, cnt2,  2/16);
AFp5h  = contour_fraction(cnt1, cnt2,  3/16);
AFp3   = contour_fraction(cnt1, cnt2,  4/16);
AFp3h  = contour_fraction(cnt1, cnt2,  5/16);
AFp1   = contour_fraction(cnt1, cnt2,  6/16);
AFp1h  = contour_fraction(cnt1, cnt2,  7/16);
AFp2h  = contour_fraction(cnt1, cnt2,  9/16);
AFp2   = contour_fraction(cnt1, cnt2, 10/16);
AFp4h  = contour_fraction(cnt1, cnt2, 11/16);
AFp4   = contour_fraction(cnt1, cnt2, 12/16);
AFp6h  = contour_fraction(cnt1, cnt2, 13/16);
AFp6   = contour_fraction(cnt1, cnt2, 14/16);
AFp8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('AF contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, AF7, AFz, AF8);
AF7h  = contour_fraction(cnt1, cnt2,  1/16);
AF5   = contour_fraction(cnt1, cnt2,  2/16);
AF5h  = contour_fraction(cnt1, cnt2,  3/16);
AF3   = contour_fraction(cnt1, cnt2,  4/16);
AF3h  = contour_fraction(cnt1, cnt2,  5/16);
AF1   = contour_fraction(cnt1, cnt2,  6/16);
AF1h  = contour_fraction(cnt1, cnt2,  7/16);
AF2h  = contour_fraction(cnt1, cnt2,  9/16);
AF2   = contour_fraction(cnt1, cnt2, 10/16);
AF4h  = contour_fraction(cnt1, cnt2, 11/16);
AF4   = contour_fraction(cnt1, cnt2, 12/16);
AF6h  = contour_fraction(cnt1, cnt2, 13/16);
AF6   = contour_fraction(cnt1, cnt2, 14/16);
AF8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('AFF contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, AFF7, AFFz, AFF8);
AFF7h  = contour_fraction(cnt1, cnt2,  1/16);
AFF5   = contour_fraction(cnt1, cnt2,  2/16);
AFF5h  = contour_fraction(cnt1, cnt2,  3/16);
AFF3   = contour_fraction(cnt1, cnt2,  4/16);
AFF3h  = contour_fraction(cnt1, cnt2,  5/16);
AFF1   = contour_fraction(cnt1, cnt2,  6/16);
AFF1h  = contour_fraction(cnt1, cnt2,  7/16);
AFF2h  = contour_fraction(cnt1, cnt2,  9/16);
AFF2   = contour_fraction(cnt1, cnt2, 10/16);
AFF4h  = contour_fraction(cnt1, cnt2, 11/16);
AFF4   = contour_fraction(cnt1, cnt2, 12/16);
AFF6h  = contour_fraction(cnt1, cnt2, 13/16);
AFF6   = contour_fraction(cnt1, cnt2, 14/16);
AFF8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('F contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, F7, Fz, F8);
F7h  = contour_fraction(cnt1, cnt2,  1/16);
F5   = contour_fraction(cnt1, cnt2,  2/16);
F5h  = contour_fraction(cnt1, cnt2,  3/16);
F3   = contour_fraction(cnt1, cnt2,  4/16);
F3h  = contour_fraction(cnt1, cnt2,  5/16);
F1   = contour_fraction(cnt1, cnt2,  6/16);
F1h  = contour_fraction(cnt1, cnt2,  7/16);
F2h  = contour_fraction(cnt1, cnt2,  9/16);
F2   = contour_fraction(cnt1, cnt2, 10/16);
F4h  = contour_fraction(cnt1, cnt2, 11/16);
F4   = contour_fraction(cnt1, cnt2, 12/16);
F6h  = contour_fraction(cnt1, cnt2, 13/16);
F6   = contour_fraction(cnt1, cnt2, 14/16);
F8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('FFC contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, FFT7, FFCz, FFT8);
FFT7h  = contour_fraction(cnt1, cnt2,  1/16);
FFC5   = contour_fraction(cnt1, cnt2,  2/16);
FFC5h  = contour_fraction(cnt1, cnt2,  3/16);
FFC3   = contour_fraction(cnt1, cnt2,  4/16);
FFC3h  = contour_fraction(cnt1, cnt2,  5/16);
FFC1   = contour_fraction(cnt1, cnt2,  6/16);
FFC1h  = contour_fraction(cnt1, cnt2,  7/16);
FFC2h  = contour_fraction(cnt1, cnt2,  9/16);
FFC2   = contour_fraction(cnt1, cnt2, 10/16);
FFC4h  = contour_fraction(cnt1, cnt2, 11/16);
FFC4   = contour_fraction(cnt1, cnt2, 12/16);
FFC6h  = contour_fraction(cnt1, cnt2, 13/16);
FFC6   = contour_fraction(cnt1, cnt2, 14/16);
FFT8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('FC contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, FT7, FCz, FT8);
FT7h  = contour_fraction(cnt1, cnt2,  1/16);
FC5   = contour_fraction(cnt1, cnt2,  2/16);
FC5h  = contour_fraction(cnt1, cnt2,  3/16);
FC3   = contour_fraction(cnt1, cnt2,  4/16);
FC3h  = contour_fraction(cnt1, cnt2,  5/16);
FC1   = contour_fraction(cnt1, cnt2,  6/16);
FC1h  = contour_fraction(cnt1, cnt2,  7/16);
FC2h  = contour_fraction(cnt1, cnt2,  9/16);
FC2   = contour_fraction(cnt1, cnt2, 10/16);
FC4h  = contour_fraction(cnt1, cnt2, 11/16);
FC4   = contour_fraction(cnt1, cnt2, 12/16);
FC6h  = contour_fraction(cnt1, cnt2, 13/16);
FC6   = contour_fraction(cnt1, cnt2, 14/16);
FT8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('FCC contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, FTT7, FCCz, FTT8);
FTT7h  = contour_fraction(cnt1, cnt2,  1/16);
FCC5   = contour_fraction(cnt1, cnt2,  2/16);
FCC5h  = contour_fraction(cnt1, cnt2,  3/16);
FCC3   = contour_fraction(cnt1, cnt2,  4/16);
FCC3h  = contour_fraction(cnt1, cnt2,  5/16);
FCC1   = contour_fraction(cnt1, cnt2,  6/16);
FCC1h  = contour_fraction(cnt1, cnt2,  7/16);
FCC2h  = contour_fraction(cnt1, cnt2,  9/16);
FCC2   = contour_fraction(cnt1, cnt2, 10/16);
FCC4h  = contour_fraction(cnt1, cnt2, 11/16);
FCC4   = contour_fraction(cnt1, cnt2, 12/16);
FCC6h  = contour_fraction(cnt1, cnt2, 13/16);
FCC6   = contour_fraction(cnt1, cnt2, 14/16);
FTT8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('CCP contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, TTP7, CCPz, TTP8);
TTP7h  = contour_fraction(cnt1, cnt2,  1/16);
CCP5   = contour_fraction(cnt1, cnt2,  2/16);
CCP5h  = contour_fraction(cnt1, cnt2,  3/16);
CCP3   = contour_fraction(cnt1, cnt2,  4/16);
CCP3h  = contour_fraction(cnt1, cnt2,  5/16);
CCP1   = contour_fraction(cnt1, cnt2,  6/16);
CCP1h  = contour_fraction(cnt1, cnt2,  7/16);
CCP2h  = contour_fraction(cnt1, cnt2,  9/16);
CCP2   = contour_fraction(cnt1, cnt2, 10/16);
CCP4h  = contour_fraction(cnt1, cnt2, 11/16);
CCP4   = contour_fraction(cnt1, cnt2, 12/16);
CCP6h  = contour_fraction(cnt1, cnt2, 13/16);
CCP6   = contour_fraction(cnt1, cnt2, 14/16);
TTP8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('CP contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, TP7, CPz, TP8);
TP7h  = contour_fraction(cnt1, cnt2,  1/16);
CP5   = contour_fraction(cnt1, cnt2,  2/16);
CP5h  = contour_fraction(cnt1, cnt2,  3/16);
CP3   = contour_fraction(cnt1, cnt2,  4/16);
CP3h  = contour_fraction(cnt1, cnt2,  5/16);
CP1   = contour_fraction(cnt1, cnt2,  6/16);
CP1h  = contour_fraction(cnt1, cnt2,  7/16);
CP2h  = contour_fraction(cnt1, cnt2,  9/16);
CP2   = contour_fraction(cnt1, cnt2, 10/16);
CP4h  = contour_fraction(cnt1, cnt2, 11/16);
CP4   = contour_fraction(cnt1, cnt2, 12/16);
CP6h  = contour_fraction(cnt1, cnt2, 13/16);
CP6   = contour_fraction(cnt1, cnt2, 14/16);
TP8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('CPP contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, TPP7, CPPz, TPP8);
TPP7h  = contour_fraction(cnt1, cnt2,  1/16);
CPP5   = contour_fraction(cnt1, cnt2,  2/16);
CPP5h  = contour_fraction(cnt1, cnt2,  3/16);
CPP3   = contour_fraction(cnt1, cnt2,  4/16);
CPP3h  = contour_fraction(cnt1, cnt2,  5/16);
CPP1   = contour_fraction(cnt1, cnt2,  6/16);
CPP1h  = contour_fraction(cnt1, cnt2,  7/16);
CPP2h  = contour_fraction(cnt1, cnt2,  9/16);
CPP2   = contour_fraction(cnt1, cnt2, 10/16);
CPP4h  = contour_fraction(cnt1, cnt2, 11/16);
CPP4   = contour_fraction(cnt1, cnt2, 12/16);
CPP6h  = contour_fraction(cnt1, cnt2, 13/16);
CPP6   = contour_fraction(cnt1, cnt2, 14/16);
TPP8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('P contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, P7, Pz, P8);
P7h  = contour_fraction(cnt1, cnt2,  1/16);
P5   = contour_fraction(cnt1, cnt2,  2/16);
P5h  = contour_fraction(cnt1, cnt2,  3/16);
P3   = contour_fraction(cnt1, cnt2,  4/16);
P3h  = contour_fraction(cnt1, cnt2,  5/16);
P1   = contour_fraction(cnt1, cnt2,  6/16);
P1h  = contour_fraction(cnt1, cnt2,  7/16);
P2h  = contour_fraction(cnt1, cnt2,  9/16);
P2   = contour_fraction(cnt1, cnt2, 10/16);
P4h  = contour_fraction(cnt1, cnt2, 11/16);
P4   = contour_fraction(cnt1, cnt2, 12/16);
P6h  = contour_fraction(cnt1, cnt2, 13/16);
P6   = contour_fraction(cnt1, cnt2, 14/16);
P8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('PPO contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, PPO7, PPOz, PPO8);
PPO7h  = contour_fraction(cnt1, cnt2,  1/16);
PPO5   = contour_fraction(cnt1, cnt2,  2/16);
PPO5h  = contour_fraction(cnt1, cnt2,  3/16);
PPO3   = contour_fraction(cnt1, cnt2,  4/16);
PPO3h  = contour_fraction(cnt1, cnt2,  5/16);
PPO1   = contour_fraction(cnt1, cnt2,  6/16);
PPO1h  = contour_fraction(cnt1, cnt2,  7/16);
PPO2h  = contour_fraction(cnt1, cnt2,  9/16);
PPO2   = contour_fraction(cnt1, cnt2, 10/16);
PPO4h  = contour_fraction(cnt1, cnt2, 11/16);
PPO4   = contour_fraction(cnt1, cnt2, 12/16);
PPO6h  = contour_fraction(cnt1, cnt2, 13/16);
PPO6   = contour_fraction(cnt1, cnt2, 14/16);
PPO8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('PO contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, PO7, POz, PO8);
PO7h  = contour_fraction(cnt1, cnt2,  1/16);
PO5   = contour_fraction(cnt1, cnt2,  2/16);
PO5h  = contour_fraction(cnt1, cnt2,  3/16);
PO3   = contour_fraction(cnt1, cnt2,  4/16);
PO3h  = contour_fraction(cnt1, cnt2,  5/16);
PO1   = contour_fraction(cnt1, cnt2,  6/16);
PO1h  = contour_fraction(cnt1, cnt2,  7/16);
PO2h  = contour_fraction(cnt1, cnt2,  9/16);
PO2   = contour_fraction(cnt1, cnt2, 10/16);
PO4h  = contour_fraction(cnt1, cnt2, 11/16);
PO4   = contour_fraction(cnt1, cnt2, 12/16);
PO6h  = contour_fraction(cnt1, cnt2, 13/16);
PO6   = contour_fraction(cnt1, cnt2, 14/16);
PO8h  = contour_fraction(cnt1, cnt2, 15/16);

fprintf('POO contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, POO7, POOz, POO8);
POO7h  = contour_fraction(cnt1, cnt2,  1/16);
POO5   = contour_fraction(cnt1, cnt2,  2/16);
POO5h  = contour_fraction(cnt1, cnt2,  3/16);
POO3   = contour_fraction(cnt1, cnt2,  4/16);
POO3h  = contour_fraction(cnt1, cnt2,  5/16);
POO1   = contour_fraction(cnt1, cnt2,  6/16);
POO1h  = contour_fraction(cnt1, cnt2,  7/16);
POO2h  = contour_fraction(cnt1, cnt2,  9/16);
POO2   = contour_fraction(cnt1, cnt2, 10/16);
POO4h  = contour_fraction(cnt1, cnt2, 11/16);
POO4   = contour_fraction(cnt1, cnt2, 12/16);
POO6h  = contour_fraction(cnt1, cnt2, 13/16);
POO6   = contour_fraction(cnt1, cnt2, 14/16);
POO8h  = contour_fraction(cnt1, cnt2, 15/16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start contouring the low electrode locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low horizontal ant-post through T9
fprintf('low horizontal left contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, Nz, T9, Iz);
AFp9 = contour_fraction(cnt1, cnt2,  3/20);
AF9  = contour_fraction(cnt1, cnt2,  4/20);
AFF9 = contour_fraction(cnt1, cnt2,  5/20);
F9   = contour_fraction(cnt1, cnt2,  6/20);
FFT9 = contour_fraction(cnt1, cnt2,  7/20);
FT9  = contour_fraction(cnt1, cnt2,  8/20);
FTT9 = contour_fraction(cnt1, cnt2,  9/20);
T9   = contour_fraction(cnt1, cnt2, 10/20);
TTP9 = contour_fraction(cnt1, cnt2, 11/20);
TP9  = contour_fraction(cnt1, cnt2, 12/20);
TPP9 = contour_fraction(cnt1, cnt2, 13/20);
P9   = contour_fraction(cnt1, cnt2, 14/20);
PPO9 = contour_fraction(cnt1, cnt2, 15/20);
PO9  = contour_fraction(cnt1, cnt2, 16/20);
POO9 = contour_fraction(cnt1, cnt2, 17/20);
I1   = contour_fraction(cnt1, cnt2, 18/20);
I1h  = contour_fraction(cnt1, cnt2, 19/20);

[cnt1, cnt2] = contour_follow(pnt, dhk, NFpz, T9h, OIz);
AFp9h = contour_fraction(cnt1, cnt2,  3/20);
AF9h  = contour_fraction(cnt1, cnt2,  4/20);
AFF9h = contour_fraction(cnt1, cnt2,  5/20);
F9h   = contour_fraction(cnt1, cnt2,  6/20);
FFT9h = contour_fraction(cnt1, cnt2,  7/20);
FT9h  = contour_fraction(cnt1, cnt2,  8/20);
FTT9h = contour_fraction(cnt1, cnt2,  9/20);
T9h   = contour_fraction(cnt1, cnt2, 10/20);
TTP9h = contour_fraction(cnt1, cnt2, 11/20);
TP9h  = contour_fraction(cnt1, cnt2, 12/20);
TPP9h = contour_fraction(cnt1, cnt2, 13/20);
P9h   = contour_fraction(cnt1, cnt2, 14/20);
PPO9h = contour_fraction(cnt1, cnt2, 15/20);
PO9h  = contour_fraction(cnt1, cnt2, 16/20);
POO9h = contour_fraction(cnt1, cnt2, 17/20);
OI1   = contour_fraction(cnt1, cnt2, 18/20);
OI1h  = contour_fraction(cnt1, cnt2, 19/20);

% low horizontal ant-post through T10
fprintf('low horizontal right contour\n');
[cnt1, cnt2] = contour_follow(pnt, dhk, Nz, T10, Iz);
AFp10 = contour_fraction(cnt1, cnt2,  3/20);
AF10  = contour_fraction(cnt1, cnt2,  4/20);
AFF10 = contour_fraction(cnt1, cnt2,  5/20);
F10   = contour_fraction(cnt1, cnt2,  6/20);
FFT10 = contour_fraction(cnt1, cnt2,  7/20);
FT10  = contour_fraction(cnt1, cnt2,  8/20);
FTT10 = contour_fraction(cnt1, cnt2,  9/20);
T10   = contour_fraction(cnt1, cnt2, 10/20);
TTP10 = contour_fraction(cnt1, cnt2, 11/20);
TP10  = contour_fraction(cnt1, cnt2, 12/20);
TPP10 = contour_fraction(cnt1, cnt2, 13/20);
P10   = contour_fraction(cnt1, cnt2, 14/20);
PPO10 = contour_fraction(cnt1, cnt2, 15/20);
PO10  = contour_fraction(cnt1, cnt2, 16/20);
POO10 = contour_fraction(cnt1, cnt2, 17/20);
I2    = contour_fraction(cnt1, cnt2, 18/20);
I2h   = contour_fraction(cnt1, cnt2, 19/20);

[cnt1, cnt2] = contour_follow(pnt, dhk, NFpz, T10h, OIz);
AFp10h = contour_fraction(cnt1, cnt2,  3/20);
AF10h  = contour_fraction(cnt1, cnt2,  4/20);
AFF10h = contour_fraction(cnt1, cnt2,  5/20);
F10h   = contour_fraction(cnt1, cnt2,  6/20);
FFT10h = contour_fraction(cnt1, cnt2,  7/20);
FT10h  = contour_fraction(cnt1, cnt2,  8/20);
FTT10h = contour_fraction(cnt1, cnt2,  9/20);
T10h   = contour_fraction(cnt1, cnt2, 10/20);
TTP10h = contour_fraction(cnt1, cnt2, 11/20);
TP10h  = contour_fraction(cnt1, cnt2, 12/20);
TPP10h = contour_fraction(cnt1, cnt2, 13/20);
P10h   = contour_fraction(cnt1, cnt2, 14/20);
PPO10h = contour_fraction(cnt1, cnt2, 15/20);
PO10h  = contour_fraction(cnt1, cnt2, 16/20);
POO10h = contour_fraction(cnt1, cnt2, 17/20);
OI2    = contour_fraction(cnt1, cnt2, 18/20);
OI2h   = contour_fraction(cnt1, cnt2, 19/20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect all the computed electrode locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lab = {
'LPA',
'RPA',
'NAS',
'INI',
'Nz',
'Fp1',
'Fpz',
'Fp2',
'AF9',
'AF7',
'AF5',
'AF3',
'AF1',
'AFz',
'AF2',
'AF4',
'AF6',
'AF8',
'AF10',
'F9',
'F7',
'F5',
'F3',
'F1',
'Fz',
'F2',
'F4',
'F6',
'F8',
'F10',
'FT9',
'FT7',
'FC5',
'FC3',
'FC1',
'FCz',
'FC2',
'FC4',
'FC6',
'FT8',
'FT10',
'T9',
'T7',
'C5',
'C3',
'C1',
'Cz',
'C2',
'C4',
'C6',
'T8',
'T10',
'TP9',
'TP7',
'CP5',
'CP3',
'CP1',
'CPz',
'CP2',
'CP4',
'CP6',
'TP8',
'TP10',
'P9',
'P7',
'P5',
'P3',
'P1',
'Pz',
'P2',
'P4',
'P6',
'P8',
'P10',
'PO9',
'PO7',
'PO5',
'PO3',
'PO1',
'POz',
'PO2',
'PO4',
'PO6',
'PO8',
'PO10',
'O1',
'Oz',
'O2',
'I1',
'Iz',
'I2',
'AFp9h',
'AFp7h',
'AFp5h',
'AFp3h',
'AFp1h',
'AFp2h',
'AFp4h',
'AFp6h',
'AFp8h',
'AFp10h',
'AFF9h',
'AFF7h',
'AFF5h',
'AFF3h',
'AFF1h',
'AFF2h',
'AFF4h',
'AFF6h',
'AFF8h',
'AFF10h',
'FFT9h',
'FFT7h',
'FFC5h',
'FFC3h',
'FFC1h',
'FFC2h',
'FFC4h',
'FFC6h',
'FFT8h',
'FFT10h',
'FTT9h',
'FTT7h',
'FCC5h',
'FCC3h',
'FCC1h',
'FCC2h',
'FCC4h',
'FCC6h',
'FTT8h',
'FTT10h',
'TTP9h',
'TTP7h',
'CCP5h',
'CCP3h',
'CCP1h',
'CCP2h',
'CCP4h',
'CCP6h',
'TTP8h',
'TTP10h',
'TPP9h',
'TPP7h',
'CPP5h',
'CPP3h',
'CPP1h',
'CPP2h',
'CPP4h',
'CPP6h',
'TPP8h',
'TPP10h',
'PPO9h',
'PPO7h',
'PPO5h',
'PPO3h',
'PPO1h',
'PPO2h',
'PPO4h',
'PPO6h',
'PPO8h',
'PPO10h',
'POO9h',
'POO7h',
'POO5h',
'POO3h',
'POO1h',
'POO2h',
'POO4h',
'POO6h',
'POO8h',
'POO10h',
'OI1h',
'OI2h',
'Fp1h',
'Fp2h',
'AF9h',
'AF7h',
'AF5h',
'AF3h',
'AF1h',
'AF2h',
'AF4h',
'AF6h',
'AF8h',
'AF10h',
'F9h',
'F7h',
'F5h',
'F3h',
'F1h',
'F2h',
'F4h',
'F6h',
'F8h',
'F10h',
'FT9h',
'FT7h',
'FC5h',
'FC3h',
'FC1h',
'FC2h',
'FC4h',
'FC6h',
'FT8h',
'FT10h',
'T9h',
'T7h',
'C5h',
'C3h',
'C1h',
'C2h',
'C4h',
'C6h',
'T8h',
'T10h',
'TP9h',
'TP7h',
'CP5h',
'CP3h',
'CP1h',
'CP2h',
'CP4h',
'CP6h',
'TP8h',
'TP10h',
'P9h',
'P7h',
'P5h',
'P3h',
'P1h',
'P2h',
'P4h',
'P6h',
'P8h',
'P10h',
'PO9h',
'PO7h',
'PO5h',
'PO3h',
'PO1h',
'PO2h',
'PO4h',
'PO6h',
'PO8h',
'PO10h',
'O1h',
'O2h',
'I1h',
'I2h',
'AFp9',
'AFp7',
'AFp5',
'AFp3',
'AFp1',
'AFpz',
'AFp2',
'AFp4',
'AFp6',
'AFp8',
'AFp10',
'AFF9',
'AFF7',
'AFF5',
'AFF3',
'AFF1',
'AFFz',
'AFF2',
'AFF4',
'AFF6',
'AFF8',
'AFF10',
'FFT9',
'FFT7',
'FFC5',
'FFC3',
'FFC1',
'FFCz',
'FFC2',
'FFC4',
'FFC6',
'FFT8',
'FFT10',
'FTT9',
'FTT7',
'FCC5',
'FCC3',
'FCC1',
'FCCz',
'FCC2',
'FCC4',
'FCC6',
'FTT8',
'FTT10',
'TTP9',
'TTP7',
'CCP5',
'CCP3',
'CCP1',
'CCPz',
'CCP2',
'CCP4',
'CCP6',
'TTP8',
'TTP10',
'TPP9',
'TPP7',
'CPP5',
'CPP3',
'CPP1',
'CPPz',
'CPP2',
'CPP4',
'CPP6',
'TPP8',
'TPP10',
'PPO9',
'PPO7',
'PPO5',
'PPO3',
'PPO1',
'PPOz',
'PPO2',
'PPO4',
'PPO6',
'PPO8',
'PPO10',
'POO9',
'POO7',
'POO5',
'POO3',
'POO1',
'POOz',
'POO2',
'POO4',
'POO6',
'POO8',
'POO10',
'OI1',
'OIz',
'OI2',
'T3',
'T5',
'T4',
'T6',
'M1',
'M2',
'A1',
'A2'};

% remove trailing blanks from the electrode labels
nlab = length(lab);
for i=1:nlab
  lab(i) = deblank2(lab(i));
end

% do not forget these
NAS = nas;
INI = ini;
LPA = lpa;
RPA = rpa;

% assign the known local electrode positions to the output matrix
elc = ones(nlab,3) * nan;
for i=1:nlab
  if exist(char(lab(i)))
    eval(sprintf('elc(%d,:) = %s;', i, char(lab(i))));
  end
end

% remove unknown electrode positions
sel = ~isnan(elc(:,1));
elc = elc(find(sel), :);
lab = lab(find(sel));

