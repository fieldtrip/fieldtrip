function SQIscore = SQI(OD1,OD2,oxy,dxy,Fs)
%%% This function reads in fNIRS data and returns a scalar value representing
%%% the signal quality score.
%%%
%%% Use as SQIscore = SQI(OD1,OD2,oxy,dxy,Fs)
%%% where
%%%     OD1, OD2    -   1xN vectors containing optical density signals from 
%%%                     each wavelength. N is the number of samples. 
%%%     oxy         -   1xN vector containing O2Hb signal (concentration
%%%                     changes in oxygenated hemoglobin)
%%%     dxy         -   1xN vector containing HHb signal (concentration
%%%                     changes in deoxygenated hemoglobin)
%%%     Fs          -   scalar, sampling rate in Hz
%%%
%%%     SQIscore    -   scalar, signal quality score ranging from 1 (very 
%%%                     low quality) to 5 (very high quality)
%%%
%%% All input arguments are necessary.
%%%
%%% The performance of this function has been tested on 10-second interval 
%%% signal segments recorded with Artinis devices (OxyMon, OctaMon, Brite23, Brite24)
%%%
%%% This script makes use of the following functions:
%%% detrend                         - from Matlab Data Import and Analysis toolbox
%%% ft_preprproc_bandpassfilter     - from Fieldtrip toolbox (https://github.com/fieldtrip/fieldtrip)
%%% xcorr                           - from Matlab Signal Processing toolbox
%%%
%%% Version 0.1, copyright (c) by Artinis Medical Systems http://www.artinis.com. 
%%% Last modified on 08-24-2020
%%% Authors: Sofia Sappia (sofia@artinis.com) and Naser Hakimi (naser@artinis.com) 
%%% Cite as: Sappia, Hakimi, Colier, Horschig (Submitted for peer review): 
%%% "Signal Quality Index: an algorithm for quantitative assessment of
%%% functional near infrared spectroscopy signal quality"
%%%
%%% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
%%% Based on a work at https://github.com/Artinis-Medical-Systems-B-V/SignalQualityIndex.
%%% Permissions beyond the scope of this license may be available upon request at science@artinis.com.
%%%-----------------------------------------------------------------------------------------------

%%%% Setting algorithm parameters
% Thresholds for features in rating stages one and two
thrUp_intensity     = 2.5;
thrLow_intensity    = 0.04;
thr_sumHbratio      = 1.95;
thr_acorrDiffODs    = 0.025;
% Slope and intercept for score conversion in rating stage three
slope       = 1.796;
intercept   = 0.846;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RATING STAGE ONE: Identifying very low quality signals %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First feature in rating stage one: counts are outside of linear range
if (any(OD1<thrLow_intensity) || any(OD1>thrUp_intensity) || any(OD2<thrLow_intensity) || any(OD2>thrUp_intensity))
    SQIscore=1;
    return;
end

% Second feature in rating stage one: at least one of the optical density signals is a flat line
if (std(OD1)==0 || std(OD2)==0)
    SQIscore=1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Filtering the signals %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero padding
OD1 = [zeros(1,2*Fs) detrend(OD1) zeros(1,2*Fs)];
OD2 = [zeros(1,2*Fs) detrend(OD2) zeros(1,2*Fs)];
oxy = [zeros(1,2*Fs) detrend(oxy) zeros(1,2*Fs)];
dxy = [zeros(1,2*Fs) detrend(dxy) zeros(1,2*Fs)];
[OD1_filt, ~, ~] = ft_preproc_bandpassfilter((OD1),Fs,[0.4,3],[],'firws','onepass-zerophase',[],[],[],[],[],[]);
[OD2_filt, ~, ~] = ft_preproc_bandpassfilter((OD2),Fs,[0.4,3],[],'firws','onepass-zerophase',[],[],[],[],[],[]);
[oxy_filt, ~, ~] = ft_preproc_bandpassfilter((oxy),Fs,[0.4,3],[],'firws','onepass-zerophase',[],[],[],[],[],[]);
[dxy_filt, ~, ~] = ft_preproc_bandpassfilter((dxy),Fs,[0.4,3],[],'firws','onepass-zerophase',[],[],[],[],[],[]);
% Deleting the padded data
OD1_filt = OD1_filt(2*Fs+1:end-2*Fs);
OD2_filt = OD2_filt(2*Fs+1:end-2*Fs);
oxy_filt = oxy_filt(2*Fs+1:end-2*Fs);
dxy_filt = dxy_filt(2*Fs+1:end-2*Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third feature in rating stage one: Different Scale of Oxy and Deoxy
if (sum(abs(oxy_filt))/sum(abs(dxy_filt)) < thr_sumHbratio)
    SQIscore=1;
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RATING STAGE TWO: Identifying very high quality signals %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature in rating stage two: The difference between the auto-correlation signals of filtered OD signals should be low
[autocorr_od1_filt,~] = xcorr(OD1_filt, OD1_filt, 'coeff');
[autocorr_od2_filt,~] = xcorr(OD2_filt, OD2_filt, 'coeff');
if ((std(autocorr_od1_filt - autocorr_od2_filt))<thr_acorrDiffODs)
    SQIscore=5;
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  RATING STAGE THREE: Signal quality rating  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature in rating stage three: the standard deviation of O2Hb is higher than 
% the standard deviation of HHb for high quality NIRS signals
logStdHb = log( std(oxy_filt)/std(dxy_filt) );
SQIscore = logStdHb*slope + intercept;

% forcing the score to be within 1 (very low quality) and 5 (very high
% quality)
if SQIscore < 1
    SQIscore = 1;
elseif SQIscore > 5
    SQIscore = 5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
