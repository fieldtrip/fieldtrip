function [sig,mixedsig]=demosig();
%
% function [sig,mixedsig]=demosig();
% 
% Returns artificially generated test signals, sig, and mixed
% signals, mixedsig. Signals are row vectors of
% matrices. Input mixedsig to FastICA to see how it works.

% @(#)$Id: demosig.m,v 1.2 2003/04/05 14:23:57 jarmo Exp $

%create source signals (independent components)
N=500; %data size

v=[0:N-1];
sig=[];
sig(1,:)=sin(v/2); %sinusoid
sig(2,:)=((rem(v,23)-11)/9).^5; %funny curve
sig(3,:)=((rem(v,27)-13)/9); %saw-tooth
sig(4,:)=((rand(1,N)<.5)*2-1).*log(rand(1,N)); %impulsive noise

for t=1:4
sig(t,:)=sig(t,:)/std(sig(t,:));
end

%remove mean (not really necessary)

[sig mean]=remmean(sig);

%create mixtures

Aorig=rand(size(sig,1));
mixedsig=(Aorig*sig);
