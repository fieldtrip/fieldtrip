function GFPRV = readBESAgfprv(filename)

% readBESAgfprv reads all information from a .dat file that contains the
% global field power (GFP) and the residual variance (RV) as exported from 
% BESA
%
% Use as
%   GFPRV = readBESAgfprv(filename)
%
% The output is a structure containing the following fields:
%   Npts: number of sample points
%   TSB: latency of the first sample
%   DI: time interval between two consecutive sample points
%   Total RV: the total residual variance in percent
%   MinRV: the minimum residual variance
%   MaxGFP: The absolut value of the maximum of the global field power
%   RV: the residual variance in percent [Npts x 1]
%   GFP: the global field power in percent [Npts x 1
%
% Last modified November 10, 2006 Karsten Hoechstetter


if isempty(findstr(filename,'.'))
  filename = [filename,'.dat'];
end
fp = fopen(filename);

header = fgetl(fp);
headerinfo = sscanf(header,'Npts= %f TSB= %f DI= %f TotalRV= %f%% MinRV= %f%% MaxGFP= %f');

GFPRV.Npts = headerinfo(1);
GFPRV.TSB = headerinfo(2);
GFPRV.DI = headerinfo(3);
GFPRV.TotalRV = headerinfo(4);
GFPRV.MinRV = headerinfo(5);
GFPRV.MaxGFP = headerinfo(6);

RVline = fgetl(fp);
GFPRV.RV = sscanf(RVline(6:end),'%f');

GFPline = fgetl(fp);
GFPRV.GFP = sscanf(GFPline(6:end),'%f');

