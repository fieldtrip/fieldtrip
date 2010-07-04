function [h,m,s] = sectotime (ts)
%[h,m,s] = sectotime (ts)
%
% ts : time in seconds
% [h,m,s] : hour, minute, second representation of ts
%
%				Ziad Saad Oct 22 97
%
h = fix(ts./3600);
m = fix(rem(ts,3600)./60);
s = rem(rem(ts,3600),60);
