function [y1]=rs(y1,T,f2)
% [y2] = rs(y1, T) resamples y1 to the target sampling rate y2 using T
% [y2] = rs(y1, f1, f2) resamples y1 with f1 to the target sampling rate f2 
%
% RS does not require overlap data. 
%
% see also: SOPEN, SREAD, SCLOSE, MAT2SEL, SAVE2TXT, SAVE2BKR
%
% Reference(s):

%	$Revision$
%	$Id$
%	Copyright (C) 1997-2004 by Alois Schloegl 
%	a.schloegl@ieee.org	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

if nargin==3,
        f1=T;
        if f1==f2
                return;
        elseif f1>f2
                D=f1/f2;
                [yr,yc]=size(y1);
                LEN=yr/D;
		y2=zeros(yr*f2/f1,yc);
                for k=0:LEN-1
                        y2(k+1,:)=sum(y1(k*D+(1:D),:),1)/D;
                end;
		y1=y2;
        else %f1<f2
		y1=y1(ceil((1:size(y1,1)*f2/f1)/f2*f1),:);                
        end;
        
elseif nargin==2,
        [f1,f2]=size(T);
        if f1==f2,
                return;
        end;
        [yr,yc]=size(y1);
        LEN=yr/f1;
	y2=zeros(yr*f2/f1,yc);
        for k=0:LEN-1
                y2(k*f2+(1:f2),:)=T'*y1(k*f1+(1:f1),:);
        end;
	y1=y2;
end;
