% topo2sph() - convert a topoplot() style 2-D polar-coordinate
%              channel locations file to a 3-D spherical-angle
%              file for use with headplot()
% Usage: 
%   >> [c h] = topo2sph('eloc_file','eloc_outfile', method, unshrink);
%   >> [c h] = topo2sph( topoarray, method, unshrink );
%
% Inputs:
%   'eloc_file'    = filename of polar 2-D electrode locations file used by 
%                    topoplot(). See >> topoplot example or cart2topo()
%   'eloc_outfile' = output file of 3-D electrode locations in spherical angle 
%                    coords. for use in headplot().
%   topoarray      = polar array of 2-D electrode locations, with polar angle
%                    in the first column and radius in the second one.
%   method         = [1|2] 1 is for Besa compatibility, 2 is for
%                    compatibility with Matlab function cart2sph(). {default: 2}
%   unshrink       = [0<real<1] unshrink factor. Enter a shrink factor used
%                    to convert spherical to topo (see sph2topo()). Only 
%                    implemented for 'method' 1 (above). Electrode 'shrink' 
%                    is now deprecated. See >> help topoplot
% Outputs:
%   c = coronal rotation
%   h = horizontal rotation
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: sph2topo(), cart2topo()

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 3-16-00 changed name to topo2sph() for compatibility with cart2topo() -sm
% 01-25-02 reformated help & license -ad 
% 03-22-02 complete remodeling for returning arguments and taking arrays -ad 

function [c, h] = topo2sph(eloc_locs,eloc_angles, method, unshrink)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end
if nargin > 1 && ~ischar(eloc_angles)
	if nargin > 2
		unshrink = method;
    end
	method = eloc_angles;
else
	method = 2;
end

if isstr(eloc_locs)
	fid = fopen(eloc_locs);
	if fid<1
	    fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
	    return
	end
	E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
	E = E';
	fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end
    
if nargin > 1 && ischar(eloc_angles)
	if exist(eloc_angles)==2
	   fprintf('topo2sph: eloc_angles file (%s) already exists and will be erased.\n',eloc_angles);
	end

	fid = fopen(eloc_angles,'a');
	if fid<1
	    fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
	    return
	end
end

if method == 2
	t = E(:,2); % theta
	r = E(:,3); % radius
	h = -t;  % horizontal rotation
	c = (0.5-r)*180;
else
	for e=1:size(E,1)
		% (t,r) -> (c,h)
		
		t = E(e,2); % theta
		r = E(e,3); % radius
		r = r*unshrink;
		if t>=0
			h(e) = 90-t; % horizontal rotation
		else
			h(e) = -(90+t);
		end
		if t~=0
			c(e) = sign(t)*180*r; % coronal rotation
		else
			c(e) = 180*r;
		end
    end
	t = t';
	r = r';
end

for e=1:size(E,1)
   if nargin > 1 && ischar(eloc_angles)
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
   end 
end

