function err = figarng (hvect,xll,yll,xur,yur,nfx)
%	err = figarng (hvect,xll,yll,xur,yur,nfx)
% arranges the figures with handles present in hvect.
% The figures are placed in the region bound at the lower left corner
% by (xll,yll) and the upper right corner by (xur,yur)
% nfx is the number of figures placed in the x direction
% to cover all the screen, use :
% xll,yll = 0 , 0
% xur,yur = 1280, 1024  (approx values).
%           3100, 1020  (dual disp)
% err = 0 if all is fine, 1 otherwise
%
% Ziad Saad			Sept 27 97
%

err = 1;
nf = length(hvect);

monpos = get(0, 'monitorposition');

mx_x = monpos(1,3)+monpos(2,3);
mx_y = monpos(1,4)+monpos(2,4);

if (nf == 0),
	fprintf (2,'figarng : No figure handles specified\n');
	return;
end

if (xll < 0 | yll < 0 | xur > mx_x | yur > mx_y),
	fprintf (2,'figarng : display boundaries out of legal range\n');
	return;
end

frmx = 13;	% frame size (approximate extra size of window decorations
frmy = 100;

nfy = (nf./nfx);

xwdth = xur - xll ;
ywdth = yur - yll ;

fwdx = (xwdth)./nfx - frmx;
fwdy = (ywdth)./nfy - frmy;

xlli = xll;
ylli = yll;

for i=1:nf,
	figure (i);
	h = setfpos ( xlli , ylli );
	setfsize (fwdx,fwdy);
	xlli = xlli + fwdx + frmx;
	if (xlli + fwdx > xur),
		xlli = xll;
		ylli = ylli + fwdy + frmy;
		if (ylli + fwdy > mx_y),
      fprintf(1,'Warning figarng: Overlapping figures!\n');
		ylli = yll;
		end
	end
end
