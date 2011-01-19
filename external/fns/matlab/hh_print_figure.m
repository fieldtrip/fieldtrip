function hh_print_figure(fname,PaperSize,options)
% This routine print the picture to the image file
% Syntax -> hh_print_figure(outfile,PaperSize,option), where
%   1. fname: Output file name. Do not need to specifu the file extension!
%   2. PaperSize: Output image sizes
%   3. options: Output image format
%	1). option = 1 -> Output to a PDF file
%	2). option = 2 -> Output to a PS file
%	3). option = 3 -> Output to a EPS file
%	4). option = 4 -> Output to a EMF file
%	5). option = 5 -> Output to a PNG file
%
% $Date: Updated on July 23, 2008$
% $author: Copyright (C) by Hung Dang$
% $email: hungptit@gmail.com$

%% Main code
mfigset = struct('PaperSize',[6,3],'margin',[0,0]);
mfigset.PaperSize = PaperSize;

% Setup figure orientation
set(gcf,'PaperOrientation','portrait')

% Paper units
set(gcf,'PaperUnits','inches');

% Paper size
set(gcf,'PaperSize',mfigset.PaperSize);

% Paper position
set(gcf,'PaperPosition',[mfigset.margin(1),mfigset.margin(2),...
    mfigset.PaperSize(1) - 2 * mfigset.margin(1),mfigset.PaperSize(2) - 2 * mfigset.margin(2)]);

% Change the figure corlor
% set(gcf,'Color',[0,0,0]);

% Print to file based on options
switch options
    case 1
        str = strcat(fname,'.pdf');
        % Print to a PDF file
        print('-dpdf',str);
    case 2
        str = strcat(fname,'.ps');
        % Print to a black and white PS file
        print('-dpsc',str);
    case 3
        str = strcat(fname,'.eps');
        % Print to a color EPS file
        print('-depsc2',str);
    case 4
        str = strcat(fname,'.emf');
        % print to a metal file
        print('-dmeta',str);
    case 5
        str = strcat(fname,'.png');
        % print to a metal file
        print('-dpng',str);
    otherwise
        error('Invalid options!');
end

% Display the debug message
fprintf('Output image: %s\n',str);

%% EOF
