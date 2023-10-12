Scalable Vector Graphics (SVG) Export of Figures

Converts 2D & 3D Matlab plots to the scalable vector format (SVG). This format is specified by W3C (http://www.w3.org) and can be viewed and printed with internet browsers.

Added preliminary support of filter, clipping, and tickmark extensions that go beyond the Matlab functionality. SVG filters are a great tool to create stylish plots! Try it out! Before you start using this new features have a look at the tutorial. More information and examples can be found on my blog http://www.zhinst.com/blogs/schwizer/.

Tested browsers and editors for basic SVG support (no filters, no animation):
  Opera 9.64, 10.50, 10.63  -> yes
  Firefox 3.5, 3.6, 12.0, 15.0 -> yes
  Inkscape 0.46, 0.47, 0.48 -> yes
  Chrome 8.0, 18.0, 21.0 -> yes
  Internet Explorer 9.0 beta -> yes
  Internet Explorer 8.0 -> no
  Internet Explorer + RENESIS -> yes

Tested browsers and editors for SVG filters:
  Opera 9.64, 10.50, 10.63  -> yes
  Firefox 3.5, 3.6, 12.0, 15.0 -> yes
  Inkscape 0.46, 0.47, 0.48 -> yes (some limitations)
  Chrome 8.0, 18.0, 21.0 -> yes
  Internet Explorer 8.0, 9.0 beta -> no
  Internet Explorer + RENESIS -> no

Editors for the SVG file format can be found at http://www.inkscape.org.

Usage:
> plot2svg   % opens a file dialog to plot the active figure
    or
> plot2svg('myfile.svg', figure handle, pixelfiletype)
    
  pixelfiletype = 'png' (default), 'jpg'         

See http://www.zhinst.com/blogs/schwizer/ to get more informations

Supported Features
- line, patch, contour, contourf, quiver, surf, ...
- markers
- image (saved as linked png pictures)
- grouping of elements
- alpha values for patches
- subplot
- colorbar
- legend
- zoom
- reverse axes
- controls are saved as png pictures
- log axis scaling
- axis scaling factors (10^x)
- labels that contain Latex commands are interpreted (with some limitations):
\alpha, \Alpha, \beta, \Beta, ... \infity, \pm, \approx
{\it.....} for italic text
{\bf.....} for bold text
^{...} for superscript
_{...} for subscript

How to use SVG files in HTML code
<object type="image/svg+xml" data="./mySVGfile.svg" width="140" height="100"></object>

Changes in Version 22-May-2005
- bugfix line color
- bugfix path of linked jpeg figures
- improved patch handling (interpolation and texture still missing, preliminary depth sorting)
- support of pcolor plots
- preliminary: surface plots are projected on the xy-plane (use 'rotate' command)

Changes in Version 12-Dec-2005
- bugfix viewBox
- improvement of the axis scaling (many thanks to Bill Denney)
- improvement handling of exponents for log-plots
- default pixel format png instead of jpeg (many thanks to Bill Denney)
- bugfix axindex
- bugfix cell array cells (many thanks to Bill Denney)
- improved handling of pixel images (many thanks to Bill Denney)
- to save original figure background use set(gcf,'InvertHardcopy','off')
- improved markers

Changes in Version 8-Jan-2006
- axes handling fully reworked (3D axes)
- rework of axes scaling (3D axes)
- clipping enabled (Use carefully, as all figure data is written to file -> may get large)
- minor grid lines are now supported for linear and log plots
- linear color interpolation on patches (The interploation needs to be emulated as SVG does not support a linear interpolation of colors between three points. This is done by combination of different patches with linear alpha gradients. See limitation for Firefox 1.5.)

Changes in Version 20-Jun-2009
- Bugfix '°','±','µ','²','³','¼''½','¾','©''®'
- Bugfix 'projection' in hggroup and hgtransform
- Added Octave functionality (thanks to Jakob Malm)
  Bugfixe cdatamapping (thanks to Tom)
  Bugfix image data writing (thanks to Tom)
  Patches includes now markers as well (needed for 'scatter'
  plots (thanks to Phil)
- Bugfix markers for Octave (thanks to Jakob Malm)
- Bugfix image scaling and orientation
  Bugfix correct backslash (thanks to Jason Merril)
- Improvment of image handling (still some remaining issues)
  Fix for -. line style (thanks to Ritesh Sood)

Changes in Version 28-Jun-2009
- Improved depth sorting for patches and surface
- Bugfix patches
- Bugfix 3D axis handling

Changes in Version 11-Jul-2009
- Support of FontWeight and FontAngle properties
- Improved markers (polygon instead of polyline for closed markers)
- Added character encoding entry to be fully SVG 1.1 conform

Changes in Version 13-Jul-2009
- Support of rectangle for 2D
- Added preliminary support for SVG filters
- Added preliminary support for clipping with pathes
- Added preliminary support for turning axis tickmarks

Changes in Version 18-Jul-2009
- Line style scaling with line width (will not match with png
  output)
- Small optimizations for the text base line
- Bugfix text rotation versus shift
- Added more SVG filters
- Added checks for filter strings

Changes in Version 21-Jul-2009
- Improved bounding box calculation for filters
- Bugfixes for text size / line distance
- Support of background box for text
- Correct bounding box for text objects

Changes in Version 06-Mar-2010
- Improved support of filters
- Experimental support of animations
- Argument checks for filters
- Rework of tex string handling
- 'sub' and 'super' workaround for Firefox and Inkscape
- Bugfix for log axes (missing minor grid for some special
  cases)
- Bugfix nomy line #1102 (thanks to Pooya Jannaty)
- Bugfix minor tickmarks for log axis scaling (thanks to
  Harke Pera)
- Added more lex symbols
- Automatic correction of illegal axis scalings by the user
  (thanks to Juergen)
- Renamed plot2svg_beta to plot2svg

Changes in Version 12-04-2010
 - Improved Octave compatibility

Changes in Version 05-05-2010
- Bugfix for ticklabels outside of the axis limits (thanks to
  Ben Scandella)

Changes in Version 30-10-2010
- Improved handling of empty cells for labels (thanks to 
  Constantine)
- Improved HTML character coding (thanks to David Mack)
- Bugfix for last ')' (thanks to Jonathon Harding and Benjamin)
- Enabled scatter plots using hggroups
- Closing patches if they do not contain NaNs

Changes in Version 10-11-2010
- Support of the 'Layer' keyword to but the grid on top of 
  of the other axis content using 'top' (Many thanks to Justin
  Ashmall)
- Tiny optimization of the grid display at axis borders

Changes in Version 25-08-2011
- Fix for degree character (thanks to Manes Recheis)
- Fix for problems with dash-arrays in Inkscape (thanks to
  Rüdiger Stirnberg)
- Modified shape of driangles (thanks to Rüdiger Stirnberg)

Changes in Version 22-10-2011
- Removed versn as return value of function fileparts (thanks
  to Andrew Scott)
- Fix for images (thanks to Roeland)

Changes in Version 20-05-2012
- Added some security checks for empty data
- Fixed rotation for multiline text

Changes in Version 25-08-2012
- Special handling of 1xn char arrays for tick labels
  (thanks to David Plavcan)
- Fix for 'Index exceeds matrix dimensions' of axis labels
  (thanks to Aslak Grinsted)
- Fix for another axis label problem (thanks to Ben Mitch)

Changes in Version 15-09-2012
- Fix for linestyle none of rectangles (thanks to Andrew)
- Enabled scatter plot functionality

Limitations:
- axis scaling factors for 3D axes
- 3D plot functionality limited (depth sorting, light)

Example of a SVG file is included to the zip file.

Reports of bugs highly welcome. 
