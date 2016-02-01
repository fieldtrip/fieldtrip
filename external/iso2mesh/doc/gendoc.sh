#!/bin/sh

# commands to update the document pages from homepage

lynx -dont_wrap_pre -dump "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?keywords=Download&embed=1" > Download_and_License.txt
lynx -dont_wrap_pre -dump "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?keywords=Doc/Installation&embed=1" > INSTALL.txt
lynx -dont_wrap_pre -dump "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?keywords=Doc/Basics&embed=1" > Get_Started.txt
lynx -dont_wrap_pre -dump "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?keywords=Doc/FAQ&embed=1" > FAQ.txt
lynx -dont_wrap_pre -dump "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?keywords=Advanced&embed=1" > Advanced_Features.txt

wget http://iso2mesh.sourceforge.net/upload/iso2mesh_workflow_v2.png -Oiso2mesh_workflow.png
