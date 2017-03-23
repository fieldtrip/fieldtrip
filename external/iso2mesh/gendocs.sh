#!/bin/sh

#============================================================
#  iso2mesh inline documentation to wiki convertor
#
#  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
#============================================================

print_help()
{
   awk '/^%/ {dp=1} /-- this function is part of iso2mesh/ {exit} \
        /-- this function is part of JSONLab/ {exit} \
        /^function/ {dp=1} /./ {if(dp==1) print;}' $1 \
     | grep -v 'Qianqian' | grep -v 'date:' | grep -v '^%\s*$'| \
     sed -e 's/^%//g' -e 's/^function\(.*$\)/\n==== function\1 ====/g'
}
print_group()
{
   for fun in $@
   do 
      print_help $fun.m
   done
   echo ''
}

func_shortcut="v2m v2s s2m s2v m2v sms"
func_mainfun="vol2mesh vol2surf surf2mesh surf2vol mesh2vol"
func_backend="binsurface cgalv2m cgals2m vol2restrictedtri surf2volz mesh2mask"
func_primitive="meshabox meshasphere meshanellip meshunitsphere meshacylinder 
                meshgrid5 meshgrid6 latticegrid"
func_inquery="finddisconnsurf surfedge volface extractloops meshconn  
                meshcentroid nodevolume elemvolume neighborelem 
		faceneighbors edgeneighbors maxsurf flatsegment orderloopedge  
		mesheuler bbxflatsegment surfplane surfinterior surfpart
                surfseeds meshquality meshedge surfacenorm uniqedges advancefront"
func_meshfix="meshcheckrepair meshreorient removedupelem 
                removedupnodes removeisolatednode removeisolatedsurf
                surfaceclean getintersecttri delendelem surfreorient"
func_remesh="meshresample remeshsurf smoothsurf sortmesh mergemesh 
                meshrefine mergesurf surfboolean"
func_fileio="saveasc savedxf savestl savebinstl saveinr saveoff 
                savesmf savesurfpoly savegts readgts
                savevrml readasc readinr readmedit readoff readsmf
	        readtetgen deletemeshfile mcpath mwpath savemedit
		savejson loadjson saveubjson loadubjson savejmesh
                savemphtxt savetetgenele savetetgennode saveabaqus"
func_binimage="bwislands fillholes3d deislands2d deislands3d 
                imedge3d internalpoint smoothbinvol 
		thickenbinvol thinbinvol"
func_plotting="plotmesh plotsurf plottetra plotedges qmeshcut"
func_misc="surfdiffuse volmap2mesh isoctavemesh getvarfrom raytrace
		getplanefrom3pt getexeext fallbackexeext iso2meshver
                raysurf getoptkey rotatevec3d rotmat2vec varargin2struct
                jsonopt mergestruct orthdisk struct2jdata"

echo === "#" Streamlined mesh generation - shortcuts ===
print_group $func_shortcut

echo === "#" Streamlined mesh generation ===
print_group $func_mainfun

echo === "#" iso2mesh main function backend ===
print_group $func_backend

echo === "#" iso2mesh primitive meshing functions ===
print_group $func_primitive

echo === "#" Mesh decomposition and query ===
print_group $func_inquery

echo === "#" Mesh processing and reparing ===
print_group $func_meshfix

echo === "#" Mesh resampling and optimization ===
print_group $func_remesh

echo === "#" File I/O ===
print_group $func_fileio

echo === "#" Volumetric image pre-processing ===
print_group $func_binimage

echo === "#" Mesh plotting ===
print_group $func_plotting

echo === "#" Miscellaneous functions ===
print_group $func_misc

