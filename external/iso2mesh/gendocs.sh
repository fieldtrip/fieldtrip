#!/bin/sh

#============================================================
#  iso2mesh inline documentation to wiki convertor
#
#  Author: Qianqian Fang <q.fang at neu.edu>
#============================================================

printhelp=$1
if [ "$#" -ge 1 ]; then
    echo "iso2mesh >> Iso2Mesh"
fi
print_help()
{
   if [ -z "$printhelp" ]; then
     awk '/^%/ {dp=1} /-- this function is part of iso2mesh/ {exit} \
        / this function is part of / {exit} \
        / this file is part of / {exit} \
        /^function/ {dp=1} /./ {if(dp==1) print;}' $1 \
       | grep -v 'Qianqian' | grep -v 'date:' | #grep -v '^%\s*$'| \
       sed -e 's/^%//g' -e 's/^function\(.*$\)/\n==== function\1 ====/g'
    else
	echo " $1" | sed -e 's/\.m$//g'
    fi
}
print_group()
{
   for fun in $@
   do 
      print_help $fun.m
   done
}
print_title ()
{
   if [ -z "$printhelp" ]; then
       echo "=== # $@ ==="
   else
       echo $@
   fi
}

func_shortcut="v2m v2s s2m s2v m2v sms i2m"
func_mainfun="vol2mesh vol2surf surf2mesh surf2vol mesh2vol img2mesh"
func_backend="binsurface cgalv2m cgals2m vol2restrictedtri surf2volz mesh2mask"
func_primitive="meshabox meshasphere meshanellip meshunitsphere meshacylinder 
                meshgrid5 meshgrid6 latticegrid extrudecurve meshcylinders"
func_inquery="finddisconnsurf surfedge volface extractloops meshconn  
                meshcentroid nodevolume elemvolume neighborelem layersurf 
		faceneighbors edgeneighbors maxsurf flatsegment orderloopedge  
		mesheuler bbxflatsegment surfplane surfinterior surfpart
                surfseeds meshquality meshedge meshface surfacenorm nodesurfnorm
                uniqedges uniqfaces advancefront innersurf outersurf surfvolume
		insurface"
func_meshfix="meshcheckrepair meshreorient removedupelem 
                removedupnodes removeisolatednode removeisolatedsurf
                surfaceclean getintersecttri delendelem surfreorient"
func_metch="proj2mesh dist2surf regpt2surf affinemap metchgui metchgui_one"
func_remesh="meshresample remeshsurf smoothsurf sortmesh mergemesh 
                meshrefine mergesurf surfboolean fillsurf highordertet
		elemfacecenter barydualmesh meshinterp meshremap extrudesurf"
func_polyline="slicesurf slicesurf3 polylinelen polylinesimplify polylineinterp closestnode"
func_fileio="saveasc savedxf savestl savebinstl saveinr saveoff 
                savesmf savesurfpoly savegts readgts savemsh
                savevrml readasc readinr readmedit readoff readsmf
	        readtetgen deletemeshfile mcpath mwpath savemedit
		savejson loadjson saveubjson loadubjson loadmsgpack savemsgpack
                savebj loadbj savemphtxt savetetgenele savetetgennode saveabaqus
		savenirfast readnirfast readnifti readmptiff"
func_jdata="savejmesh loadjnifti savejnifti loadnifti savenifti jdataencode 
                jdatadecode jload jsave decodevarname encodevarname jnifticreate
                nifticreate nii2jnii jnii2nii niicodemap niiformat savebnii savejnii"
func_compression="zlibencode zlibdecode gzipencode gzipdecode lzmaencode 
                lzmadecode lzipencode lzipdecode lz4encode lz4decode lz4hcencode 
		lz4hcdecode base64decode base64encode"
func_binimage="bwislands fillholes3d deislands2d deislands3d ndgaussian ndimfilter
                imedge3d internalpoint smoothbinvol 
		thickenbinvol thinbinvol maskdist"
func_plotting="plotmesh plotsurf plottetra plotedges qmeshcut"
func_misc="surfdiffuse volmap2mesh isoctavemesh getvarfrom raytrace linextriangle
		getplanefrom3pt getexeext fallbackexeext iso2meshver
                raysurf getoptkey rotatevec3d rotmat2vec varargin2struct
                jsonopt mergestruct orthdisk nestbracket2dim fast_match_bracket 
		match_bracket memmapstream"

print_title Streamlined mesh generation - shortcuts 
print_group $func_shortcut

print_title Streamlined mesh generation 
print_group $func_mainfun

print_title Iso2mesh main function backend 
print_group $func_backend

print_title Iso2mesh primitive meshing functions 
print_group $func_primitive

print_title Mesh decomposition and query 
print_group $func_inquery

print_title Mesh processing and reparing 
print_group $func_meshfix

print_title Mesh registration - Metch Toolbox 
print_group $func_metch

print_title Polyline handling 
print_group $func_polyline

print_title Mesh resampling and optimization 
print_group $func_remesh

print_title File I/O 
print_group $func_fileio

print_title JData functions 
print_group $func_jdata

print_title Data compression 
print_group $func_compression

print_title Volumetric image pre-processing 
print_group $func_binimage

print_title Mesh plotting 
print_group $func_plotting

print_title Miscellaneous functions 
print_group $func_misc

