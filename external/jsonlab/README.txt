===============================================================================
=                                 JSONLab                                     =
= compact, portable, robust JSON/binary-JSON encoder/decoder for MATLAB/Octave=
===============================================================================

* Copyright (c) 2011-2022  Qianqian Fang <q.fang at neu.edu>
* License: BSD or GNU General Public License version 3 (GPL v3), see License*.txt
* Version: 2.9.8 (Micronus Prime - Beta)
* URL: https://neurojson.org/jsonlab
* JData Specification Version: V1 Draft-3 (http://github.com/NeuroJSON/jdata)
* Binary JData Specification Version: V1 Draft-2 (http://github.com/NeuroJSON/bjdata)
* Compatibility: MATLAB R2008 or newer, GNU Octave 3.8 or newer
* Acknowledgement: This project is supported by US National Institute of Health (NIH) \
  grant [https://reporter.nih.gov/project-details/10308329 U24-NS124027]
-------------------------------------------------------------------------------

Table of Content:

0.   What's New
I.   Introduction
II.  Installation
III. Using JSONLab
IV.  Using `jsave/jload` to share workspace
V.   Sharing JSONLab created data files in Python
VI.  Known Issues and TODOs
VII. Contribution and feedback
VIII.Acknowledgement

-------------------------------------------------------------------------------

0. What's New

We are excited to announce that the JSONLab project, as the official reference library
for both the [https://neurojson.org/jdata/draft3 JData] and [https://neurojson.org/bjdata/draft2 BJData]
specifications, is funded by the US National Institute of Health (NIH) as
part of the NeuroJSON project (https://neurojson.org).
The goal of the NeuroJSON project is to develop human-readable, scalable and 
future-proof neuroimaging data standards and data sharing services. All data
produced from the NeuroJSON project will be using JSON/Binary JData formats as the
underlying serialization standards and use the lightweight JData specification as
language-independent data annotation standard, all of which have been evolved from
the over a decade development of JSONLab.

JSONLab v2.9.8 - code named "Micronus - beta" - is the beta-release of the next milestone -
JSONLab v3.0 - containing a number of key feature enhancement and bug fixes. The major 
new features include

# exporting JSON Memory-Map for rapid disk-map like JSON/binary JSON reading and writing, 
# supporting JSONPath query to MATLAB data and JSON/binary JSON file and streams, 
# ('''breaking''') upgrading the supported BJData spec to [https://neurojson.org/bjdata/draft2 V1 Draft 2] \
  where the default numerical data byte order changed from Big-Endian to '''Little-Endian'''
# adding initial support to JData [https://github.com/NeuroJSON/jdata/blob/master/JData_specification.md#data-referencing-and-links _DataLink_] \
  decoding to link multiple JSON/binary JSON files

There have been many major updates added to this release since the previous 
release v2.0 in June 2020. A list of the major changes are summarized below
(with key features marked by *), including the support to BJData Draft-2 specification,
new interface functions `savejd/loadjd`, and options to use MATLAB/Octave built-in
`jsonencode/jsondecode` functions. The `octave-jsonlab` package has also been
included in the official distributions of Debian Bullseye and Ubuntu 21.04 or newer.

* 2022-04-19*[2278bb1] stop escaping / to \/ in JSON string, see https://mondotondo.com/2010/12/29/the-solidus-issue/
* 2022-04-01*[fb711bb] add loadjd and savejd as the unified JSON/binary JSON file interface
* 2022-03-30 [4433a21] improve datalink uri handling to consider : inside uri
* 2022-03-30 [6368409] make datalink URL query more robust
* 2022-03-29 [dd9e9c6] when file suffix is missing, assume JSON feed
* 2022-03-29*[07c58f3] initial support for _DataLink_ of online/local file with JSONPath ref
* 2022-03-29 [897b7ba] fix test for older octave
* 2022-03-20 [bf03eff] force msgpack to use big-endian
* 2022-03-13 [46bbfa9] support empty name key, which is valid in JSON, fix #79
* 2022-03-12 [9ab040a] increase default float number digits from 10 to 16, fix #78
* 2022-03-11 [485ea29] update error message on the valid root-level markers
* 2022-02-23 [aa3913e] disable TFN marker in optimized header due to security risk and low benefit
* 2022-02-23 [f2c3223] support SCH{[ markers in optimized container type
* 2022-02-14 [540f95c] add optional preceding whitespace, explain format
* 2022-02-13 [3dfa904] debugged and tested mmap, add mmapinclude and mmapexclude options
* 2022-02-10*[6150ae1] handle uncompressed raw data (only base64 encoded) in jdatadecode
* 2022-02-10 [88a59eb] give a warning when jdatadecode fails, but still return the raw data
* 2022-02-03*[05edb7a] fast reading and writing json data record using mmap and jsonpath
* 2022-02-02*[b0f0ebd] return disk-map or memory-map table in loadjson
* 2022-02-01 [0888218] correct typos and add additional descriptions in README
* 2022-02-01*[03133c7] fix row-major ('formatversion',1.8) ND array storage order, update demo outputs
* 2022-02-01 [5998c70] revert variable name encoding to support unicode strings
* 2022-01-31 [16454e7] test flexible whitespaces in 1D/2D arrays, test mixed array from string
* 2022-01-31*[5c1ef15] accelerate fastarrayparser by 200%! jsonlab_speedtest cuts from 11s to 5.8s
* 2022-01-30 [9b25e20] fix octave 3.8 error on travis, it does not support single
* 2022-01-30 [5898f6e] add octave 5.2 to travis
* 2022-01-30*[2e3344c] [bjdata:breaking] Upgrade savebj/loadbj to BJData v1-draft 2, use little-endian by default
* 2022-01-30*[2e3344c] [bjdata:breaking] Fix optimized ND array element order (previously used column-major)
* 2022-01-30*[2e3344c] optimize loadjson and loadbj speed
* 2022-01-30*[2e3344c] add 'BuiltinJSON' option for savejson/loadjson to call jsonencode/jsondecode
* 2022-01-30*[2e3344c] more robust tests on ND array when parsing JSON numerical array construct
* 2021-06-23 [632531f] fix inconsistency between singlet integer and float values, close #70
* 2021-06-23 [f7d8226] prevent function calls when parsing array strings using eval, fix #75
* 2021-06-23 [b1ae5fa] fix #73 as a regression to #22
* 2021-11-22*[       ] octave-jsonlab is officially in Debian Testing/Bullseye
* 2020-09-29 [d0cb3b8] Fix for loading objects.
* 2020-07-26 [d0fb684] Add travis badge
* 2020-07-25 [708c36c] drop octave 3.2
* 2020-07-25 [436d84e] debug octave 3.2
* 2020-07-25 [0ce96ec] remove windows and osx targets from travis-ci
* 2020-07-25 [0d8baa4] fix ruby does not support error on windows
* 2020-07-25*[faa7921] enable travis-ci for jsonlab
* 2020-07-08 [321ab1a] add Debian and Ubuntu installation commands
* 2020-07-08 [e686828] update author info
* 2020-07-08*[ce40fdf] supports ND cell array, fix #66
* 2020-07-07 [6a8ce93] fix string encoding over 399 characters, close #65
* 2020-06-14 [5a58faf] fix DESCRIPTION date bug
* 2020-06-14 [9d7e94c] match octave description file and upstream version number
* 2020-06-14 [a5b6170] fix warning about lz4encode file name


Please note that the `savejson/loadjson` in both JSONLab v2.0-v3.0 are
compliant with JData Spec Draft 3; the savebj/loadbj` in JSONLab v3.0 is
compatible to BJData spec Draft 2, which contains breaking feature changes
compared to those in JSONLab v2.0.

The BJData spec was derived from UBJSON spec Draft 12, with the 
following breaking differences:

* BJData adds 4 new numeric data types: `uint16 [u]`, `uint32 [m]`, `uint64 [M]` \
  and `float16 [h]` (supported in JSONLab v2.0 or newer)
* BJData supports an optimized ND array container (supported in JSONLab since 2013)
* BJData does not convert `NaN/Inf/-Inf` to `null` (supported in JSONLab since 2013)
* BJData Draft 2 changes the default byte order to Little-Endian instead of Big-Endian (JSONLab 3.0 or later)
* BJData only permits non-zero-fixed-length data types as the optimized array type, i.e. only `UiuImlMLhdDC` are allowed

To avoid using the new features, one should attach `'UBJSON',1` and `'Endian','B'`
in the `savebj` command as

   savebj('',data,'FileName','myfile.bjd','UBJSON',1, 'Endian','B');

To read BJData data files generated by JSONLab v2.0, you should call

   data=loadbj('my_old_data_file.bjd','Endian','B')

You are strongly encouraged to convert all pre-v2.9 JSONLab generated BJD or .jamm
files using the new format.


-------------------------------------------------------------------------------

I.  Introduction

JSONLab is a free and open-source JSON/UBJSON/MessagePack encoder and 
decoder written in the native MATLAB language. It can be used to convert a MATLAB 
data structure (array, struct, cell, struct array, cell array, and objects) into 
JSON/UBJSON/MessagePack formatted strings and files, or to parse a 
JSON/UBJSON/MessagePack file into MATLAB data structure. JSONLab supports both 
MATLAB and [http://www.gnu.org/software/octave GNU Octave] (a free MATLAB clone).

Compared to other MATLAB/Octave JSON parsers, JSONLab is uniquely lightweight, 
ultra-portable, producing dependable outputs across a wide-range of MATLAB 
(tested on R2008) and Octave (tested on v3.8) versions. It also uniquely supports 
BinaryJData/UBJSON/MessagePack data files as binary-JSON-like formats, designed 
for efficiency and flexibility with loss-less binary storage. As a parser written
completely with the native MATLAB language, it is surprisingly fast when reading 
small-to-moderate sized JSON files (1-2 MB) with simple hierarchical structures,
and is heavily optimized for reading JSON files containing large N-D arrays
(known as the "fast array parser" in `loadjson`).

JSON ([http://www.json.org/ JavaScript Object Notation]) is a highly portable, 
human-readable and [http://en.wikipedia.org/wiki/JSON "fat-free"] text format 
to represent complex and hierarchical data, widely used for data-exchange in applications.
UBJSON ([http://ubjson.org/ Universal Binary JSON]) is a binary JSON format,  
specifically designed to specifically address the limitations of JSON, permitting 
efficient storage of binary data with strongly typed data records, resulting in smaller
file sizes and fast encoding and decoding. MessagePack is another binary
JSON-like data format widely used in data exchange in web/native applications.
It is slightly more compact than UBJSON, but is not directly readable compared
to UBJSON.

We envision that both JSON and its binary counterparts will play important 
rules not only for light-weight data storage, but also for storage and interchange
of scientific data. It has both the flexibility and generality as in other general-purpose 
file specifications, such as [http://www.hdfgroup.org/HDF5/whatishdf5.html HDF5] 
but has significantly reduced complexity and excellent readability.

Towards this goal, we have developed the JData Specification (http://github.com/NeuroJSON/jdata) 
to standardize serializations of complex scientific data structures, such as
N-D arrays, sparse/complex-valued arrays, trees, maps, tables and graphs using
JSON/binary JSON constructs. The text and binary formatted JData files are
syntactically compatible with JSON/UBJSON formats, and can be readily parsed 
using existing JSON and UBJSON parsers. JSONLab is not just a parser and writer 
of JSON/UBJSON data files, but one that systematically converts complex scientific
data structures into human-readable and universally supported JSON forms using the
standardized JData data annotations.

-------------------------------------------------------------------------------

II. Installation

The installation of JSONLab is no different from installing any other
MATLAB toolbox. You only need to download/unzip the JSONLab package
to a folder, and add the folder's path to MATLAB/Octave's path list
by using the following command:

    addpath('/path/to/jsonlab');

If you want to add this path permanently, you need to type `pathtool`, 
browse to the root folder of JSONLab and add to the list, then click "Save".
Then, run `rehash` in MATLAB, and type "which savejson", if you see an 
output, that means JSONLab is installed for MATLAB/Octave.

If you use MATLAB in a shared environment such as a Linux server, the
best way to add path is to type 

   mkdir ~/matlab/
   nano ~/matlab/startup.m

and type `addpath('/path/to/jsonlab')` in this file, save and exit the editor.
MATLAB will execute this file every time it starts. For Octave, the file
you need to edit is `~/.octaverc` , where `"~"` represents your home directory.

To use the data compression features, please download the ZMat toolbox from
https://github.com/fangq/zmat/releases/latest and follow the instruction to
install ZMat first. The ZMat toolbox is required when compression is used on 
MATLAB running in the `-nojvm` mode or GNU Octave, or 'lzma/lzip/lz4/lz4hc' 
compression methods are specified. ZMat can also compress large arrays that 
MATLAB's Java-based compression API does not support.


=== Install JSONLab on Fedora 24 or later ===

JSONLab has been available as an official Fedora package since 2015. You may
install it directly using the below command

   sudo dnf install octave-jsonlab

To enable data compression/decompression, you need to install `octave-zmat` using

   sudo dnf install octave-zmat


=== Install JSONLab on Debian ===

JSONLab is currently available on Debian Bullseye. To install, you may run

   sudo apt-get install octave-jsonlab

One can alternatively install `matlab-jsonlab` if MATLAB is available.

=== Install JSONLab on Ubuntu ===

JSONLab is currently available on Ubuntu 21.04 or newer as package
`octave-jsonlab`. To install, you may run

   sudo apt-get install octave-jsonlab

For older Ubuntu releases, one can add the below PPA
https://launchpad.net/~fangq/+archive/ubuntu/ppa

To install, please run

   sudo add-apt-repository ppa:fangq/ppa
   sudo apt-get update

to add this PPA, and then use

   sudo apt-get install octave-jsonlab

to install the toolbox. `octave-zmat` will be automatically installed.

=== Install JSONLab on Arch Linux ===

JSONLab is also available on Arch Linux. You may install it using the below command

   sudo pikaur -S jsonlab

-------------------------------------------------------------------------------

III.Using JSONLab

JSONLab provides a pair of functions, `loadjson` -- a JSON parser, and `savejson` -- 
a MATLAB-to-JSON encoder, to read/write text-based JSON; it also provides
three equivallent pairs -- `loadbj/savebj` for binary JData, `loadubjson/saveubjson`
for UBJSON and `loadmsgpack/savemsgpack` for MessagePack. The `load*` functions 
for the 3 supported data formats share almost the same input parameter format,
similarly for the 3 `save*` functions (`savejson/saveubjson/savemsgpack`).
These encoders and decoders are capable of processing/sharing almost all
data structures supported by MATLAB, thanks to `jdataencode/jdatadecode` - 
a pair of in-memory data converters translating complex data structures
to their easy-to-serialized forms according to the JData specifications.
The detailed help information can be found in the `Contents.m` file. 

In JSONLab 2.9.8 and later versions, a unified file loading and saving interface
is provided for JSON, binary JSON and HDF5, including `loadjd` and `savejd`
for reading and writing below files types:

* JSON based files: `.json`, `.jdt` (text JData file), `.jmsh` (text JMesh file), \
  `.jnii` (text JNIfTI file), `.jnirs` (text JSNIRF file)
* BJData based files: `.bjd`, `.jdb` (binary JData file), `.bmsh` (binary JMesh file), \
  `.bnii` (binary JNIfTI file), `.bnirs` (binary JSNIRF file), `.jamm` (MATLAB session file)
* UBJSON based files: `.ubj`
* MessagePack based files: `.msgpack`
* HDF5 based files: `.h5`, `.hdf5`, `.snirf` (SNIRF fNIRS data files) - require [https://github.com/fangq/easyh5 EasyH5 toolbox]

In the below section, we provide a few examples on how to us each of the 
core functions for encoding/decoding JSON/UBJSON/MessagePack data.

=== savejson.m ===

       jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 
                'MeshElem',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],...
                'MeshSurf',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;...
                           2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],...
                'MeshCreator','FangQ','MeshTitle','T6 Cube',...
                'SpecialData',[nan, inf, -inf]);
       savejson(jsonmesh)
       savejson('jmesh',jsonmesh)
       savejson('',jsonmesh,'compact',1)
       savejson('jmesh',jsonmesh,'outputfile.json')
       savejson('',jsonmesh,'ArrayIndent',0,'FloatFormat','\t%.5g','FileName','outputfile2.json')
       savejson('cpxrand',eye(5)+1i*magic(5))
       savejson('ziparray',eye(10),'Compression','zlib','CompressArraySize',1)
       savejson('',jsonmesh,'ArrayToStruct',1)
       savejson('',eye(10),'UseArrayShape',1)

=== loadjson.m ===

       loadjson('{}')
       dat=loadjson('{"obj":{"string":"value","array":[1,2,3]}}')
       dat=loadjson(['examples' filesep 'example1.json'])
       dat=loadjson(['examples' filesep 'example1.json'],'SimplifyCell',0)

=== savebj.m (saveubjson.m as an alias) ===

       a={single(rand(2)), struct('va',1,'vb','string'), 1+2i};
       savebj(a)
       savebj('rootname',a,'testdata.ubj')
       savebj('zeros',zeros(100),'Compression','gzip')

=== loadbj.m (loadubjson.m as an alias) ===

       obj=struct('string','value','array',single([1 2 3]),'empty',[],'magic',uint8(magic(5)));
       ubjdata=savebj('obj',obj);
       dat=loadbj(ubjdata)
       class(dat.obj.array)
       isequaln(obj,dat.obj)
       dat=loadbj(savebj('',eye(10),'Compression','zlib','CompressArraySize',1))

=== jdataencode.m ===

      jd=jdataencode(struct('a',rand(5)+1i*rand(5),'b',[],'c',sparse(5,5)))
      savejson('',jd)

=== jdatadecode.m ===

      rawdata=struct('a',rand(5)+1i*rand(5),'b',[],'c',sparse(5,5));
      jd=jdataencode(rawdata)
      newjd=jdatadecode(jd)
      isequaln(newjd,rawdata)


=== examples ===

Under the `examples` folder, you can find several scripts to demonstrate the
basic utilities of JSONLab. Running the `demo_jsonlab_basic.m` script, you 
will see the conversions from MATLAB data structure to JSON text and backward.
In `jsonlab_selftest.m`, we load complex JSON files downloaded from the Internet
and validate the loadjson/savejson functions for regression testing purposes.
Similarly, a `demo_ubjson_basic.m` script is provided to test the `saveubjson`
and `loadubjson` functions for various matlab data structures, and
`demo_msgpack_basic.m` is for testing `savemsgpack` and `loadmsgpack`.

Please run these examples and understand how JSONLab works before you use
it to process your data.


== unit testing ===

Under the `test` folder, you can find a script to test individual data types and
inputs using various encoders and decoders. This unit testing script also serves as
a '''specification validator''' to the JSONLab functions and ensure that the outputs
are compliant to the underlying specifications.


-------------------------------------------------------------------------------

IV.Using `jsave/jload` to share workspace

Starting from JSONLab v2.0, we provide a pair of functions, `jsave/jload` to store
and retrieve variables from the current workspace, similar to the `save/load` 
functions in MATLAB and Octave. The files that `jsave/jload` reads/writes is by  
default a binary JData file with a suffix `.jamm`. The file size is comparable
(can be smaller if use `lzma` compression) to `.mat` files. This feature
is currently experimental.

The main benefits of using .jamm file to share matlab variables include

* a `.jamm` file can be 50% smaller than a `.mat` file when using \
  `jsave(..., "compression","lzma")`; the only drawback is longer saving time.
* a `.jamm` file can be readily read/opened among many programming environments, including \
  Python, JavaScript, Go, Java etc, where `.mat` file support is not generally available. \
  Parsers of `.jamm` is largely compatible with UBJSON's parsers available at \
  http://ubjson.org/?page_id=48
* a `.jamm` file is quasi-human-readable, one can see the internal data fields \
  even in a command line, for example using `strings -n 2 file.jamm | astyle`, \
  making the binary data easy to be understood, shared and reused. 
* `jsave/jload` can also use MessagePack and JSON formats as the underlying \
  data storage format, addressing needs from diverse applications. \
  MessagePack parsers are readily available at https://msgpack.org/


=== jsave.m ===

      jsave    % save the current workspace to jamdata.jamm
      jsave mydata.jamm
      jsave('mydata.jamm','vars',{'var1','var2'})
      jsave('mydata.jamm','compression','lzma')
      jsave('mydata.json','compression','gzip')

=== jload.m ===

      jload    % load variables from jamdata.jamm to the current workspace
      jload mydata.jamm   % load variables from mydata.jamm
      vars=jload('mydata.jamm','vars',{'var1','var2'}) % return vars.var1, vars.var2
      jload('mydata.jamm','simplifycell',0)
      jload('mydata.json')

-------------------------------------------------------------------------------

V. Sharing JSONLab created data files in Python

Despite the use of portable data annotation defined by the JData Specification, 
the output JSON files created by JSONLab are 100% JSON compatible (with
the exception that long strings may be broken into multiple lines for better
readability). Therefore, JSONLab-created JSON files (`.json, .jnii, .jnirs` etc) 
can be readily read and written by nearly all existing JSON parsers, including
the built-in `json` module parser in Python.

However, we strongly recommend one to use a lightweight `jdata` module, 
developed by the same author, to perform the extra JData encoding and decoding
and convert JSON data directly to convenient Python/Numpy data structures.
The `jdata` module can also directly read/write UBJSON/Binary JData outputs
from JSONLab (`.bjd, .ubj, .bnii, .bnirs, .jamm` etc). Using binary JData
files are exptected to produce much smaller file sizes and faster parsing,
while maintainining excellent portability and generality.

In short, to conveniently read/write data files created by JSONLab into Python,
whether they are JSON based or binary JData/UBJSON based, one just need to download
the below two light-weight python modules:

* `jdata`: PyPi: https://pypi.org/project/jdata/  ; Github: https://github.com/NeuroJSON/pyjdata
* `bjdata` PyPi: https://pypi.org/project/bjdata/ ; Github: https://github.com/NeuroJSON/pybj

To install these modules on Python 2.x, please first check if your system has
`pip` and `numpy`, if not, please install it by running (using Ubuntu/Debian as example)

      sudo apt-get install python-pip python3-pip python-numpy python3-numpy

After the installation is done, one can then install the `jdata` and `bjdata` modules by

      pip install jdata --user
      pip install bjdata --user

To install these modules for Python 3.x, please replace `pip` by `pip3`.
If one prefers to install these modules globally for all users, simply
execute the above commands using `sudo` and remove the `--user` flag.

The above modules require built-in Python modules `json` and NumPy (`numpy`).

Once the necessary modules are installed, one can type `python` (or `python3`), and run

      import jdata as jd
      import numpy as np
      from collections import OrderedDict

      data1=jd.loadt('myfile.json',object_pairs_hook=OrderedDict);
      data2=jd.loadb('myfile.ubj',object_pairs_hook=OrderedDict);
      data3=jd.loadb('myfile.jamm',object_pairs_hook=OrderedDict);

where `jd.loadt()` function loads a text-based JSON file, performs
JData decoding and converts the enclosed data into Python `dict`, `list` 
and `numpy` objects. Similarly, `jd.loadb()` function loads a binary 
JData/UBJSON file and performs similar conversions. One can directly call
`jd.load()` to open JSONLab (and derived toolboxes such as '''jnifti''': 
https://github.com/NeuroJSON/jnifti or '''jsnirf''': https://github.com/NeuroJSON/jsnirf) 
generated files based on their respective file suffix.

Similarly, the `jd.savet()`, `jd.saveb()` and `jd.save` functions
can revert the direction and convert a Python/Numpy object into JData encoded
data structure and store as text-, binary- and suffix-determined output files,
respectively.

-------------------------------------------------------------------------------

VI. Known Issues and TODOs

JSONLab has several known limitations. We are striving to make it more general
and robust. Hopefully in a few future releases, the limitations become less.

Here are the known issues:

# 3D or higher dimensional cell/struct-arrays will be converted to 2D arrays;
# When processing names containing multi-byte characters, Octave and MATLAB \
can give different field-names; you can use feature('DefaultCharacterSet','latin1') \
in MATLAB to get consistant results
# `savejson` can only export the properties from MATLAB classes, but not the methods
# `saveubjson` converts a logical array into a uint8 ([U]) array
# a special N-D array format, as defined in the JData specification, is implemented in \
`saveubjson`. You may use `saveubjson(...,'NestArray',1)` to create UBJSON \
Draft-12 compliant files 
# loadubjson can not parse all UBJSON Specification (Draft 9) compliant \
files, however, it can parse all UBJSON files produced by `saveubjson`.

-------------------------------------------------------------------------------

VII. Contribution and feedback

JSONLab is an open-source project. This means you can not only use it and modify
it as you wish, but also you can contribute your changes back to JSONLab so
that everyone else can enjoy the improvement. For anyone who want to contribute,
please download JSONLab source code from its source code repositories by using the
following command:

      git clone https://github.com/fangq/jsonlab.git jsonlab

or browsing the github site at

      https://github.com/fangq/jsonlab

Please report any bugs or issues to the below URL:

      https://github.com/fangq/jsonlab/issues

Sometimes, you may find it is necessary to modify JSONLab to achieve your 
goals, or attempt to modify JSONLab functions to fix a bug that you have 
encountered. If you are happy with your changes and willing to share those
changes to the upstream author, you are recommended to create a pull-request
on github. 

To create a pull-request, you first need to "fork" jsonlab on Github by 
clicking on the "fork" button on top-right of jsonlab's github page. Once you forked
jsonlab to your own directory, you should then implement the changes in your
own fork. After thoroughly testing it and you are confident the modification 
is complete and effective, you can then click on the "New pull request" 
button, and on the left, select fangq/jsonlab as the "base". Then type
in the description of the changes. You are responsible to format the code
updates using the same convention (tab-width: 8, indentation: 4 spaces) as
the upstream code.

We appreciate any suggestions and feedbacks from you. Please use the following
mailing list to report any questions you may have regarding JSONLab:

      https://github.com/fangq/jsonlab/issues

(Subscription to the mailing list is needed in order to post messages).

-------------------------------------------------------------------------------

VIII.  Acknowledgement

This toolbox contains modified functions from the below toolboxes:

=== loadjson.m ===

The `loadjson.m` function was significantly modified from the earlier parsers 
(BSD 3-clause licensed) written by the below authors

* Nedialko Krouchev: http://www.mathworks.com/matlabcentral/fileexchange/25713
    created on 2009/11/02
* François Glineur: http://www.mathworks.com/matlabcentral/fileexchange/23393
    created on  2009/03/22
* Joel Feenstra:
    http://www.mathworks.com/matlabcentral/fileexchange/20565
    created on 2008/07/03


=== loadmsgpack.m ===

* Author: Bastian Bechtold
* URL: https://github.com/bastibe/matlab-msgpack/blob/master/parsemsgpack.m
* License: BSD 3-clause license

Copyright (c) 2014,2016 Bastian Bechtold
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this 
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, 
  this list of conditions and the following disclaimer in the documentation 
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its contributors 
  may be used to endorse or promote products derived from this software without 
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

== zlibdecode.m, zlibencode.m, gzipencode.m, gzipdecode.m, base64encode.m, base64decode.m ==

* Author: Kota Yamaguchi
* URL:https://www.mathworks.com/matlabcentral/fileexchange/39526-byte-encoding-utilities
* License: BSD License, see below

Copyright (c) 2012, Kota Yamaguchi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
