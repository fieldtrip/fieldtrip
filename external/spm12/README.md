```
   ___  ____  __  __
  / __)(  _ \(  \/  )  
  \__ \ )___/ )    (   Statistical Parametric Mapping
  (___/(__)  (_/\/\_)  SPM - https://www.fil.ion.ucl.ac.uk/spm/
```
 
This README gives a brief introduction to the SPM software. Full details can be found
on the [SPM website](https://www.fil.ion.ucl.ac.uk/spm/).
 
See also `Contents.m`, `AUTHORS.txt` and `LICENCE.txt`.
  
# SPM
 
Statistical Parametric Mapping refers to the construction and assessment of spatially
extended statistical process used to test hypotheses about functional imaging data.
These ideas have been instantiated in software that is called SPM.  The SPM software
package has been designed for the analysis of brain imaging data sequences.  The
sequences can be a series of images from different cohorts, or time-series from the
same subject.  The current release is designed for the analysis of fMRI, PET, SPECT,
EEG and MEG.

Please refer to this version as "**SPM12**" in papers and communications.

SPM was written to organise and interpret our data (at The Wellcome Centre for Human
Neuroimaging).  The distributed version is the same as that we use ourselves.

SPM is made freely available to the [neuro]imaging community, to promote
collaboration and a common analysis scheme across laboratories.

# Software

The SPM software is a suite of MATLAB functions, scripts and data files, with some
externally compiled C routines, implementing Statistical Parametric Mapping.  MATLAB,
a commercial engineering mathematics package, is required to use SPM.  MATLAB is
produced by [MathWorks, Natick, MA, USA](https://www.mathworks.com/).

SPM requires only core MATLAB to run (no special toolboxes are required).

SPM12 is written for MATLAB version 7.4 (R2007a) onwards under Windows, Linux and Mac
(SPM12 will not work with versions of MATLAB prior to 7.4).  Binaries of the external
C-MEX routines are provided for Windows, Linux and Mac.  The source code is supplied
and can be compiled with a C compiler (Makefile provided).

See https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ for details.

Later versions of MATLAB (released after SPM12), will probably need additional
patches in order to run.  Once developed, these will be made available from:
https://www.fil.ion.ucl.ac.uk/spm/download/spm12_updates/

Although SPM12 will read image files from previous versions of SPM, there are
differences in the algorithms, templates and models used.  Therefore, we recommend
you use a single SPM version for any given project.

The SPM12 Release Notes can be found online: 
https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

# File format

SPM12 uses the NIFTI-1 data format as standard.  Take a look at
https://nifti.nimh.nih.gov/ for more information on the NIFTI-1 file format.

The old SPM2 version of Analyze format can be read straight into SPM12, but results
will be written out as NIFTI-1.  If you still use this format, then it is important
that you ensure that `spm_flip_analyze_images` has been set appropriately for your
data.

The MINC and ECAT7 formats can not be read straight into SPM12, although conversion
utilities have been provided.  Similarly, a number of DICOM flavours can also be
converted to NIFTI-1 using tools in SPM12.

# Resources

The SPM website is the central repository for SPM resources:
https://www.fil.ion.ucl.ac.uk/spm/

Introductory material, installation details, documentation, course details and
patches are published on the site.

There is an SPM email discussion list, hosted at <spm@jiscmail.ac.uk>.  The list is
monitored by the authors, and discusses theoretical, methodological and practical
issues of Statistical Parametric Mapping and SPM.  The SPM website has further
details:
https://www.fil.ion.ucl.ac.uk/spm/support/

Please report bugs to the authors at <fil.spm@ucl.ac.uk>.

Peculiarities may actually be features, and should be raised on the SPM  email
discussion list, <spm@jiscmail.ac.uk>.

# Authors

SPM is developed under the auspices of Functional Imaging Laboratory (FIL), The
Wellcome Centre for Human NeuroImaging, in the Queen Square Institute of Neurology at
University College London (UCL), UK.

SPM94 was written primarily by Karl Friston in the first half of 1994, with
assistance from John Ashburner (MRC-CU), Jon Heather (WDoIN), and Andrew Holmes
(Department of Statistics, University of Glasgow).  Subsequent development, under the
direction of Prof. Karl Friston at the Wellcome Department of Imaging Neuroscience,
has benefited from substantial input (technical and theoretical) from: John Ashburner
(WDoIN), Andrew Holmes (WDoIN & Robertson Centre for Biostatistics, University of
Glasgow, Scotland), Jean-Baptiste Poline (WDoIN & CEA/DRM/SHFJ, Orsay, France),
Christian Buechel (WDoIN), Matthew Brett (MRC-CBU, Cambridge, England), Chloe Hutton
(WDoIN) and Keith Worsley (Department of Statistics, McGill University, Montreal
Canada).

See `AUTHORS.txt` for a complete list of SPM co-authors.

We would like to thank everyone who has provided feedback on SPM.

# Disclaimer, copyright & licencing

SPM (being the collection of files given in the manifest in the `Contents.m` file) is
free but copyright software, distributed under the terms of the GNU General Public
Licence as published by the Free Software Foundation (either version 2, as given in
file `LICENCE.txt`, or at your option, any later version).  Further details on
"copyleft" can be found at https://www.gnu.org/copyleft/.  In particular, SPM is
supplied as is.  No formal support or maintenance is provided or implied.

```
Copyright (C) 1991,1994-2020 Wellcome Centre for Human Neuroimaging
$Id: README.md 7765 2020-01-02 16:29:48Z spm $
```
