# MRtrix

[![Build Status](https://travis-ci.org/MRtrix3/mrtrix3.svg?branch=updated_syntax)](https://travis-ci.org/MRtrix3/mrtrix3)

Please visit the [official website for MRtrix](http://www.mrtrix.org) to access the [documentation for MRtrix3](http://mrtrix.readthedocs.org/), including detailed installation instructions.

## Getting help

Support and general discussion is hosted on the [MRtrix3 Community Forum](http://community.mrtrix.org/). Please
address all MRtrix3-related queries there. You can use you GitHub or Google login to post questions.

## Quick install

1. Install dependencies by whichever means your system uses. 
   These include: Python (>=2.6), a C++ compiler with full C++11 support (`g++` 4.9 or later, `clang++`), 
   Eigen (>=3.2.8), zlib, OpenGL (>=3.3), and Qt (>=4.8, or at least 5.1 on MacOSX).

2. Clone Git repository and compile:

        $ git clone https://github.com/MRtrix3/mrtrix3.git
        $ cd mrtrix3/
        $ ./configure
        $ ./build

3. Set the `PATH`:

    * Bash shell:

      run the `set_path` script provided:

            $ ./set_path

      or edit the startup `~/.bashrc` or `/etc/bash.bashrc` file manually by adding this line:

            $ export PATH=/<edit as appropriate>/mrtrix3/bin:$PATH

    * C shell:

      edit the startup `~/.cshrc` or `/etc/csh.cshrc` file manually by adding this line:

            $ setenv PATH /<edit as appropriate>/mrtrix3/bin:$PATH

4. Test installation:

    Command-line:

        $ mrconvert

    GUI:

        $ mrview

## Keeping MRtrix3 up to date

1. You can update your installation at any time by opening a terminal in the mrtrix3 folder, and typing:

        git pull
        ./build

2. If this doesn't work immediately, it may be that you need to re-run the configure script:

        ./configure

    and re-run step 1 again.

## Building a specific release of MRtrix3

You can build a particular release of MRtrix3 by checking out the corresponding _tag_, and using the same procedure as above to build it:

    git checkout 3.0_RC3
    ./configure
    ./build
