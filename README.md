# Platelet analysis

This code implements the methods proposed in the following paper:

[1] Johanna C. Clauser, Judith Maas, Jutta Arens Thomas Schmitz-Rode, Ulrich Steinseifer, and Benjamin Berkels. Automation of Hemocompatibility Analysis Using Image Segmentation and Supervised Classification. *Engineering Applications of Artificial Intelligence*, 97, January 2021. [[DOI](https://dx.doi.org/10.1016/j.engappai.2020.104009) | [arXiv](https://arxiv.org/abs/2010.06245)]

It is based on the [QuocMesh software library](https://archive.ins.uni-bonn.de/numod.ins.uni-bonn.de/software/quocmesh/index.html), AG Rumpf, Institute for Numerical Simulation, University of Bonn, and distributed under the terms of the [Common Development and Distribution License](LICENSE.txt). Particularly, you can **download and use** this software **free of cost**.

Some subdirectories contain code from other software projects, which is redistributed under the respective licenses:

* `quocmesh/external/bz2/`

  Code from the [bzip2 library by Julian R Seward](http://www.bzip.org/), License under `external/bz2/LICENSE`

* `quocmesh/external/dirent/dirent.h`

  Code from the [Dirent API by Toni Ronkko](https://github.com/tronkko/dirent), License at the beginning of external/dirent/dirent.h

* `quocmesh/cmake/FindSUITESPARSE.cmake`

  CMake module from [OpenFlipper](https://www.openflipper.org/), GNU Lesser General Public License

We appreciate any feedback on your experience with our methods. In case you encounter any problems when using this software, please do not hesitate to contact us: <berkels@aices.rwth-aachen.de>

## Prerequisites
To build the software, cmake and a compiler like GCC or clang need to be installed. Furthermore, a few external libraries are needed. Under Linux, usually all of the dependencies can be installed with the package manager that comes with the Linux distribution. Under OS X/macOS, the installation of a third party package manager like [MacPorts](https://www.macports.org/) or [Homebrew](https://brew.sh/) is recommended. For instance, if MacPorts is installed, one can install cmake with:

    sudo port install cmake

## Compiling

First create a clone of our git repository e.g. by

    git clone https://github.com/berkels/platelet-analysis

In the clone directory, there should be two directories, `quocmesh` and `quocGCC`. The source code is in `quocmesh`, the compiled binaries will go into `quocGCC`.

To compile the code, go into the directory `quocGCC` and call the script called `goLinux.sh` (or `goWindows.bat` if compiling under Windows with MinGW). This invokes cmake with all necessary options and will create makefiles.

After the go script has been run, the source can be compiled with `make` like any other project with makefiles. There is also a test target to check whether the basics of the code actually work.

In short, to get started the following steps should be enough: Decompress the archive and change into the directory `quocGCC`. In there, call

    ./goLinux.sh
    make
    make test

Under Windows with MinGW call

    ./goWindows.bat
    mingw32-make
    mingw32-make test

instead.

If all dependencies are installed correctly, this should work without any errors (although there may be warnings during the first two steps).

To speed up the compilation, one can have make use as many threads as desired, by calling `make -jN` (where `N` is the number of threads) instead of just `make`.

## Running
TODO

Note: So far this repository only contains the segmentation code. The remaining parts will be added soon.
