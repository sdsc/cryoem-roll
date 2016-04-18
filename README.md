# SDSC "cryoem" roll

## Overview

This roll bundles relion, eman2 and frealign

For more information about the various packages included in the cryoem roll please visit their official web pages:

- <a href="http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page" target="_blank">relion</a> is  is a stand-alone computer program that employs an empirical Bayesian approach to refinement of (multiple) 3D reconstructions or 2D class averages in electron cryo-microscopy (cryo-EM)
- <a href="http://grigoriefflab.janelia.org/frealign" target="_blank">frealign</a> is a program for high-resolution refinement of 3D reconstructions from cryo-EM images of single particles
- <a href="http://ncmi.bcm.tmc.edu/ncmi/software/software_details?selected_software=counter_222" target="_blank">eman2</a>  is a broadly based greyscale scientific image processing suite with a primary focus on processing data from transmission electron microscopes


## Requirements

To build/install this roll you must have root access to a Rocks development
machine (e.g., a frontend or development appliance).

If your Rocks development machine does *not* have Internet access you must
download the appropriate cryoem source file(s) using a machine that does
have Internet access and copy them into the `src/<package>` directories on your
Rocks development machine.


## Dependencies

The sdsc-roll must be installed on the build machine, since the build process
depends on make include files provided by that roll.

The roll sources assume that modulefiles provided by SDSC compiler and mpi
rolls are available, but it will build without them as long as the environment
variables they provide are otherwise defined.

The eman2 package depends on the boost and hdf rolls.
The relion package depends on the fftw roll.


## Building

To build the cryoem-roll, execute this on a Rocks development
machine (e.g., a frontend or development appliance):

```shell
% make 2>&1 | tee build.log
```

A successful build will create the file `cryoem-*.disk1.iso`.  If you built
the roll on a Rocks frontend, proceed to the installation step. If you built the
roll on a Rocks development appliance, you need to copy the roll to your Rocks
frontend before continuing with installation.

This roll source supports building with different compilers and for different
MPI flavors.  The `ROLLCOMPILER` and `ROLLMPI` make variables can be used to
specify the names of compiler and MPI modulefiles to use for building the
software, e.g.,

```shell
make ROLLCOMPILER=intel ROLLMPI=mvapich2_ib 2>&1 | tee build.log
```

The build process recognizes "gnu", "intel" or "pgi" as the value for the
`ROLLCOMPILER` variable; any MPI modulefile name may be used as the value of
the `ROLLMPI` variable.  The default values are "gnu" and "rocks-openmpi".

Note that eman2 is always built with the gnu compiler.
<ENDIF>

## Installation

To install, execute these instructions on a Rocks frontend:

```shell
% rocks add roll *.iso
% rocks enable roll cryoem
% cd /export/rocks/install
% rocks create distro
% rocks run roll cryoem | bash
```

In addition to the software itself, the roll installs cryoem environment
module files in:

```shell
/opt/modulefiles/applications/cryoem
```


## Testing

The cryoem-roll includes a test script which can be run to verify proper
installation of the roll documentation, binaries and module files. To
run the test scripts execute the following command(s):

```shell
% /root/rolltests/cryoem.t 
```


