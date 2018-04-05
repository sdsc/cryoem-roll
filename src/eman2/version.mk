override ROLLCOMPILER = gnu
COMPILERNAME := $(firstword $(subst /, ,$(ROLLCOMPILER)))

NAME                = sdsc-eman2
VERSION             = 2.21a
RELEASE             = 0
PKGROOT             = /opt/eman2

SRC_SUBDIR          = eman2

SOURCE_NAME         = eman2
SOURCE_SUFFIX       = zip
SOURCE_VERSION      = $(VERSION)
SOURCE_PKG          = $(SOURCE_NAME)-$(SOURCE_VERSION).$(SOURCE_SUFFIX)
SOURCE_DIR          = $(SOURCE_PKG:%.$(SOURCE_SUFFIX)=%)

MINICONDA_NAME      = miniconda
MINICONDA_SUFFIX    = tar.gz
MINICONDA_VERSION   = 4.4.10
MINICONDA_PKG       = $(MINICONDA_NAME)-$(MINICONDA_VERSION).$(MINICONDA_SUFFIX)
MINICONDA_DIR       = $(MINICONDA_PKG:%.$(MINICONDA_SUFFIX)=%)

CONDADOWNLOADS      = backports-1.0-py27h63c9359_1 backports.shutil_get_terminal_size-1.0.0-py27h5bc021e_2 binutils_impl_linux-64-2.28.1-h04c84fa_2 binutils_linux-64-7.2.0-25 boost-1.63.0-py27_5 bzip2-1.0.6-h9a117a8_4 boost-cpp-1.63.0-2 bsddb-1.0-py27_1 ca-certificates-2018.03.07-0 cairo-1.12.18-1 decorator-4.2.1-py27_0 db-5.3.28-1 eman-deps-7.0-np19py27_0 fftw-mpi-3.3.6-2 freetype-2.5.2-1 ftgl-2.1.3-1 gcc_impl_linux-64-7.2.0-hc5ce805_2 gcc_linux-64-7.2.0-25 gsl-2.2.1-h0c605f7_3 gxx_impl_linux-64-7.2.0-hd3faf3d_2 gxx_linux-64-7.2.0-25 hdf5-1.10.1-h9caa474_1 icu-58.2-h9c2bf20_1 intel-openmp-2018.0.0-hc7b2577_8 ipython-5.4.1-py27_2 ipython_genutils-0.2.0-py27h89fb69b_0 jpeg-9b-h024ee3a_2 libgcc-7.2.0-h69d50b8_2 libgfortran-ng-7.2.0-h9f7466a_2 libgpuarray-0.7.5-h14c3975_0 libpng-1.5.13-1 libtiff-4.0.9-h28f6b97_0 mako-1.0.7-py27h3d58d4b_0 markupsafe-1.0-py27h97b2822_1 matplotlib-1.4.3-np19py27_1 mkl-2018.0.1-h19d6760_4 mkl-service-1.1.2-py27hb2d42c5_4 numpy-1.9.3-py27h7e35acb_3 openmpi-2.0.2-0 pathlib2-2.3.0-py27h6e9d198_0 openssl-1.0.2o-h20670df_0 pexpect-4.4.0-py27_0 pickleshare-0.7.4-py27h09770e1_0 pixman-0.26.2-0 prompt_toolkit-1.0.15-py27h1b593e1_0 ptyprocess-0.5.2-py27h4ccb14c_0 pydusa-1.15-np19_9 pygments-2.2.0-py27h4a8b6f5_0 pyopengl-3.1.0-py27_0 pygpu-0.7.5-py27h14c3975_0 pyparsing-2.0.3-py27_0 pyqt-4.10.4-py27_0 py2cairo-1.10.0-py27_2 python-dateutil-2.6.1-py27h4ca5741_1 pytz-2018.3-py27_0 qt-4.8.5-0 scandir-1.7-py27h14c3975_0 scikit-learn-0.19.1-py27h445a80a_0 scipy-1.0.0-py27hf5f0f52_0 simplegeneric-0.8.1-py27_2 sip-4.15.5-py27_0 theano-1.0.1-py27h6bb024c_0 traitlets-4.3.2-py27hd6ce930_0 wcwidth-0.1.7-py27h9e3e1ab_0 xz-5.2.3-h55aa19d_2

CONDADOWNLOAD_PKGS = $(addsuffix .tar.bz2,$(CONDADOWNLOADS) )

ZIP_PKGS            = $(SOURCE_PKG)
TAR_GZ_PKGS         = $(MINICONDA_PKG)

RPM.EXTRAS     = AutoReq:No\nAutoProv:No\n%global __os_install_post %(echo '%{__os_install_post}' | sed -e 's!/usr/lib[^[:space:]]*/brp-python-bytecompile[[:space:]].*$!!g')
RPM.PREFIX     = $(PKGROOT)
