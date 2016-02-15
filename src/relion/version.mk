ifndef ROLLCOMPILER
  ROLLCOMPILER = gnu
endif
COMPILERNAME := $(firstword $(subst /, ,$(ROLLCOMPILER)))

ifndef ROLLMPI
  ROLLMPI = rocks-openmpi
endif
MPINAME := $(firstword $(subst /, ,$(ROLLMPI)))

NAME           = sdsc-relion
VERSION        = 1.4
RELEASE        = 0
PKGROOT        = /opt/relion

SRC_SUBDIR     = relion

SOURCE_NAME    = relion
SOURCE_SUFFIX  = tar.bz2
SOURCE_VERSION = $(VERSION)
SOURCE_PKG     = $(SOURCE_NAME)-$(SOURCE_VERSION).$(SOURCE_SUFFIX)
SOURCE_DIR     = $(SOURCE_PKG:%.$(SOURCE_SUFFIX)=%)

VFLTK          = fltk-1.3.0

TAR_BZ2_PKGS       = $(SOURCE_PKG)

RPM.EXTRAS     = AutoReq:No
