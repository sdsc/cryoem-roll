override ROLLCOMPILER = gnu
COMPILERNAME := $(firstword $(subst /, ,$(ROLLCOMPILER)))

ifndef ROLLMPI
  ROLLMPI = rocks-openmpi
endif
MPINAME := $(firstword $(subst /, ,$(ROLLMPI)))

ifndef ROLLPY
  ROLLPY = python
endif

NAME           = sdsc-eman2
VERSION        = 2.1
RELEASE        = 6
PKGROOT        = /opt/eman2

SRC_SUBDIR     = eman2

SOURCE_NAME    = eman2
SOURCE_SUFFIX  = zip
SOURCE_VERSION = $(VERSION)
SOURCE_PKG     = $(SOURCE_NAME)-$(SOURCE_VERSION).$(SOURCE_SUFFIX)
SOURCE_DIR     = $(SOURCE_PKG:%.$(SOURCE_SUFFIX)=%)

FTGL_NAME    = ftgl
FTGL_SUFFIX  = tar.gz
FTGL_VERSION = 2.1.3-rc5
FTGL_PKG     = $(FTGL_NAME)-$(FTGL_VERSION).$(FTGL_SUFFIX)
FTGL_DIR     = $(FTGL_PKG:%.$(FTGL_SUFFIX)=%)

ZIP_PKGS       = $(SOURCE_PKG)
TAR_GZ_PKGS    = $(FTGL_PKG)

RPM.EXTRAS     = AutoReq:No
