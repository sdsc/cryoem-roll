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
RELEASE        = 9
PKGROOT        = /opt/eman2

SRC_SUBDIR     = eman2

SOURCE_NAME    = eman2
SOURCE_SUFFIX  = zip
SOURCE_VERSION = $(VERSION)
SOURCE_PKG     = $(SOURCE_NAME)-$(SOURCE_VERSION).$(SOURCE_SUFFIX)
SOURCE_DIR     = $(SOURCE_PKG:%.$(SOURCE_SUFFIX)=%)

BOOST_NAME     = boost
BOOST_SUFFIX   = tar.gz
BOOST_VERSION  = 1.55.0
BOOST_PKG      = $(BOOST_NAME)-$(BOOST_VERSION).$(BOOST_SUFFIX)
BOOST_DIR      = $(BOOST_PKG:%.$(BOOST_SUFFIX)=%)

FTGL_NAME      = ftgl
FTGL_SUFFIX    = tar.gz
FTGL_VERSION   = 2.1.3-rc5
FTGL_PKG       = $(FTGL_NAME)-$(FTGL_VERSION).$(FTGL_SUFFIX)
FTGL_DIR       = $(subst -rc,~rc,$(FTGL_PKG:%.$(FTGL_SUFFIX)=%))

ZIP_PKGS       = $(SOURCE_PKG)
TAR_GZ_PKGS    = $(BOOST_PKG) $(FTGL_PKG)

RPM.EXTRAS     = AutoReq:No\nAutoProv:No\n%define __os_install_post /usr/lib/rpm/brp-compress
RPM.PREFIX     = $(PKGROOT)
