PACKAGE     = eman2

NAME        = sdsc-eman2
VERSION     = 2.21a
RELEASE     = 0
PKGROOT     = /opt/eman2

VERSION_SRC = $(REDHAT.ROOT)/src/$(PACKAGE)/version.mk
VERSION_INC = version.inc
include $(VERSION_INC)

RPM.EXTRAS  = AutoReq:No\nAutoProv:No\n%global __os_install_post %(echo '%{__os_install_post}' | sed -e 's!/usr/lib[^[:space:]]*/brp-python-bytecompile[[:space:]].*$!!g')
RPM.PREFIX  = $(PKGROOT)
