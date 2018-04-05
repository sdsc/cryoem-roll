PACKAGE     = eman2

NAME        = sdsc-eman2
VERSION     = 2.21a
RELEASE     = 0
PKGROOT     = /opt/eman2

VERSION_SRC = $(REDHAT.ROOT)/src/$(PACKAGE)/version.mk
VERSION_INC = version.inc
include $(VERSION_INC)

RPM.EXTRAS  = AutoReq:No\nAutoProv:No
RPM.PREFIX  = $(PKGROOT)
