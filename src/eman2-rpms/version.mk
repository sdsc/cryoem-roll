PACKAGE     = eman2

NAME        = sdsc-eman2
VERSION     = 2.21a
RELEASE     = 0
PKGROOT     = /opt/eman2

VERSION_SRC = $(REDHAT.ROOT)/src/$(PACKAGE)/version.mk
VERSION_INC = version.inc
include $(VERSION_INC)

RPM.EXTRAS  = AutoReq:No\nAutoProv:No\n%define __os_install_post /usr/lib/rpm/brp-compress\n%define         __prelink_undo_cmd     %{nil}\n%define __os_install_post /usr/lib/rpm/brp-python-bytecompile
RPM.PREFIX  = $(PKGROOT)
