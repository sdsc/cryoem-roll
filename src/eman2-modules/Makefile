REDHAT.ROOT = $(CURDIR)/../../

-include $(ROCKSROOT)/etc/Rules.mk
include Rules.mk
include $(SDSCDEVEL)/Rules.mk

ifndef PYVERSION
  PYVERSION = 2.6
endif

$(VERSION_INC): $(VERSION_SRC)
	/bin/grep 'VERSION.*=' $(VERSION_SRC) > $@
 
build:

install:: typical-modulefile-install

clean::
	rm -f $(VERSION_INC)
