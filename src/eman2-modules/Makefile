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
	module load $(ROLLCOMPILER) $(ROLLMPI) python boost; \
	perl -pi -e "s#BOOSTHOME#$${BOOSTHOME}#;" \
	         -e "s#PYVERS#$(PYVERSION)#;" \
	         -e "s#PYTHONROOT#$${PYTHONROOT}#" $(ROOT)/$(PKGROOT)/$(VERSION)

clean::
	rm -f $(VERSION_INC)
