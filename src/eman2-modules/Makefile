REDHAT.ROOT = $(CURDIR)/../../

-include $(ROCKSROOT)/etc/Rules.mk
include Rules.mk
include $(SDSCDEVEL)/Rules.mk

$(VERSION_INC): $(VERSION_SRC)
	/bin/grep 'VERSION.*=' $(VERSION_SRC) > $@
 
build:

install:: typical-modulefile-install
	module load $(ROLLCOMPILER) $(ROLLMPI) boost $(ROLLPY); \
	pyvers=`python --version 2>&1 | grep -o '[0-9][0-9]*\.[0-9][0-9]*'`; \
	perl -pi -e "s#BOOSTHOME#$${BOOSTHOME}#;" \
	         -e "s#PYVERS#$${pyvers}#;" \
	         -e "s#PYTHONROOT#$${PYTHONROOT}#" $(ROOT)/$(PKGROOT)/$(VERSION)

clean::
	rm -f $(VERSION_INC)
