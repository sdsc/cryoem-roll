#
# @Copyright@
#
#                               Rocks(r)
#                        www.rocksclusters.org
#                        version 7.0 (Manzanita)
#
# Copyright (c) 2000 - 2017 The Regents of the University of California.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice unmodified and in its entirety, this list of conditions and the
# following disclaimer in the documentation and/or other materials provided
# with the distribution.
#
# 3. All advertising and press materials, printed or electronic, mentioning
# features or use of this software must display the following acknowledgement:
#
#       "This product includes software developed by the Rocks(r)
#       Cluster Group at the San Diego Supercomputer Center at the
#       University of California, San Diego and its contributors."
#
# 4. Except as permitted for the purposes of acknowledgment in paragraph 3,
# neither the name or logo of this software nor the names of its
# authors may be used to endorse or promote products derived from this
# software without specific prior written permission.  The name of the
# software includes the following terms, and any derivatives thereof:
# "Rocks", "Rocks Clusters", and "Avalanche Installer".  For licensing of
# the associated name, interested parties should contact Technology
# Transfer & Intellectual Property Services, University of California,
# San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
# Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu
#
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @Copyright@
#
# $Log$
#

.phony: $(CONDADOWNLOADS)

REDHAT.ROOT = $(CURDIR)/../../

-include $(ROCKSROOT)/etc/Rules.mk
include Rules.mk
include $(SDSCDEVEL)/Pull.mk
include $(SDSCDEVEL)/Rules.mk

NODE_FILE = $(REDHAT.ROOT)/nodes/cryoem-common.xml
TEMP_DIRS = $(REDHAT.ROOT)/src


$(NAME).spec: DESCRIPTION LICENSE

DESCRIPTION:
	$(MODULE_LOAD_CC); \
	( \
	  $(DESCRIBE_PKG); \
	  $(DESCRIBE_CC); \
	) > $@

LICENSE: $(SRC_DIRS)
	$(CHECK_LICENSE_FILES)
	echo "BSD/GNU" > $@

build: $(SRC_DIRS) $(CONDADOWNLOAD_PKGS)
	( \
          patch $(SOURCE_DIR)/examples/fftspeed.py patch-files/fftspeed.py.patch; \
          cd $(MINICONDA_DIR); \
          bash Miniconda2-latest-Linux-x86_64.sh -b -p $(PWD)/miniconda2; \
          export PYTHONPATH=$(CURDIR)/miniconda2; \
          export PATH=$(CURDIR)/miniconda2/bin:$${PATH}; \
          for package in $(CONDADOWNLOADS); do \
	    conda install --offline $(CURDIR)/$${package}.tar.bz2; \
          done; \
	  $(MODULE_LOAD_CC); \
	  $(MODULE_LOAD_CMAKE); \
          export LDFLAGS="-lGL -lGLU"; \
          export PATH=$(CURDIR)/miniconda2/bin:$${PATH}; \
          export PYTHON_ROOT=$${PYTHONROOT}; \
          export PYTHON_VERSION=`python -c "from __future__ import print_function;import sys; print(sys.version[:3])"`; \
          cd ../$(SOURCE_DIR); \
          mkdir build; \
          cd build; \
          cmake -DENABLE_OPTIMIZE_MACHINE=ON ..; \
          make -j 16; \
          make install; \
          sed -i '1s/.*/#!$(PKGROOT_ESCAPE)\/bin\/python/' $(PWD)/miniconda2/bin/conda; \
          packageCount=1; \
          mkdir  -p $(TEMP_DIRS)/build-eman2_$${packageCount}; \
          cd  $(CURDIR)/miniconda2; \
          rm -rf pkgs/*; \
          find . ! -type d -exec ls -l {} \; >$(TEMP_DIRS)/$$$$.tmp; \
          sum=0; \
          total=`du -b -s|awk '{print $$1}'`; \
          maxPackageSize=$$(( $$total / ( $(PACKAGES) ) )); \
          echo "MAX SIZE=$$maxPackageSize"; \
          exec <$(TEMP_DIRS)/$$$$.tmp; \
          cat /dev/null > $(TEMP_DIRS)/.files; \
          while read line; do \
               fileSize=`echo $$line |awk '{print $$5}'`; \
               fileName=`echo $$line |awk '{print $$9}'`; \
               sum=$$(( $$sum + $$fileSize )); \
               if test $$sum -ge $$maxPackageSize; then \
                  sum=$$(( $$sum - $$fileSize )); \
                  echo "NEW PACKAGE: OLD SUM $$sum"; \
                  tar cvf $(TEMP_DIRS)/package.tar `cat $(TEMP_DIRS)/.files`; \
                  cd $(TEMP_DIRS)/build-eman2_$$packageCount; \
                  tar xvf $(TEMP_DIRS)/package.tar; \
                  cd $(CURDIR)/miniconda2; \
                  sum=$$fileSize; \
                  packageCount=$$(( $$packageCount + 1 )); \
                  mkdir  -p $(TEMP_DIRS)/build-eman2_$$packageCount; \
                  cat /dev/null > $(TEMP_DIRS)/.files; \
               fi; \
               echo $$fileName >> $(TEMP_DIRS)/.files; \
          done; \
          echo $$packageCount > $(TEMP_DIRS)/packageCount; \
          tar cvf $(TEMP_DIRS)/package.tar `cat $(TEMP_DIRS)/.files`; \
          cd $(TEMP_DIRS)/build-eman2_$$packageCount; \
          tar xvf $(TEMP_DIRS)/package.tar; \
          cd $(CURDIR)/miniconda2; \
          rm -f $(TEMP_DIRS)/$$$$.tmp; \
          rm -f $(TEMP_DIRS)/.files; \
          rm -f $(TEMP_DIRS)/package.tar \
        )

$(CONDADOWNLOAD_PKGS):
	@echo "::: Downloading $(DL.SERVER)/$(DL.PATH)/$@ :::"
	@ if [ -f $@ ]; then \
		echo "::: $@ exists. Skipping... :::" ; \
	else \
		$(DL.CMD) $(DL.OPTS) $(DL.SERVER)/$(DL.PATH)/$@ ; \
	        echo "::: Verifying size of $@ :::" ; \
	        $(GREP.CMD) `$(STAT.CMD) --printf="%s" $@` $(VERIFY.HASHES) ; \
	        echo "::: Verifying hash of $@ :::" ; \
        	$(GREP.CMD) `$(VERIFY.CMD) $(VERIFY.OPTS) $@ ` $(VERIFY.HASHES) ; \
        	echo "" ; \
	fi

install::
	mkdir -p -m 755 $(ROOT)/$(PKGROOT)
	( \
          packageCount=`cat $(TEMP_DIRS)/packageCount`; \
          for PACKAGE in `seq 1 $${packageCount}`; do \
               perl -pi -e "print and s/$(NAME)_PACKAGE/$(NAME)_$${PACKAGE}/ if m/$(NAME)_PACKAGE/" $(NODE_FILE); \
          done; \
          perl -pi -e '$$_ = "" if m/$(NAME)_PACKAGE/' $(NODE_FILE) \
        )
	$(INSTALL_LICENSE_FILES)
	@echo "::: Package $(NAME) install complete :::"
	

clean::
	rm -rf $(SOURCE_DIR) $(MINICONDA_DIR) $(CONDADOWNLOAD_PKGS) DESCRIPTION LICENSE
