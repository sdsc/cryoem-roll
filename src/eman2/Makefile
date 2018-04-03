#
# $Id$
#
# @Copyright@
# 
# 				Rocks(r)
# 		         www.rocksclusters.org
# 		       version 6.1.1 (Sand Boa)
# 
# Copyright (c) 2000 - 2014 The Regents of the University of California.
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
# 	"This product includes software developed by the Rocks(r)
# 	Cluster Group at the San Diego Supercomputer Center at the
# 	University of California, San Diego and its contributors."
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

REDHAT.ROOT = $(CURDIR)/../../

-include $(ROCKSROOT)/etc/Rules.mk
include Rules.mk
include $(SDSCDEVEL)/Pull.mk
include $(SDSCDEVEL)/Rules.mk

DESCRIBE_GSL = echo built with gsl $(call GET_MODULE_VERSION, gsl)
MODULE_LOAD_GSL = $(call MODULE_LOAD_PACKAGE, gsl, GSLHOME)

CACHE_DIR=$(REDHAT.ROOT)/cache

CUDABUILD = no
CUDAOPTS=
ifneq ("$(ROLLOPTS)", "$(subst cuda,,$(ROLLOPTS))")
  CUDABUILD = yes
  CUDAOPTS=-DENABLE_SPARX_CUDA=yes -DENABLE_EMAN_CUDA=yes
endif

BOOST_TOOLSET = gcc
CFLAGS = 
ifeq ("$(COMPILERNAME)", "intel")
  BOOST_TOOLSET = intel-linux
else ifeq ("$(COMPILERNAME)", "pgi")
  BOOST_TOOLSET = pgi
endif

$(NAME).spec: DESCRIPTION LICENSE

DESCRIPTION:
	$(MODULE_LOAD_CC); \
	( \
	  $(DESCRIBE_PKG); \
	  $(DESCRIBE_CC); \
	  $(DESCRIBE_MPI); \
	  $(DESCRIBE_MKL); \
	  $(DESCRIBE_FFTW); \
	  $(DESCRIBE_HDF5); \
	  $(DESCRIBE_SCIPY); \
	  $(DESCRIBE_GSL) \
	) > $@
	$(MODULE_LOAD_PYTHON); \
	$(DESCRIBE_PYTHON) >> $@; \
	if test "$(CUDABUILD)" = "yes"; then \
	  $(DESCRIBE_CUDA) >> $@; \
	fi

LICENSE: $(SRC_DIRS)
	$(CHECK_LICENSE_FILES)
	echo "BSD/GNU" > $@

build: $(SRC_DIRS)
	patch $(SOURCE_DIR)/examples/fftspeed.py patch-files/fftspeed.py.patch
	( \
	  $(MODULE_LOAD_CC); \
	  $(MODULE_LOAD_MPI); \
	  $(MODULE_LOAD_CMAKE); \
	  $(MODULE_LOAD_MKL); \
	  $(MODULE_LOAD_FFTW); \
	  $(MODULE_LOAD_HDF5); \
	  $(MODULE_LOAD_SCIPY); \
	  $(MODULE_LOAD_GSL); \
	  $(MODULE_LOAD_PYTHON); \
	  if test "$(CUDABUILD)" = "yes"; then \
	    $(MODULE_LOAD_CUDA); \
	  fi; \
	  pyvers=`python --version 2>&1 | grep -o '[0-9][0-9]*\.[0-9][0-9]*'`; \
          export LDFLAGS="-lGL -lGLU"; \
          cd $(FTGL_DIR); \
          ./configure CC=$(CC) CXX=$(CXX) --prefix=$(CACHE_DIR)/ftgl; \
          make; \
          make install; \
	  cd ../$(BOOST_DIR); \
	  cp tools/build/v2/user-config.jam .; \
	  ./bootstrap.sh --prefix=$(CACHE_DIR)/boost --with-toolset=$(BOOST_TOOLSET); \
	  ./b2 cxxflags="-fPIC" cflags="-fPIC" -a install; \
          export FFTWDIR=$${FFTWHOME}; \
          export GSLDIR=$${GSLHOME}; \
          export HDF5DIR=$${HDF5HOME}; \
          export PYTHON_ROOT=$${PYTHONROOT}; \
          export PYTHON_VERSION=`python -c "from __future__ import print_function;import sys; print(sys.version[:3])"`; \
          cd ..; \
          mkdir $(SOURCE_DIR)/build; \
          cd $(SOURCE_DIR)/build; \
	  numpyIncludePath=`find $${PYTHONROOT}/lib/python$${pyvers} -name arrayobject.h | sed 's/numpy.arrayobject.h//'`; \
          cmake -DEMAN_INSTALL_PREFIX=$(ROOT)/$(PKGROOT) \
           -DCMAKE_C_COMPILER=mpicc \
           -DCMAKE_CXX_COMPILER=mpicxx \
           -DFTGL_INCLUDE_PATH=$(CACHE_DIR)/ftgl/include \
           -DFTGL_LIBRARY=$(CACHE_DIR)/ftgl/lib \
           -DBOOST_INCLUDE_PATH=$(CACHE_DIR)/boost/include \
           -DFREETYPE_INCLUDE_PATH=/usr/include/freetype2 \
           -DBOOST_LIBRARY=$(CACHE_DIR)/boost/lib $(CUDAOPTS) \
           -DNUMPY_INCLUDE_PATH=$$numpyIncludePath \
           -DPYTHON_LIBRARY=$${PYTHONROOT}/lib  .. \
	)
	@echo "::: Package $(NAME) build complete :::"
	
install::
	mkdir -p $(ROOT)/$(PKGROOT)
	( \
	  cd $(SOURCE_DIR)/build; \
	  $(MODULE_LOAD_CC); \
          $(MODULE_LOAD_MPI); \
          $(MODULE_LOAD_CMAKE); \
          $(MODULE_LOAD_MKL); \
          $(MODULE_LOAD_HDF5); \
          if test "$(CUDABUILD)" = "yes"; then \
             $(MODULE_LOAD_CUDA); \
          fi; \
	  $(MODULE_LOAD_SCIPY); \
	  $(MODULE_LOAD_GSL); \
          make install; \
	  cp -r $(CACHE_DIR)/boost/lib/* $(ROOT)/$(PKGROOT)/lib/ \
	)
	$(INSTALL_LICENSE_FILES)
	@echo "::: Package $(NAME) install complete :::"

clean::
	rm -rf $(SOURCE_DIR) DESCRIPTION LICENSE
