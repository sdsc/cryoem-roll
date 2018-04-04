SRCDIRS = `find * -prune\
	  -type d 	\
	  ! -name CVS	\
            -not -name eman2-rpms \
            -not -name build-\* \
	  ! -name .` eman2-rpms
