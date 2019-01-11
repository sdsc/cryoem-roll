ifndef ROLLCOMPILER
  ROLLCOMPILER = gnu
endif
COMPILERNAME := $(firstword $(subst /, ,$(ROLLCOMPILER)))

ifndef ROLLMPI
  ROLLMPI = rocks-openmpi
endif
MPINAME := $(firstword $(subst /, ,$(ROLLMPI)))

CUDAVERSION=cuda
ifneq ("$(ROLLOPTS)", "$(subst cuda=,,$(ROLLOPTS))")
  CUDAVERSION = $(subst cuda=,,$(filter cuda=%,$(ROLLOPTS)))
endif


NAME           = sdsc-relion
VERSION        = 2.1
RELEASE        = 0
PKGROOT        = /opt/relion

SRC_SUBDIR     = relion

SOURCE_NAME    = relion
SOURCE_SUFFIX  = tar.gz
SOURCE_VERSION = $(VERSION)
SOURCE_PKG     = $(SOURCE_NAME)-$(SOURCE_VERSION).$(SOURCE_SUFFIX)
SOURCE_DIR     = $(SOURCE_PKG:%.$(SOURCE_SUFFIX)=%)

RELIONEXES     = relion relion_autopick relion_autopick_mpi relion_ctf_toolbox relion_display relion_find_tiltpairs relion_helix_toolbox relion_image_handler relion_localsym relion_localsym_mpi relion_maingui relion_manualpick relion_mask_create relion_particle_polish relion_particle_polish_mpi relion_particle_reposition relion_particle_sort relion_particle_sort_mpi relion_particle_symmetry_expand relion_pipeliner relion_postprocess relion_postprocess_mpi relion_prepare_subtomo relion_preprocess relion_preprocess_mpi relion_project relion_reconstruct relion_refine relion_refine_mpi relion_run_ctffind relion_run_ctffind_mpi relion_run_motioncorr relion_run_motioncorr_mpi relion_stack_create relion_star_combine relion_star_compare relion_tiltpair_plot

TAR_GZ_PKGS    = $(SOURCE_PKG)

RPM.EXTRAS     = AutoReq:No\nAutoProv:No
RPM.PREFIX     = $(PKGROOT)
