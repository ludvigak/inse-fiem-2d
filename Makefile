all: fmm

fmm:
	cd external/mbh2dfmm; ./configure_makefile.sh; make lib
