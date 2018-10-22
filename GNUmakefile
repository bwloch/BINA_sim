# $Id: GNUmakefile,v 1.2 2000/10/19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := bina_breakup2.11_dTarget
G4TARGET := $(name)
G4EXLIB := true

#ifndef G4INSTALL
  G4INSTALL := /home/menorzinho/GEANT4/geant4-install/share/Geant4-10.4.2/geant4make
#endif
G4WORKDIR := /home/menorzinho/GEANT4/work_dir
G4SYSTEM := Linux-g++
G4INCLUDE := /home/menorzinho/GEANT4/geant4-install/include/Geant4
G4LIB := /home/menorzinho/GEANT4/geant4-install/lib
G4BASE := /home/menorzinho/GEANT4/geant4.10.04-build/source
#G4LIB_USE_GRANULAR := 0
.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

TEMP := ./bin

visclean:
	rm -f $(TEMP)/g4*.prim $(TEMP)/g4*.eps $(TEMP)/g4*.wrl
	rm -f $(TEMP)/.DAWN_*

