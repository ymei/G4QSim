# --------------------------------------------------------------
# GNUmakefile
# --------------------------------------------------------------

name := G4QSim
G4TARGET := $(name)
G4EXLIB := true

G4DEBUG := 0

ifndef G4INSTALL
  G4INSTALL = /usr/share/Geant4-9.6.0/geant4make
endif

QTFLAGS   += -I/usr/include/qt4 -I/usr/include/qt4/Qt
QTLIBS    += -L/usr/lib/qt4

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS)
EXTRALIBS +=$(ROOTGLIBS)
EXTRALIBS +=-L/opt/lib

CPPFLAGS += $(ROOTCFLAGS)

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
