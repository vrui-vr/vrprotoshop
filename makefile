########################################################################
# Makefile for VR ProtoShop.
# Copyright (c) 2010-2025 Oliver Kreylos
#
# This file is part of the WhyTools Build Environment.
# 
# The WhyTools Build Environment is free software; you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
# 
# The WhyTools Build Environment is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the WhyTools Build Environment; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307 USA
########################################################################

# Directory containing the Vrui build system. The directory below
# matches the default Vrui installation; if Vrui's installation
# directory was changed during Vrui's installation, the directory below
# must be adapted.
VRUI_MAKEDIR := /usr/local/share/Vrui-12.3/make

# Base installation directory for VR ProtoShop. If this is set to the
# default of $(PWD), VR ProtoShop does not have to be installed to be
# run. Created executables, configuration files, and resources will be
# installed in the bin, etc, and share directories under the given base
# directory, respectively. Important note: Do not use ~ as an
# abbreviation for the user's home directory here; use $(HOME) instead.
INSTALLDIR := $(shell pwd)

########################################################################
# Everything below here should not have to be changed
########################################################################

# Version number for installation subdirectories. This is used to keep
# subsequent release versions of VR ProtoShop from clobbering each
# other. The value should be identical to the major.minor version
# number found in VERSION in the root package directory.
PACKAGE_VERSION = 4.4
PACKAGE_NAME = VRProtoShop-$(PACKAGE_VERSION)

# Set up the source directory structure
PACKAGEROOT := $(shell pwd)
CONFIGDIR = etc
RESOURCEDIR = share

# Include definitions for the system environment and system-provided
# packages
include $(VRUI_MAKEDIR)/SystemDefinitions
include $(VRUI_MAKEDIR)/Packages.System
include $(VRUI_MAKEDIR)/Configuration.Vrui
include $(VRUI_MAKEDIR)/Packages.Vrui

# Check if the Vrui Collaboration Infrastructure is installed
-include $(VRUI_MAKEDIR)/Configuration.Collaboration
-include $(VRUI_MAKEDIR)/Packages.Collaboration
ifdef COLLABORATION_VERSION
  HAVE_COLLABORATION = 1
else
  HAVE_COLLABORATION = 0
endif

# Set up installation directory structure:
EXECUTABLEINSTALLDIR = $(INSTALLDIR)/$(EXEDIR)
ifeq ($(INSTALLDIR),$(PACKAGEROOT))
  ETCINSTALLDIR = $(INSTALLDIR)/$(CONFIGDIR)
  SHAREINSTALLDIR = $(INSTALLDIR)/$(RESOURCEDIR)
else ifneq ($(findstring $(PACKAGE_NAME),$(INSTALLDIR)),)
  ETCINSTALLDIR = $(INSTALLDIR)/$(CONFIGDIR)
  SHAREINSTALLDIR = $(INSTALLDIR)/$(RESOURCEDIR)
else
  ETCINSTALLDIR = $(INSTALLDIR)/$(CONFIGDIR)/$(PACKAGE_NAME)
  SHAREINSTALLDIR = $(INSTALLDIR)/$(RESOURCEDIR)/$(PACKAGE_NAME)
endif

########################################################################
# Specify additional compiler and linker flags
########################################################################

# Add base directory to include path:
EXTRACINCLUDEFLAGS += -I$(PACKAGEROOT)

CFLAGS += -Wall -pedantic

########################################################################
# List common packages used by all components of this project
# (Supported packages can be found in $(VRUI_MAKEDIR)/Packages.*)
########################################################################

PACKAGES = MYGLGEOMETRY MYGLSUPPORT MYGLWRAPPERS MYGEOMETRY MYMATH MYCOMM MYTHREADS MYMISC GL

########################################################################
# Specify all final targets
########################################################################

EXECUTABLES = 
COLLABORATIONPLUGINS = 

EXECUTABLES += $(EXEDIR)/VRProtoShop

ifneq ($(HAVE_COLLABORATION),0)
  # Build the ProtoShop server-side collaboration plug-in
  PROTOSHOP_NAME = ProtoShop
  PROTOSHOP_VERSION = 1
  COLLABORATIONPLUGINS += $(call COLLABORATIONPLUGIN_SERVER_TARGET,PROTOSHOP)
endif

ALL = $(EXECUTABLES) $(COLLABORATIONPLUGINS)

.PHONY: all
all: $(ALL)

########################################################################
# Pseudo-target to print configuration options and configure the package
########################################################################

.PHONY: config config-invalidate
config: config-invalidate $(DEPDIR)/config

config-invalidate:
	@mkdir -p $(DEPDIR)
	@touch $(DEPDIR)/Configure-Begin

$(DEPDIR)/Configure-Begin:
	@mkdir -p $(DEPDIR)
	@echo "---- VR ProtoShop configuration options: ----"
	@echo "Collaborative visualization enabled"
	@touch $(DEPDIR)/Configure-Begin

$(DEPDIR)/Configure-Install: $(DEPDIR)/Configure-Begin
	@echo "---- VR ProtoShop installation configuration ----"
	@echo "Installation directory: $(INSTALLDIR)"
	@echo "Executable directory: $(EXECUTABLEINSTALLDIR)"
	@echo "Configuration directory: $(ETCINSTALLDIR)"
	@echo "Resource directory: $(SHAREINSTALLDIR)"
	@echo "Collaboration plug-in directory: $(COLLABORATIONPLUGINS_LIBDIR)"
	@touch $(DEPDIR)/Configure-Install

$(DEPDIR)/Configure-End: $(DEPDIR)/Configure-Install
	@echo "---- End of VR ProtoShop configuration options ----"
	@touch $(DEPDIR)/Configure-End

$(DEPDIR)/config: $(DEPDIR)/Configure-End
	@touch $(DEPDIR)/config

########################################################################
# Specify other actions to be performed on a `make clean'
########################################################################

.PHONY: extraclean
extraclean:

.PHONY: extrasqueakyclean
extrasqueakyclean:
	-rm -rf $(LIBEXT)

# Include basic makefile
include $(VRUI_MAKEDIR)/BasicMakefile

########################################################################
# Specify build rules for executables
########################################################################

# The VR ProtoShop application
VRPROTOSHOP_SOURCES = DragBox.cpp \
                      Atom.cpp \
                      Protein.cpp \
                      ParsePdbFile.cpp \
                      ProteinRenderer.cpp \
                      ReadStandards.cpp \
                      ReadPredictionFile.cpp \
                      SetDihedrals.cpp \
                      CreateProtein.cpp \
                      VolumeRenderer.cpp \
                      PaletteRenderer.cpp \
                      DenseMatrix.cpp \
                      ProtoShopProtocol.cpp \
                      ProtoShopClient.cpp \
                      IK.cpp \
                      UndoBuffer.cpp \
                      ProteinInteractor.cpp \
                      VRProtoShop.cpp

$(OBJDIR)/VRProtoShop.o: CFLAGS += -DVRPROTOSHOP_CONFIGFILENAME='"$(ETCINSTALLDIR)/ProtoShop.cfg"' -DVRPROTOSHOP_STANDARDSDIR='"$(SHAREINSTALLDIR)/Standards/"'

$(VRPROTOSHOP_SOURCES:%.cpp=$(OBJDIR)/%.o): | $(DEPDIR)/config

$(EXEDIR)/VRProtoShop: PACKAGES += MYCOLLABORATION2CLIENT MYVRUI MYGLMOTIF MYCLUSTER MYREALTIME
$(EXEDIR)/VRProtoShop: $(VRPROTOSHOP_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: VRProtoShop
VRProtoShop: $(EXEDIR)/VRProtoShop

# The collaborative ProtoShop server plugin
PROTOSHOPSERVER_SOURCES = Atom.cpp \
                          Protein.cpp \
                          ProtoShopProtocol.cpp \
                          ParsePdbFile.cpp \
                          ProtoShopServer.cpp

$(call PLUGINOBJNAMES,$(PROTOSHOPSERVER_SOURCES)): | $(DEPDIR)/config

$(call COLLABORATIONPLUGIN_SERVER_TARGET,PROTOSHOP): PACKAGES += MYGEOMETRY MYMATH
$(call COLLABORATIONPLUGIN_SERVER_TARGET,PROTOSHOP): $(call PLUGINOBJNAMES,$(PROTOSHOPSERVER_SOURCES))
.PHONY: ProtoShopServer
ProtoShopServer: $(call COLLABORATIONPLUGIN_SERVER_TARGET,PROTOSHOP)

########################################################################
# Specify installation rules for header files, libraries, executables,
# configuration files, and shared files.
########################################################################

install: $(ALL)
	@echo Installing VR ProtoShop in $(INSTALLDIR)...
	@install -d $(INSTALLDIR)
	@install -d $(EXECUTABLEINSTALLDIR)
	@install $(EXECUTABLES) $(EXECUTABLEINSTALLDIR)
	@install $(COLLABORATIONPLUGINS) $(COLLABORATIONPLUGINS_LIBDIR)
	@install -d $(ETCINSTALLDIR)
	@install etc/ProtoShop.cfg $(ETCINSTALLDIR)
	@install -d $(SHAREINSTALLDIR)/Standards
	@install share/VRProtoShop/Standards/* $(SHAREINSTALLDIR)/Standards
