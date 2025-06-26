
###############################################################################
################### MECH Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project 
# MODULE_DIR       - Location of the MOOSE modules directory
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
###############################################################################
MECH_DIR           ?= $(shell pwd)
MOOSE_DIR          ?= $(shell dirname $(MECH_DIR))/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
SRC_DIR	           ?= $(MECH_DIR)/src
###############################################################################
# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
# dep apps
APPLICATION_DIR    := $(MECH_DIR)
APPLICATION_NAME   := mech
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here

#clean::
#	@./auxiliaryCleaner.sh;

