###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# HERD_TRUNK_DIR   - Location of the HERD repository
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
WORKSPACE_DIR         ?= $(shell dirname `pwd`)
MOOSE_DIR             ?= $(shell dirname $(MODULE_DIR))/../moose
FRAMEWORK_DIR         ?= $(MOOSE_DIR)/framework
###############################################################################
# Additional special case targets should be added here


include $(UTOPIA_FE_DIR)/config/utopia_fe_config.makefile

ADDITIONAL_INCLUDES+=$(UTOPIA_FE_INCLUDES)

EXTERNAL_FLAGS=$(UTOPIA_FE_LIBRARIES)
