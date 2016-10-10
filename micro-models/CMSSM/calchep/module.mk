DIR		:= calchep
MODNAME		:= calchep

CALCHEP_MK	:= \
		$(DIR)/module.mk

CALCHEP_MODELS_DIR  := $(DIR)/models
CALCHEP_RESULTS_DIR := $(DIR)/results
CALCHEP_TEMP_DIR    := $(DIR)/tmp

CALCHEP_SH          := $(DIR)/calchep
CALCHEP_INI         := $(DIR)/calchep.ini

CALCHEP_MODEL_NUM   := 1

CALCHEP_MODEL	:= \
		$(CALCHEP_MODELS_DIR)/extlib$(CALCHEP_MODEL_NUM).mdl \
		$(CALCHEP_MODELS_DIR)/func$(CALCHEP_MODEL_NUM).mdl \
		$(CALCHEP_MODELS_DIR)/lgrng$(CALCHEP_MODEL_NUM).mdl \
		$(CALCHEP_MODELS_DIR)/prtcls$(CALCHEP_MODEL_NUM).mdl \
		$(CALCHEP_MODELS_DIR)/vars$(CALCHEP_MODEL_NUM).mdl

.PHONY: all-$(MODNAME) clean-$(MODNAME) mostly-clean-$(MODNAME) \
	distclean-$(MODNAME)

all-$(MODNAME):

# BEGIN: NOT EXPORTED ##########################################
CALCHEP_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

install-micro-model::
		install -d $(CALCHEP_INSTALL_DIR)
		install -d $(CALCHEP_INSTALL_DIR)/results
		install -d $(CALCHEP_INSTALL_DIR)/tmp
		install -d $(CALCHEP_INSTALL_DIR)/models
		$(INSTALL_STRIPPED) $(CALCHEP_MK) $(CALCHEP_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CALCHEP_MODEL) $(CALCHEP_INSTALL_DIR)/models
		install -m u=rwx,g=r,o=r $(CALCHEP_SH) $(CALCHEP_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CALCHEP_INI) $(CALCHEP_INSTALL_DIR)
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME):
		-rm -rf $(CALCHEP_TEMP_DIR)/*
		-rm -rf $(CALCHEP_RESULTS_DIR)/*

mostly-clean-$(MODNAME): clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)

clean::		clean-$(MODNAME)

mostly-clean::	mostly-clean-$(MODNAME)

distclean::	distclean-$(MODNAME)
