DIR		:= work
MODNAME		:= work

LIBWORK_MK	:= \
		$(DIR)/module.mk

LIBWORK_AUTO	    := $(DIR)/autoprot.h

LIBWORK_CHEP        := $(DIR)/calchep

LIBWORK_PATH_TMPL   := $(DIR)/path.c.in
LIBWORK_PATH        := $(DIR)/path.c
WORK_PATH	    := $(ABSBASEDIR)/$(DIR)

LIBWORK_VERTEX      := $(DIR)/VandP.c

LIBWORK_MODELS_DIR  := $(DIR)/models
LIBWORK_PROCESS_DIR := $(DIR)/so_generated
LIBWORK_RESULTS_DIR := $(DIR)/results
LIBWORK_TEMP_DIR    := $(DIR)/tmp

LIBWORK_LOCKS	    := \
			$(DIR)/lock_ \
			$(DIR)/LOCK

LIBWORK_SRC := \
		$(LIBWORK_VERTEX) \
		$(LIBWORK_PATH)

LIBWORK_OBJ := \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBWORK_SRC)))

LIBWORK_DEP := \
		$(LIBWORK_OBJ:.o=.d)

LIBWORK	:= $(DIR)/work_aux.a

LIBWORK_MODEL_NUM := 1

LIBWORK_MODEL	:= \
		$(LIBWORK_MODELS_DIR)/extlib$(LIBWORK_MODEL_NUM).mdl \
		$(LIBWORK_MODELS_DIR)/func$(LIBWORK_MODEL_NUM).mdl \
		$(LIBWORK_MODELS_DIR)/lgrng$(LIBWORK_MODEL_NUM).mdl \
		$(LIBWORK_MODELS_DIR)/prtcls$(LIBWORK_MODEL_NUM).mdl \
		$(LIBWORK_MODELS_DIR)/vars$(LIBWORK_MODEL_NUM).mdl

.PHONY:	all-$(MODNAME) clean-$(MODNAME) mostly-clean-$(MODNAME) \
	distclean-$(MODNAME)

all-$(MODNAME): $(LIBWORK)

# BEGIN: NOT EXPORTED ##########################################
WORK_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

install-micro-model::
		install -d $(WORK_INSTALL_DIR)
		install -d $(WORK_INSTALL_DIR)/so_generated
		install -d $(WORK_INSTALL_DIR)/results
		install -d $(WORK_INSTALL_DIR)/tmp
		install -d $(WORK_INSTALL_DIR)/models
		$(INSTALL_STRIPPED) $(LIBWORK_MK) $(WORK_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(LIBWORK_PATH_TMPL) $(WORK_INSTALL_DIR)
		install -m u=rwx,g=r,o=r $(LIBWORK_CHEP) $(WORK_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBWORK_MODEL) $(WORK_INSTALL_DIR)/models
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME)-dep:
		-rm -f $(LIBWORK_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBWORK_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBWORK)
		-rm -f $(LIBWORK_PROCESS_DIR)/*
		-rm -rf $(LIBWORK_RESULTS_DIR)/*
		-rm -f $(LIBWORK_TEMP_DIR)/*
		-rm -f $(LIBWORK_LOCKS)
		-rm -rf $(WORK_PATH)/_*_*_

mostly-clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBWORK)
		-rm -rf $(LIBWORK_RESULTS_DIR)/*
		-rm -f $(LIBWORK_TEMP_DIR)/*
		-rm -f $(LIBWORK_LOCKS)
		-rm -rf $(WORK_PATH)/_*_*_ 

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(LIBWORK_AUTO)
		-rm -f $(LIBWORK_SRC)

clean::		clean-$(MODNAME)

mostly-clean::	mostly-clean-$(MODNAME)

distclean::	distclean-$(MODNAME)

$(LIBWORK_PATH):
		sed -e "s|@WORKDIR@|$(WORK_PATH)|" \
		< $(LIBWORK_PATH_TMPL) > $(LIBWORK_PATH)

$(LIBWORK_VERTEX): $(LIBWORK_MODEL)
		$(CALCHEP_DIR)/bin/make_VandP $(LIBWORK_MODELS_DIR) $(LIBWORK_MODEL_NUM)
		mv $(ABSBASEDIR)/VandP.c $(LIBWORK_VERTEX)
		mv $(ABSBASEDIR)/autoprot.h $(LIBWORK_AUTO)

$(LIBWORK): $(LIBWORK_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBWORK_DEP)
ALLSRC += $(LIBWORK_SRC)
ALLLIB += $(LIBWORK)
