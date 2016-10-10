DIR              := src
MODNAME	         := CMSSM

LIBCMSSM_MK	 := \
		$(DIR)/module.mk

LIBCMSSM_SRC := \
		$(DIR)/CMSSM_utilities.c \
		$(DIR)/read_slha.c

LIBCMSSM_HDR := \
		$(DIR)/CMSSM_utilities.h \
		$(DIR)/read_slha.h

EXECMSSM_SRC := \
		$(DIR)/run_slha_file.c

LIBCMSSM_OBJ := \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBCMSSM_SRC)))

EXECMSSM_OBJ := \
		$(patsubst %.c, %.o, $(filter %.c, $(EXECMSSM_SRC)))

EXECMSSM_EXE := \
		$(patsubst %.c, %.x, $(filter %.c, $(EXECMSSM_SRC)))

LIBCMSSM_DEP := \
		$(LIBCMSSM_OBJ:.o=.d)

EXECMSSM_DEP := \
		$(EXECMSSM_OBJ:.o=.d)

LIBCMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:		all-$(MODNAME) clean-$(MODNAME) mostly-clean-$(MODNAME) \
		distclean-$(MODNAME)

all-$(MODNAME):	$(LIBCMSSM)

# BEGIN: NOT EXPORTED ##########################################
SRC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

install-micro-model::
		install -d $(SRC_INSTALL_DIR)

		$(INSTALL_STRIPPED) $(LIBCMSSM_MK) $(SRC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(LIBCMSSM_HDR) $(SRC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCMSSM_SRC) $(EXECMSSM_SRC) $(SRC_INSTALL_DIR)
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME)-dep:
		-rm -f $(LIBCMSSM_DEP)
		-rm -f $(EXECMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCMSSM_OBJ)
		-rm -f $(EXECMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCMSSM)
		-rm -f $(EXECMSSM_EXE)

mostly-clean-$(MODNAME): clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)

clean::		clean-$(MODNAME)

mostly-clean::	mostly-clean-$(MODNAME)

distclean::	distclean-$(MODNAME)

$(LIBCMSSM_DEP) $(EXECMSSM_DEP) $(LIBCMSSM_OBJ) $(EXECMSSM_OBJ): CPPFLAGS += $(MICROMEGASFLAGS) $(CALCHEPFLAGS)

$(LIBCMSSM): $(LIBCMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCMSSM) $(LIBWORK)
		$(CC) -o $@ $(call abspathx,$<) $(call abspathx,$(LIBCMSSM)) $(MICROMEGASLIBS) $(call abspathx,$(LIBWORK)) $(CALCHEPMODELLIBS) $(CALCHEPBASELIBS) $(X11LIBS) $(FLIBS) $(DYNAMICLIBS) $(MATHLIBS)

ALLDEP += $(LIBCMSSM_DEP) $(EXECMSSM_DEP)
ALLSRC += $(LIBCMSSM_SRC) $(EXECMSSM_SRC)
ALLLIB += $(LIBCMSSM)
ALLEXE += $(EXECMSSM_EXE)
