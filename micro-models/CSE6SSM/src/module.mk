DIR              := src
MODNAME	         := CSE6SSM

LIBCSE6SSM_MK	 := \
		$(DIR)/module.mk

LIBCSE6SSM_SRC := \
		$(DIR)/CSE6SSM_utilities.c \
		$(DIR)/read_slha.c

LIBCSE6SSM_HDR := \
		$(DIR)/CSE6SSM_utilities.h \
		$(DIR)/read_slha.h

EXECSE6SSM_SRC := \
		$(DIR)/run_slha_file.c

LIBCSE6SSM_OBJ := \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBCSE6SSM_SRC)))

EXECSE6SSM_OBJ := \
		$(patsubst %.c, %.o, $(filter %.c, $(EXECSE6SSM_SRC)))

EXECSE6SSM_EXE := \
		$(patsubst %.c, %.x, $(filter %.c, $(EXECSE6SSM_SRC)))

LIBCSE6SSM_DEP := \
		$(LIBCSE6SSM_OBJ:.o=.d)

EXECSE6SSM_DEP := \
		$(EXECSE6SSM_OBJ:.o=.d)

LIBCSE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:		all-$(MODNAME) clean-$(MODNAME) mostly-clean-$(MODNAME) \
		distclean-$(MODNAME)

all-$(MODNAME):	$(LIBCSE6SSM)

# BEGIN: NOT EXPORTED ##########################################
SRC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

install-micro-model::
		install -d $(SRC_INSTALL_DIR)

		$(INSTALL_STRIPPED) $(LIBCSE6SSM_MK) $(SRC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(LIBCSE6SSM_HDR) $(SRC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCSE6SSM_SRC) $(EXECSE6SSM_SRC) $(SRC_INSTALL_DIR)
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME)-dep:
		-rm -f $(LIBCSE6SSM_DEP)
		-rm -f $(EXECSE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCSE6SSM_OBJ)
		-rm -f $(EXECSE6SSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCSE6SSM)
		-rm -f $(EXECSE6SSM_EXE)

mostly-clean-$(MODNAME): clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)

clean::		clean-$(MODNAME)

mostly-clean::	mostly-clean-$(MODNAME)

distclean::	distclean-$(MODNAME)

$(LIBCSE6SSM_DEP) $(EXECSE6SSM_DEP) $(LIBCSE6SSM_OBJ) $(EXECSE6SSM_OBJ): CPPFLAGS += $(MICROMEGASFLAGS) $(CALCHEPFLAGS)

$(LIBCSE6SSM): $(LIBCSE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCSE6SSM) $(LIBWORK)
		$(CC) -o $@ $(call abspathx,$<) $(call abspathx,$(LIBCSE6SSM)) $(MICROMEGASLIBS) $(call abspathx,$(LIBWORK)) $(CALCHEPMODELLIBS) $(CALCHEPBASELIBS) $(X11LIBS) $(FLIBS) $(DYNAMICLIBS) $(MATHLIBS)

ALLDEP += $(LIBCSE6SSM_DEP) $(EXECSE6SSM_DEP)
ALLSRC += $(LIBCSE6SSM_SRC) $(EXECSE6SSM_SRC)
ALLLIB += $(LIBCSE6SSM)
ALLEXE += $(EXECSE6SSM_EXE)
