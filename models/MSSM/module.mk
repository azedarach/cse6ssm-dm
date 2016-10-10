DIR          := models/MSSM
MODNAME      := MSSM
SARAH_MODEL  := MSSM

MSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSM_MK     := \
		$(DIR)/module.mk

MSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSM_BETAS_MK := \
		$(MSSM_SUSY_BETAS_MK) \
		$(MSSM_SOFT_BETAS_MK)

CMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CMSSM

lowMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowMSSM

CMSSM_GNUPLOT := \
		$(DIR)/CMSSM_plot_rgflow.gnuplot \
		$(DIR)/CMSSM_plot_spectrum.gnuplot

lowMSSM_GNUPLOT := \
		$(DIR)/lowMSSM_plot_rgflow.gnuplot \
		$(DIR)/lowMSSM_plot_spectrum.gnuplot

MSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSM_SRC := \
		$(DIR)/CMSSM_slha_io.cpp \
		$(DIR)/CMSSM_scan_utilities.cpp \
		$(DIR)/lowMSSM_slha_io.cpp \
		$(DIR)/MSSM_info.cpp \
		$(DIR)/MSSM_mass_eigenstates.cpp \
		$(DIR)/MSSM_physical.cpp \
		$(DIR)/MSSM_soft_parameters.cpp \
		$(DIR)/MSSM_susy_parameters.cpp \
		$(DIR)/MSSM_utilities.cpp

EXEMSSM_SRC :=

LIBMSSM_HDR := \
		$(DIR)/CMSSM_convergence_tester.hpp \
		$(DIR)/CMSSM_high_scale_constraint.hpp \
		$(DIR)/CMSSM_initial_guesser.hpp \
		$(DIR)/CMSSM_input_parameters.hpp \
		$(DIR)/CMSSM_low_scale_constraint.hpp \
		$(DIR)/CMSSM_model.hpp \
		$(DIR)/CMSSM_model_slha.hpp \
		$(DIR)/CMSSM_slha_io.hpp \
		$(DIR)/CMSSM_susy_scale_constraint.hpp \
		$(DIR)/CMSSM_scan_utilities.hpp \
		$(DIR)/lowMSSM_convergence_tester.hpp \
		$(DIR)/lowMSSM_high_scale_constraint.hpp \
		$(DIR)/lowMSSM_initial_guesser.hpp \
		$(DIR)/lowMSSM_input_parameters.hpp \
		$(DIR)/lowMSSM_low_scale_constraint.hpp \
		$(DIR)/lowMSSM_model.hpp \
		$(DIR)/lowMSSM_model_slha.hpp \
		$(DIR)/lowMSSM_slha_io.hpp \
		$(DIR)/lowMSSM_susy_scale_constraint.hpp \
		$(DIR)/MSSM_info.hpp \
		$(DIR)/MSSM_mass_eigenstates.hpp \
		$(DIR)/MSSM_physical.hpp \
		$(DIR)/MSSM_soft_parameters.hpp \
		$(DIR)/MSSM_susy_parameters.hpp \
		$(DIR)/MSSM_utilities.hpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMSSM_SRC += \
		$(DIR)/CMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/CMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/CMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/CMSSM_two_scale_input_parameters.cpp \
		$(DIR)/CMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/CMSSM_two_scale_model.cpp \
		$(DIR)/CMSSM_two_scale_model_slha.cpp \
		$(DIR)/CMSSM_two_scale_susy_scale_constraint.cpp \
		$(DIR)/lowMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/lowMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/lowMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/lowMSSM_two_scale_input_parameters.cpp \
		$(DIR)/lowMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/lowMSSM_two_scale_model.cpp \
		$(DIR)/lowMSSM_two_scale_model_slha.cpp \
		$(DIR)/lowMSSM_two_scale_susy_scale_constraint.cpp
EXEMSSM_SRC += \
		$(DIR)/run_CMSSM.cpp \
		$(DIR)/run_cmd_line_CMSSM.cpp \
		$(DIR)/gridscan_CMSSM.cpp \
		$(DIR)/scan_CMSSM.cpp \
		$(DIR)/run_lowMSSM.cpp \
		$(DIR)/run_cmd_line_lowMSSM.cpp \
		$(DIR)/scan_lowMSSM.cpp
LIBMSSM_HDR += \
		$(DIR)/CMSSM_spectrum_generator.hpp \
		$(DIR)/CMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/CMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/CMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/CMSSM_two_scale_input_parameters.hpp \
		$(DIR)/CMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/CMSSM_two_scale_model.hpp \
		$(DIR)/CMSSM_two_scale_model_slha.hpp \
		$(DIR)/CMSSM_two_scale_susy_scale_constraint.hpp \
		$(DIR)/lowMSSM_spectrum_generator.hpp \
		$(DIR)/lowMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/lowMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/lowMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/lowMSSM_two_scale_input_parameters.hpp \
		$(DIR)/lowMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/lowMSSM_two_scale_model.hpp \
		$(DIR)/lowMSSM_two_scale_model_slha.hpp \
		$(DIR)/lowMSSM_two_scale_susy_scale_constraint.hpp
endif

ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
LIBMSSM_SRC += \
		$(DIR)/CMSSM_semi_two_scale_constraint_handler.cpp \
		$(DIR)/CMSSM_semi_two_scale_convergence_tester.cpp \
		$(DIR)/CMSSM_semi_two_scale_high_scale_constraint.cpp \
		$(DIR)/CMSSM_semi_two_scale_initial_guesser.cpp \
		$(DIR)/CMSSM_semi_two_scale_input_parameters.cpp \
		$(DIR)/CMSSM_semi_two_scale_low_scale_constraint.cpp \
		$(DIR)/CMSSM_semi_two_scale_model.cpp \
		$(DIR)/CMSSM_semi_two_scale_model_slha.cpp \
		$(DIR)/CMSSM_semi_two_scale_susy_scale_constraint.cpp \
		$(DIR)/CMSSM_susy_two_scale_convergence_tester.cpp

EXEMSSM_SRC += \
		$(DIR)/run_database_scan.cpp \
		$(DIR)/run_semianalytic_CMSSM.cpp \
		$(DIR)/gridscan_semianalytic_CMSSM.cpp

ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
EXEMSSM_SRC += \
		$(DIR)/get_susyhd_higgs_mass.cpp
endif

LIBMSSM_HDR += \
		$(DIR)/CMSSM_semi_constraint_handler.hpp \
		$(DIR)/CMSSM_semi_convergence_tester.hpp \
		$(DIR)/CMSSM_semi_high_scale_constraint.hpp \
		$(DIR)/CMSSM_semi_initial_guesser.hpp \
		$(DIR)/CMSSM_semi_input_parameters.hpp \
		$(DIR)/CMSSM_semi_low_scale_constraint.hpp \
		$(DIR)/CMSSM_semi_model.hpp \
		$(DIR)/CMSSM_semi_model_slha.hpp \
		$(DIR)/CMSSM_semi_susy_scale_constraint.hpp \
		$(DIR)/CMSSM_semianalytic_spectrum_generator.hpp \
		$(DIR)/CMSSM_semi_two_scale_constraint_handler.hpp \
		$(DIR)/CMSSM_semi_two_scale_convergence_tester.hpp \
		$(DIR)/CMSSM_semi_two_scale_high_scale_constraint.hpp \
		$(DIR)/CMSSM_semi_two_scale_initial_guesser.hpp \
		$(DIR)/CMSSM_semi_two_scale_input_parameters.hpp \
		$(DIR)/CMSSM_semi_two_scale_low_scale_constraint.hpp \
		$(DIR)/CMSSM_semi_two_scale_model.hpp \
		$(DIR)/CMSSM_semi_two_scale_model_slha.hpp \
		$(DIR)/CMSSM_semi_two_scale_susy_scale_constraint.hpp \
		$(DIR)/CMSSM_susy_convergence_tester.hpp \
		$(DIR)/CMSSM_susy_two_scale_convergence_tester.hpp
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MSSM_SUSY_BETAS_MK)
-include $(MSSM_SOFT_BETAS_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all algorithms are used
LIBMSSM_SRC := $(sort $(LIBMSSM_SRC))
EXEMSSM_SRC := $(sort $(EXEMSSM_SRC))

LIBMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSM_SRC)))

EXEMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSM_SRC)))

LIBMSSM_DEP := \
		$(LIBMSSM_OBJ:.o=.d)

EXEMSSM_DEP := \
		$(EXEMSSM_OBJ:.o=.d)

LIBMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_CMSSM_OBJ := $(DIR)/run_CMSSM.o
RUN_CMSSM_EXE := $(DIR)/run_CMSSM.x

RUN_CMD_LINE_CMSSM_OBJ := $(DIR)/run_cmd_line_CMSSM.o
RUN_CMD_LINE_CMSSM_EXE := $(DIR)/run_cmd_line_CMSSM.x

GRIDSCAN_CMSSM_OBJ := $(DIR)/gridscan_CMSSM.o
GRIDSCAN_CMSSM_EXE := $(DIR)/gridscan_CMSSM.x

SCAN_CMSSM_OBJ := $(DIR)/scan_CMSSM.o
SCAN_CMSSM_EXE := $(DIR)/scan_CMSSM.x

RUN_DATABASE_OBJ := $(DIR)/run_database_scan.o
RUN_DATABASE_EXE := $(DIR)/run_database_scan.x

RUN_SEMI_CMSSM_OBJ := $(DIR)/run_semianalytic_CMSSM.o
RUN_SEMI_CMSSM_EXE := $(DIR)/run_semianalytic_CMSSM.x

GRIDSCAN_SEMI_CMSSM_OBJ := $(DIR)/gridscan_semianalytic_CMSSM.o
GRIDSCAN_SEMI_CMSSM_EXE := $(DIR)/gridscan_semianalytic_CMSSM.x

GET_HIGGS_MASS_OBJ := $(DIR)/get_susyhd_higgs_mass.o
GET_HIGGS_MASS_EXE := $(DIR)/get_susyhd_higgs_mass.x

RUN_lowMSSM_OBJ := $(DIR)/run_lowMSSM.o
RUN_lowMSSM_EXE := $(DIR)/run_lowMSSM.x

RUN_CMD_LINE_lowMSSM_OBJ := $(DIR)/run_cmd_line_lowMSSM.o
RUN_CMD_LINE_lowMSSM_EXE := $(DIR)/run_cmd_line_lowMSSM.x

SCAN_lowMSSM_OBJ := $(DIR)/scan_lowMSSM.o
SCAN_lowMSSM_EXE := $(DIR)/scan_lowMSSM.x

METACODE_STAMP_CMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSM_SRC) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSM_HDR) $(MSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSM_SRC) $(MSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSM_MK) $(MSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSM_BETAS_MK) $(MSSM_INSTALL_DIR)
ifneq ($(CMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CMSSM_SLHA_INPUT) $(MSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CMSSM_GNUPLOT) $(MSSM_INSTALL_DIR)
ifneq ($(lowMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(lowMSSM_SLHA_INPUT) $(MSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(lowMSSM_GNUPLOT) $(MSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSM_DEP)
		-rm -f $(EXEMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSM_OBJ)
		-rm -f $(EXEMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBMSSM)
		-rm -f $(RUN_CMSSM_EXE)
		-rm -f $(RUN_CMD_LINE_CMSSM_EXE)
		-rm -f $(GRIDSCAN_CMSSM_EXE)
		-rm -f $(SCAN_CMSSM_EXE)
		-rm -f $(RUN_DATABASE_EXE)
		-rm -f $(RUN_SEMI_CMSSM_EXE)
		-rm -f $(GRIDSCAN_SEMI_CMSSM_EXE)
		-rm -f $(GET_HIGGS_MASS_EXE)
		-rm -f $(RUN_lowMSSM_EXE)
		-rm -f $(RUN_CMD_LINE_lowMSSM_EXE)
		-rm -f $(SCAN_lowMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSM_TARBALL) \
		$(LIBMSSM_SRC) $(LIBMSSM_HDR) \
		$(EXEMSSM_SRC) \
		$(MSSM_MK) $(MSSM_BETAS_MK) \
		$(CMSSM_SLHA_INPUT) $(CMSSM_GNUPLOT) \
		$(lowMSSM_SLHA_INPUT) $(lowMSSM_GNUPLOT)

$(LIBMSSM_SRC) $(LIBMSSM_HDR) $(EXEMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CMSSM)"
		@echo "Note: to regenerate CMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CMSSM):
		@true
endif

$(LIBMSSM_DEP) $(EXEMSSM_DEP) $(LIBMSSM_OBJ) $(EXEMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(MLINKFLAGS) $(SQLITEFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSM_DEP) $(EXEMSSM_DEP) $(LIBMSSM_OBJ) $(EXEMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBMSSM): $(LIBMSSM_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBMSSM): $(LIBMSSM_OBJ)
		$(MAKELIB) $@ $^ $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)
endif

$(RUN_CMSSM_EXE): $(RUN_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_SEMI_CMSSM_EXE): $(RUN_SEMI_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_DATABASE_EXE): $(RUN_DATABASE_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(SQLITELIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CMD_LINE_CMSSM_EXE): $(RUN_CMD_LINE_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GRIDSCAN_CMSSM_EXE): $(GRIDSCAN_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_CMSSM_EXE): $(SCAN_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GRIDSCAN_SEMI_CMSSM_EXE): $(GRIDSCAN_SEMI_CMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GET_HIGGS_MASS_EXE): $(GET_HIGGS_MASS_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBMATHLINK) $(LIBSUSYHDLINK)
		$(CXX) -Wl,-no-as-needed -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(MLINKLIBS) $(EXTRA_MLINK_LIBS)

$(RUN_lowMSSM_EXE): $(RUN_lowMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CMD_LINE_lowMSSM_EXE): $(RUN_CMD_LINE_lowMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_lowMSSM_EXE): $(SCAN_lowMSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ALLDEP += $(LIBMSSM_DEP) $(EXEMSSM_DEP)
ALLSRC += $(LIBMSSM_SRC) $(EXEMSSM_SRC)
ALLLIB += $(LIBMSSM)
ifneq ($(findstring two_scale,$(ALGORITHMS)),)
ALLEXE += $(RUN_CMSSM_EXE) $(RUN_CMD_LINE_CMSSM_EXE) $(GRIDSCAN_CMSSM_EXE) $(SCAN_CMSSM_EXE) $(RUN_lowMSSM_EXE) $(RUN_CMD_LINE_lowMSSM_EXE) $(SCAN_lowMSSM_EXE)
endif
ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
ALLEXE += $(RUN_SEMI_CMSSM_EXE) $(GRIDSCAN_SEMI_CMSSM_EXE) $(RUN_DATABASE_EXE)
ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
ALLEXE += $(GET_HIGGS_MASS_EXE)
endif
endif
