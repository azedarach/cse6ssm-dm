DIR          := models/CSE6SSM
MODNAME      := CSE6SSM
SARAH_MODEL  := SE6SSM

CSE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CSE6SSM_MK     := \
		$(DIR)/module.mk

CSE6SSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CSE6SSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CSE6SSM_BETAS_MK := \
		$(CSE6SSM_SUSY_BETAS_MK) \
		$(CSE6SSM_SOFT_BETAS_MK)

CSE6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CSE6SSM

CSE6SSM_GNUPLOT := \
		$(DIR)/CSE6SSM_plot_rgflow.gnuplot \
		$(DIR)/CSE6SSM_plot_spectrum.gnuplot

CSE6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCSE6SSM_SRC := \
		$(DIR)/CSE6SSM_higgs_upper_bound.cpp \
		$(DIR)/CSE6SSM_info.cpp \
		$(DIR)/CSE6SSM_mass_eigenstates.cpp \
		$(DIR)/CSE6SSM_physical.cpp \
		$(DIR)/CSE6SSM_slha_io.cpp \
		$(DIR)/CSE6SSM_soft_parameters.cpp \
		$(DIR)/CSE6SSM_susy_parameters.cpp \
		$(DIR)/CSE6SSM_utilities.cpp \
		$(DIR)/CSE6SSM_scan_utilities.cpp

EXECSE6SSM_SRC :=

LIBCSE6SSM_HDR := \
		$(DIR)/CSE6SSM_constraint_handler.hpp \
		$(DIR)/CSE6SSM_convergence_tester.hpp \
		$(DIR)/CSE6SSM_higgs_upper_bound.hpp \
		$(DIR)/CSE6SSM_high_scale_constraint.hpp \
		$(DIR)/CSE6SSM_info.hpp \
		$(DIR)/CSE6SSM_initial_guesser.hpp \
		$(DIR)/CSE6SSM_input_parameters.hpp \
		$(DIR)/CSE6SSM_low_scale_constraint.hpp \
		$(DIR)/CSE6SSM_mass_eigenstates.hpp \
		$(DIR)/CSE6SSM_model.hpp \
		$(DIR)/CSE6SSM_model_slha.hpp \
		$(DIR)/CSE6SSM_physical.hpp \
		$(DIR)/CSE6SSM_slha_io.hpp \
		$(DIR)/CSE6SSM_soft_parameters.hpp \
		$(DIR)/CSE6SSM_susy_parameters.hpp \
		$(DIR)/CSE6SSM_susy_scale_constraint.hpp \
		$(DIR)/CSE6SSM_utilities.hpp \
		$(DIR)/CSE6SSM_scan_utilities.hpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBCSE6SSM_SRC += \
		$(DIR)/CSE6SSM_scan_parameters.cpp \
		$(DIR)/CSE6SSM_two_scale_constraint_handler.cpp \
		$(DIR)/CSE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/CSE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/CSE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/CSE6SSM_two_scale_input_parameters.cpp \
		$(DIR)/CSE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/CSE6SSM_two_scale_model.cpp \
		$(DIR)/CSE6SSM_two_scale_model_slha.cpp \
		$(DIR)/CSE6SSM_two_scale_susy_scale_constraint.cpp
EXECSE6SSM_SRC += \
		$(DIR)/gridscan_CSE6SSM.cpp \
		$(DIR)/gridscan_rge_coeffs_CSE6SSM.cpp \
		$(DIR)/run_CSE6SSM.cpp \
		$(DIR)/run_cmd_line_CSE6SSM.cpp \
		$(DIR)/scan_CSE6SSM.cpp
LIBCSE6SSM_HDR += \
		$(DIR)/CSE6SSM_scan_parameters.hpp \
		$(DIR)/CSE6SSM_spectrum_generator.hpp \
		$(DIR)/CSE6SSM_two_scale_constraint_handler.hpp \
		$(DIR)/CSE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/CSE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/CSE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/CSE6SSM_two_scale_input_parameters.hpp \
		$(DIR)/CSE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/CSE6SSM_two_scale_model.hpp \
		$(DIR)/CSE6SSM_two_scale_model_slha.hpp \
		$(DIR)/CSE6SSM_two_scale_susy_scale_constraint.hpp
endif

ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
LIBCSE6SSM_SRC += \
		$(DIR)/CSE6SSM_semi_two_scale_constraint_handler.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_convergence_tester.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_high_scale_constraint.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_initial_guesser.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_input_parameters.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_low_scale_constraint.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_model.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_model_slha.cpp \
		$(DIR)/CSE6SSM_semi_two_scale_susy_scale_constraint.cpp \
		$(DIR)/CSE6SSM_susy_two_scale_convergence_tester.cpp

EXECSE6SSM_SRC += \
		$(DIR)/run_semianalytic_CSE6SSM.cpp \
		$(DIR)/gridscan_semianalytic_CSE6SSM.cpp \
		$(DIR)/scan_semianalytic_CSE6SSM.cpp \
		$(DIR)/generate_semianalytic_points.cpp

ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
EXECSE6SSM_SRC += \
		$(DIR)/run_susyhd_CSE6SSM.cpp \
		$(DIR)/get_susyhd_higgs_mass.cpp
endif

LIBCSE6SSM_HDR += \
		$(DIR)/CSE6SSM_semi_constraint_handler.hpp \
		$(DIR)/CSE6SSM_semi_convergence_tester.hpp \
		$(DIR)/CSE6SSM_semi_high_scale_constraint.hpp \
		$(DIR)/CSE6SSM_semi_initial_guesser.hpp \
		$(DIR)/CSE6SSM_semi_input_parameters.hpp \
		$(DIR)/CSE6SSM_semi_low_scale_constraint.hpp \
		$(DIR)/CSE6SSM_semi_model.hpp \
		$(DIR)/CSE6SSM_semi_model_slha.hpp \
		$(DIR)/CSE6SSM_semi_susy_scale_constraint.hpp \
		$(DIR)/CSE6SSM_semianalytic_spectrum_generator.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_constraint_handler.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_convergence_tester.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_high_scale_constraint.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_initial_guesser.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_input_parameters.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_low_scale_constraint.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_model.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_model_slha.hpp \
		$(DIR)/CSE6SSM_semi_two_scale_susy_scale_constraint.hpp \
		$(DIR)/CSE6SSM_susy_convergence_tester.hpp \
		$(DIR)/CSE6SSM_susy_two_scale_convergence_tester.hpp
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(CSE6SSM_SUSY_BETAS_MK)
-include $(CSE6SSM_SOFT_BETAS_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CSE6SSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CSE6SSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
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
LIBCSE6SSM_SRC := $(sort $(LIBCSE6SSM_SRC))
EXECSE6SSM_SRC := $(sort $(EXECSE6SSM_SRC))

LIBCSE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCSE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCSE6SSM_SRC)))

EXECSE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECSE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECSE6SSM_SRC)))

LIBCSE6SSM_DEP := \
		$(LIBCSE6SSM_OBJ:.o=.d)

EXECSE6SSM_DEP := \
		$(EXECSE6SSM_OBJ:.o=.d)

LIBCSE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

GRIDSCAN_CSE6SSM_OBJ := $(DIR)/gridscan_CSE6SSM.o
GRIDSCAN_CSE6SSM_EXE := $(DIR)/gridscan_CSE6SSM.x

GRIDSCAN_RGE_CSE6SSM_OBJ := $(DIR)/gridscan_rge_coeffs_CSE6SSM.o
GRIDSCAN_RGE_CSE6SSM_EXE := $(DIR)/gridscan_rge_coeffs_CSE6SSM.x

RUN_CSE6SSM_OBJ := $(DIR)/run_CSE6SSM.o
RUN_CSE6SSM_EXE := $(DIR)/run_CSE6SSM.x

RUN_CMD_LINE_CSE6SSM_OBJ := $(DIR)/run_cmd_line_CSE6SSM.o
RUN_CMD_LINE_CSE6SSM_EXE := $(DIR)/run_cmd_line_CSE6SSM.x

SCAN_CSE6SSM_OBJ := $(DIR)/scan_CSE6SSM.o
SCAN_CSE6SSM_EXE := $(DIR)/scan_CSE6SSM.x

RUN_SEMI_CSE6SSM_OBJ := $(DIR)/run_semianalytic_CSE6SSM.o
RUN_SEMI_CSE6SSM_EXE := $(DIR)/run_semianalytic_CSE6SSM.x

RUN_SUSYHD_CSE6SSM_OBJ := $(DIR)/run_susyhd_CSE6SSM.o
RUN_SUSYHD_CSE6SSM_EXE := $(DIR)/run_susyhd_CSE6SSM.x

GRIDSCAN_SEMI_CSE6SSM_OBJ := $(DIR)/gridscan_semianalytic_CSE6SSM.o
GRIDSCAN_SEMI_CSE6SSM_EXE := $(DIR)/gridscan_semianalytic_CSE6SSM.x

SCAN_SEMI_CSE6SSM_OBJ := $(DIR)/scan_semianalytic_CSE6SSM.o
SCAN_SEMI_CSE6SSM_EXE := $(DIR)/scan_semianalytic_CSE6SSM.x

GET_HIGGS_MASS_OBJ := $(DIR)/get_susyhd_higgs_mass.o
GET_HIGGS_MASS_EXE := $(DIR)/get_susyhd_higgs_mass.x

GENERATE_POINTS_OBJ := $(DIR)/generate_semianalytic_points.o
GENERATE_POINTS_EXE := $(DIR)/generate_semianalytic_points.x

METACODE_STAMP_CSE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CSE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCSE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CSE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCSE6SSM_SRC) $(CSE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCSE6SSM_HDR) $(CSE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECSE6SSM_SRC) $(CSE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CSE6SSM_MK) $(CSE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CSE6SSM_BETAS_MK) $(CSE6SSM_INSTALL_DIR)
ifneq ($(CSE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CSE6SSM_SLHA_INPUT) $(CSE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CSE6SSM_GNUPLOT) $(CSE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCSE6SSM_DEP)
		-rm -f $(EXECSE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCSE6SSM_OBJ)
		-rm -f $(EXECSE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCSE6SSM)
		-rm -f $(RUN_CSE6SSM_EXE)
		-rm -f $(RUN_CMD_LINE_CSE6SSM_EXE)
		-rm -f $(SCAN_CSE6SSM_EXE)
		-rm -f $(GRIDSCAN_CSE6SSM_EXE)
		-rm -f $(GRIDSCAN_RGE_CSE6SSM_EXE)
		-rm -f $(RUN_SEMI_CSE6SSM_EXE)
		-rm -f $(RUN_SUSYHD_CSE6SSM_EXE)
		-rm -f $(GRIDSCAN_SEMI_CSE6SSM_EXE)
		-rm -f $(SCAN_SEMI_CSE6SSM_EXE)
		-rm -f $(GET_HIGGS_MASS_EXE)
		-rm -f $(GENERATE_POINTS_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CSE6SSM_TARBALL) \
		$(LIBCSE6SSM_SRC) $(LIBCSE6SSM_HDR) \
		$(EXECSE6SSM_SRC) \
		$(CSE6SSM_MK) $(CSE6SSM_BETAS_MK) \
		$(CSE6SSM_SLHA_INPUT) $(CSE6SSM_GNUPLOT)

$(LIBCSE6SSM_SRC) $(LIBCSE6SSM_HDR) $(EXECSE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CSE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CSE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CSE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CSE6SSM)"
		@echo "Note: to regenerate CSE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CSE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CSE6SSM):
		@true
endif

$(LIBCSE6SSM_DEP) $(EXECSE6SSM_DEP) $(LIBCSE6SSM_OBJ) $(EXECSE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(MLINKFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCSE6SSM_DEP) $(EXECSE6SSM_DEP) $(LIBCSE6SSM_OBJ) $(EXECSE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBCSE6SSM): $(LIBCSE6SSM_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBCSE6SSM): $(LIBCSE6SSM_OBJ)
		$(MAKELIB) $@ $^ $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)
endif

$(GRIDSCAN_CSE6SSM_EXE): $(GRIDSCAN_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GRIDSCAN_RGE_CSE6SSM_EXE): $(GRIDSCAN_RGE_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CSE6SSM_EXE): $(RUN_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CMD_LINE_CSE6SSM_EXE): $(RUN_CMD_LINE_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_CSE6SSM_EXE): $(SCAN_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_SEMI_CSE6SSM_EXE): $(RUN_SEMI_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_SUSYHD_CSE6SSM_EXE): $(RUN_SUSYHD_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBMATHLINK) $(LIBSUSYHDLINK)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(MLINKLIBS) $(EXTRA_MLINK_LIBS)

$(GRIDSCAN_SEMI_CSE6SSM_EXE): $(GRIDSCAN_SEMI_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_SEMI_CSE6SSM_EXE): $(SCAN_SEMI_CSE6SSM_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GET_HIGGS_MASS_EXE): $(GET_HIGGS_MASS_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBMATHLINK) $(LIBSUSYHDLINK)
		$(CXX) -Wl,-no-as-needed -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(MLINKLIBS) $(EXTRA_MLINK_LIBS)

$(GENERATE_POINTS_EXE): $(GENERATE_POINTS_OBJ) $(LIBCSE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ALLDEP += $(LIBCSE6SSM_DEP) $(EXECSE6SSM_DEP)
ALLSRC += $(LIBCSE6SSM_SRC) $(EXECSE6SSM_SRC)
ALLLIB += $(LIBCSE6SSM)
ifneq ($(findstring two_scale,$(ALGORITHMS)),)
ALLEXE += $(GRIDSCAN_CSE6SSM_EXE) $(GRIDSCAN_RGE_CSE6SSM_EXE) $(RUN_CSE6SSM_EXE) $(RUN_CMD_LINE_CSE6SSM_EXE) $(SCAN_CSE6SSM_EXE)
endif
ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
ALLEXE += $(RUN_SEMI_CSE6SSM_EXE) $(GRIDSCAN_SEMI_CSE6SSM_EXE) $(SCAN_SEMI_CSE6SSM_EXE) $(GENERATE_POINTS_EXE)
ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
ALLEXE += $(RUN_SUSYHD_CSE6SSM_EXE) $(GET_HIGGS_MASS_EXE)
endif
endif
