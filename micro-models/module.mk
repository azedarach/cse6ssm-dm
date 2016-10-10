DIR     := micro-models
MODNAME := micro

MICROMEGAS_MODEL_DIRS := \
		$(DIR)/CMSSM \
		$(DIR)/CSE6SSM

INSTALL_TARGETS := \
		install-micro-CSE6SSM \
		install-micro-CMSSM

.PHONY:		clean-$(MODNAME) distclean-$(MODNAME)

$(INSTALL_TARGETS): install-micro-%:
		cd $(DIR)/$* && make install-micro-model

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		-@for d in $(MICROMEGAS_MODEL_DIRS); do \
			(cd $$d && make distclean); \
		done

clean::		clean-$(MODNAME)

distclean::	distclean-$(MODNAME)
