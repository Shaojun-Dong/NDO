

include make.inc

lib: 
	( cd $(TNSDIR); $(MAKE) )
	( cd $(EIG); $(MAKE) )
	( cd $(OPTI); $(MAKE); )
	(mkdir -p $(outputName);)
	(cp $(TNSDIR)/*.mod $(outputName);)
	(cp $(EIG)/*.mod $(outputName);)
	(cp $(OPTI)/*.mod $(outputName);)
	( ar rc $(outputName)/$(outputName).a $(TNSDIR)/*.o  $(OPTI)/*.o $(EIG)/*.o)
	make echo

old: 
	( cd $(oldTNSDIR); $(MAKE) )
	( cd $(oldSYMM); $(MAKE) )
	( cd $(oldEIG); $(MAKE) )
	( cd $(oldOPTI); $(MAKE); )
	(mkdir -p $(oldoutputName);)
	(cp $(oldTNSDIR)/*.mod $(oldoutputName);)
	(cp $(oldSYMM)/*.mod $(oldoutputName);)
	(cp $(oldEIG)/*.mod $(oldoutputName);)
	(cp $(oldOPTI)/*.mod $(oldoutputName);)
	( ar rc $(oldoutputName)/$(oldoutputName).a $(oldTNSDIR)/*.o $(oldSYMM)/*.o $(oldEIG)/*.o $(oldOPTI)/*.o )
	make echoold

clean:
	( find ./src/ -name "*.mod" -exec rm  {} \; )
	( find ./src/ -name "*.o" -exec rm  {} \;)
	

cleanAll:
	( find . -name "*.mod" -exec rm  {} \; )
	( find . -name "*.o" -exec rm  {} \;)
	(rm  -rf libTensor-*)

cleanlib:
	(rm  -rf $(outputName))

cleanAlllib:
	(rm  -rf libTensor-*)

cleanAllMod:
	( rm */*/*.mod )
	( rm */*/*.o )
echoold:
	@$(ECHO) " "
	@$(ECHO) " "
	@$(ECHO) "Write out the info of the package to $(outputName)/package.info" 
	@$(ECHO) " "
	@$(ECHO) " ========== info of the package  ================"
	@$(ECHO) "   Package Name              : $(outputName)"
	@$(ECHO) "   Tensor version            : $(version)"
	@$(ECHO) "   Symmetry Tensor version   : $(Symmetry_version)"
	@$(ECHO) "   Eigen value version       : $(eigvalue_version)"
	@$(ECHO) "   optimization Tool version : $(optimization_version)"
	@$(ECHO) "   Compile optionl           : $(openmp) $(Optionl)"
	@$(ECHO) " "
	@$(ECHO) " $(FC) version:"
	@$(ECHO) " "
	@$(FC) --version 



echo:
	@$(ECHO) " "
	@$(ECHO) " "
	@$(ECHO) "Write out the info of the package to $(outputName)/package.info" 
	@$(ECHO) " "
	@$(ECHO) " ========== info of the package  ================"
	@$(ECHO) "   Package Name              : $(outputName)"
	@$(ECHO) "   Tensor version            : $(version)"
	@$(ECHO) "   optimization Tool version : $(optimization_version)"
	@$(ECHO) "   Compile optionl           : $(openmp) $(Optionl)"
	@$(ECHO) " "
	@$(ECHO) " $(FC) version:"
	@$(ECHO) " "
	@$(FC) --version 
	@$(ECHO) " ">$(outputName)/package.info
	@$(ECHO) " ========== info of the package  ================">>$(outputName)/package.info
	@$(ECHO) "   Package Name              : $(outputName)">>$(outputName)/package.info
	@$(ECHO) "   Tensor version            : $(version)">>$(outputName)/package.info
	@$(ECHO) "   optimization Tool version : $(optimization_version)">>$(outputName)/package.info
	@$(ECHO) "   Compile optionl           : $(openmp) $(Optionl)">>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
	@$(ECHO) " $(FC) version:">>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
	@$(FC) --version >>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
	@$(ECHO) " "
	@$(ECHO) " Report any bugs to Shaojun Dong via ">>$(outputName)/package.info
	@$(ECHO) "    sj.dong@outlook.com ">>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
	@$(ECHO) " ">>$(outputName)/package.info
