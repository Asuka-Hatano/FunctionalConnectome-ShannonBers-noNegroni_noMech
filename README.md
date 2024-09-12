FunctionalConnectome_recalcForTau.m --- main program, with comments. need two input files: Metrics***.csv and pyfinal***.csv  Just run.

parfor_calling.m ---------------------- parallelized version of the above code. 
FunctionalConnectome_recalcForTau_forParallel.m

ChangeFilteringRule_ReFilterCSV.m ----- Metrics thresholds have arbitrality. Read output afterwards and refilter with different threshods.

TauCalcFunc.m ------------------------- In main routine, code for tau calculation is long... I made subroutine for that, but I haven't merge this with main code.

