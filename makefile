




include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

laplace_test: global_mod.o laplace_test.o chkopts
	        -${FLINKER} -o laplace_test laplace_test.o global_mod.o ${PETSC_DM_LIB}
		        ${RM} -f laplace_test.o global_mod.o global_mod global_mod.mod


laplace: laplace.o  chkopts
		-${FLINKER} -o laplace laplace.o  ${PETSC_KSP_LIB}
			${RM} laplace.o

runlaplace_test:
	        -@${MPIEXEC} -n 4 ./laplace_test  > laplace_test.tmp 2>&1;   \
			           if (${DIFF} output/laplace_test.out laplace_test.tmp) then true; \
				              else printf "${PWD}\nPossible problem with laplace_test, diffs above\n=========================================\n"; fi; \
					                 ${RM} -f laplace_test.tmp


runlaplace:
		-@${MPIEXEC} -n 1 ./laplace -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > laplace_1.tmp 2>&1;	  \
				   if (${DIFF} output/laplace_1.out laplace_1.tmp) then true; \
				   	   else printf "${PWD}\nPossible problem with laplace_1, diffs above\n=========================================\n"; fi; \
					   	   ${RM} -f laplace_1.tmp



