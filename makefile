




include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules



laplace: laplace.o  chkopts
		-${FLINKER} -o laplace laplace.o  ${PETSC_KSP_LIB}
			${RM} laplace.o




runlaplace:
		-@${MPIEXEC} -n 1 ./laplace -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > laplace_1.tmp 2>&1;	  \
				   if (${DIFF} output/laplace_1.out laplace_1.tmp) then true; \
				   	   else printf "${PWD}\nPossible problem with laplace_1, diffs above\n=========================================\n"; fi; \
					   	   ${RM} -f laplace_1.tmp



