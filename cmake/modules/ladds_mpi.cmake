option(LADDS_MPI "Activates distributed memory parallelization via MPI" ON)

# actual action
if (LADDS_MPI)
    # set flag that will be used to include the macro with the same name
    set(LADDS_INCLUDE_MPI true)
    # Activate MPI in AutoPas. No further linking required as this is handled by the AutoPas target.
    set(AUTOPAS_INCLUDE_MPI true)
endif ()
