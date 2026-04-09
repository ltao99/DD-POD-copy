subroutine solver_MUMPS_multiple_solves(b, sol, r)
    USE commons_3d
    USE sizes
    USE globalvars
    USE mumps_handler
    USE mpi
    implicit none
    integer, intent(in) :: r
    complex(DP), intent(in) :: b(1:r)   ! this rank's RHS
    complex(DP), intent(out) :: sol(1:r)
    integer :: nnz_local

    ! ---------------- Force per-process MUMPS context ----------------
    mumps_par%COMM = MPI_COMM_SELF
    if (use_pod) then
        mumps_par%SYM = 0
    else
        mumps_par%SYM = 1
    end if

    ! Each rank assembles its own local COO (no root/centralization)
    nnz_local = size(reduced_stiffnes_coo(sub, mnodos)%apj_rb)

    ! ---------- One-time init + analysis + factorization ----------
    if (.NOT. is_factorized) then
        start_aux = MPI_WTIME()
        mumps_par%PAR = 1
        mumps_par%JOB = -1
        CALL ZMUMPS(mumps_par)

        mumps_par%N = r
        mumps_par%NNZ = nnz_local
        if (nnz_local > 0) then
            ALLOCATE(mumps_par%IRN(nnz_local), mumps_par%JCN(nnz_local), mumps_par%A(nnz_local))
            mumps_par%IRN = reduced_stiffnes_coo(sub, mnodos)%ip_rb
            mumps_par%JCN = reduced_stiffnes_coo(sub, mnodos)%ij_rb
            mumps_par%A   = reduced_stiffnes_coo(sub, mnodos)%apj_rb
        end if

        ! Quiet
        mumps_par%ICNTL(1:4) = 0
        mumps_par%ICNTL(7)   = 0

        mumps_par%JOB = 1   ! analysis
        CALL ZMUMPS(mumps_par)
        mumps_par%JOB = 2   ! factorization
        CALL ZMUMPS(mumps_par)

        end_aux = MPI_WTIME()
        factorization_time = factorization_time + (end_aux - start_aux)
        is_factorized = .TRUE.
    end if

    ! ------------------------- Solve (single RHS) -------------------------
    start_aux = MPI_WTIME()
    if (ASSOCIATED(mumps_par%RHS)) DEALLOCATE(mumps_par%RHS)
    ALLOCATE(mumps_par%RHS(r))
    mumps_par%RHS = b

    mumps_par%JOB = 3
    CALL ZMUMPS(mumps_par)

    sol = mumps_par%RHS
    DEALLOCATE(mumps_par%RHS)

    end_aux = MPI_WTIME()
    solve_time = solve_time + (end_aux - start_aux)
end subroutine solver_MUMPS_multiple_solves


subroutine finalize_solver_MUMPS()
    USE commons_3d
    USE mumps_handler

    if (is_factorized) then
        mumps_par%JOB = -2
        CALL ZMUMPS(mumps_par)

        If (ASSOCIATED(mumps_par%IRN)) DEALLOCATE(mumps_par%IRN)
        If (ASSOCIATED(mumps_par%JCN)) DEALLOCATE(mumps_par%JCN)
        If (ASSOCIATED(mumps_par%A))   DEALLOCATE(mumps_par%A)

        is_factorized = .FALSE.
    end if
end subroutine finalize_solver_MUMPS







