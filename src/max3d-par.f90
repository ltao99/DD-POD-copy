PROGRAM maxwell3D
   ! Mixed hybrid finite element domain decomposition method to solve 3D-magnetotellurics
   ! With relaxation: facbeta = 1.0 , factor = .5

   USE sizes
   USE globalvars
   USE commons_3d
   USE initialization
   USE mpi

   INTEGER :: status(MPI_STATUS_SIZE)
   REAL(DP) :: timer, start, finish
   INTEGER :: jk, I, k, l, j, cond1, cond2, cond3, init
   REAL(DP) :: start1, finish1, aux1, aux2, l2val, l2val_old
   REAL(DP) :: start_assembly, end_assembly, assembly_time
   REAL(DP) :: start_source, end_source, source_time
   REAL(DP) :: start_nodes, finish_nodes, nodes_time
   REAL(DP) :: start_basis, finish_basis, basis_time
   REAL(DP) :: start_init, finish_init, init_time
   REAL(DP) :: start_coef, finish_coef,coef_time
   REAL(DP) :: start_rhs, finish_rhs, rhs_time
   REAL(DP) :: start_xproj, finish_xproj, xproj_time
   REAL(DP) :: start_map, finish_map, map_time
   REAL(DP) :: start_calc, finish_calc, calc_time
   REAL(DP) :: start_sol, finish_sol, sol_time
   REAL(DP) :: start_comm, end_comm, comm_time
   REAL(DP) :: start_gather, end_gather, gather_time
   REAL(DP) :: start_read, end_read, read_time
   REAL(DP) :: start_output, end_output, output_time, total_time, steps_time
   CHARACTER(40) ::  prefix
   CHARACTER(200) :: logname
   COMPLEX(DP) :: ep, hp

   !---------------- MPI Initialization ----------------

  CALL MPI_INIT(mpierr)
  IF (mpierr .NE. 0) THEN
      WRITE(6,*) ' MPI error !!! (1)'
      STOP
  END IF

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
  IF (mpierr .NE. 0) THEN
      WRITE(6,*) ' MPI error !!! (2)'
      STOP
  END IF
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, procid, mpierr)
  IF (mpierr .NE. 0) THEN
      WRITE(6,*) ' MPI error !!! (2)'
      STOP
  END IF

  IF (MOD(nprocs,2) /= 0) THEN
     IF (procid == 0) WRITE(*,*) 'ERROR: world_size must be even!'
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, mpierr)
  END IF

  ! --- Map global rank to (sub_id, mode_id) ---
  mode_rank = MOD(procid, 2)    ! 0 or 1
  blockid  = procid / 2        ! 0..nsub-1

  ! --- Split communicators ---
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, blockid,  mode_rank, COMM_SUBD, mpierr)
  CALL MPI_COMM_SIZE(COMM_SUBD, mode_size, ierr)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, mode_rank, blockid,  COMM_MODE, mpierr)

   start = MPI_WTIME() 
 
   IF (mode_rank .EQ. 0) THEN
        mnodos = 1
        modo = "TE"
   ELSEIF (mode_rank .EQ. 1) THEN
        mnodos = 2
        modo = "TM"
   ENDIF

   !---------------- Read Mesh and Init ----------------
   start_read = MPI_WTIME()
   CALL readata                   ! Read mesh and frequencies
   CALL set_subdomain_coordinates ! Read mesh and frequencies

   CALL subtype                   ! Subdomain type
   CALL model                     ! Read Conductivity model
   end_read = MPI_WTIME()
   read_time = end_read - start_read

   !---------------- Get Connectivity Matrix and Interface NOD ----------------
   start_nodes = MPI_WTIME()
   CALL set_nodes_loc
   CALL set_nodes_interface
   finish_nodes = MPI_WTIME()
   nodes_time = finish_nodes - start_nodes


   !---------------- Frequency Loop ----------------
   DO jk = 1, numfreq
      start1 = MPI_WTIME()

      freq = frequency(jk)
      omega = twopi * freq
      coef_i = aim / (amu * omega)

      IF (blockid .EQ. 0 .AND. mode_rank .EQ. 0) THEN
         WRITE(6,*) 'frequency ', jk, ' from ', numfreq
         WRITE(6,104) freq
104      FORMAT(/1x,' frequency = ',g20.8,' Hz ')
         WRITE(6,108) omega
108      FORMAT(/1x,'omega = ',g20.8,' Hz ')
      END IF

      !---------------- POD ----------------
      ! Prepare POD bases for each subdomain and polarization mode (TE/TM).
      ! Computes and saves SVD bases, and registers them if POD is enabled.
      start_basis = MPI_WTIME()
      CALL prepare_pod_bases
      CALL mpi_Barrier(MPI_COMM_WORLD, mpierr)
      finish_basis = MPI_WTIME()
      basis_time = finish_basis - start_basis  

      !---------------Compute Parameters and Operators
      start_source = MPI_WTIME()
      CALL msource            !Compute Source from the background
      end_source = MPI_WTIME()
      source_time = end_source - start_source
      start_coef = MPI_WTIME()
      
      ! Parameters
      CALL inbeta
      CALL setalpha
      CALL setdelta
      end_coef = MPI_WTIME()
      coef_time = end_coef - start_coef

      !---------------- Assemble Stiffness Matrix (Full/Reduced) ----------------
      start_assembly = MPI_WTIME()
      ! Stiffness Matrices
      DO sub = 1, nsub_proc
         CALL mat_coef_global
         IF (use_pod) THEN
            ! Convert full stiffness matrix to CSR and compute/store reduced matrices K_rb
            ! in COO format for each subdomain, for each each MODE
            CALL assemble_Krb_from_coo
         ENDIF
         CALL mpi_Barrier(COMM_MODE, mpierr)
       end_assembly = MPI_WTIME()
         assembly_time = end_assembly - start_assembly
      END DO

      !------------- SOLVING ----------------
      !Set initial guess of the unknowns
      CALL inguess

      !Set Times
        rank = basis_map(1,mnodos)
        mumps_time = 0.0_dp
        solve_time = 0.0_dp
        factorization_time = 0.0_dp
        comm_time = 0.0_dp
        rhs_time = 0.0_dp
        gather_time = 0.0_dp
        xproj_time = 0.0_dp
        calc_time = 0.0_dp
        map_time = 0.0_dp
    
      ! Assemble Constant component of RHS
      start_rhs = MPI_WTIME()
      DO sub = 1,nsub_proc
         CALL rhs_glob_constant(rank)
      END DO
      finish_rhs = MPI_WTIME()
      rhs_time = rhs_time + (finish_rhs - start_rhs)

      !-----------Iterative Scheme-------------------
         DO i = 1, maxiter
            !Communication between subdomains
            start_comm = MPI_WTIME()
            IF (nprocs .GT. 1) CALL setbuff
            CALL mpi_Barrier(COMM_MODE, mpierr)
            end_comm = MPI_WTIME()
            comm_time = comm_time + (end_comm - start_comm)

            DO sub = 1, nsub_proc
               ! Set Size of Reduced LSE
               rank = basis_map(sub,mnodos)
               ! Compute RHS and Project RHS into reduced space if POD active
               start_rhs = MPI_WTIME()
               CALL rhs_glob_dd(rank)
               CALL prepare_rhs_local_space(rank)
               CALL mpi_Barrier(COMM_MODE, mpierr)
	            finish_rhs = MPI_WTIME()
	            rhs_time = rhs_time + (finish_rhs - start_rhs)
               ! Solve LSE
               CALL solver_MUMPS_multiple_solves(bglobal_rb(1:rank), sol_MUMPS_rb(1:rank), rank)
               CALL mpi_Barrier(COMM_MODE, mpierr)
              ! Project reduced solution into full space if POD active
               start_xproj = MPI_WTIME()
               CALL recover_full_solution(rank)
               CALL mpi_Barrier(COMM_MODE, mpierr)
               finish_xproj = MPI_WTIME()
	            xproj_time = xproj_time + (finish_xproj - start_xproj)
               ! Store solution as epsilon^n+1 lambda^n+1
               start_map = MPI_WTIME()
               CALL change_sol(sol_MUMPS(1:total_l_nod))
               CALL mpi_Barrier(COMM_MODE, mpierr)
               finish_map = MPI_WTIME()
            ENDDO

            ! Apply relaxation between epsilon^n+1 lambda^n+1 and epsilon^n lambda^n
            start_map = MPI_WTIME()
            CALL update_e
            CALL update_l
            CALL mpi_Barrier(COMM_MODE, mpierr)
            finish_map = MPI_WTIME()
	         map_time = map_time + (finish_map - start_map)
            start_calc = MPI_WTIME()
            !Compute relative L2 NORM between epsilon^n+1 and epsilon^n
            CALL calc_error(errn, den)
            CALL mpi_Barrier(COMM_MODE, mpierr)
            finish_calc = MPI_WTIME()
	         calc_time = calc_time + (finish_calc - start_calc)
            start_map = MPI_WTIME()
            ! Update epsilon^n+1 lambda^n+1 as epsilon^n lambda^n
            CALL update_eo_lo
            CALL mpi_Barrier(COMM_MODE, mpierr)
            finish_map = MPI_WTIME()
	         map_time = map_time + (finish_map - start_map)

            !Sum relative L2 norm from all subdomains 
            IF (i .GT. 1) THEN
               start_gather = MPI_WTIME()
               suml(1) = errn
               suml(2) = den
               CALL MPI_ALLREDUCE(suml, swork, 2, MPI_DOUBLE_PRECISION, MPI_SUM, COMM_MODE, ierr)
               errn = swork(1)
               den = swork(2)
               end_gather = MPI_WTIME()
               gather_time = gather_time + (end_gather - start_gather)
               IF (den .LT. 1.d-36) den = 1.d0
               errn = 2.d0 * dsqrt(errn / den)
               IF (blockid .EQ. ncontrol .AND. MOD(i,100) .EQ. 0) THEN
                  WRITE(6,200) i, errn
200               FORMAT(/1x,'  inner iteration  = ',i5,/,1x,'  relative error   = ',g20.8)
               END IF

               IF (errn .LT. convergence_tol) THEN
                  ! Check Convergence
                  IF (blockid .EQ. ncontrol) THEN
                     ! Finish Iterative Scheme
                     WRITE(6,66291) i, errn
66291                FORMAT(/1x,' convergence achieved in ',i5,' inner iterations ',/, &
                                 ' relative error value = ',e20.8/)
                  END IF
                  GOTO 1001
               END IF
            END IF
         END DO ! Iterations Mode 

1001     CONTINUE

         !Close MUMPS
         start_solve = MPI_WTIME()
         CALL finalize_solver_MUMPS
         CALL mpi_Barrier(COMM_MODE, mpierr)
         finish_solve = MPI_WTIME()
	 solve_time = solve_time + (finish_solve - start_solve)

         ! Calculate QoIs and output them
         start_output = MPI_WTIME()
         CALL output
         end_output = MPI_WTIME()
         output_time = end_output - start_output
         CALL MPI_Barrier(MPI_COMM_WORLD, mpierr)

         !Write Time log
         IF (blockid .EQ. 0) THEN
            ! Select filename prefix depending on use_pod
            IF (use_pod) THEN
               prefix = 'svdtimings'
            ELSE
               prefix = 'timings'
            END IF

            ! Construct full log filename
            WRITE(logname, '(A, A, A, I0, A, ES8.1, A, F5.3, A, F5.3, A, I0, A, ES8.1, A)') &
    	   'outputs/', TRIM(prefix), &
           '_nsnap', Nsnapshots, &
           '_svdtol', tol_svd, &
           '_relx', relx, &
           '_pnl', penal, &
           '_miter', maxiter, &
           '_ctol', convergence_tol, &
           '.txt'

            OPEN(UNIT=22+mnodos, FILE=TRIM(logname), STATUS='REPLACE')

            steps_time = read_time + nodes_time + basis_time + source_time &
            + assembly_time + coef_time + factorization_time + solve_time &
            + comm_time + rhs_time + xproj_time + map_time + calc_time + gather_time + output_time

            finish = MPI_WTIME()
            total_time = finish - start - basis_time

            WRITE(22+mnodos,'(A,F10.4)') 'Frequency (Hz):          ', freq
            WRITE(22+mnodos,'(A,I6)')    'Processes:               ', nprocs
            WRITE(22+mnodos,'(A,A)')     'Mode:                    ', TRIM(modo)
            WRITE(22+mnodos,'(A,3(I6,1X))') 'Mesh size (ngx, ngy, ngz): ', ngx, ngy, ngz
            WRITE(22+mnodos,'(A,2(I6,1X))') 'Subdomains (nsy, nsz):     ', nsy, nsz
            WRITE(22+mnodos,'(A,I4)')    'Iterations:                   ', i
            WRITE(22+mnodos,'(A,F10.4)') 'Read Model (ONCE) (s):        ', read_time
            WRITE(22+mnodos,'(A,F10.4)') 'Nodes Time (ONCE) (s):        ', nodes_time
            WRITE(22+mnodos,'(A,F10.4)') 'Basis Time (ONCE) (s):        ', basis_time
            WRITE(22+mnodos,'(A,F10.4)') 'Source Time (ONCE) (s):       ', source_time
            WRITE(22+mnodos,'(A,F10.4)') 'Assembly Time (ONCE) (s):     ', assembly_time
            WRITE(22+mnodos,'(A,F10.4)') 'Compute Coef Time (ONCE) (s): ', coef_time
            WRITE(22+mnodos,'(A,F10.4)') 'Analysis and Factor (s):      ', factorization_time
            WRITE(22+mnodos,'(A,F10.4)') 'Solve Time (s):               ', solve_time
            WRITE(22+mnodos,'(A,F10.4)') 'Communication Time (s):       ', comm_time
            WRITE(22+mnodos,'(A,F10.4)') 'RHS Time (s):                 ', rhs_time
            WRITE(22+mnodos,'(A,F10.4)') 'x Projection Time (s):        ', xproj_time
            WRITE(22+mnodos,'(A,F10.4)') 'Mapping Time (s):             ', map_time
            WRITE(22+mnodos,'(A,F10.4)') 'Calc Error Time (s):          ', calc_time
            WRITE(22+mnodos,'(A,F10.4)') 'Gather Time (s):              ', gather_time
            WRITE(22+mnodos,'(A,F10.4)') 'Output Time (s):              ', output_time
            WRITE(22+mnodos,'(A,F10.4)') 'Overhead Time (s):            ', total_time - steps_time + basis_time
            WRITE(22+mnodos,'(A,F10.4)') 'Total Time (s):               ', total_time
            CLOSE(22+mnodos)
         END IF
   END DO     ! Frequency loop


CALL MPI_Barrier(MPI_COMM_WORLD, mpierr)



! Finalize MPI
CALL MPI_FINALIZE(mpierr)
IF (mpierr .NE. 0) THEN
   WRITE(6,*) ' MPI error !!! (6)'
   STOP
END IF

END PROGRAM maxwell3D

 !---------------- Timer Function ----------------
 !REAL(DP) FUNCTION timer()
 !  USE sizes
 !  INTEGER :: mclock, time
 !  INTEGER :: time0 = 0
! 
!   time = mclock()
!   timer = (time - time0) / 100.d0
!   time0 = time
!   RETURN
! END FUNCTION timer



 