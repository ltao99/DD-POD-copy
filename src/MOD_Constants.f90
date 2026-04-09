!######################################################################
module globalvars
   USE sizes
   implicit none
   !Constants
   REAL(DP), parameter :: twopi= 6.28318530717958_dp ! 2pi
   COMPLEX(DP), parameter :: aim= (0.d0,1.d0)           ! Imag Unit
   REAL(DP), parameter  :: amu= 2.0*twopi*1e-7_dp     ! Magnetic Permeability

   !Numerical Parameters 
   REAL(DP), parameter :: ACCURMACHINE = 1.e-13      ! Compare with 0
   REAL(DP) :: f1= 781.d0/5040.d0         ! FEM Coefficient
   REAL(DP) :: f2= 1036.d0/784.d0         ! FEM Coefficient
   REAL(DP) :: f3= -59.d0/5040.d0         ! FEM Coefficient
   REAL(DP) :: f4= 269.d0/5040.d0         ! FEM Coefficient
   
   !POD
   CHARACTER(400) :: pod_path = '/scratch/dasilva/100x100x60_10Hz/'
   INTEGER :: Nsnapshots = 1024
   LOGICAL :: use_pod = .TRUE.
   REAL(DP) :: tol_svd = 1.0e-2

   !Simulations
   REAL(DP) :: relx = 0.05D0
   REAL(DP) :: penal = 0.05D0
   INTEGER :: maxiter = 999
   REAL :: convergence_tol = 5.0e-4

   end module globalvars

   
    