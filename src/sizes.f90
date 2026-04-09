MODULE sizes

  USE iso_fortran_env, ONLY: DP => REAL64
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: nsy= 1, nsz= 9

  INTEGER, PARAMETER :: nprocs_input = 18
  
  INTEGER, PARAMETER :: nsub = nsy * nsz

  INTEGER, PARAMETER :: nsub_proc = 2*nsub/nprocs_input

  INTEGER, PARAMETER :: ngx= 40, ngy= 40, ngz=36, mnmat= ngx * ngy * ngz, iter= 1000

  ! EN SERIE

  !INTEGER, PARAMETER :: mnx= ngx, mny= ngy, mnz= ngz ! mny= ngy/nsy,!mnz= ngz/nsz if nsub_proc > 1

  ! EN PARALLEL

  INTEGER, PARAMETER :: mnx= ngx, mny= ngy/nsy, mnz= ngz/nsz ! mny= ngy/nsy,!mnz= ngz/nsz if nsub_proc = 1

  INTEGER, PARAMETER :: mnmax= max(mny, mnz)
    
  INTEGER, PARAMETER :: nxs= 10*mnx+2, bw= 16 !bw es el ancho de banda de la matriz
  !                                           de coeficientes de cada tira, el 1
  !                                           tuve que agregarlo para tener menos lios
  !                                           al pasar de matrices locales a globales 
END MODULE sizes
  
  
