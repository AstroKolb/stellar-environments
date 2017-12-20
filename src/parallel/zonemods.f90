module gridsize

! integer, parameter :: imax = 512, jmax = 128, kmax = 384   ! Memory dimensions
 integer, parameter :: imax = 196, jmax = 48,  kmax = 144
! integer, parameter :: imax = 128, jmax = 32,  kmax = 96
! integer, parameter :: imax = 96,  jmax = 24,  kmax = 72
 integer, parameter :: pey =  4, pez =  4                  ! number of MPI tasks

end module gridsize

module zone 
!=======================================================================
! (formerly zone.h) global (3D) data arrays
!======================================================================= 

 use gridsize
 
 integer :: isy, isz, js, ks
 integer :: npe, npeyy, npey, npez               ! # of pes 
 integer :: mype, mypeyy, mypey, mypez           ! local pe number
 integer :: ngeomx, ngeomy, ngeomz               ! XYZ Geometry flag
 integer :: nleftx, nlefty, nleftz               ! XYZ Lower Boundary Condition
 integer :: nrightx,nrighty,nrightz              ! XYZ Upper Boundary Condition
 integer :: mpi_comm_row, mpi_comm_col, mpi_comm_yy
 integer :: vh1_datatype
 integer :: krow, jcol
 
 integer, dimension(jmax) :: knmin, knmax
 integer, dimension(kmax) :: jnmin, jnmax
 integer, dimension(jmax,kmax) :: ongrid
 
 real, dimension(imax,jmax/pey,kmax/pez) :: zro, zpr, zux, zuy, zuz, zfl, zcl
 real, dimension(6,   jmax/pey,kmax/pez) :: rin, pin, cin, uin, vin, win
 real, dimension(imax) :: zxa, zdx, zxc
 real, dimension(jmax) :: zya, zdy, zyc
 real, dimension(kmax) :: zza, zdz, zzc
 real, dimension(imax) :: kap
  
 real, dimension(12,kmax/pez,imax/pey) :: ybr, ybp, ybu, ybv, ybw, ybf, ybc
 real, dimension(12,jmax/pey,imax/pez) :: zbr, zbp, zbu, zbv, zbw, zbf, zbc

 real, dimension(7,kmax/pez,jmax/pey,imax) :: send1
 real, dimension(7,kmax/pez,imax/pey,jmax) :: send2
 real, dimension(7,jmax/pey,kmax/pez,imax) :: send3
 real, dimension(7,jmax/pey,imax/pez,kmax) :: send4

 real, dimension(7,kmax/pez,jmax/pey,imax/pey,pey) :: recv1
 real, dimension(7,kmax/pez,imax/pey,jmax/pey,pey) :: recv2
 real, dimension(7,jmax/pey,kmax/pez,imax/pez,pez) :: recv3
 real, dimension(7,jmax/pey,imax/pez,kmax/pez,pez) :: recv4
 
 equivalence ( send1(1,1,1,1),   send2(1,1,1,1),   send3(1,1,1,1),   send4(1,1,1,1) )
 equivalence ( recv1(1,1,1,1,1), recv2(1,1,1,1,1), recv3(1,1,1,1,1), recv4(1,1,1,1,1) )


end module zone

module yinyang
!=======================================================================
! Arrays needed to pass boundary information between Yin and Yang grids
!======================================================================= 
 
 use gridsize
 
 integer :: mdy, mdz, syc_tot, szc_tot
 integer, dimension(20000) :: jysource, kysource, kzsource, jzsource
 integer, dimension(pez) :: sycount, rycount, syoffset, ryoffset
 integer, dimension(pey) :: szcount, rzcount, szoffset, rzoffset
 integer, dimension(48*kmax/pez) :: jytarget, kytarget, nytarget
 integer, dimension(48*jmax/pey) :: jztarget, kztarget, nztarget
 real, dimension(12,kmax/pez) :: cy1, cy2, cy3, cy4, cy5, cy6
 real, dimension(jmax/pey,12) :: cz1, cz2, cz3, cz4, cz5, cz6
 real, dimension(:,:,:), allocatable :: sybuff, szbuff, sybufg, szbufg
 real, dimension(7,imax/pey,48*kmax/pez) :: rybuff  ! 48 = (6+6)*4 neighbors
 real, dimension(7,imax/pez,48*jmax/pey) :: rzbuff
 
end module yinyang
