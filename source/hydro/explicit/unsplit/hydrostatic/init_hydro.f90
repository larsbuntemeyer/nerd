!
!----------------------------------------------------------------------------------------------
!
subroutine init_hydro
   !
   use mo_dynamics
   use mo_database, only: ps, hydrop, hydrodp, hydromp
   use mo_grid, only: ke=>nz, ak, bk, akh, bkh, &
                      dak, dbk, ptop
   !
   implicit none
   !
   integer :: k
   real :: zpno,zpnu
   real, external :: getp
   !
   call allocate_dynamics
   !
   DO K = 1,KE
      ZPNO   = GETP(AK(K  ), BK(K  ), 1.0E5, PTOP)
      ZPNU   = GETP(AK(K+1), BK(K+1), 1.0E5, PTOP)
      IF(850.0E2 .GE. ZPNO .AND. 850.0E2 .LT. ZPNU) KFL850 = K
      IF(800.0E2 .GE. ZPNO .AND. 800.0E2 .LT. ZPNU) KFL800 = K
      IF(500.0E2 .GE. ZPNO .AND. 500.0E2 .LT. ZPNU) KFL500 = K
      IF(400.0E2 .GE. ZPNO .AND. 400.0E2 .LT. ZPNU) KFL400 = K
      IF(300.0E2 .GE. ZPNO .AND. 300.0E2 .LT. ZPNU) KFL300 = K
   ENDDO
   !
   PS = 101325.
   !
   CALL PIBER(PS, 1, AK, BK, AKH, BKH, DAK, DBK,         &
        HYDROP, HYDRODP, HYDROMP)
   CALL PIBER(PS, 2, AK, BK, AKH, BKH, DAK, DBK,         &
        HYDROP, HYDRODP, HYDROMP)
        !
end subroutine init_hydro
