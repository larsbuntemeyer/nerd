
   !
   SUBROUTINE READ_AKBK
     use mo_grid, only: nz,nz1,ak,bk,dak,dbk,ptop,lvl,lvl2
     use mo_namelist
     implicit none
     integer :: k,ke,ke1
     ke=nz
     ke1=nz1
     ! read aks and bks from file
     OPEN(UNIT=110,FILE=TRIM(AKBK_FILE))
     K = 1
     DO WHILE (K <= KE1)
       READ(110,*) AK(K), BK(K)        
       K=K+1                          
     END DO
     CLOSE(UNIT=110)                 
     PTOP = AK(1)
     LVL2(1) = KE1                  
     DO K = 1,KE 
        DAK(K) = -AK(K) + AK(K+1)  
        DBK(K) = -BK(K) + BK(K+1) 
        LVL(K) = KE-K+1
        LVL2(K+1) = KE1-K        
     END DO
   END SUBROUTINE READ_AKBK     
   !
