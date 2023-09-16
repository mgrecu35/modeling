!  SFM  04/06/2013  Code module added for WOlson
!
! -------------------------- fband.f90 -------------------------------
!       Translated from C language by J.-P. Moreau, March, 1999.    
!   Reference: Numerical Algotithms with C, Springer, 1996 [BIBLI 11]

!  Subroutine band(mode,n,ld,ud,pmat,b,perm,signd,code)
Subroutine band(mode,n,ww,bb,code)
  implicit none
  integer maxlyr, ld, n, mode
  INTEGER code,ud,signd,perm(0:n-1)
  REAL*8 pmat(0:n-1,0:n-1),b(0:n-1)
  PARAMETER    ( MAXLYR = 250 )
  real *8 BB(2*MAXLYR), WW(2*MAXLYR,2*MAXLYR) !WSO 04/07/2013
  integer:: i, k
  INTEGER rc
  ld=2
  ud=2
  pmat=0
  do i=0,n-1
     b(i)=bb(i+1)
     pmat(i,2)=ww(i+1,i+1)
     do k=-2,2
        if(i+1+k>0 .and. i+1+k<=n) then
           pmat(i,2+k)=ww(i+1,i+1+k)
        endif
     enddo
  enddo
  
!======================================================================
!*                                                                    *
!*  The subroutine band solves a linear banded system:  pmat * x = b  *
!*  Here pmat is a nonsingular n x n matrix in condensed form, i.e.   *
!*  represented in an ld+ud+1 x n matrix for its ld lower and ud upper*
!*  co-diagonals. b denotes the right hand side of the system, and x  *
!*  is the solution.                                                  *
!*                                                                    *
!*  band uses the Gauss algorithm with column pivot search.           *
!*  The result of pivoting are min( ud, ld) additional columns,       *
!*  so that pmat needs all in all a n x (ld+1+ud+min(ld,ud)) matrix.  *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Applications:                                                    *
!*   -------------                                                    *
!*      Solve linear systems with nonsingular banded system matrices. *
!*      Particularly useful for large sparse and banded matrices with *
!*      n >> ld+1+ud.                                                 *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Control parameter:                                               *
!*   ------------------                                               *
!*      mode     integer mode : calling modus for band.               *
!*       = 0     factor matrix and solve linear system                *
!*       = 1     factor only, store factors in pmat                   *
!*       = 2     solve linear system only; for this call, the factors *
!*               must already be available in pmat, such as when      *
!*               many systems are solved for differing right hand     *
!*               sides and the same system matrix.                    *
!*                                                                    *
!*   Input parameters:                                                *
!*   -----------------                                                *
!*      n        integer n  ( n > 2 )                                 *
!*               Dimension of pmat, size of b                         *
!*      ld       integer ld; ( ld >= 0 )                              *
!*               number of lower co-diagonals                         *
!*      ud       integer ud; ( ud >= 0 )                              *
!*               number of upper co-diagonals                         *
!*      pmat     REAL*8 pmat(0:n-1,0:n-1)                             *
!*               mode = 0, 1:                                         *
!*               Matrix of the system in condensed form.              *
!*               Each row has length at least ld + 1 + ud + min(ld,ud)*
!*               where the columns 0, .., ld-1 denote the lower       *
!*               co-diagonals, column ld stores the diagonal and the  *
!*               columns ld+1, .., ld+ud contain the upper            *
!*               co-diagonals.                                        *
!*               If A is the original uncondensed band matrix, then : *
!*               A[i][k] = pmat[i][ld+k-i],                           *
!*                           for k,i inside the band                  *
!*               mode = 2                                             *
!*               LU factors in condensed form.                        *
!*      b        REAL*8 b(0:n-1)           ( for mode = 0, 2 )        *
!*               Right hand side                                      *
!*      perm     integer perm(0:n-1)       ( for mode = 2 )           *
!*               row permutation vector                               *
!*      signd    integer signd             ( for mode= 2 )            *
!*               sign of perm. The determinant of A can be computed   *
!*               as the product of the diagonal entries times signd.  *
!*                                                                    *
!*   Output parameters:                                               *
!*   ------------------                                               *
!*      pmat     REAL*8 pmat(0:n-1,0:n-1)  ( for mode= 0, 1 )         *
!*               LU factorization in condensed form                   *
!*      perm     integer perm[n]           ( for mode = 0, 1 )        *
!*               row permutation vector                               *
!*      b        REAL*8 b(0:n-1)           ( for mode = 0, 2 )        *
!*               solution vector for the system                       *
!*      signd    integer signd;            ( for mode = 0, 1 )        *
!*               sign of perm              ( value:  1 or -1 )        *
!*                                                                    *
!*   Return value code                                                *
!*   -----------------                                                *
!*      = 0      all ok                                               *
!*      = 1      n < 3 or other incorrect input parameter             *
!*      = 2      not used here                                        *
!*      = 3      Matrix is numerically singular                       *
!*      = 4      wrong calling modus                                  *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Subroutines used:                                                *
!*   ----------------                                                 *
!*      banddec : Factor matrix                                       *
!*      bandsol : solve linear system                                 *
!*                                                                    *
!======================================================================
 

  if (n < 1.or.ld < 0.or.ud < 0) then
    code=1
    return
  end if

  SELECT CASE (mode)
    CASE (0)  ! factor and solve system ............................... 
      call banddec (n, ld, ud, pmat, perm, signd,rc)
      if (rc.eq.0) then
        call bandsol (n, ld, ud, pmat, b, perm, code)
      else
        code=rc
      end if
    CASE (1)   ! factor only ........................................... 
       call banddec (n, ld, ud, pmat, perm, signd, code)
    CASE (2)   ! solve only ............................................ 
       call bandsol (n, ld, ud, pmat, b, perm, code)
    CASE DEFAULT
       code = 4
  END SELECT

  do i=0,n-1
     bb(i+1)=b(i)
  enddo
  return
  end

  Subroutine banddec(n,ld,ud,pmat,perm,signd,code)
  PARAMETER(MACH_EPS=1.D-15)
  INTEGER code,ud,signd,perm(0:n-1)
  REAL*8 pmat(0:n-1,0:n-1)
!======================================================================
!*                                                                    *
!*  The subroutine banddec factors a condensed banded matrix pmat.    *
!*  banddec uses the Gauss algorithm with column pivot search.        *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Input parameters:                                                *
!*   ----------------                                                 *
!*      n        integer n   ( n > 2 )                                *
!*               Dimension of pmat, size of b                         *
!*      ld       integer ld  ( ld >= 0 )                              *
!*               number of lower co-diagonals                         *
!*      ud       integer ud  ( ud >= 0 )                              *
!*               number of upper co-diagonals                         *
!*      pmat     REAL*8 pmat(0:n-1,0:n-1)                             *
!*               Matrix of the system in comndensed form.             *
!*               Each row has length at least ld + 1 + ud + min(ld,ud)*
!*               where the columns 0, .., ld-1 denote the lower       *
!*               co-diagonals, column ld stores the diagonal and the  *
!*               columns ld+1, .., ld+ud contain the upper            *
!*               co-diagonals.                                        *
!*      perm     integer perm(0:n-1) : row permutation vector         *
!*      signd    integer signd (value: 1 or -1)                       *
!*               sign of perm. The determinant of A can be computed   *
!*               as the product of the diagonal entries times signd.  *
!*                                                                    *
!*   Output parameters:                                               *
!*   ------------------                                               *
!*      pmat     REAL*8 pmat(0:n-1,0:n-1)                             *
!*               LU factorization in condensed form                   *
!*      perm     integer perm(0:n-1)                                  *
!*               row permutation vector                               *
!*      signd    integer signd : sign of perm                         *
!*                                                                    *
!*   Return value :                                                   *
!*   -------------                                                    *
!*      = 0      all ok                                               *
!*      = 1      n < 3 or other incorrect input parameter             *
!*      = 2      not used here                                        *
!*      = 3      Matrix is numerically singular                       *
!*                                                                    *
!*   Function used : user defined SWAP()                              *
!*   -------------                                                    *
!======================================================================
  INTEGER j0, mm, up, istart, iend, step, kstart,kend, kjend, km, jm, jk
  INTEGER k, j, i, ii
  REAL*8 piv

  if (ld < 0.or.ud < 0.or.n < 1) then
    code=1                              ! Invalid input parameter     
    return
  end if

  mm = ld + 1 + ud + min (ld, ud);

  if (ld <= ud) then                   ! up = 0 ==> transform integero    
    up=1                               ! lower triangular matrix      
  else
    up=0
  end if

  do i = 0, n-1
    do k = ld + ud + 1, mm-1           ! initialize needed extra
      pmat(i,k) = 0.D0                 ! columns    
    end do
  end do

  signd = 1                               ! initialize signd        
  if (up.eq.1) then 
    istart = 0; iend = n-1; step = 1      ! Find start, end and     
    kstart = 1                            ! directions of loops                                             ! depending of up         
  else
    istart = n-1; iend = 0; step = -1
    kstart = -1
  end if
  
  i = istart
  do while (i.ne.iend)
    !kend = (up ? min (ld+1, n-i) : max (-i-1, -ud-1));  
    if (up.eq.1) then 
	  kend = min (ld+1, n-i)
    else 
	  kend = max (-i-1, -ud-1)
    end if

    j0 = 0
    piv = DABS (pmat(i,ld))                ! choose pivot           

    k=kstart
    do while (k.ne.kend)
      if (DABS(pmat(k+i,ld-k)) > piv)  then
        piv = DABS(pmat(k+i,ld-k))
        j0 = k
      end if
      k = k + step
    end do

    if (piv < MACH_EPS) then
      code = 3                            ! If piv = 0, matrix is   
      return                              ! singular                
    end if

    perm(i) = j0
    !kjend =  (up ? min (j0+ud+1, n-i) : max ( -i-1, j0-ld-1));  
    if (up.eq.1) then 
	  kjend = min (j0+ud+1, n-i)
    else 
	  kjend = max ( -i-1, j0-ld-1)
    end if

    if (j0.ne.0) then
      signd = - signd                     ! swap rows                
      k = 0
      do while (k.ne.kjend)
        km = k + ld
        if (km < 0) then 
		  km = km + mm
        end if

        call SWAP(pmat(i,km), pmat(i+j0,k+ld-j0))  !see basis_r.f90

        k = k + step
      end do
    end if

    k = kstart
    do while (k.ne.kend)
                                                   ! below row i        
      pmat(k+i,ld-k) = pmat(k+i,ld-k) / pmat(i,ld)
      j=kstart
      do while (j.ne.kjend)
                                                   ! loop over all        
        jk = j + ld - k                            ! columns right of i   
        jm = j + ld
                                      ! additional columns from pivoting  
        if (jk < 0) then              ! are stored to right of column     
		  jk = jk +mm
        end if
        if (jm < 0) then 
		  jm = jm + mm                !  ud+ld+1                          
        end if 

        pmat(k+i,jk) = pmat(k+i,jk) - pmat(k+i,ld-k) * pmat(i,jm)

        j = j + step
      end do  ! while j
	  k = k + step  
    end do  ! do while k
	i = i + step  
  end do  ! do while i  

  if (up.eq.1) then 
    ii = n-1 
  else 
    ii=0
  end if	

  piv = DABS (pmat(ii,ld))                     ! choose pivot             
                                               ! If piv = 0, matrix is    
  if (piv < MACH_EPS) then                     ! singular                 
    code=3
    return
  end if

  perm(iend)=0.d0
  code = 0
  return
  end


  Subroutine bandsol(n,ld,ud,pmat,b,perm,code)
  PARAMETER(MACH_EPS=1.D-15)
  INTEGER code,ud,perm(0:n-1)
  REAL*8 pmat(0:n-1,0:n-1),b(0:n-1)
!======================================================================
!*                                                                    *
!*  The subroutine bandsol solves a factored linear banded system in  *
!*  condensed form using banddec.                                     *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Input parameters:                                                *
!*   -----------------                                                *
!*      n        integer n;  ( n > 2 )                                *
!*               Dimension of pmat, size of b                         *
!*      ld       integer ld; ( ld >= 0 )                              *
!*               number of lower co-diagonals                         *
!*      ud       integer ud; ( ud >= 0 )                              *
!*               number of upper co-diagonals                         *
!*      pmat     REAL*8  pmat(0:n-1,0:n-1)                            *
!*               Matrices for the factored system in comndensed form. *
!*      perm     integer perm(0:n-1)                                  *
!*               row permutation vector                               *
!*                                                                    *
!*   Output parameters:                                               *
!*   ------------------                                               *
!*      b        REAL*8 b(0:n-1)                                      *
!*               solution vector of linear system                     *
!*                                                                    *
!*   Return value :                                                   *
!*   -------------                                                    *
!*      = 0      all ok                                               *
!*      = 1      n < 3 or other incorrect input parameter             *
!*                                                                    *
!====================================================================== 
  INTEGER i, k, s, mm, up, istart, iend, step, kstart, kend, km          
		                               
  if (ld < 0.or.ud < 0.or.n < 1) then    ! invalid input            
    code = 1
    return
  end if

  mm = ld + ud + 1 + min (ld, ud)        ! mm = max. column number     

  if (ld <= ud) then                     ! up = 0 ==> transform into    
    up=1                                 ! lower triangular matrix      
  else
    up=0
  end if

  if (up.eq.1) then
    istart = 0; iend = n-1; step = 1     ! determine bounds and       
    kstart = 1; s = -1                   ! direction of loop          
  else                                   ! depending on up            
    istart = n-1; iend = 0; step = -1
    kstart = -1; s = 1
  end if

  i=istart
  do while (i.ne.iend)
    if (perm(i).ne.0) then
      call SWAP( b(i), b(i+perm(i)))
    end if

    !kend = (up ? min (ld+1, n-i) : max (-i-1, -ud-1));  
    if (up.eq.1) then 
	  kend = min(ld+1, n-i)
    else 
	  kend = max(-i-1, -ud-1)
    end if

    k=kstart
    do while (k.ne.kend)
      b(k+i) = b(k+i) - pmat(k+i,ld-k) * b(i)
      k=k+step
    end do
    i=i+step
  end do

  i=iend
  do while (i.ne.istart+s)
    !kend =  (up ? min (ld+ud+1, n-i) : max (-i-1, -ud-ld));  
    if (up.eq.1) then 
	  kend = min(ld+ud+1, n-i)
    else 
	  kend = max(-i-1, -ud-ld)
    end if

    k=kstart
    do while (k.ne.kend)
      km = k + ld                      ! update and                   
      if (km < 0) km = km + mm         ! back substitute              
      b(i) = b(i) - pmat(i,km) * b(i+k)
      k=k+step
    end do
    if (DABS(pmat(i,ld)) > MACH_EPS) then
      b(i) = b(i) / pmat(i,ld)
    else
      b(i)=0.d0
    end if 
    i=i-step
  end do

  code = 0
  return
  end  !bandsol 

!***********************
! Swap two real values
!***********************
  Subroutine Swap(a,b)
  REAL*8 a,b,temp
  temp=a
  a=b
  b=temp
  return
  end subroutine Swap


! --------------------------- END fband.f90 -------------------------  
