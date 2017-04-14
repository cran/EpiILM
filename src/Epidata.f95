
    module  subprograms
    use ISO_C_BINDING
    implicit none
    public :: dataxy,dataxysir,datacon,dataconsir

    contains
!####################
    subroutine initrandomseed()
    !The seed for the random number generation method random_number() has been reset

    implicit none

    integer :: i
    integer :: n
    integer :: clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)

    deallocate(seed)
    end subroutine initrandomseed

!####################
    subroutine dataxy (x,y,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,tau) bind(C, name="dataxy_")

    implicit none

    integer (C_INT),intent(in):: n,tmax,ns,ni,tmin
    real (C_DOUBLE), intent(in)  :: alpha(ns),beta(ni),spark
    real (C_DOUBLE), intent(in) ::covmat(n,ns)
    real (C_DOUBLE), intent(in) :: x(n),y(n)
    integer (C_INT),intent(inout) :: tau(n)

    integer::i,j,t
    double precision :: u, dx, p
    double precision :: eu (n,n), Somega(n)


    !Calculate the distance matrix
    do i=1,n
        do j=i,n
            eu(i,j)= sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i)=eu(i,j)
        end do
    end do

    call initrandomseed()

    Somega= matmul(covmat,alpha)

    !Calculate the probabilities of susceptible being exposed and update tau
    do t=tmin,tmax
        do i=1,n
            if (tau(i)==0) then
            dx=0.0d0
            do j=1,n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                dx=dx+(eu(i,j)**(-beta(ni)))
                end if
            end do
            p=1.0d0 - exp(-(Somega(i)*dx + spark))
            call random_number(u)
            if (p>u) then
            tau(i)=t+1
            end if
            end if
        end do
    end do


    end subroutine dataxy
!####################

    subroutine dataxysir (n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,lambda,x,y,tau,remt) bind(C, name="dataxysir_")
    implicit none

    integer (C_INT),intent(in) :: n,tmax,ns,ni,tmin
    integer (C_INT), intent(inout) :: lambda(n)
    real (C_DOUBLE), intent(in) :: alpha(ns),beta(ni),spark
    real (C_DOUBLE), intent(in)::covmat(n,ns)
    real (C_DOUBLE), intent(in) :: x(n),y(n)
    integer (C_INT), intent(inout)  ::  remt(n),tau(n)

    integer::i,j,t
    double precision :: u, dx, p
    double precision :: eu (n,n), Somega(n)

    call initrandomseed()

    !Calculate the distance matrix
    do i=1,n
        do j=i,n
            eu(i,j)= sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i)=eu(i,j)
        end do
    end do


    Somega= matmul(covmat,alpha)

    !Calculate the probabilities of susceptible being exposed and update tau
    do t=tmin,tmax
        do i=1,n
            if (tau(i)==0) then
            dx=0.0d0
            do j=1,n
                if ((tau(j) .NE. 0) ) then
                if ((tau(j) .LE. t) .and. ((tau(j)+lambda(j)) .GT. t)) then
                dx=dx+(eu(i,j)**(-beta(ni)))
                end if
                end if
            end do
            p=1.0d0- exp(-(Somega(i)*dx+spark))
            call random_number(u)
            if (p>u) then
            tau(i)=t+1
            remt(i)=t+1+lambda(i)
            end if
            end if
        end do
    end do

    end subroutine dataxysir
!####################

    subroutine datacon(n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,network,&
                         & tau) bind(C, name="datacon_")
    implicit none

    integer (C_INT),intent(in):: n,ns,ni,tmax,tmin
    real (C_DOUBLE), intent(in)  :: alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in) :: network(n,n,ni),covmat(n,ns)
    integer (C_INT),intent(inout)::  tau(n)

    integer::i,j,t,k
    double precision :: u, dx,p,Somega(n)


    call initrandomseed()

    Somega= matmul(covmat,alpha)

    do t=tmin,tmax
        do i=1,n
            if (tau(i)==0) then
            dx=0.0d0
            do j=1,n
                if ((tau(j) .NE. 0) .and. (tau(j) .LE. t)) then
                    do k=1,ni
                    dx=dx+ beta(k)*network(i,j,k)
                    end do
                end if
            end do
            p=1.0d0 -exp(-(Somega(i)*dx+spark))
            call random_number(u)
            if (p>u) then
            tau(i)=t+1
            end if
            end if
        end do
    end do


    end subroutine datacon

!###############################

    subroutine dataconsir(n,tmin,tmax,ns,ni,lambda,alpha,beta,spark,covmat,network,&
                & tau,remt) bind(C, name="dataconsir_")

    implicit none

    integer (C_INT),intent(in):: n,tmax,ns,ni,tmin
    integer (C_INT),intent(in):: lambda(n)
    real (C_DOUBLE), intent(in)  ::  alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in)  ::  network(n,n,ni),covmat(n,ns)
    integer (C_INT),intent(inout)  :: tau(n)
    integer (C_INT), intent(inout)  ::  remt(n)

    integer::i,j,t,k
    double precision :: u, dx, p, Somega(n)


    call initrandomseed()

    Somega= matmul(covmat,alpha)


    do t=tmin,tmax
        do i=1,n
            if (tau(i)==0) then
            dx=0.0d0
            do j=1,n
                if (tau(j) .NE. 0)  then
                  if ((tau(j) .LE. t).and. ((tau(j)+lambda(j)) .GT. t))then
                    do k=1,ni
                    dx=dx+ beta(k)*network(i,j,k)
                    end do
                  end if
                end if
            end do
            p=1.0d0 - exp(-(Somega(i)*dx+spark))
            call random_number(u)
            if (p>u) then
            tau(i)=t+1
            remt(i)=t+1+lambda(i)
            end if
            end if
        end do
    end do

   

 end subroutine dataconsir


!###############################

    end module subprograms







