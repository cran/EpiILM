    
    module  subprogramr0
    use ISO_C_BINDING
    implicit none
    public :: rxysir,rconsir

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

    subroutine rxysir (n,tmax,ns,ni,alpha,beta,spark,covmat,lambda,x,y,&
                    & sim,val,countinf) bind(C, name="rxysir_")
    implicit none

    integer (C_INT),intent(in) :: n,tmax,ns,ni,sim
    integer (C_INT), intent(in) ::  lambda(n)

    real (C_DOUBLE), intent(in) :: alpha(ns),beta(ni),spark
    real (C_DOUBLE), intent(in)::covmat(n,ns)
    real (C_DOUBLE), intent(in) :: x(n),y(n)

    real (C_DOUBLE), intent(inout) :: val
    integer (C_INT), intent(out):: countinf(sim)


    integer::i,j,t,tau(n),A
    double precision :: u, dx, p
    double precision :: eu(n,n), Somega(n)

    call initrandomseed()

    !Calculate the distance matrix
    do i=1,n
        do j=i,n
            eu(i,j)= sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
            eu(j,i)=eu(i,j)
        end do
    end do


    Somega= matmul(covmat,alpha)

    do j=1,sim  !simulation starts here

    call random_number(u)
    A=int((u*tmax)+1)

    do i=1,n
    tau(i)=0
    end do
    tau(A)=1

    !Calculate the initial probability of susceptible being exposed and update newtau
    do t=1,tmax
        do i=1,n
        dx=0.0d0
        if (tau(i)==0) then
            if ((tau(A) .LE. t) .and. ((tau(A)+lambda(A)) .GT. t)) then
                dx=eu(i,A)**(-beta(ni))
            p=1.0d0 - exp(-(Somega(i)*dx+spark))
            call random_number(u)
            if (p .GT. u) then
                tau(i)=t+1
            end if
        end if
        end if
        end do
    end do

    countinf(j)=0
    do i=1,n
    if (tau(i) .NE. 0) then
        countinf(j)=countinf(j)+1
    end if
    end do

    countinf(j)=countinf(j)-1

    end do ! simulation ends


    val = sum(countinf)/dble(sim)


    end subroutine rxysir

!###############################

    subroutine rconsir(n,tmax,ns,ni,lambda,alpha,beta,spark,covmat,network,&
                & sim,val,countinf) bind(C, name="rconsir_")

    implicit none

    integer (C_INT),intent(in):: n,tmax,ns,ni,sim
    integer (C_INT),intent(in):: lambda(n)
    real (C_DOUBLE), intent(in)  ::  alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in)  ::  network(n,n,ni),covmat(n,ns)

    real (C_DOUBLE), intent(inout) :: val
    integer (C_INT),intent(out) :: countinf(sim)


    integer::i,j,t,k,tau(n),A
    double precision :: u, dx, p, Somega(n)


    call initrandomseed()

    Somega= matmul(covmat,alpha)

    do j=1,sim  ! simulation starts

    call random_number(u)
    A=int((u*tmax)+1)

    do i=1,n
    tau(i)=0
    end do
    tau(A)=1

    do t=1,tmax
        do i=1,n
        dx=0.0d0
        if (tau(i)==0) then
           if ((tau(A) .LE. t) .and. ((tau(A)+lambda(A)) .GT. t)) then
                    do k=1,ni
                    dx=dx+ beta(k)*network(i,A,k)
                    end do
            p=1.0d0 - (exp(-(Somega(i)*dx+spark)))
            call random_number(u)
            if (p .GT. u) then
            tau(i)=t+1
            end if
           end if
        end if
        end do
    end do

    countinf(j)=0
    do i=1,n
    if (tau(i) .NE. 0) then
    countinf(j)=countinf(j)+1
    end if
    end do

    countinf(j)=countinf(j)-1

    end do ! simulation ends


    val = sum(countinf)/dble(sim)

 end subroutine rconsir


!###############################

    end module subprogramr0







