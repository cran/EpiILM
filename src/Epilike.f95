    module  likesubprograms
    use ISO_C_BINDING
    implicit none
    public :: like, likesir,likecon,likeconsir

    contains
    !####################

    !subroutine for SI likelihood

    subroutine like(x,y,tau,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,val) bind(C, name="like_")
    implicit none

    !Declarations
    integer (C_INT), intent(in) :: n,tau(n),tmax,ns,ni,tmin
    real (C_DOUBLE), intent(in) :: alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in) :: x(n), y(n)
    real (C_DOUBLE), intent(in) ::covmat(n,ns)
    real (C_DOUBLE), intent(out):: val

    integer::i,j,t
    double precision :: eu(n,n),Somega(n)
    double precision  :: dx1,dx2, p1,p2



    !Calculate the distance matrix
    do i=1,n
    do j=i,n
    eu(i,j)= sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
    eu(j,i)=eu(i,j)
    end do
    end do


    Somega= matmul(covmat,alpha)

    ! Calculate likelihood
    val=0.0d0
    do t=tmin,(tmax-1)
    do i=1,n
    !infectious period
    if (tau(i)==(t+1)) then
    dx1=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
    dx1=dx1+(eu(i,j)**(-beta(ni)))
  end if
    end do
    p1=1.0d0-exp(-(Somega(i)*dx1+spark))
    val=val+log(p1)
  end if
    !susceptible period
    if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
    dx2=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
    dx2=dx2+(eu(i,j)**(-beta(ni)))
  end if
    end do
    p2=exp(-(Somega(i)*dx2+spark))
    val=val+log(p2)
  end if
    end do
    end do

    end subroutine like
    !####################

    !subroutine for SIR likelihood

   subroutine likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,val) &
       & bind(C, name="likesir_")

    implicit none

    !Declarations
    integer (C_INT), intent(in) :: n,tau(n),lambda(n),tmax,ns,ni,tmin
    real (C_DOUBLE), intent(in) :: alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in) :: x(n), y(n)
    real (C_DOUBLE), intent(in) ::covmat(n,ns)
    real (C_DOUBLE), intent(out):: val

    integer::i,j,t
    double precision :: eu(n,n),Somega(n)
    double precision  :: dx1,dx2, p1,p2



    !Calculate the distance matrix
    do i=1,n
    do j=i,n
    eu(i,j)= sqrt(((x(i)-x(j))**2)+((y(i)-y(j))**2))
    eu(j,i)=eu(i,j)
    end do
    end do


    Somega= matmul(covmat,alpha)

    ! Calculate likelihood
    val=0.0d0
    do t=tmin,(tmax-1)
    do i=1,n
    !infectious period
    if (tau(i)==(t+1)) then
    dx1=0.0d0
    do j=1,n
    if (tau(j) .NE. 0) then
    if ((tau(j) .LT. (t+1)) .AND. (tau(j)+lambda(j) .GE. (t+1))) then
    dx1=dx1+(eu(i,j)**(-beta(ni)))
    end if
    end if
    end do
    p1=1.0d0-exp(-(Somega(i)*dx1+spark))
    val=val+log(p1)
    end if
    !susceptible period
    if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
    dx2=0.0d0
    do j=1,n
    if (tau(j) .NE. 0) then
    if ((tau(j) .LT. (t+1)) .AND. (tau(j)+lambda(j) .GE. (t+1)))then
    dx2=dx2+(eu(i,j)**(-beta(ni)))
    end if
    end if
    end do
    p2=exp(-(Somega(i)*dx2+spark))
    val=val+log(p2)
    end if
    end do
    end do

    end subroutine likesir


  !###############
    subroutine likecon(tau,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,val) &
    & bind(C, name="likecon_")

    implicit none

    !Declarations
    integer (C_INT), intent(in) :: n,tmax,ns,ni,tmin
    integer (C_INT), intent(in) :: tau(n)
    real (C_DOUBLE), intent(in) :: alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in) :: covmat(n,ns)
    real (C_DOUBLE), intent(in) :: network(n,n,ni)
    real (C_DOUBLE), intent(out) ::  val

    integer::i,j,t,k
    double precision  :: dx1,dx2, p1,p2,Somega(n)


    Somega= matmul(covmat,alpha)
    val=0.0d0


    do t=tmin,(tmax-1)
    do i=1,n
    !infectious period
    if (tau(i)==(t+1)) then
    dx1=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0)) then
    do k=1,ni
    dx1=dx1+ beta(k)*network(i,j,k)
    end do
  end if
    end do
    p1=1.0d0- exp(-(Somega(i)*dx1+spark))
    val=val+log(p1)
  end if
    !susceptible period
    if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
    dx2=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0))then
    do k=1,ni
    dx2=dx2+ beta(k)*network(i,j,k)
    end do
  end if
    end do
    p2=exp(-(Somega(i)*dx2+spark))
    val=val+log(p2)
  end if

    end do
    end do



    end subroutine likecon


    !#################################

    subroutine likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,&
         & val) bind(C, name="likeconsir_")

    implicit none


    !Declarations
    integer (C_INT), intent(in) :: n,tmax,ns,ni,tmin
    integer (C_INT), intent(in) :: tau(n),lambda(n)
    real (C_DOUBLE), intent(in) :: alpha(ns), beta(ni),spark
    real (C_DOUBLE), intent(in) ::  covmat(n,ns)
    real (C_DOUBLE), intent(in) :: network(n,n,ni)
    real (C_DOUBLE), intent(out) ::  val

    integer::i,j,t,k
    double precision  :: dx1,dx2, p1,p2,somega(n)

    somega= matmul(covmat,alpha)

    val=0.0d0



    do t=tmin,(tmax-1)
    do i=1,n
    !infectious period
    if (tau(i)==(t+1)) then
    dx1=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0) .AND. &
    & (tau(j)+lambda(j) .GE. (t+1))) then
    do k=1,ni
    dx1=dx1+ beta(k)*network(i,j,k)
    end do
  end if
    end do
    p1=1.0d0- exp(-(somega(i)*dx1+spark))
    val=val+log(p1)
  end if
    !susceptible period
    if ((tau(i) .GT. (t+1)) .OR. (tau(i)== 0))then
    dx2=0.0d0
    do j=1,n
    if ((tau(j) .LT. (t+1)) .AND. (tau(j) .NE. 0) .AND.&
    & (tau(j)+lambda(j) .GE. (t+1)))then
    do k=1,ni
    dx2=dx2+ beta(k)*network(i,j,k)
    end do
  end if
    end do
    p2=exp(-(somega(i)*dx2+spark))
    val=val+log(p2)
  end if
    end do
    end do

  

    end subroutine likeconsir




    !#######################################

    end module likesubprograms









