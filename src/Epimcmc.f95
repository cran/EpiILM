   !##################################### module begins ##############################

    module mcmcdata
    use ISO_C_BINDING
    implicit none
    public:: mcmc, conmcmc

    contains

!############################## mcmc ##############################


    subroutine mcmc(tnum,x,y,tau,n,lambda,tmin,tmax,nsim,aalpha,ns,ni,bbeta,covmat,&
        & prostda,prostdb,anum,bnum,halfvar,unifmin,unifmax,gshape,gscale,&
        & halfvarb,unifminb,unifmaxb,gshapeb,gscaleb,simalpha,simbeta, &
        & sspark,flag,prostdsp,snum,halfvarsp,unifminsp,unifmaxsp,gshapesp,&
        & gscalesp,simspark,llikeval) bind(C, name="mcmc_")
    implicit none

        !Declarations
    integer (C_INT) , intent(in) :: n, nsim, tmax, tau(n),ns,ni,tmin
    integer (C_INT) , intent(in) :: tnum,anum,bnum,lambda(n)
    real (C_DOUBLE), intent(in) :: aalpha(ns), bbeta(ni)
    real (C_DOUBLE), intent(in) ::  prostda(ns),prostdb
    real (C_DOUBLE), intent(in) :: halfvar(ns),gshape(ns),gscale(ns)
    real (C_DOUBLE), intent(in) :: unifmin(ns),unifmax(ns)
    real (C_DOUBLE), intent(in) :: halfvarb,gshapeb,gscaleb,unifminb,unifmaxb
    real (C_DOUBLE), intent(in) :: x(n),y(n), covmat(n,ns)
    real (C_DOUBLE), intent(out) :: simalpha(nsim,ns),simbeta(nsim,ni),simspark(nsim)
    real (C_DOUBLE), intent(out) :: llikeval(nsim)

    integer (C_INT) , intent(in) :: flag,snum
    real (C_DOUBLE), intent(in) :: sspark, prostdsp
    real (C_DOUBLE), intent(in) :: halfvarsp,unifminsp,unifmaxsp,gshapesp,gscalesp

    double precision :: alpha(ns),beta(ni)
    double precision:: z_alpha, z_beta, alpha_n(ns), beta_n(ni)
    double precision :: ratio,ratio_b,value,value_ini
    double precision :: u,psi,psib
    double precision :: value_b, value_ini_b,value_fb,value_fb_ini
    double precision :: value_f,value_f_ini
    integer :: m,i,check

    double precision :: spark,z_sp,spark_n,ratio_sp,psisp
    double precision :: value_sp,value_ini_sp, value_fsp,value_fsp_ini

!############# initial values of alpha and beta #############
    call initrandomseed()

    psi=0.d0
    psib=0.0d0
    psisp=0.0d0


    do i=1,ns
    alpha(i)= aalpha(i)
    end do
    beta(ni)=bbeta(ni)
    spark=sspark

    do i=1,ns
    simalpha(1,i)=alpha(i)
    end do
    simbeta(1,ni)=beta(ni)
    simspark(1)=spark

    do i=1,ns
    alpha_n(i)=alpha(i)
    end do

!############# likelihood for initial values  #############
    SELECT CASE (tnum)
    CASE (1)
    call like(x,y,tau,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,value_ini)
    CASE (2)
    call likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,value_ini)
    END SELECT


    !############# ############# loglikelihood  #############

    llikeval(1)= value_ini

    value_ini_b=0.0d0
    value_ini_sp=0.0d0

!############# ############# simulation starts  #############


    do m=1,(nsim-1)

!############# start do loop for checking each alpha value #############

        do i=1,ns

            z_alpha = rand_normal(0.0d0,prostda(i)) !generate random numbers from proposal

            alpha_n(i)= alpha(i) + z_alpha  ! update alpha

            if (alpha_n(i) .GT. 0.0d0) then

        !############# likelihood  #############

            SELECT CASE (tnum)
            CASE (1)
            call like(x,y,tau,n,tmin,tmax,ns,ni,alpha_n,beta,spark,covmat,value)
            CASE (2)
            call likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha_n,beta,spark,covmat,value)
            END SELECT
        !############# prior distribution selection #############

            SELECT CASE (anum)
            CASE (1)
            value_f= gamma_density(alpha_n(i),gshape(i),gscale(i))
            value_f_ini= gamma_density(alpha(i),gshape(i),gscale(i))
            !calcualte AP
            ratio= (value-value_ini)+(value_f-value_f_ini)
            psi= min(1.0d0,exp(ratio))

            CASE (2)
            value_f= half_normal(alpha_n(i),halfvar(i))
            value_f_ini= half_normal(alpha(i),halfvar(i))
            !calcualte AP
            ratio= (value-value_ini)+(value_f-value_f_ini)
            psi= min(1.0d0,exp(ratio))

            CASE (3)
            if ((unifmin(i) .LT. alpha_n(i)) .AND. (alpha_n(i).LT. unifmax(i))) then
            !calcualte AP
            ratio= (value-value_ini)
            psi= min(1.0d0,exp(ratio))
            else
            psi=-1.0d0
            end if
            END SELECT

        !#############  Acceptance probbaility check for alpha #############

            call random_number(u)
            if (psi>=u) then
            simalpha(m+1,i)=alpha_n(i)
            else
            simalpha(m+1,i)=alpha(i)
            end if

            else
            simalpha(m+1,i)=alpha(i)
            end if

            alpha_n(i)= simalpha(m+1,i)

    !############# end do loop for (ns): checking each alpha value #############
        end do


    !#############  beta parameter update #############

        z_beta= rand_normal(0.0d0,prostdb) !generate from proposal

        beta_n(ni)= beta(ni) + z_beta

    !############# start do loop for checking beta value #############

        if (beta_n(ni) >0.0d0) then

    !############# likelihood  #############
        SELECT CASE (tnum)
        CASE (1)
        call like(x,y,tau,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark,covmat,value_b)
        check=0
        do i=1,ns
            if (alpha_n(i) .EQ. alpha(i)) then
            check=check+1
            end if
        end do
        if (check .EQ. ns) then
        value_ini_b= value_ini
        else
        value_ini_b=value
        end if

        CASE (2)
        call likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark,covmat,value_b)
        check=0
        do i=1,ns
            if (alpha_n(i) .EQ. alpha(i)) then
            check=check+1
            end if
        end do
        if (check .EQ. ns) then
        value_ini_b= value_ini
        else
        value_ini_b=value
        end if
        END SELECT

    !############# prior distribution selection #############

        SELECT CASE (bnum)
        CASE (1)
        value_fb= gamma_density(beta_n(ni),gshapeb,gscaleb)
        value_fb_ini= gamma_density(beta(ni),gshapeb,gscaleb)
        !calcualte AP
        ratio_b= (value_b-value_ini_b)+(value_fb-value_fb_ini)
        psib= min(1.0d0,exp(ratio_b))

        CASE (2)
        value_fb= half_normal(beta_n(ni),halfvarb)
        value_fb_ini= half_normal(beta(ni),halfvarb)
        !calcualte AP
        ratio_b= (value_b-value_ini_b)+(value_fb-value_fb_ini)
        psib= min(1.0d0,exp(ratio_b))

        CASE (3)
        if ((unifminb .LT. beta_n(ni)).AND. (beta_n(ni) .LT. unifmaxb)) then
        !calcualte AP
        ratio_b= value_b-value_ini_b
        psib= min(1.0d0,exp(ratio_b))
        else
        psib=-1.0d0
        end if
        END SELECT

    !#############  Acceptance probbaility check for beta #############

        call random_number(u)
        if (psib >=u) then
        simbeta(m+1,ni)=beta_n(ni)
        else
        simbeta(m+1,ni)=beta(ni)
        end if

        else
        simbeta(m+1,ni)=beta(ni)
        end if

        beta_n(ni)= simbeta(m+1,ni)


!#############  spark parameter estimation #############


    if (flag .EQ. 1) then

    z_sp= rand_normal(0.0d0,prostdsp) !generate from proposal
    spark_n = spark + z_sp

    if (spark_n >0.0d0) then

    SELECT CASE (tnum)
    CASE (1)
    call like(x,y,tau,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark_n,covmat,value_sp)
    if (beta_n(ni) .EQ. beta(ni))  then
    value_ini_sp= value_ini_b
    else
    value_ini_sp=value_b
    end if
    CASE (2)
    call likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark_n,covmat,value_sp)
    if (beta_n(ni) .EQ. beta(ni))  then
    value_ini_sp= value_ini_b
    else
    value_ini_sp=value_b
    end if
    END SELECT
!############# prior distribution selection #############

    SELECT CASE (snum)
    CASE (1)
    value_fsp= gamma_density(spark_n,gshapesp,gscalesp)
    value_fsp_ini= gamma_density(spark,gshapesp,gscalesp)
    !calcualte AP
    ratio_sp= (value_sp-value_ini_sp)+(value_fsp-value_fsp_ini)
    psisp= min(1.0d0,exp(ratio_sp))

    CASE (2)
    value_fsp= half_normal(spark_n,halfvarsp)
    value_fsp_ini= half_normal(spark,halfvarsp)
    !calcualte AP
    ratio_sp= (value_sp-value_ini_sp)+(value_fsp-value_fsp_ini)
    psisp= min(1.0d0,exp(ratio_sp))

    CASE (3)
    if ((unifminsp .LT. spark_n).AND. (spark_n .LT. unifmaxsp)) then
    !calcualte AP
    ratio_sp= value_sp-value_ini_sp
    psisp= min(1.0d0,exp(ratio_sp))
    else
    psisp=-1.0d0
    end if
    END SELECT

!#############  Acceptance probbaility check for spark #############

    call random_number(u)
    if (psisp >=u) then
    simspark(m+1)=spark_n
    else
    simspark(m+1)=spark
    end if

    else
    simspark(m+1)=spark
    end if

    spark_n= simspark(m+1)
    else
    spark=0.0d0
    spark_n=0.0d0
    end if


    SELECT CASE (tnum)
    CASE (1)
    call like(x,y,tau,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark_n,covmat,value_ini)
    llikeval(m+1)= value_ini
    CASE (2)
    call likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha_n,beta_n,spark_n,covmat,value_ini)
    llikeval(m+1)= value_ini
    END SELECT

   do i=1,ns
   alpha(i)=alpha_n(i)
   end do
   beta(ni)=beta_n(ni)
   spark= spark_n


!############# simulation ends here #############
    end do

    end subroutine mcmc



!############# contact network #############


    subroutine conmcmc(tnum,tau,n,lambda,tmin,tmax,nsim,aalpha,ns,ni,bbeta,covmat,&
   		& network,prostda,prostdb,anum,bnum,halfvar,unifmin,unifmax,gshape,gscale,&
    	& halfvarb,unifminb,unifmaxb,gshapeb,gscaleb,simalpha,simbeta,&
        & sspark,flag,prostdsp,snum,halfvarsp,unifminsp,unifmaxsp,gshapesp,&
        & gscalesp,simspark,llikeval) bind(C, name="conmcmc_")
    implicit none

    !Declarations
    integer (C_INT) , intent(in) :: n, nsim, tmax,ns,ni,tmin
    integer (C_INT) , intent(in) :: tnum,anum,bnum,lambda(n),tau(n)
    real (C_DOUBLE), intent(in) :: aalpha(ns), bbeta(ni)
    real (C_DOUBLE), intent(in) ::  prostda(ns),prostdb(ni)
    real (C_DOUBLE), intent(in) :: halfvar(ns),gshape(ns),gscale(ns)
    real (C_DOUBLE), intent(in) :: unifmin(ns),unifmax(ns)
    real (C_DOUBLE), intent(in) :: halfvarb(ni),gshapeb(ni),gscaleb(ni)
    real (C_DOUBLE), intent(in) :: unifminb(ni),unifmaxb(ni)
    real (C_DOUBLE), intent(in) :: covmat(n,ns),network(n,n,ni)
    real (C_DOUBLE), intent(out) :: simalpha(nsim,ns),simbeta(nsim,ni),simspark(nsim)
    real (C_DOUBLE), intent(out) :: llikeval(nsim)

    integer (C_INT) , intent(in) :: flag,snum
    real (C_DOUBLE), intent(in) :: sspark, prostdsp
    real (C_DOUBLE), intent(in) :: halfvarsp,unifminsp,unifmaxsp,gshapesp,gscalesp

    double precision :: alpha(ns),beta(ni)
    double precision:: z_alpha, z_beta, alpha_n(ns), beta_n(ni)
    double precision :: ratio,ratio_b,value,value_ini
    double precision :: u,psi,psib
    double precision :: value_b,value_ini_b,value_fb,value_fb_ini
    double precision :: value_f,value_f_ini
    integer :: m,i,j,check

    double precision :: spark,z_sp,spark_n,ratio_sp,psisp
    double precision :: value_sp,value_ini_sp,value_fsp,value_fsp_ini

    !############# initial values of alpha and beta #############
    call initrandomseed()

    psi=0.d0
    psib=0.0d0
    psisp=0.0d0


    do i=1,ns
    alpha(i)= aalpha(i)
    simalpha(1,i)=alpha(i)
    alpha_n(i)=alpha(i)
    end do

    do i=1,ni
    beta(i)=bbeta(i)
    simbeta(1,i)=beta(i)
    beta_n(i)=beta(i)
    end do

    spark=sspark
    simspark(1)=spark


    !############# likelihood for initial values  #############
    SELECT CASE (tnum)
    CASE (1)
    call likecon(tau,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,value_ini)
    CASE (2)
    call likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,value_ini)
    END SELECT

  !############# ############# loglikelihood  #############

     llikeval(1)= value_ini

     value_ini_b=0.0d0
     value_ini_sp=0.0d0

  !############# ############# simulation starts  #############

    do m=1,(nsim-1)

        !############# start do loop for checking each alpha value #############

        do i=1,ns
            z_alpha = rand_normal(0.0d0,prostda(i)) !generate random numbers from proposal

            alpha_n(i)= alpha(i) + z_alpha  ! update alpha

            if (alpha_n(i) .GT. 0.0d0) then

        SELECT CASE (tnum)
        CASE (1)
        call likecon(tau,n,ns,ni,tmin,tmax,alpha_n,beta,spark,covmat,network,value)
        CASE (2)
        call likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha_n,beta,spark,covmat,network,value)
        END SELECT
            !############# prior distribution selection #############

            SELECT CASE (anum)
            CASE (1)
            value_f= gamma_density(alpha_n(i),gshape(i),gscale(i))
            value_f_ini= gamma_density(alpha(i),gshape(i),gscale(i))
            !calcualte AP
            ratio= (value-value_ini)+(value_f-value_f_ini)
            psi= min(1.0d0,exp(ratio))

            CASE (2)
            value_f= half_normal(alpha_n(i),halfvar(i))
            value_f_ini= half_normal(alpha(i),halfvar(i))
            !calcualte AP
            ratio= (value-value_ini)+(value_f-value_f_ini)
            psi= min(1.0d0,exp(ratio))

            CASE (3)
            if ((unifmin(i) .LT. alpha_n(i)).AND. (alpha_n(i) .LT. unifmax(i))) then
            !calcualte AP
            ratio= (value-value_ini)
            psi= min(1.0d0,exp(ratio))
            else
            psi=-1.0d0
            end if
            END SELECT

            !#############  Acceptance probbaility check for alpha #############

            call random_number(u)
            if (psi>=u) then
            simalpha(m+1,i)=alpha_n(i)
            else
            simalpha(m+1,i)=alpha(i)
            end if

            else
            simalpha(m+1,i)=alpha(i)
            end if

            alpha_n(i)= simalpha(m+1,i)

        !############# end do loop for (ns): checking each alpha value #############
        end do


        !#############  beta parameter update #############
        do i=1,ni

            z_beta= rand_normal(0.0d0,prostdb(i)) !generate from proposal

            beta_n(i)= beta(i) + z_beta   !update beta

            !############# start do loop for checking beta value #############

            if (beta_n(i) >0.0d0) then

        SELECT CASE (tnum)
        CASE (1)
            call likecon(tau,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark,covmat,network,value_b)
            check=0
            do j=1,ns
            if (alpha_n(j) .EQ. alpha(j)) then
            check=check+1
            end if
            end do
            if (check .EQ. ns) then
            value_ini_b= value_ini
            else
            value_ini_b=value
            end if
        CASE (2)
            call likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark,covmat,network,value_b)
            check=0
            do j=1,ns
            if (alpha_n(j) .EQ. alpha(j)) then
            check=check+1
            end if
            end do
            if (check .EQ. ns) then
            value_ini_b= value_ini
            else
            value_ini_b=value
            end if

        END SELECT

            !############# prior distribution selection #############

            SELECT CASE (bnum)
            CASE (1)
            value_fb= gamma_density(beta_n(i),gshapeb(i),gscaleb(i))
            value_fb_ini= gamma_density(beta(i),gshapeb(i),gscaleb(i))
            !calcualte AP
            ratio_b= (value_b-value_ini_b)+(value_fb-value_fb_ini)
            psib= min(1.0d0,exp(ratio_b))

            CASE (2)
            value_fb= half_normal(beta_n(i),halfvarb(i))
            value_fb_ini= half_normal(beta(i),halfvarb(i))
            !calcualte AP
            ratio_b= (value_b-value_ini_b)+(value_fb-value_fb_ini)
            psib= min(1.0d0,exp(ratio_b))

            CASE (3)
            if ((unifminb(i) .LT. beta_n(i)).AND. (beta_n(i).LT. unifmaxb(i))) then
            !calcualte AP
            ratio_b= value_b-value_ini_b
            psib= min(1.0d0,exp(ratio_b))
            else
            psib=-1.0d0
            end if
            END SELECT

            !#############  Acceptance probbaility check for beta #############

            call random_number(u)
            if (psib >=u) then
            simbeta(m+1,i)=beta_n(i)
            else
            simbeta(m+1,i)=beta(i)
            end if

            else
            simbeta(m+1,i)=beta(i)
            end if

            beta_n(i)= simbeta(m+1,i)

!############# end do loop for (ni): checking each beta value #############
        end do


    !#############  spark parameter update #############
    if (flag .EQ. 1) then

    z_sp= rand_normal(0.0d0,prostdsp) !generate from proposal
    spark_n = spark + z_sp

    if (spark_n > 0.0d0) then

    SELECT CASE (tnum)
    CASE (1)
    call likecon(tau,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark_n,covmat,network,value_sp)
    check=0
    do i=1,ni
    if (beta_n(i) .EQ. beta(i)) then
    check=check+1
    end if
    end do
    if (check .EQ. ni) then
    value_ini_sp= value_ini_b
    else
    value_ini_sp=value_b
    end if

    CASE (2)
    call likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark_n,covmat,network,value_sp)
    check=0
    do i=1,ni
    if (beta_n(i) .EQ. beta(i)) then
    check=check+1
    end if
    end do
    if (check .EQ. ni) then
    value_ini_sp= value_ini_b
    else
    value_ini_sp=value_b
    end if
    END SELECT
    !############# prior distribution selection #############

    SELECT CASE (snum)
    CASE (1)
    value_fsp= gamma_density(spark_n,gshapesp,gscalesp)
    value_fsp_ini= gamma_density(spark,gshapesp,gscalesp)
    !calcualte AP
    ratio_sp= (value_sp-value_ini_sp)+(value_fsp-value_fsp_ini)
    psisp= min(1.0d0,exp(ratio_sp))

    CASE (2)
    value_fsp= half_normal(spark_n,halfvarsp)
    value_fsp_ini= half_normal(spark,halfvarsp)
    !calcualte AP
    ratio_sp= (value_sp-value_ini_sp)+(value_fsp-value_fsp_ini)
    psisp= min(1.0d0,exp(ratio_sp))

    CASE (3)
    if ((unifminsp .LT. spark_n).AND. (spark_n .LT. unifmaxsp)) then
    !calcualte AP
    ratio_sp= value_sp-value_ini_sp
    psisp= min(1.0d0,exp(ratio_sp))
    else
    psisp=-1.0d0
    end if
    END SELECT

    !#############  Acceptance probbaility check for spark #############

    call random_number(u)
    if (psisp >=u) then
    simspark(m+1)=spark_n
    else
    simspark(m+1)=spark
    end if

    else
    simspark(m+1)=spark
    end if

    spark_n= simspark(m+1)
    else
    spark=0.0d0
    spark_n=0.0d0
    end if


    !############# likelihood for final values  #############
    SELECT CASE (tnum)
    CASE (1)
    call likecon(tau,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark_n,covmat,network,value_ini)
     llikeval(m+1)= value_ini
    CASE (2)
    call likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha_n,beta_n,spark_n,covmat,network,value_ini)
     llikeval(m+1)= value_ini
    END SELECT

    do i=1,ns
    alpha(i)=alpha_n(i)
    end do
    do i=1,ni
    beta(i)=beta_n(i)
    end do
    spark= spark_n



    !############# simulation ends here #############
    end do
    end subroutine conmcmc

!############################## Random seed initialization ##############################

    subroutine initrandomseed()
    !The seed for the random number generation method random_number() has been reset
    implicit none

    integer :: i,n,clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)

    deallocate(seed)
    end subroutine initrandomseed

!#################### Prior: HALFNORMAL distribution ######################

    FUNCTION half_normal(alpha,b) RESULT(pdf)
    implicit none
    double precision, parameter:: pi = 3.141592653589793D+00
    double precision :: alpha,b
    double precision :: val,pdf

    val=  sqrt (2 / pi*b) * exp ( - 0.5D+00 * ((alpha)**2 / b) )
    pdf=log(val)

    END FUNCTION  ! returns log value

!#################### Prior: GAMMA distribution ######################

    FUNCTION gamma_density(x,a,b) RESULT(pdf)
    implicit none
    double precision :: x,a,b
    double precision :: dn,pdf

    dn = (x**(a-1)) * exp(- (x*b))
    pdf= log(dn)

    END FUNCTION ! returns log value

!#################### Proposal: NORMAL distribution ######################


    FUNCTION rand_normal(mean,stdev) RESULT(c)
    implicit none
    double precision:: mean,stdev,c,temp(2),r,theta
    double precision, parameter:: pi = 3.141592653589793D+00

    call random_number(temp)
    r=(-2.0d0*log(temp(1)))**0.5
    theta = 2.0d0*pi*temp(2)
    c= mean+stdev*r*sin(theta)

    END FUNCTION


!############################## Likelihood calculation

    subroutine like(x,y,tau,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,val)
    implicit none

    !Declarations
    integer, intent(in) :: n,tau(n),tmax,ns,ni,tmin
    double precision, intent(in) :: alpha(ns), beta(ni),spark
    double precision, intent(in) :: x(n), y(n)
    double precision, intent(in) ::covmat(n,ns)
    double precision, intent(out):: val

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
    p1=1.0d0 -exp(-(Somega(i)*dx1+spark))
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

    subroutine likesir(x,y,tau,lambda,n,tmin,tmax,ns,ni,alpha,beta,spark,covmat,val)

    implicit none

    !Declarations
    integer, intent(in) :: n,tau(n),lambda(n),tmax,ns,ni,tmin
    double precision, intent(in) :: alpha(ns), beta(ni),spark
    double precision, intent(in) :: x(n), y(n)
    double precision, intent(in) ::covmat(n,ns)
    double precision, intent(out):: val

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
        p1=1.0d0 -exp(-(Somega(i)*dx1+spark))
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
!####################

    subroutine likecon(tau,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,val)


    implicit none

    !Declarations
    integer, intent(in) :: n,tmax,ns,ni,tmin
    integer, intent(in) :: tau(n)
    double precision, intent(in) :: alpha(ns), beta(ni),spark
    double precision, intent(in) :: covmat(n,ns)
    double precision, intent(in) :: network(n,n,ni)
    double precision, intent(out) ::  val

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

    subroutine likeconsir(tau,lambda,n,ns,ni,tmin,tmax,alpha,beta,spark,covmat,network,val)
    implicit none


    !Declarations
    integer, intent(in) :: n,tmax,ns,ni,tmin
    integer, intent(in) :: tau(n),lambda(n)
    double precision, intent(in) :: alpha(ns), beta(ni),spark
    double precision, intent(in) ::  covmat(n,ns)
    double precision, intent(in) :: network(n,n,ni)
    double precision, intent(out) ::  val

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
    p1= 1.0d0 - exp(-(somega(i)*dx1+spark))
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


!#####################################

    end module mcmcdata

!##################################### module ends ##############################

