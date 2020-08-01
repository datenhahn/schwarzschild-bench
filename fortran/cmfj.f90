module m_cmfj
    implicit none


    integer,parameter :: dp = selected_real_kind(15,307) !  kind(1.d0)
    real(dp),parameter :: const_g       = 9.81    ! graviational const [m/s**2]
    real(dp),parameter :: const_cp      = 1.006e3 ! specific heat capac, dry air [J/kg/K]
    real(dp),parameter :: base_pressure = 1013e2  ! base pressure at surface [Pa]

    type t_atmosphere
        integer :: ke, ke1 ! number of layers/ interfaces
        double precision,allocatable,dimension(:) :: T ! Temperature in layers [K]
        double precision,allocatable,dimension(:) :: pm ! pressure in layers [Pa]
        double precision,allocatable,dimension(:) :: pt ! pressure on interfaces [Pa]
    end type

contains

    elemental function T_2_potentialTemp(T,p,p0)
        double precision,intent(in) :: T,p,p0
        double precision :: T_2_potentialTemp
        double precision,parameter :: kappa = 2./7.

        T_2_potentialTemp = (p0/p)**kappa * T
    end function

    elemental function potentialTemp_2_T(theta,p,p0)
        double precision,intent(in) :: theta,p,p0
        double precision :: potentialTemp_2_T
        double precision,parameter :: kappa = 2./7.

        potentialTemp_2_T = theta / (p0/p)**kappa 
    end function

    function derivative(arr)
        double precision,intent(in) :: arr(:)
        double precision :: derivative(size(arr)-1)

        derivative = arr(1:size(arr)-1) - arr(2:size(arr))
    end function

    function mean_profile(arr)
        double precision,intent(in) :: arr(:)
        double precision :: mean_profile(size(arr)-1)

        mean_profile = (arr(1:size(arr)-1) + arr(2:size(arr)))/2
    end function


    function init_atmosphere(Nlayer) result(atm)
        integer,intent(in) :: Nlayer
        type(t_atmosphere) :: atm

        integer :: k

        atm%ke  = Nlayer
        atm%ke1 = atm%ke+1
        allocate(atm%T (atm%ke) )
        allocate(atm%pm(atm%ke) )
        allocate(atm%pt(atm%ke1))

        atm%pt = (/ ( base_pressure/atm%ke *(k-1) , k=1,atm%ke1) /)
        atm%pm = mean_profile(atm%pt)

        atm%T = 288
        !call print_atm(atm)
    end function

    subroutine driver(atm)
        type(t_atmosphere),intent(inout) :: atm
        double precision :: divE(atm%ke) ! absorbed radiative energy
        double precision :: Tinc(atm%ke) ! Temperature tendency [K/s]

        double precision :: time, dt ! adaptive timestep [s]

        time = 0
        timeloop: do 
            call radiative_transfer(atm, divE)
            !divE = 0
            !divE(atm%ke) = 235-5.67e-8*atm%T(atm%ke)**4

            Tinc = divE *const_g / const_cp / abs( derivative(atm%pt) )

            dt = min(3600._dp, 3_dp / max(epsilon(Tinc), maxval(Tinc))) ! restrict timestep to max. 3K step
            atm%T = atm%T + Tinc*dt

            time = time+dt
            !        print *,
            !        print *,time/3600/24,'[d]    dt',dt,'[s]'
            !        print *,

            call convection(atm)

            if(check_equilibrium(atm)) exit timeloop
        enddo timeloop
        call print_atm(atm)
        print *,'model time:',time/3600/24,'[days]'
    end subroutine

    subroutine radiative_transfer(atm, divE)
        type(t_atmosphere),intent(in) :: atm
        double precision,intent(out) :: divE(:)

        real(dp),dimension(atm%ke)  :: dtau
        real(dp),dimension(atm%ke) :: planck
        real(dp),dimension(atm%ke)  :: Ediv_thermal, Ediv_solar
        real(dp),dimension(atm%ke1) :: Edn, Eup
        real(dp),dimension(atm%ke1) :: Edn_int, Eup_int
        real(dp),dimension(atm%ke)  :: Ediv_thermal_int

        integer :: k
        logical, parameter :: lgray=.True.

        Ediv_solar = 0
        Ediv_solar(atm%ke) = 235

        if(lgray) then
            dtau = 1._dp/atm%ke 
            !planck = planck_int(atm%T,1e-6_dp,100e-6_dp)
            planck = 5.67e-8*atm%T**4 / 3.1415

            call schwarzschild(dtau, 0._dp, planck, Edn_int, Eup_int, Ediv_thermal_int)
        else
            Eup_int = 0
            Edn_int = 0
            Ediv_thermal_int = 0

            dtau = 10._dp/atm%ke 
            planck = planck_int(atm%T,1e-6_dp,8e-6_dp)
            call schwarzschild(dtau, 0._dp, planck, Edn, Eup, Ediv_thermal)

            Eup_int = Eup_int + Eup
            Edn_int = Edn_int + Edn
            Ediv_thermal_int = Ediv_thermal_int + Ediv_thermal

            dtau = .5_dp/atm%ke 
            planck = planck_int(atm%T,8e-6_dp,12e-6_dp)
            call schwarzschild(dtau, 0._dp, planck, Edn, Eup, Ediv_thermal)

            Eup_int = Eup_int + Eup
            Edn_int = Edn_int + Edn
            Ediv_thermal_int = Ediv_thermal_int + Ediv_thermal

            dtau = 1._dp/atm%ke 
            planck = planck_int(atm%T,12e-6_dp,100e-6_dp)
            call schwarzschild(dtau, 0._dp, planck, Edn, Eup, Ediv_thermal)

            Eup_int = Eup_int + Eup
            Edn_int = Edn_int + Edn
            Ediv_thermal_int = Ediv_thermal_int + Ediv_thermal

        endif

        divE = Ediv_solar + Ediv_thermal_int

    contains

        elemental function planck_int(T, l1, l2)
            real(dp), intent(in) :: T, l1, l2
            real(dp) :: planck_int
            real(dp),parameter :: delta_lambda = 1e-6_dp
            integer :: i,N 
            real(dp) :: l,dl 

            !if(l2.gt.l1) stop 'You may not call planck integration with l2.gt.l1'

            N = (l2-l1)/delta_lambda+.5
            dl = (l2-l1)/N

            planck_int = 0
            do i=1,N
                l = l1 + (i-.5)*dl
                planck_int = planck_int + plk(T,l)*dl
            enddo
        end function

        elemental function plk(T,l)
            real(dp),intent(in) :: T,l
            real(dp) :: plk
            real(dp),parameter :: c=299792458, k_B=1.3806e-23_dp, h=6.626e-34_dp
            plk = 2*h*c**2 / l**5 / (exp(h*c/(l*k_B*T))-1);
        end function

        subroutine schwarzschild(dtau, albedo, planck, Edn, Eup, Ediv)
            real(dp),intent(in),dimension(:) :: dtau
            real(dp),intent(in) :: albedo
            real(dp),dimension(:),intent(in) :: planck
            real(dp),dimension(:),intent(out):: Edn, Eup, Ediv

            integer :: imu,k,ke,ke1

            integer,parameter :: Nmu = 5
            real(dp),parameter  :: one=1, zero=0, dmu = one/Nmu, pi=3.1415


            real(dp) :: T(size(dtau)) ! Transmission coefficients
            real(dp) :: Lup, Ldn
            real(dp) :: mu

            Edn=zero
            Eup=zero

            ke = size(dtau)
            ke1 = ke+1

            ! Transmission coefficients
            do imu=1,Nmu
                mu = (imu-.5_dp)*dmu  
                T = exp(- dtau/mu)

                ! Boundary conditions at surface
                Lup = planck(ke)
                Eup(ke1) = Eup(ke1) + Lup*mu

                ! zero incoming radiation at TOA
                Ldn = zero
                Edn(1) = Edn(1) + Ldn*mu

                do k=ke,1,-1
                    Lup = Lup * T(k) + planck(k)*(one-T(k))
                    Eup(k) = Eup(k) + Lup*mu
                enddo
                do k=1,ke
                    Ldn = Ldn * T(k) + planck(k)*(one-T(k))
                    Edn(k+1) = Edn(k+1) + Ldn*mu
                enddo

            enddo

            Eup = Eup*2*pi*dmu
            Edn = Edn*2*pi*dmu

            Ediv = Eup(2:ke1) + Edn(1:ke) - Eup(1:ke) - Edn(2:ke1)
            Ediv(ke) = Ediv(ke) + Edn(ke1) - Eup(ke1)

        end subroutine

    end subroutine


    subroutine convection(atm)
        type(t_atmosphere),intent(inout) :: atm

        logical,parameter :: lconvection_mixing = .False.
        double precision :: theta(atm%ke)

        if(lconvection_mixing) then

            call mixing(atm)
        else
            theta = T_2_potentialTemp(atm%T, atm%pm, base_pressure)
            call bubble_sort(theta)
            atm%T = potentialTemp_2_T(theta, atm%pm, base_pressure)
        endif

    contains
        subroutine mixing(atm)
            type(t_atmosphere),intent(inout) :: atm

            double precision :: theta(atm%ke-1), dT
            integer :: iter,k
            integer,parameter :: maxiter=1e6
            logical :: lmixed

            theta = T_2_potentialTemp(atm%T(1:atm%ke-1), atm%pm(1:atm%ke-1), atm%pm(2:atm%ke) ) ! Temperature when brought one layer down

            repeat: do iter=1,maxiter
                lmixed=.False.

                do k=atm%ke-1,1,-1
                    dT = atm%T(k+1) - theta(k)
                    if(dT.gt.1e-4) then ! convection if lower temp is smaller than potT of layer above
                        lmixed=.True.
                        atm%T(k+1) = atm%T(k+1)-dT/2
                        atm%T(k)   = atm%T(k)  +dT/2
                    endif
                enddo

                if(.not.lmixed) exit repeat
            enddo repeat
            print *,'number of convective mixing:',iter,'imix',lmixed

        end subroutine

        subroutine bubble_sort(arr)
            double precision, intent(inout) :: arr(:)
            double precision :: temp
            integer :: i, j
            logical :: swapped

            do j = size(arr)-1, 1, -1
                swapped = .FALSE.
                do i = 1, j
                    if (arr(i) < arr(i+1)) then
                        temp = arr(i)
                        arr(i) = arr(i+1)
                        arr(i+1) = temp
                        swapped = .TRUE.
                    endif 
                enddo 
                if (.not. swapped) exit
            enddo 
        end subroutine            
    end subroutine

    function check_equilibrium(atm) result(equilibrium)
        type(t_atmosphere) :: atm
        type(t_atmosphere),allocatable,save :: last_atm
        logical :: equilibrium
        real(dp) :: residual

        equilibrium = .False.

        if(.not.allocated(last_atm)) then
            allocate(last_atm, source=atm)
            return
        endif

        residual = rmse(last_atm%T, atm%T)
        !            print *,'Temperature residual',residual

        ! check if temperature profile changed in last layer:
        if (residual.lt.1e-4_dp) equilibrium = .True.


        ! overwrite last_atm with current one
        last_atm = atm

    contains
        double precision function rmse(a,b)
            double precision,intent(in) :: a(:),b(:)
            rmse = sqrt( sum((a-b)**2) / size(a) )
        end function
    end function



    subroutine print_atm(atm)
        type(t_atmosphere),intent(in) :: atm
        integer :: k
        double precision :: theta(atm%ke)

        theta = T_2_potentialTemp(atm%T, atm%pm, base_pressure)

        do k=1,atm%ke
            print *,k,'pm',atm%pm(k)/1e2,'[hPa] ::  T',atm%T(k),'[K] :: theta',theta(k) !,T2(k)
        enddo
    end subroutine

end module

program main
    use m_cmfj
    type(t_atmosphere) :: atm

    atm = init_atmosphere(100)
    call driver(atm)

end program
