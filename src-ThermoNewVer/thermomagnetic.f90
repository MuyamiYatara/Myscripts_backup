!---------------------------------------------------------!
!> Calculate electrical and thermal transport with        !
!> R.G.Chambers's formula based on Boltzmann transport    !
!> The main subroutine and electrical part is written by  !
!> QuanSheng Wu (wuquansheng@gmail.com), and the thermal  !
!> part is added by Hanqi Pi (hqpi1999@gmail.com)         !       
!> Thanks Yi Liu for helpful discussions                  !
!>                                                        !
!> References :                                           !
!> [1] Electrons in metals and semiconductors,            !
!>     R.G. Chambers,                                     !
!> [2] Ab initio investigation of magnetic transport      !
!>     properties by Wannier interpolation, PHYSICAL      !
!>     REVIEW B 79, 245123 (2009), Yi Liu, Hai-Jun Zhang, !
!>     and Yugui Yao                                      !
!> [3] Magnetoresistance from Fermi surface topology,     !
!>     ShengNan Zhang, QuanSheng Wu, Yi Liu, and          !
!>     Oleg V. Yazyev,Phys. Rev. B 99, 035142 (2019)      !
!>                                                        !
!> Implemented on Oct. 07 2017                            !
!> uploaded on Sep. 05. 2019                              !
!---------------------------------------------------------!

subroutine thermomag
!-----------------------------------------------------------!
!> This is the main routine to calculate the electrical and !
!> thermal transport coefficients in the presence/absence of!
!> magnetic field based on Boltzmann transport equation.    !
!>                                                          !
!> We calculate the conductivity and resistivity under      !
!> the band-resolved constant relaxation time.              !
!>                                                          !
!> This version impose symmetry constrain to get the        !
!> transport coefficients under magnetic field.             !
!>                                                          !
!> It produces 4 types of output files,                     !
!> every type has OmegaNum files at different chemical      !
!> potentials
!>                                                          ! 
!> 1. rhotau_bands_mu_***eV.dat :                           !
!>                                                          !
!> The first two files are output on a grid of              !
!> (Nband_Fermi_Level, NumT, NBTau), while the last two     !
!> files are output on a grid of (NumT, NBTau). NumT is the !
!> number of temperature points, NBTau is the number of Btau! 
!> and Nband_Fermi_Level is the number of bands.            !
!>                                                          !
!> The output rho*tau is in the unit of Ohm*m*s.            !
!-----------------------------------------------------------!
   use wmpi
   use para
   implicit none

   integer :: iLn, icoord, iBTau, imu, iT, iband, i, j, ikt, ie,ierr

   !> Bands crossing Fermi level
   integer :: Nband_Fermi_Level
   integer, allocatable :: bands_fermi_level_temp(:)
   integer, allocatable :: bands_fermi_level(:)

   real(dp) :: KBT 
   ! K_Boltzmann*Temperature in eV
   ! real(dp) :: mu ! chemical potential relative to Fermi level in eV

   !> In this method, we can't treat the magnetic field and relaxation time individually
   !> They always come together as Btau. Omega= eB/m*
   !> BTau is in units of Tesla*ps where ps means 10^-12 second.
   !> BTau is in order of 1
   !> the relaxation time tau is in order of 1
   !> For Si, at zero temperature, tau=1ps
   !> For Ge, tau= 0.26ps
   !> For GaAs, tau= 0.48ps
   !> For InAs, tau= 0.08ps 
   !> reference  http://www.iue.tuwien.ac.at/phd/palankovski/node51.html
   real(dp), allocatable :: BTau_array(:), mu_array(:), KBT_array(:)

   real(dp) :: time_start, time_end

   !> Ln_tau_bands(3, coordinates, Btau, mu, KBT, Nbands)
   real(dp), allocatable :: Ln_tau_bands(:, :, :, :, :, :)
   !> L(0)/tau=Ln_tau(1,:,:,:,:),  L(1)/tau=Ln_tau(2,:,:,:,:), L(2)/tau=Ln_tau(3,:,:,:,:)
   real(dp), allocatable :: Ln_tau(:,:,:,:,:)

   !> Lxx_tau (coordinate, coordinate)
   !> Lxx_tau = Lxx/tau
   real(dp), allocatable :: L11_tau(:, :) !> L11_tau=L(0)/tau
   real(dp), allocatable :: L11_tau_inv(:, :) !> L11_tau_inv=(L11_tau)^-1
   real(dp), allocatable :: L12_tau(:, :) !> L12_tau=L(1)/tau/electron/Temperature
   real(dp), allocatable :: L21_tau(:, :) !> L21_tau=L(1)/tau/electron
   real(dp), allocatable :: L22_tau(:, :) !> L22_tau=L(2)/tau/electron**2/Temperature
  
   real(dp), allocatable :: numerator_local(:, :) 
   real(dp), allocatable :: denominator_local(:, :) 

   !@@@@@ velocity test and Umkehr effect
   !integer :: velocity_bar_fileindex

   real(dp), allocatable :: rhotau(:,:), seebeck(:,:), etatau(:,:), thermalhall(:,:)

   !> file name
   character(50) :: Lnfname, muname, rhotaufname, seebeckfname, etataufname, thermalhallfname

   !> file index
   integer :: rhotaufindex, seebeckfindex, etataufindex, thermalhallfindex

   real(dp), external :: det3
 
   !> Nband_Fermi_Level and bands_fermi_level will be updated after
   !> get_bands_cross_fermilevel, we hope 1000 is quite enough
   Nband_Fermi_Level= 1000
   allocate(bands_fermi_level_temp(Nband_Fermi_Level))
   ! allocate(myfileindex(Nband_Fermi_Level))

   if (NumberofSelectedBands/=0) then
      Nband_Fermi_Level= NumberofSelectedBands
      allocate(bands_fermi_level(Nband_Fermi_Level))
      do i=1, Nband_Fermi_Level
         bands_fermi_level(i)= Selected_band_index(i) 
      enddo
   else
      !> First we calculate how many and which bands cross the Fermi level
      call get_bands_cross_fermilevel(Nband_Fermi_Level, bands_fermi_level_temp)
      allocate(bands_fermi_level(Nband_Fermi_Level))
      bands_fermi_level= bands_fermi_level_temp(1:Nband_Fermi_Level)
   endif


   !> set Btau array
   allocate(BTau_array(NBTau))
   BTau_array= 0d0
   if (NBTau>1) then
      do i=1, NBTau
         BTau_array(i)= (i-1.0d0)/(NBTau-1)*BTauMax
      enddo
   else
      BTau_array = BTauMax
   endif

   !> set chemical potential array
   allocate(mu_array(OmegaNum))
   mu_array= 0d0
   if (OmegaNum>1) then
      do i=1, OmegaNum
         mu_array(i)= OmegaMin+ (i-1.0d0)/(OmegaNum-1.0d0)*(OmegaMax-OmegaMin)
      enddo
   else
      mu_array = OmegaMin
   endif
   !write(*,*) mu_array
   !> set temperature array
   allocate(KBT_array(NumT))
   KBT_array= 0d0
   if (NumT>1) then
      do i=1, NumT
         KBT_array(i)= Tmin+ (i-1.0d0)/(NumT-1.0d0)*(Tmax-Tmin)
      enddo
   else
      KBT_array = Tmin
   endif

   if (cpuid.eq.0) then
      write(stdout, *) ' '
      write(stdout, *)' KBT array in the calculation in unit of Kelvin'
      write(stdout, '(10f8.2)') KBT_array
      write(stdout, *) ' '
   endif

   !> transform from Kelvin to eV
   !> The SI unit of temperature is the kelvin (K), but using the above relation the electron temperature is often expressed in
   !> terms of the energy unit electronvolt (eV). Each kelvin (1 K) corresponds to 8.6173324(78)×10−5 eV; this factor is the ratio
   !> of the Boltzmann constant to the elementary charge. After version 2.6, we 
   !> adopt the atomic unit
   KBT_array= KBT_array*8.6173324E-5*eV2Hartree

   !>> calculate the band resolved conductivity tensor
   !> The tensor is like
   !> xx  xy  xz
   !> yx  yy  yz
   !> zx  zy  zz
   !> tensor1=xx, tensor2=xy, tensor3=xz
   !> tensor4=yx, tensor5=yy, tensor6=yz
   !> tensor7=zx, tensor8=zy, tensor9=zz

   allocate(Ln_tau_bands(3, 9, NBTau, OmegaNum, NumT, Nband_Fermi_Level))
   allocate(Ln_tau(3, 9, NBTau, OmegaNum, NumT))
   allocate(L11_tau(3, 3))
   allocate(L11_tau_inv(3, 3))
   allocate(L12_tau(3, 3))
   allocate(L21_tau(3, 3))
   allocate(L22_tau(3, 3))
   allocate(rhotau(3, 3))
   allocate(seebeck(3, 3))
   allocate(etatau(3, 3))
   allocate(thermalhall(3, 3))
   Ln_tau_bands = 0d0
   Ln_tau  = 0d0
   L11_tau = 0d0
   L11_tau_inv = 0d0
   L12_tau = 0d0
   L21_tau = 0d0
   L22_tau = 0d0
   rhotau = 0d0
   seebeck = 0d0
   etatau = 0d0
   thermalhall = 0d0

   time_start= 0d0
   time_end= 0d0
   !if ( cpuid .eq. 0 ) then
   !   !@@@@@
   !   velocity_bar_fileindex = outfileindex+1
   !   open(unit=velocity_bar_fileindex, file="velocity_bar_compare")
   !   outfileindex = outfileindex+2
   !end if
 
   !> We obtain Ln/tau from 'Ln_calc'
   call Ln_tau_calc_symm(mu_array, KBT_array, BTau_array, Nband_Fermi_Level, bands_fermi_level, Ln_tau_bands)
   
   !> We sum over all bands to get the total Ln/tau
   do iLn = 1, 3
      do icoord = 1, 9
         do iBTau = 1, NBTau
            do imu = 1, OmegaNum
               do iT = 1, NumT
                  Ln_tau(iLn, icoord, iBTau, imu, iT) = sum(Ln_tau_bands(iLn, icoord, iBTau, imu, iT, :))
               enddo
            enddo
         enddo
      enddo 
   enddo

   if(cpuid.eq.0)write(stdout, *)' '
   if(cpuid.eq.0)write(stdout, *)'>> Start to calculate the transport properties in magnetic field'

   if (cpuid.eq.0) then
      do ie=1, OmegaNum
         write(muname, '(f12.3)')mu_array(ie)/eV2Hartree
         
         rhotaufindex      = outfileindex+1
         seebeckfindex     = outfileindex+2
         etataufindex      = outfileindex+3
         thermalhallfindex = outfileindex+4
         outfileindex      = outfileindex+4


      
         !> write rhotau into file
         write(rhotaufname, '(6a)')'rhotau_total','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=rhotaufindex, file=rhotaufname)
         write(rhotaufindex, '(a20, (a, f16.4, a))')'# rhotau tensor ', &
            ' with chemical potential set to ', mu_array(ie)/eV2Hartree, ' eV'

         !> write seebeck coefficient into file
         write(seebeckfname, '(6a)')'seebeck_total','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=seebeckfindex, file=seebeckfname)
         write(seebeckfindex, '(a20, (a, f16.4, a))')'# seebeck tensor ', &
            ' with chemical potential set to ', mu_array(ie)/eV2Hartree, ' eV'

         !> write Ettinghausen coefficient into file
         write(etataufname, '(6a)')'etatau_total','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=etataufindex, file=etataufname)
         write(etataufindex, '(a20, (a, f16.4, a))')'# Ettinghausen tensor ', &
               ' with chemical potential set to ', mu_array(ie)/eV2Hartree, ' eV'

         !> write Thermal Hall coefficient into file
         write(thermalhallfname, '(6a)')'thermalhall_total','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=thermalhallfindex, file=thermalhallfname)
         write(thermalhallfindex, '(a20, (a, f16.4, a))')'# ThermalHall tensor ', &
            ' with chemical potential set to ', mu_array(ie)/eV2Hartree, ' eV'
         
         do ikt=1, NumT
            KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree

            write(rhotaufindex, '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
            write(seebeckfindex, '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
            write(etataufindex, '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
            write(thermalhallfindex, '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'

            write(rhotaufindex, '("#",11a16)')'BTau (T.ps)', 'OmegaTau (eV.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
            write(seebeckfindex, '("#",11a16)')'BTau (T.ps)', 'OmegaTau (eV.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
            write(etataufindex, '("#",11a16)')'BTau (T.ps)', 'OmegaTau (eV.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
            write(thermalhallfindex, '("#",11a16)')'BTau (T.ps)', 'OmegaTau (eV.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
            
            do i=1, NBTau
               !> L11_tau=L(0)/tau=Ln_tau(1,:,:,:,:)
               !> It has the units of CV^-1m^-1s^-2
               L11_tau(1, 1:3)=Ln_tau(1, 1:3, i, ie, ikt)
               L11_tau(2, 1:3)=Ln_tau(1, 4:6, i, ie, ikt)
               L11_tau(3, 1:3)=Ln_tau(1, 7:9, i, ie, ikt)
               !> L11_tau_inv=1/L11_tau
               !> It has the units of C^-1Vms^2
               L11_tau_inv = L11_tau
               !> L12_tau=L(1)/tau/electron/Temperature=Ln_tau(2,:,:,:,:)/electron/Temperature
               !> It has the units of Cm^-1s^-2K^-1
               L12_tau(1, 1:3)=-Ln_tau(2, 1:3, i, ie, ikt)/Echarge/KBT
               L12_tau(2, 1:3)=-Ln_tau(2, 4:6, i, ie, ikt)/Echarge/KBT
               L12_tau(3, 1:3)=-Ln_tau(2, 7:9, i, ie, ikt)/Echarge/KBT
               !> L21_tau=L(1)/tau/electron=Ln_tau(2,:,:,:,:)/electron
               !> It has the units of Cm^-1s^-2
               L21_tau(1, 1:3)=-Ln_tau(2, 1:3, i, ie, ikt)/Echarge
               L21_tau(2, 1:3)=-Ln_tau(2, 4:6, i, ie, ikt)/Echarge
               L21_tau(3, 1:3)=-Ln_tau(2, 7:9, i, ie, ikt)/Echarge
               !> L22_tau=L(2)/tau/electron**2/Temperature=Ln_tau(3,:,:,:,:)/electron**2/Temperature
               !> It has the units of CVm^-1s^-2K^-1
               L22_tau(1, 1:3)=Ln_tau(3, 1:3, i, ie, ikt)/Echarge**2/KBT
               L22_tau(2, 1:3)=Ln_tau(3, 4:6, i, ie, ikt)/Echarge**2/KBT
               L22_tau(3, 1:3)=Ln_tau(3, 7:9, i, ie, ikt)/Echarge**2/KBT
               if (abs(det3(L11_tau))>eps6) then
                  call inv_r(3, L11_tau_inv)
                  !> It has the units of Vms^2/C = Omega*m*s
                  rhotau = L11_tau_inv
                  !> write rhotau into files
                  write(rhotaufindex, '(11E16.6)')&
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau*0.175874356d0, &
                     rhotau(1,:), rhotau(2,:), rhotau(3,:)

                  !> It has the units of V/K
                  seebeck = matmul(L11_tau_inv, L12_tau)
                  !> write seebeck into files
                  write(seebeckfindex, '(11E16.6)')&
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau*0.175874356d0, &
                  seebeck(1,:), seebeck(2,:), seebeck(3,:)

                  !> to obtain Ettinghausen coefficient
                  numerator_local = matmul(L21_tau, L11_tau_inv)
                  denominator_local = L22_tau - matmul(L21_tau, matmul(L11_tau_inv, L12_tau))
                  do j = 1, 3
                    if (abs(denominator_local(j, j))>eps6) then
                        !> It has the units of  msK/A
                        etatau(j, 1:3) = numerator_local(j, 1:3)/denominator_local(j, j)

                        !> to obtain Thermal Hall coefficient
                        thermalhall(j, 1:3) = -denominator_local(j, 1:3)/denominator_local(j, j)
                    else
                        write(etataufindex, '(a)')'# error: etatau is zero since no k points contribute to the calculations of MR'
                        write(thermalhallfindex, '(a)')'# error: thermalhall is zero since no k points contribute to the calculations of MR'
                    endif
                  enddo ! j
                  !> write etatau into files
                  write(etataufindex, '(11E16.6)')&
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau*0.175874356d0, &
                  etatau(1,:), etatau(2,:), etatau(3,:)
                  
                  !> write thermalhall into files 
                  write(thermalhallfindex, '(11E16.6)')&
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                  BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau*0.175874356d0, &
                  thermalhall(1,:),thermalhall(2,:), thermalhall(3,:)
               else
                  write(rhotaufindex, '(a)')'# error: rhotau is zero since no k points contribute to the calculations of MR'
                  write(seebeckfindex,'(a)')'# error: seebeck is zero since no k points contribute to the calculations of MR'
                  write(etataufindex, '(a)')'# error: etatau is zero since no k points contribute to the calculations of MR'
                  write(thermalhallfindex, '(a)')'# error: thermalhall is zero since no k points contribute to the calculations of MR'
               endif
            enddo ! i, NBTau

            write(rhotaufindex,'(a)') ''
            write(seebeckfindex,'(a)') ''
            write(etataufindex,'(a)') ''
            write(thermalhallfindex,'(a)') ''

         enddo  ! ikt=1, numT
         close(rhotaufindex)
         close(seebeckfindex)
         close(etataufindex)
         close(thermalhallfindex)
      enddo  ! ie=1, OmegaNum
   endif ! cpuid=0

   if (cpuid.eq.0) write(stdout, '(a)') ' <<< The thermomagnetic calculation finished'
   if (cpuid.eq.0) write(stdout, '(a)') ' '

   return
end subroutine thermomag


subroutine Ln_tau_calc_symm(mu_array, KBT_array, BTau_array, Nband_Fermi_Level, bands_fermi_level, Ln_tau_bands)
!-----------------------------------------------------------!
!> This subroutine aims to calculate Ln/tau and is organized!
!> as followings                                            !
!>                                                          !
!> 1. Distribute the kpoints in the reduced BZ into MPI     !
!>    threads
!> transport coefficients under magnetic field.             !
!>                                                          !
!> It produces 4 types of output files,                     !
!> every type has OmegaNum files at different chemical      !
!>                                                          ! 
!> 1. rhotau_bands_mu_***eV.dat :                           !
!>                                                          !
!> The first two files are output on a grid of              !
!> (Nband_Fermi_Level, NumT, NBTau), while the last two     !
!> files are output on a grid of (NumT, NBTau). NumT is the !
!> number of temperature points, NBTau is the number of Btau! 
!> and Nband_Fermi_Level is the number of bands.            !
!>                                                          !
!> The output rho*tau is in the unit of Ohm*m*s.            !
!-----------------------------------------------------------!
   use wmpi
   use para
   implicit none


   integer, intent(inout) :: Nband_Fermi_Level
   integer, intent(inout) :: bands_fermi_level(Nband_Fermi_Level)
   real(dp), intent(in) :: KBT_array(NumT) ! K_Boltzmann*Temperature in eV
   real(dp), intent(in) :: mu_array(OmegaNum) ! chemical potential relative to Fermi level in eV
   real(dp), intent(in) :: BTau_array(NBTau) ! omega*tau without units
   real(dp), intent(inout) :: Ln_tau_bands(3, 9, NBTau, OmegaNum, NumT, Nband_Fermi_Level) 

   real(dp) :: coeff, mu, BTau, KBT, exponent_max
   integer :: ie, ibtau, ikt, n

   integer :: Nk_total, Nk_current, Nk_start, Nk_end
   integer  :: knv3_left, knv3_left_mod, knv3, knv3_mod
   integer :: ik, ib, ik1, ik2, ik3, ik_temp
   integer :: ierr, it, i, ix, j1, j2, j
   integer :: nrecevs

   real(dp) :: v_t(3), v_k(3)
   real(dp) :: k(3), k_start(3), magnetic_field(3)
   real(dp) :: Ln_tau_symm_t(9)
   
   real(dp) :: time_start, time_end
   real(dp) :: time_start0, time_end0
   integer :: NSlice_Btau_inuse

   !> Btau slices for Runge-Kutta integration
   real(dp) :: Btau_start, Btau_final, DeltaBtau
   real(dp), allocatable :: exp_factors(:)
   logical :: fail
   integer :: icycle


   !> energy bands 
   real(dp) :: EE
   real(dp), allocatable :: Ek(:)  ! in eV
   real(dp), allocatable :: Enk(:, :)

   !> minus fermi derivation
   real(dp) :: minusdfde

   !> 3-component velocity for each band and each k point 
   real(dp), allocatable :: velocity_k(:, :)
   real(dp), allocatable :: velocity_bar_k(:)
   real(dp), allocatable :: velocity(:, :, :)
   
   !> 3-component velocity for each band and each k point 
   !> Eq.(3) in PRB 79, 245123(2009)
   real(dp), allocatable :: velocity_bar(:, :, :)

   type(kcube_type) :: KCube3D_total
   type(kcube_type) :: KCube3D_left(Nband_Fermi_Level)

   !> some arrays for mpi
   integer, allocatable :: info(:, :) !> howmany integers to be sent on each cpu 
   integer, allocatable :: Displs(:, :)

   !> number of steps used in the Runge-Kutta integration
   integer :: NSlice_Btau
   integer :: NSlice_Btau_local
   real(dp), allocatable :: kout(:, :)

   !@@@@@ velocity test and Umkehr effect
   !integer, intent(in) :: velocity_bar_fileindex
   real(dp) :: M_be_inverse(2,3,3)
   real(dp) :: magnetic_field_antisym_tensor(3,3)
   real(dp) :: tau_anistropic_alpha(2,3,3)
   real(dp) :: tau_T_c_60K
   real(dp) :: Lax_gamma_prime(2), Lax_gamma(2)
   real(dp) :: Identity3(3,3) = (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/) 
   real(dp) :: velocity_bar_k_article(3)

   !@@@@@ failed kpoints test file
   integer :: failed_k_fileindex
   character(50) :: failed_k_filename

   
   !> inverse of group operator
   real(dp) :: Tmat(3, 3)
   real(dp), allocatable :: pgop_cart_inverse(:, :, :)

   character(50) :: Lnfname, muname


   !> define some arrays for different bands. Since there are different number
   !> of k points left for different bands.
   type klist_iband_type
      !> dim=3*NSlice_Btau
      real(dp), allocatable :: klist_rkfs(:, :)
      !> dim=3*NSlice_Btau
      real(dp), allocatable :: velocity_k(:, :)

      !> calculate -df(e)/de, where f(e) is the Fermi-Dirac distribution
      !> dim=(NumT, OmegaNum))
      real(dp), allocatable :: minusdfde(:, :)
   end type klist_iband_type
   type(klist_iband_type) :: klist_iband(Nband_Fermi_Level)

   type Ln_tau_iband_type
         ! Ln_tau_k(3, coordinate, NBTau, OmegaNum, NumT) 
         real(dp), allocatable :: Ln_tau_k_mpi(:, :, :, :,  :)
         real(dp), allocatable :: Ln_tau_k(:, :, :, :, :)
         ! Ln_tensor_kz(9, NBTau, OmegaNum, NumT, Nk3) 
         ! real(dp), allocatable :: Ln_tensor_kz(:, :, :, :, :)

         !> time cost for k point
         real(dp), allocatable :: time_cost(:)
         real(dp), allocatable :: time_cost_mpi(:)
   end type Ln_tau_iband_type
   type(Ln_tau_iband_type) :: Ln_tau_iband_k(Nband_Fermi_Level)

   !> default value
   icycle = 1
   fail = .False.

   !> STEP ONE
   !> Distribute k kpoints in the 3DCube into different MPI threads
   !> knv3 and Kcube3D_total%Nk_total are the totalnumber of kpoints in reduced BZ
   allocate(pgop_cart_inverse(3, 3, number_group_operators))
   do i= 1, number_group_operators
      Tmat= pgop_cart(:, :, i) 
      call inv_r(3, Tmat)
      pgop_cart_inverse(:, :, i)= Tmat
   enddo
   
   knv3= KCube3D_symm%Nk_total_symm
   !write(*,*) knv3
   KCube3D_total%Nk_total= knv3
   knv3_mod= mod(knv3, num_cpu)
   !> Nk_current is the number of k points for every MPI thread
   if (knv3_mod==0) then  !> perfect divided
      KCube3D_total%Nk_current= knv3/num_cpu
      KCube3D_total%Nk_start=1+ knv3*cpuid/num_cpu
      KCube3D_total%Nk_end  =(1+cpuid)*knv3/num_cpu
   else if (knv3/num_cpu==0) then    !> Number of MPI threads is large than knv3
      KCube3D_total%Nk_current= 1 !> one k piont per MPI thread
      KCube3D_total%Nk_start= cpuid+ 1 !> one k piont per MPI thread
      KCube3D_total%Nk_end  = cpuid+ 1
      if (cpuid+1 > knv3) then
         KCube3D_total%Nk_start= 1
         KCube3D_total%Nk_end  = 0
      endif
   else
      KCube3D_total%Nk_current= knv3/num_cpu+ 1
      if (cpuid< knv3_mod) then
         KCube3D_total%Nk_start= 1+ cpuid*KCube3D_total%Nk_current
         KCube3D_total%Nk_end  = (1+cpuid)*KCube3D_total%Nk_current
      else
         KCube3D_total%Nk_start= knv3_mod*KCube3D_total%Nk_current+ &
            (cpuid-knv3_mod)*(KCube3D_total%Nk_current-1)+1
         KCube3D_total%Nk_end  = knv3_mod*KCube3D_total%Nk_current+ &
            (cpuid-knv3_mod+1)*(KCube3D_total%Nk_current-1)
      endif
   endif

   !> obtain the direct coordinates and weight of k points at every MPI thread
   allocate(KCube3D_total%k_direct(3, KCube3D_total%Nk_start:KCube3D_total%Nk_end))
   allocate(KCube3D_total%weight_k(KCube3D_total%Nk_start:KCube3D_total%Nk_end))

   do ik= KCube3D_total%Nk_start, KCube3D_total%Nk_end
      !> KCube3D_symm%ik_array_symm(ik) stores the reduced k points
      j1= KCube3D_symm%ik_array_symm(ik)
      ik1= (j1-1)/(Nk2*Nk3)+1
      ik2= ((j1-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
      ik3= (j1-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
      KCube3D_total%k_direct(1, ik)=  (ik1-1d0)/dble(Nk1)  
      KCube3D_total%k_direct(2, ik)=  (ik2-1d0)/dble(Nk2)  
      KCube3D_total%k_direct(3, ik)=  (ik3-1d0)/dble(Nk3)  
      KCube3D_total%weight_k(ik)= KCube3D_symm%weight_k(ik)
   enddo

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   !> STEP TWO
   !> assign values to important variables for every MPI thread

   !> setup NSlice_Btau
   !> NSlice_Btau should be the integer times of NBTau
   if (NBTau>1) then
      NSlice_Btau= (NBTau-1)*(Nslice_BTau_Max/(NBTau-1))
   else
      NSlice_Btau= 1
   endif
   if (cpuid.eq.0) write(stdout, *) ' NSlice_Btau :', NSlice_Btau

   Nk_total= KCube3D_total%Nk_total
   Nk_current= KCube3D_total%Nk_current
   Nk_start= KCube3D_total%Nk_start
   Nk_end= KCube3D_total%Nk_end

   allocate( Ek(Nband_Fermi_Level)) 
   allocate( Enk(Nk_start:Nk_end, Nband_Fermi_Level))
   allocate( velocity(3, Nk_start:Nk_end, Nband_Fermi_Level))
   allocate( velocity_k(3, Nband_Fermi_Level))
   allocate( velocity_bar(3, Nk_start:Nk_end, Nband_Fermi_Level))
   allocate( velocity_bar_k(3))
   Ek= 0d0
   Enk= 0d0
   velocity= 0d0
   velocity_k= 0d0
   velocity_bar= 0d0
   velocity_bar_k= 0d0
   velocity_bar_k_article= 0d0

   !> STEP THREE
   !> calculate the velocity at every k point
   time_start= 0d0
   time_end= 0d0
   do ik= Nk_start, Nk_end
      if (cpuid.eq.0.and. mod(ik, 100).eq.0) &
         write(stdout, '(a, i18, "   /", i18, a, f10.3, "s")') 'ik/NK', &
         ik-Nk_start,Nk_current, 'time left', &
         (Nk_current+Nk_start-ik)*(time_end-time_start)/num_cpu

      call now(time_start)
      k= KCube3D_total%k_direct(:, ik)

      call velocity_calc(Nband_Fermi_Level, bands_fermi_level, k, velocity_k, Ek)
      velocity(:, ik, :)= velocity_k
      Enk(ik, :)= Ek
      call now(time_end)
   enddo

   !@@@@@
   !write(*,*)EF_integral_range
   !> STEP FOUR
   !> exclude all kpoints away from Fermi level
   !write(*,*) Nband_Fermi_Level , "Nband"
   do ib= 1, Nband_Fermi_Level 
      !write(*,*) Nk_start, Nk_end, "Nk"
      !> first check howmany k points left
      it= 0
      do ik= Nk_start, Nk_end
         !> check whether v x B=0, which means the magnetic field is parallel with velocity
         v_t= velocity(:, ik, ib)
         !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
         !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
         !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
         !if (abs(Enk(ik, ib))<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
         
         if (abs(Enk(ik, ib))/eV2Hartree<EF_integral_range) then
            it = it+ 1
         endif
      enddo ! ik

      !> define KCube3D_left to store those left k points
      KCube3D_left(ib)%Nk_current= it
      !@@@@@
      !write(*,*) it
      if (it>0) then
         allocate(KCube3D_left(ib)%ik_array(it))
         allocate(KCube3D_left(ib)%Ek_local(it))
         allocate(KCube3D_left(ib)%vx_local(it))
         allocate(KCube3D_left(ib)%vy_local(it))
         allocate(KCube3D_left(ib)%vz_local(it))
         allocate(KCube3D_left(ib)%weight_k_local(it))
      else 
         allocate(KCube3D_left(ib)%ik_array(1))  !> only useful for mpi_allgatherv
         allocate(KCube3D_left(ib)%Ek_local(1))
         allocate(KCube3D_left(ib)%vx_local(1))
         allocate(KCube3D_left(ib)%vy_local(1))
         allocate(KCube3D_left(ib)%vz_local(1))
         allocate(KCube3D_left(ib)%weight_k_local(1))
      endif
      KCube3D_left(ib)%ik_array= 0
      KCube3D_left(ib)%Ek_local= 0d0
      KCube3D_left(ib)%vx_local= 0d0
      KCube3D_left(ib)%vy_local= 0d0
      KCube3D_left(ib)%vz_local= 0d0
      KCube3D_left(ib)%weight_k_local= 0d0


      it= 0
      do ik= Nk_start, Nk_end
         !> check whether v x B=0, which means the magnetic field is parallel with velocity
         v_t= velocity(:, ik, ib)
         !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
         !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
         !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
         !if (abs(Enk(ik, ib))<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
         if (abs(Enk(ik, ib))/eV2Hartree<EF_integral_range) then
            it = it+ 1
            KCube3D_left(ib)%weight_k_local(it) = KCube3D_symm%weight_k(ik)
            KCube3D_left(ib)%ik_array(it) = KCube3D_symm%ik_array_symm(ik)
            KCube3D_left(ib)%Ek_local(it) = Enk(ik, ib)
            KCube3D_left(ib)%vx_local(it) = v_t(1)
            KCube3D_left(ib)%vy_local(it) = v_t(2)
            KCube3D_left(ib)%vz_local(it) = v_t(3)
         endif
      enddo ! ik
   enddo ! ib

   !> try to get the total number of k points left for each band 
   do ib=1, Nband_Fermi_Level
#if defined (MPI)
      call mpi_allreduce(KCube3D_left(ib)%Nk_current,KCube3D_left(ib)%Nk_total,1,&
                     mpi_in,mpi_sum,mpi_cmw,ierr)
#else
      KCube3D_left(ib)%Nk_total= KCube3D_left(ib)%Nk_current
#endif
   enddo

   !> gather the number of k points left into a array info
   allocate(info(num_cpu, Nband_Fermi_Level))
   info= -1
   do ib= 1, Nband_Fermi_Level
      nrecevs= KCube3D_left(ib)%Nk_current
      if (nrecevs<0) nrecevs= 0
#if defined (MPI)
      call mpi_allgather(nrecevs, 1, mpi_in, info(:, ib), 1, mpi_in, mpi_cmw, ierr)
#else
      info(1, ib)= KCube3D_left(ib)%Nk_current
#endif
   enddo


   !> An array for mpi_allgatherv
   allocate(Displs(num_cpu+1, Nband_Fermi_Level))
   Displs= 0
   do ib=1, Nband_Fermi_Level
      Displs(1, ib)=0
      do i=2, num_cpu+1
         Displs(i, ib)=Displs(i-1, ib)+ info(i-1, ib)
      enddo
   enddo ! ib
#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   !> put all the kpoints left together
   do ib=1, Nband_Fermi_Level
      allocate(KCube3D_left(ib)%IKleft_array(KCube3D_left(ib)%Nk_total))
      allocate(KCube3D_left(ib)%Ek_total(KCube3D_left(ib)%Nk_total))
      allocate(KCube3D_left(ib)%vx_total(KCube3D_left(ib)%Nk_total))
      allocate(KCube3D_left(ib)%vy_total(KCube3D_left(ib)%Nk_total))
      allocate(KCube3D_left(ib)%vz_total(KCube3D_left(ib)%Nk_total))
      allocate(KCube3D_left(ib)%weight_k(KCube3D_left(ib)%Nk_total))
      KCube3D_left(ib)%IKleft_array = 0
      KCube3D_left(ib)%Ek_total= 0d0
      KCube3D_left(ib)%vx_total= 0d0
      KCube3D_left(ib)%vy_total= 0d0
      KCube3D_left(ib)%vz_total= 0d0
      KCube3D_left(ib)%weight_k= 0d0
   enddo  ! ib
   !> gather Enk and velocity 
#if defined (MPI)
   do ib=1, Nband_Fermi_Level
      !nrecevs= KCube3D_left(ib)%Nk_end-KCube3D_left(ib)%Nk_start+ 1
      nrecevs= KCube3D_left(ib)%Nk_current
      if (nrecevs<0) nrecevs= 0
      call mpi_allgatherv(KCube3D_left(ib)%ik_array, nrecevs, &
                           mpi_in, KCube3D_left(ib)%IKleft_array, &
                           info(:, ib), Displs(:, ib), mpi_in, mpi_cmw, ierr)
      call mpi_allgatherv(KCube3D_left(ib)%Ek_local, nrecevs, &
                           mpi_dp, KCube3D_left(ib)%Ek_total, &
                           info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
      call mpi_allgatherv(KCube3D_left(ib)%vx_local, nrecevs, &
                           mpi_dp, KCube3D_left(ib)%vx_total, &
                           info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
      call mpi_allgatherv(KCube3D_left(ib)%vy_local, nrecevs, &
                           mpi_dp, KCube3D_left(ib)%vy_total, &
                           info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
      call mpi_allgatherv(KCube3D_left(ib)%vz_local, nrecevs, &
                           mpi_dp, KCube3D_left(ib)%vz_total, &
                           info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
      call mpi_allgatherv(KCube3D_left(ib)%weight_k_local, nrecevs, &
                           mpi_dp, KCube3D_left(ib)%weight_k, &
                           info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
   enddo ! ib
#else
   do ib=1, Nband_Fermi_Level
      KCube3D_left(ib)%IKleft_array= KCube3D_left(ib)%ik_array
      KCube3D_left(ib)%Ek_total= KCube3D_left(ib)%Ek_local
      KCube3D_left(ib)%vx_total= KCube3D_left(ib)%vx_local
      KCube3D_left(ib)%vy_total= KCube3D_left(ib)%vy_local
      KCube3D_left(ib)%vz_total= KCube3D_left(ib)%vz_local
      KCube3D_left(ib)%weight_k= KCube3D_left(ib)%weight_k_local
   enddo ! ib
#endif

   !> redistribute all those k points into different cpus
   if (cpuid.eq.0) then
      write(stdout, '(a)')' '
      write(stdout, '(a, i10, a)')' There are ', KCube3D_total%Nk_total, ' k points generated by the input files' 
      do ib= 1, Nband_Fermi_Level
         write(stdout, '(a, i10, a, i10)')' However there are only ', KCube3D_left(ib)%Nk_total, &
            ' k points contribute to the Ln/tau calculations for band ', bands_fermi_level(ib)
      enddo
   endif


   !> redistribute the left kpoints into different processors
   !> for different bands, the number of left kpoints is different.
   do ib= 1, Nband_Fermi_Level
      knv3_left= KCube3D_left(ib)%Nk_total
      knv3_left_mod= mod(knv3_left, num_cpu)
      if (knv3_left_mod==0) then  !> perfect divided
         KCube3D_left(ib)%Nk_current= knv3_left/num_cpu
         KCube3D_left(ib)%Nk_start=1+ knv3_left*cpuid/num_cpu
         KCube3D_left(ib)%Nk_end  =(1+cpuid)*knv3_left/num_cpu
      else if (knv3_left/num_cpu==0) then    !> Number of MPI threads is large than knv3_left
         KCube3D_left(ib)%Nk_current= 1 !> one k piont per MPI thread
         KCube3D_left(ib)%Nk_start= cpuid+ 1 !> one k piont per MPI thread
         KCube3D_left(ib)%Nk_end  = cpuid+ 1
         if (cpuid+1 > knv3_left) then
            KCube3D_left(ib)%Nk_start= 1
            KCube3D_left(ib)%Nk_end  = 0
         endif
      else
         KCube3D_left(ib)%Nk_current= knv3_left/num_cpu+ 1
         if (cpuid< knv3_left_mod) then
            KCube3D_left(ib)%Nk_start= 1+ cpuid*KCube3D_left(ib)%Nk_current
            KCube3D_left(ib)%Nk_end  = (1+cpuid)*KCube3D_left(ib)%Nk_current
         else
            KCube3D_left(ib)%Nk_start= knv3_left_mod*KCube3D_left(ib)%Nk_current+ &
               (cpuid-knv3_left_mod)*(KCube3D_left(ib)%Nk_current-1)+1
            KCube3D_left(ib)%Nk_end  = knv3_left_mod*KCube3D_left(ib)%Nk_current+ &
               (cpuid-knv3_left_mod+1)*(KCube3D_left(ib)%Nk_current-1)
         endif
      endif


      allocate(KCube3D_left(ib)%k_direct(3, KCube3D_left(ib)%Nk_start:KCube3D_left(ib)%Nk_end))
   
      do ik= KCube3D_left(ib)%Nk_start, KCube3D_left(ib)%Nk_end
         i= KCube3D_left(ib)%IKleft_array(ik)
         ik1= (i-1)/(Nk2*Nk3)+1
         ik2= ((i-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
         ik3= (i-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
         KCube3D_left(ib)%k_direct(1, ik)= (ik1-1)/dble(Nk1) 
         KCube3D_left(ib)%k_direct(2, ik)= (ik2-1)/dble(Nk2) 
         KCube3D_left(ib)%k_direct(3, ik)= (ik3-1)/dble(Nk3) 
      enddo ! ik
      
      !> allocate array to store the Ln/tau
      allocate( Ln_tau_iband_k(ib)%Ln_tau_k(3, 9, NBTau, OmegaNum, NumT), stat= ierr)
      if (ierr>0) stop ' Error : not enough memory'
      allocate( Ln_tau_iband_k(ib)%Ln_tau_k_mpi(3, 9, NBTau, OmegaNum, NumT), stat= ierr)
      if (ierr>0) stop ' Error : not enough memory'
      allocate( Ln_tau_iband_k(ib)%time_cost(KCube3D_left(ib)%Nk_total))
      allocate( Ln_tau_iband_k(ib)%time_cost_mpi(KCube3D_left(ib)%Nk_total))
      Ln_tau_iband_k(ib)%time_cost= 0d0
      Ln_tau_iband_k(ib)%time_cost_mpi= 0d0
      Ln_tau_iband_k(ib)%Ln_tau_k= 0d0
      Ln_tau_iband_k(ib)%Ln_tau_k_mpi= 0d0
   enddo  ! ib=1, Nband_Fermi_Level


   !> gather the number of receive buffs left into a array info
   info= -1
   do ib= 1, Nband_Fermi_Level
      nrecevs= (KCube3D_left(ib)%Nk_end-KCube3D_left(ib)%Nk_start+1)*9*NBTau*OmegaNum*NumT*3
      if (nrecevs<0) nrecevs=0
#if defined (MPI)
      call mpi_allgather(nrecevs, 1, mpi_in, info(:, ib), 1, mpi_in, mpi_cmw, ierr)
#else
      info(1, ib)= nrecevs
#endif
   enddo

   if (cpuid.eq.0) write(stdout, '(a, 10i8)') ' '
   do ib=1, Nband_Fermi_Level
      if (cpuid.eq.0) write(stdout, '(a, i7)')'>> Number of k points at different CPUs at band', ib
      if (cpuid.eq.0) write(stdout, '(10i8)')info(:, ib)/9/NBTau/OmegaNum/NumT/3
   enddo
   if (cpuid.eq.0) write(stdout, '(a, 10i8)') ' '


   !> An array for mpi_allgatherv
   Displs= 0
   do ib=1, Nband_Fermi_Level
      Displs(1, ib)=0
      do i=2, num_cpu+1
         Displs(i, ib)=Displs(i-1, ib)+ info(i-1, ib)
      enddo
   enddo ! ib

   do ib=1, Nband_Fermi_Level
      allocate(klist_iband(ib)%klist_rkfs(3, NSlice_Btau))
      allocate(klist_iband(ib)%velocity_k(3, NSlice_Btau))
      klist_iband(ib)%klist_rkfs= 0d0
      klist_iband(ib)%velocity_k= 0d0
   enddo

   !> a temp array used in RKFS
   allocate(kout(3, NSlice_Btau))
   kout= 0d0

   !> now we turn to use Runge-Kutta method to get all the kpoints from (0, BTauMax)
   !> and we calculate the conductivity/Tau over different bands and different k points
   time_start= 0d0
   time_end  = 0d0

   !> we get the kpath by Btau_final=-exponent_max*BTauMax, but we only use half of them
      !>  means that we can reach the accuracy as to exp(-exponent_max)
   exponent_max= 30d0

   do ib= 1, Nband_Fermi_Level
      call now(time_start0)
   
      !> dim=(Nk_start: Nk_end, NumT, OmegaNum))
      allocate(klist_iband(ib)%minusdfde(OmegaNum, NumT))
      do ik= KCube3D_left(ib)%Nk_start, KCube3D_left(ib)%Nk_end
         if (cpuid.eq.0) &
            write(stdout, '(a, i8, a, i18, "   /", i18, a, f10.3, "s", a, f10.3, "s")') &
            'In Ln iband', ib, ' ik/NK', &
            ik-KCube3D_left(ib)%Nk_start+1,KCube3D_left(ib)%Nk_current, &
            ' time cost', time_end-time_start, &
            ' time left', &
            (KCube3D_left(ib)%Nk_current+KCube3D_left(ib)%Nk_start-ik)*(time_end-time_start)

         call now(time_start)
         EE= KCube3D_left(ib)%Ek_total(ik)
         v_k(1)= KCube3D_left(ib)%vx_total(ik)
         v_k(2)= KCube3D_left(ib)%vy_total(ik)
         v_k(3)= KCube3D_left(ib)%vz_total(ik)
         
         !> calculate df/de for each k point and each band
         do ikt=1, NumT
            KBT= KBT_array(ikt)
            do ie=1, OmegaNum
               mu= mu_array(ie)
               call minusdfde_calc_single(EE, KBT, mu,  minusdfde)
               klist_iband(ib)%minusdfde(ie, ikt)= minusdfde
            enddo ! ie
         enddo ! ikt


         !> start to get the evolution of k points under magnetic field using Runge-Kutta method
         k_start= KCube3D_left(ib)%k_direct(:, ik)
         kout= 0d0
         Btau_start= 0d0

         !@@@@@ open the failed kpoints file in each cpu 
         !failed_k_fileindex = 114514
         !write(failed_k_filename, '(a, i6, a)') 'failed_kpoints_cpuid=', cpuid, '.dat'
         !open(unit=failed_k_fileindex+cpuid, file=failed_k_filename)

         !> we get the kpath by Btau_final=-30*BTauMax, but we only use half of them
         Btau_final= -exponent_max*BTauMax   !< -15 means that we can reach the accuracy as to exp(-15d0)
         !Btau_final= -10d0*BTauMax   !< test for arXiv:2201.03292 (2022)

         !> Runge-Kutta only applied with BTauMax>0
         !> if the magnetic field is zero. 
         if (BTauMax>eps3) then
            NSlice_Btau_inuse = NSlice_Btau
            call RKF45_pack(magnetic_field, bands_fermi_level(ib),  &
                  NSlice_Btau_inuse, k_start, Btau_start, Btau_final, kout, icycle, fail)
         else
            do ibtau=1, NSlice_Btau
               kout(:, ibtau)= k_start(:)
            enddo
            NSlice_Btau_inuse = NSlice_Btau
         endif

         !> NSlice_Btau_inuse is the number of slices that we really use, COMMENT BY PHQ
         !> This can be combined with the following
         if (NSlice_Btau_inuse==1) then
            write(stdout, '(a, i6, a, i4, a, i6, a, 3f12.6)')&
               '>>> NSlice_Btau_inuse=1 at cpuid=', cpuid, ' ib=', ib, ' ik', ik, ' k', k_start
         endif
 
         if (fail) then
            write(stdout, '(a, i6, a, i4, a, i6, a, 3f12.6)')&
               '>>> Runge-Kutta integration fails at cpuid=', cpuid, ' ib=', ib, ' ik', ik, ' k', k_start
            write(stdout, *)' '

            !@@@@@ failed kpoints write out
            !write(failed_k_fileindex + cpuid, '(a, i6, a, i4, a, 3f12.6)')&
            !   'cpuid=', cpuid, ' ib=', ib, ' k', k_start
            !cycle
         endif

         if (NSlice_Btau_inuse==1)then
            call velocity_calc_iband(bands_fermi_level(ib), k_start, v_t)
            do ibtau=1, NSlice_Btau
               kout(:, ibtau)= k_start(:)
               klist_iband(ib)%velocity_k(:, ibtau)= v_t
            enddo
            NSlice_Btau_inuse = NSlice_Btau
         else
            do it= 1, icycle
               k= kout(:, it) 
               call velocity_calc_iband(bands_fermi_level(ib), k, v_t)
               klist_iband(ib)%velocity_k(:, it)= v_t
            enddo ! integrate over time step

            !> periodic kpath in the BZ can be reused
            do i=2, NSlice_Btau_inuse/icycle
               do j=1, icycle
                  klist_iband(ib)%velocity_k(:, j+(i-1)*icycle)= klist_iband(ib)%velocity_k(:, j)
               enddo
            enddo
            do i=(NSlice_Btau_inuse/icycle)*icycle+1, NSlice_Btau
               klist_iband(ib)%velocity_k(:, i)= klist_iband(ib)%velocity_k(:, i-(NSlice_Btau_inuse/icycle)*icycle)
            enddo
         endif

         !> print the icycle into stdout
         if (cpuid.eq.0) write(stdout, '(a, i6, a, i4, a, i6, a, i6)')&
            '>>> icycle=', icycle, ' at cpuid=', cpuid, ' ib=', ib, ' ik', ik

         !> calculate the L(n)
         do ikt=1, NumT
            KBT= KBT_array(ikt)
            do ie=1, OmegaNum
               mu= mu_array(ie)
               minusdfde= klist_iband(ib)%minusdfde(ie, ikt)

               do ibtau=1, NBTau
                  BTau= BTau_array(ibtau)

                  if (NBTau==1)then
                     NSlice_Btau_local= 2
                  else
                     NSlice_Btau_local= (ibtau-1)*NSlice_Btau_inuse/(NBTau-1)/2
                     if (NSlice_Btau_local==0)then
                        NSlice_Btau_local= 2
                     else
                        DeltaBtau= exponent_max/2d0/NSlice_Btau_local
                     endif
                  endif

                  if (allocated(exp_factors)) deallocate(exp_factors)
                  allocate(exp_factors(NSlice_Btau_local))

                  !> here, velocity is in the cartesian coordinate
                  !> The core of Chamber formular is to get the average of velocity 
                  !> in the relaxation time approximation
                  v_k= klist_iband(ib)%velocity_k(:, 1)
                  if (BTau>eps3.and.NSlice_Btau_local>2) then
                     velocity_bar_k= 0d0
                     ! do it=1, NSlice_Btau_local
                     !    velocity_bar_k= velocity_bar_k+ &
                     !       DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(ib)%velocity_k(:, it)
                     ! enddo
                     !> vectorization the process
                     do i = 1, NSlice_Btau_local
                        exp_factors(i) = -(i - 1.0_dp) * DeltaBtau
                     end do
                     exp_factors = exp(exp_factors)
                     velocity_bar_k = DeltaBtau * MATMUL(klist_iband(ib)%velocity_k(:, 1:NSlice_Btau_local),exp_factors)
                     velocity_bar_k= velocity_bar_k &
                     - 0.5d0*DeltaBtau*(exp(-(NSlice_Btau_local-1d0)*DeltaBtau)&
                     * klist_iband(ib)%velocity_k(:, NSlice_Btau_local)  &
                     + klist_iband(ib)%velocity_k(:, 1))
                  else
                     velocity_bar_k= v_k
                  endif

               !   !@@@@@ 
               !   magnetic_field_antisym_tensor(:,1) = (/0d0, Bdirection(3), -Bdirection(2)/)
               !   magnetic_field_antisym_tensor(:,2) = (/-Bdirection(3), 0d0, Bdirection(1)/)
               !   magnetic_field_antisym_tensor(:,3) = (/Bdirection(2), -Bdirection(1), 0d0/)
               !   magnetic_field_antisym_tensor = magnetic_field_antisym_tensor*BTau
               !   !iband = 5 
               !   M_be_inverse(1,:,1) =(/14.75d0, 0d0, 0d0 /)
               !   M_be_inverse(1,:,2) =(/0d0, 14.75d0, 0d0 /)
               !   M_be_inverse(1,:,3) =(/0d0, 0d0, 1.387d0 /)
               !   Lax_gamma(1) = EE
               !   Lax_gamma_prime(1) = 1d0
               !
               !   !iband = 7
               !   M_be_inverse(2,:,1) =(/806d0, 0d0, 0d0 /)
               !   M_be_inverse(2,:,2) =(/0d0, 7.95d0, 37.6d0/)
               !   M_be_inverse(2,:,3) =(/0d0, 37.6d0, 349d0 /)
               !   Lax_gamma(2) = EE*(1d0+EE/0.1d0)
               !   Lax_gamma_prime(2) = 1d0+2d0*EE/0.1d0                      
               !
               !   velocity_bar_k_article = matmul(v_k ,(Identity3-&
               !   (1d0/Lax_gamma_prime(2))*matmul(magnetic_field_antisym_tensor,M_be_inverse(2,:,:))) )
               !
               !   
               !   if (cpuid == 0) then
               !      write(unit=velocity_bar_fileindex, fmt='(a,a, i5)')"#####", " iband = ",ib
               !      write(unit=velocity_bar_fileindex, fmt='(a," ",f16.10," ",f16.10," ",f16.10)') "v_bar", velocity_bar_k
               !      write(unit=velocity_bar_fileindex, fmt='(a," ",f16.10," ",f16.10," ",f16.10)') "v_bar_article", velocity_bar_k_article
               !      write(unit=velocity_bar_fileindex, fmt='(a, 3f16.12)') "k", KCube3D_left(ib)%k_direct(:, ik)
               !   end if
               !
               !   velocity_bar_k =matmul(v_k ,(Identity3-&
               !   (1d0/Lax_gamma_prime(2))*matmul(magnetic_field_antisym_tensor,M_be_inverse(2,:,:))) )

                  !> calculate the Ln/tau now
                  do n = 1, 3
                     Ln_tau_symm_t= 0d0
                     !> Apply point group operations to the velocities, and average them
                     do j1=1, 3
                     do j2=1, 3
                     do j=1, number_group_operators
                     !> Ln_tau_xx
                        Ln_tau_symm_t(1)= Ln_tau_symm_t(1)+ &
                        pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(1, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_xy
                        Ln_tau_symm_t(2)= Ln_tau_symm_t(2)+ &
                        pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(2, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_xz
                        Ln_tau_symm_t(3)= Ln_tau_symm_t(3)+ &
                        pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(3, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_yx
                        Ln_tau_symm_t(4)= Ln_tau_symm_t(4)+ &
                        pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(1, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_yy
                        Ln_tau_symm_t(5)= Ln_tau_symm_t(5)+ &
                        pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(2, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_yz
                        Ln_tau_symm_t(6)= Ln_tau_symm_t(6)+ &
                        pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(3, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_zx
                        Ln_tau_symm_t(7)= Ln_tau_symm_t(7)+ &
                        pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(1, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_zy
                        Ln_tau_symm_t(8)= Ln_tau_symm_t(8)+ &
                        pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(2, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     !> Ln_tau_zz
                        Ln_tau_symm_t(9)= Ln_tau_symm_t(9)+ &
                        pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(3, j2, j) &
                        *v_k(j1)*velocity_bar_k(j2)*(EE-mu)**(n-1)
                     enddo
                     enddo
                     enddo
                     Ln_tau_iband_k(ib)%Ln_tau_k_mpi(n, :, ibtau, ie, ikt) = &
                     Ln_tau_iband_k(ib)%Ln_tau_k_mpi(n, :, ibtau, ie, ikt) + &
                     Ln_tau_symm_t/dble(number_group_operators)*minusdfde*KCube3D_left(ib)%weight_k(ik)
                  enddo ! n

               enddo ! ibtau  Btau
            enddo ! ie  mu
         enddo ! ikt KBT
         call now(time_end)
         Ln_tau_iband_k(ib)%time_cost_mpi(ik)= time_end- time_start
         if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for this loop is ', time_end- time_start, ' s'
      enddo ! ik  kpoints

      call now(time_start)
#if defined (MPI)
      call mpi_allreduce(Ln_tau_iband_k(ib)%time_cost_mpi, Ln_tau_iband_k(ib)%time_cost, &
                           size(Ln_tau_iband_k(ib)%time_cost), &
                           mpi_dp,mpi_sum,mpi_cmw,ierr)
      ! nrecevs= (KCube3D_left(ib)%Nk_end-KCube3D_left(ib)%Nk_start+1)*9*NBTau*OmegaNum*NumT*3
      call mpi_allreduce(Ln_tau_iband_k(ib)%Ln_tau_k_mpi, Ln_tau_iband_k(ib)%Ln_tau_k, &
                           size(Ln_tau_iband_k(ib)%Ln_tau_k), mpi_dp, mpi_sum, mpi_cmw, ierr)                       
      if (ierr>0) then
         write(stdout, *)'>>>Error happends in mpi_allreduce at cpuid', cpuid, ' ierr=', ierr
         stop
      endif

#else
      Ln_tau_iband_k(ib)%Ln_tau_k= Ln_tau_iband_k(ib)%Ln_tau_k_mpi
      Ln_tau_iband_k(ib)%time_cost= Ln_tau_iband_k(ib)%time_cost_mpi
#endif

      call now(time_end)
      if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for mpi_allreduce is ', time_end- time_start, ' s'
      if (cpuid.eq.0) then
         write(stdout, '(a)')' '
         write(stdout, '(a, i10)')'>> Time cost for each k point at iband ', ib
         write(stdout, '(10f16.2)')(Ln_tau_iband_k(ib)%time_cost(ik), ik= 1, KCube3D_left(ib)%Nk_total)
         write(stdout, '(a)')' '
      endif

   ! --------------------------------------------------------------------
   ! At this point Ln_tau corresponding to L^(n-1) with n=1,2,3 contains
   !
   !           (1/N) sum_k (-df/dE) v vbar (E-mu)^(n-1),
   !
   ! an approximation to
   !
   !           V_c int dk/[8(pi)^3] (-df/dE) v vbar (E-mu)^(n-1),
   !
   ! V_c is the cell volume and has the unit of Bohr_radius^(3), dk has the unit of Bohr_radius^(-3)
   ! (-df/dE) and (E-mu) have the unit of Hatree^(-1) and Hartree, respectively
   ! v and vbar have the unit of Hartree*Bohr_radius/hbar
   ! Thus the above quantity has the unit of Bohr_radius^(2)*Hartree^(n)/hbar^2
   !
   ! We want
   !
   !           Ln/tau=( e^2 ) int dk/[8(pi)^3] (-df/dE) v vbar (E-mu)^(n-1)
   !
   ! (i) Hence we need to multiply by factor1 = 1/V_c 
   !     where 'V_c' in Bohr_radius^3, 'e' is 1 in atomic units
   !
   ! (ii) We then convert atomic units to SI units, 
   !      So we multiply by factor2 = Bohr_radius^(-1)*Echarge^2*Hartree^(n)/hbar^2
   !
   ! (iii) In summary, we need to multiply the Ln_tau_bands with 
   !       factor=factor1*factor2= 2*e^n/(hbar^2*V_c)*(eV2Hartree)^(2-n)/Bohr_radius,
   !       Then Ln_tau_bands*factor has the units of C^(n)*V^(n-2)*m^-1*s^-2
   !
   ! Take sigma/tau = L1/tau = L1_tau*factor as an example, 
   ! it has units of C*V^-1*m^-1*s^-2 = Omega^-1*m^-1*s^-1
   ! --------------------------------------------------------------------

      do ikt=1, NumT
         do ie=1, OmegaNum
            do ibtau=1, NBTau
               do i=1, 9
                  do n =1, 3
                     !> different coeff for Ln/tau
                     coeff= Echarge**2/Bohr_radius*Hartree_SI**n/hbar**2/Origin_cell%CellVolume & 
                            /kCubeVolume*Origin_cell%ReciprocalCellVolume
                     !> multiply coeff
                     Ln_tau_bands(n, i, ibtau, ie, ikt, ib)=  &
                     Ln_tau_iband_k(ib)%Ln_tau_k(n, i, ibtau, ie, ikt)*coeff
                  enddo ! n
               enddo ! i
            enddo ! ibtau
         enddo ! ie
      enddo ! ikt
      
      call now(time_end0)
      if (cpuid.eq.0) write(stdout, '(a, i6, a, f16.2, a)')'>> Time cost for calculate Ln/tau at ib=', ib, &
      'is ', time_end0- time_start0, ' s'

   enddo ! ib=1, Nband_Fermi_Level

   if (cpuid.eq.0) then
      do ie=1, OmegaNum
         write(muname, '(f12.3)')mu_array(ie)/eV2Hartree

         !> write Ln/tau for every bands into file
         do n = 1, 3
            outfileindex = outfileindex+1
            write(Lnfname, '(a,I1,3a)')'L', n, '_bands_mu_',trim(adjustl(muname)),'eV.dat'
            open(unit=outfileindex, file=Lnfname)
            write(outfileindex, '(a,I1,a)')'# L', n, '/tau  for every contributing band' 
            write(outfileindex, '(a,3I6)')'# NBAND  NumT  NumBtau  =  ',Nband_Fermi_Level, NumT, NBTau 
            write(outfileindex, '(a,100I5)')'# SELECTEDBANDS  =  ', bands_fermi_level(:)
            write(outfileindex, '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
            do ib= 1, Nband_Fermi_Level
               write(outfileindex, '(2a, i5)')'# ',' iband = ', bands_fermi_level(ib)
               write(outfileindex,'(a)') ''
               do ikt = 1, NumT
                  KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
                  write(outfileindex, '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
                     write(outfileindex, '("# Column", i5, 100i16)')(i, i=1, 10)
                     write(outfileindex, '("#",19a16)')'BTau (T.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
                  !> write out the Ln/tau into file
                  do i=1, NBTau
                     write(outfileindex, '(19E16.6)') &
                        BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                        Ln_tau_bands(n, :, i, ie, ikt, ib)
                  enddo ! i, NBTau
                  write(outfileindex,'(a)') ''
               enddo
               write(outfileindex,'(a)') ''
            enddo
            close(outfileindex)
         enddo !> n = 1, 3
      enddo
   endif
  

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   return
end subroutine Ln_tau_calc_symm

subroutine single_band_seebeck()
!================================================================
! This subroutine is used to calculate the magneto-seebek effect
! of a single band with H(k) = \hbar \sum_i (1-cos(k_i))/m_i, 
! where i=x,y,z 
!================================================================
   use wmpi
   use para
   implicit none

   integer :: ie, iT, ik, ikx, iky, ikz, i, j, ib, ibtau, ikT, iband, ierr, &
              ibtau_local
   integer :: knv3, Nslice_Btau, Nslice_Btau_local
   real(dp) :: time_start, time_end, minusdfde, EE, mu, KBT, alpha_tau, hstep
   real(dp) ::  k(3), kxyz(3), Rxyz(3), vk(3), mass(3), L11_inv(3, 3), seebeck(3, 3)
   real(dp), allocatable :: mu_array(:), T_array(:), KBT_array(:), Btau_array(:),&
                            vout(:, :), kout(:, :), vkbar(:, :)
   character(len=100) :: seebeckfilename, muname, L11filename, L12filename

   real(dp), allocatable :: L11(:, :, :, :, :), L12(:, :, :, :, :), &
                            L11_mpi(:, :, :, :, :), L12_mpi(:, :, :, :, :)

   real(dp), external :: det3

   allocate( L11_mpi(3, 3, OmegaNum, NumT, NBTau))
   allocate( L12_mpi(3, 3, OmegaNum, NumT, NBTau))
   allocate( L11(3, 3, OmegaNum, NumT, NBTau))
   allocate( L12(3, 3, OmegaNum, NumT, NBTau))
   L11= 0d0
   L12= 0d0
   L11_mpi= 0d0
   L12_mpi= 0d0

   allocate(mu_array(OmegaNum))
   allocate(T_array(NumT))
   allocate(KBT_array(NumT))
   allocate(Btau_array(NBTau))
   allocate(vkbar(3, NBTau))
   vkbar = 0d0

   !> energy
   do ie=1, OmegaNum
      if (OmegaNum>1) then
         mu_array(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
      else
         mu_array= OmegaMin
      endif
   enddo ! ie

   !> temperature
   do iT=1, NumT
       if (NumT>1) then
          T_array(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
       else
          T_array= Tmin
       endif
   enddo ! iT

   !> transform from Kelvin to eV
   !> The SI unit of temperature is the kelvin (K), but using the above relation the electron temperature is often expressed in
   !> terms of the energy unit electronvolt (eV). Each kelvin (1 K) corresponds to 8.6173324(78)×10−5 eV; this factor is the ratio
   !> of the Boltzmann constant to the elementary charge. After version 2.6, we 
   !> adopt the atomic unit
   KBT_array= T_array*8.6173324E-5*eV2Hartree

   !> number of Btau is NBTau, total number of Btau points is Nslice_BTau which does not include Btau = 0
   !> number of Btau points per line is NSlice_Btau_local (not include the origin but include the end point of every line)
   do ibtau=1, NBTau
      if (NBTau>1) then
         Btau_array(ibtau)= (ibtau-1.0d0)/dble(NBTau-1)*BTauMax 
      else
         Btau_array= 0      
      endif
   enddo ! ibtau

   if (NBTau>1) then
      NSlice_Btau= (NBTau-1)*(Nslice_BTau_Max/(NBTau-1))
      Nslice_Btau_local = Nslice_BTau_Max/(NBTau-1)
   else
      NSlice_Btau= 1
      Nslice_Btau_local = 1
   endif
   if (cpuid == 0) write(stdout, '(a, i6, a, i6)')'>> NSlice_Btau=', NSlice_Btau, ' Nslice_Btau_local=', Nslice_Btau_local

   allocate(vout(3, NSlice_Btau+1))
   allocate(kout(3, NSlice_Btau+1))
   vout= 0d0
   kout= 0d0
   
   knv3= Nk1*Nk2*Nk3
   alpha_tau = 15d0
   Rxyz(1) = Origin_cell%Rua(1)
   Rxyz(2) = Origin_cell%Rub(2)
   Rxyz(3) = Origin_cell%Ruc(3)
   mass(1) = m1x
   mass(2) = m1y
   mass(3) = m1z

   call now(time_start) 
   do ik= 1+ cpuid, knv3, num_cpu
      if (cpuid.eq.0.and. mod(ik/num_cpu, 1000).eq.0) then
         call now(time_end) 
         write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
         ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/1000d0
         time_start= time_end
      endif

      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
       + K3D_vec2_cube*(iky-1)/dble(nk2)  &
       + K3D_vec3_cube*(ikz-1)/dble(nk3)

      kxyz = k(1)*Origin_cell%Kua + k(2)*Origin_cell%Kub + k(3)*Origin_cell%Kuc

      EE = (1-cos(k(1)*twopi))/m1x+(1-cos(k(2)*twopi))/m1y+(1-cos(k(3)*twopi))/m1z
      EE = EE*eV2Hartree

      do i=1,3
         vk(i) = sin(kxyz(i)*Rxyz(i))*Rxyz(i)/mass(i)*eV2Hartree
      enddo
      ! write(stdout, '(a, 4f12.6)') 'v_k=', vk(:), EE

      if (NBTau == 1) then
         kout(:, 1) = kxyz
         vout(:, 1) = vk
      else 
         hstep = alpha_tau*BTauMax/Nslice_Btau
         call ode45(kout, kxyz, hstep, Nslice_Btau+1, m1x, m1y, m1z)
         do ibtau_local = 1, Nslice_Btau+1
            do i = 1, 3
               vout(i, ibtau_local) = sin(kout(i, ibtau_local)*Rxyz(i))*Rxyz(i)/mass(i)*eV2Hartree
            enddo
         enddo
      endif
      
      vkbar(:, 1) = vk
      if (NBTau >1) then
         do ibtau = 2, NBTau 
            vkbar(:, ibtau) = alpha_tau/(Nslice_Btau_local*(ibtau-1))* &
                           (vout(:,1)+vout(:, Nslice_Btau_local*(ibtau-1)+1)*exp(-alpha_tau))/2
            ! if (ibtau == 2) then
            !       write(stdout, '(a, 6f12.6)') 'alpha_tau/(Nslice_Btau_local*(ibtau-1)) = ', alpha_tau/(Nslice_Btau_local*(ibtau-1))
            !       write(stdout, '(a, 9f12.6)') 'vk, vout(:, use)*exp(), vkbar = ', vk(:), vout(:, Nslice_Btau_local*(ibtau-1))*exp(-alpha_tau),  vkbar(:, ibtau)
            ! endif
            do ibtau_local = 2, (Nslice_Btau_local*(ibtau-1))
               vkbar(:, ibtau) = vkbar(:, ibtau) + alpha_tau/(Nslice_Btau_local*(ibtau-1))* &
                     vout(:, ibtau_local)*exp(-alpha_tau*ibtau_local/(Nslice_Btau_local*(ibtau-1)))
               ! if (ibtau == 2) then
               !    write(stdout, '(a, I5, 7f12.6)') 'ibtau_local, vkbar(i, ibtau), vout(:, ibtau_local) = ', ibtau_local, vkbar(:, ibtau), vout(:, ibtau_local), exp(-alpha_tau*ibtau_local/(Nslice_Btau_local*(ibtau-1)))
               ! endif
            enddo
         enddo
      endif
      
      ! outfileindex = outfileindex+1
      ! open(unit=outfileindex, file='kout.dat')
      ! write(outfileindex, '("#", 6a12)')'kx', 'ky', 'kz', 'vx', 'vy', 'vz'
      ! do ibtau_local = 1, Nslice_Btau+1
      !    write(outfileindex, '(6f12.6)') kout(:, ibtau_local)*Angstrom2atomic/Pi, vout(:, ibtau_local)/Angstrom2atomic
      ! enddo
      ! write(outfileindex,'(2a12)')'vxvybar','vyvxbar'
      ! write(outfileindex,'(2f12.6)') vk(1)*vkbar(2, 2), vk(2)*vkbar(1, 2) 
      ! write(stdout, '("#", 6a12)')'kx', 'ky', 'kz'
      ! write(stdout, '(6f12.6)') k(:)
      ! write(stdout,'(2f12.6)') vk(1)*vkbar(2, 2), vk(2)*vkbar(1, 2) 


      !> calculate df/de for each k point and each band
      do ibtau = 1, Nbtau
         do ikt=1, NumT
            KBT= KBT_array(ikt)
            do ie=1, OmegaNum
               mu= mu_array(ie)
               call minusdfde_calc_single(EE, KBT, mu, minusdfde)
               do i = 1, 3
                  do j = 1, 3
                     L11_mpi(i, j, ie, ikt, ibtau)= L11_mpi(i, j, ie, ikt, ibtau)+ &
                                          minusdfde*vk(i)*vkbar(j, ibtau)/dble(knv3)
                     L12_mpi(i, j, ie, ikt, ibtau)= L12_mpi(i, j, ie, ikt, ibtau)+ &
                                          minusdfde*vk(i)*vkbar(j, ibtau)*(EE-mu)/dble(knv3)
                  enddo
               enddo 
            enddo
         enddo 
      enddo 
     
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(L11_mpi, L11, size(L11), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(L12_mpi, L12, size(L12), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   L11 = L11_mpi
   L12 = L12_mpi
#endif
   L11 = L11*Echarge**2/Bohr_radius*Hartree_SI/hbar**2/Origin_cell%CellVolume & 
         /kCubeVolume*Origin_cell%ReciprocalCellVolume
   L12 = L12*Echarge**2/Bohr_radius*Hartree_SI**2/hbar**2/Origin_cell%CellVolume &
         /kCubeVolume*Origin_cell%ReciprocalCellVolume

   
   if(cpuid .eq. 0) then
      do ie=1, OmegaNum
         write(muname, '(f12.3)')mu_array(ie)/eV2Hartree
         outfileindex = outfileindex+1
         write(seebeckfilename, '(6a)')'seebeck_single','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=seebeckfilename)
         write(outfileindex, '("#",a)')'Seebeck coefficient  (V/K)'
         write(outfileindex, "('#column', i7, i13,3i16)")(i, i=1, 5)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T \cdot ps)', 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_zx', 'S_zy', 'S_zz'
         do iT = 1, NumT
            do ibtau = 1, NBTau
               L11_inv = L11(:, :, ie, iT, ibtau)
               if (abs(det3(L11_inv))>eps6) then
                  call inv_r(3, L11_inv)
                  seebeck = matmul(L11_inv, -L12(:, :, ie, iT, ibtau)/Echarge/T_array(iT))
                  write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                            seebeck(1, 1), seebeck(1, 2), seebeck(1, 3), seebeck(2, 1), seebeck(2, 2), seebeck(2, 3), &
                                                            seebeck(3, 1), seebeck(3, 2), seebeck(3, 3)
               else 
                  write(outfileindex, '(a)') 'No solution because the conductivity tensor is singular'
               endif
            enddo
         write(outfileindex,'(a)') ''
         enddo
         close(outfileindex)

         outfileindex = outfileindex+1
         write(L11filename, '(6a)')'L11_single','_mu_',trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=L11filename)
         write(outfileindex, '("#",a)')'L11  (S/m)'
         write(outfileindex, "('#column', i7, i13,3i16)")(i, i=1, 5)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T \cdot ps)', 'L11_xx', 'L11_xy', 'L11_xz', 'L11_yx', 'L11_yy', 'L11_yz', 'L11_zx', 'L11_zy', 'L11_zz'
         do iT = 1, NumT
            do ibtau = 1, NBTau
               write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                         L11(1, 1, ie, iT, ibtau), L11(1, 2, ie, iT, ibtau), L11(1, 3, ie, iT, ibtau), &
                                                         L11(2, 1, ie, iT, ibtau), L11(2, 2, ie, iT, ibtau), L11(2, 3, ie, iT, ibtau), &
                                                         L11(3, 1, ie, iT, ibtau), L11(3, 2, ie, iT, ibtau), L11(3, 3, ie, iT, ibtau)
            enddo
            write(outfileindex,'(a)') ''
         enddo 
         close(outfileindex)

         outfileindex = outfileindex+1
         write(L12filename, '(6a)')'L12_single','_mu_',trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=L12filename)
         write(outfileindex, '("#",a)')'L12  (S/m)'
         write(outfileindex, "('#column', i7, i13,3i16)")(i, i=1, 5)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T \cdot ps)', 'L12_xx', 'L12_xy', 'L12_xz', 'L12_yx', 'L12_yy', 'L12_yz', 'L12_zx', 'L12_zy', 'L12_zz'
         do iT = 1, NumT
            do ibtau = 1, NBTau
               write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                         L12(1, 1, ie, iT, ibtau), L12(1, 2, ie, iT, ibtau), L12(1, 3, ie, iT, ibtau), &
                                                         L12(2, 1, ie, iT, ibtau), L12(2, 2, ie, iT, ibtau), L12(2, 3, ie, iT, ibtau), &
                                                         L12(3, 1, ie, iT, ibtau), L12(3, 2, ie, iT, ibtau), L12(3, 3, ie, iT, ibtau)
            enddo
            write(outfileindex,'(a)') ''
         enddo
         close(outfileindex)


      enddo
   endif

end subroutine single_band_seebeck

subroutine two_band_seebeck()
!================================================================
! This subroutine is used to calculate the magneto-seebek effect
! of a two-band system with H(k) =
!================================================================
   use wmpi
   use para
   implicit none

   integer :: ie, iT, ik, ikx, iky, ikz, i, j, ib, ibtau, ikT, iband, ierr, &
              ibtau_local
   integer :: knv3, Nslice_Btau, Nslice_Btau_local
   real(dp) :: time_start, time_end, minusdfde, mu, KBT, alpha_tau, hstep
   real(dp) :: EE(2), k(3), kxyz(3), Rxyz(3), vk(2, 3), mass(2, 3), L11_inv(3, 3), L12_tmp(3,3), seebeck(3, 3)
   real(dp), allocatable :: mu_array(:), T_array(:), KBT_array(:), Btau_array(:),&
                            vout(:, :, :), kout(:, :, :), vkbar(:, :, :)
   character(len=100) :: seebeckfilename, muname, L11filename, L12filename

   real(dp), allocatable :: L11(:, :, :, :, :, :), L12(:, :, :, :, :, :), &
                            L11_mpi(:, :, :, :, :, :), L12_mpi(:, :, :, :, :, :)


   real(dp), external :: det3

   allocate( L11_mpi(2, 3, 3, OmegaNum, NumT, NBTau))
   allocate( L12_mpi(2, 3, 3, OmegaNum, NumT, NBTau))
   allocate( L11(2, 3, 3, OmegaNum, NumT, NBTau))
   allocate( L12(2, 3, 3, OmegaNum, NumT, NBTau))
   L11= 0d0
   L12= 0d0
   L11_mpi= 0d0
   L12_mpi= 0d0

   allocate(mu_array(OmegaNum))
   allocate(T_array(NumT))
   allocate(KBT_array(NumT))
   allocate(Btau_array(NBTau))
   allocate(vkbar(2, 3, NBTau))
   vkbar = 0d0

   !> energy
   do ie=1, OmegaNum
      if (OmegaNum>1) then
         mu_array(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
      else
         mu_array= OmegaMin
      endif
   enddo ! ie

   !> temperature
   do iT=1, NumT
       if (NumT>1) then
          T_array(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
       else
          T_array= Tmin
       endif
   enddo ! iT

   !> transform from Kelvin to eV
   !> The SI unit of temperature is the kelvin (K), but using the above relation the electron temperature is often expressed in
   !> terms of the energy unit electronvolt (eV). Each kelvin (1 K) corresponds to 8.6173324(78)×10−5 eV; this factor is the ratio
   !> of the Boltzmann constant to the elementary charge. After version 2.6, we 
   !> adopt the atomic unit
   KBT_array= T_array*8.6173324E-5*eV2Hartree

   !> number of Btau is NBTau, total number of Btau points is Nslice_BTau which does not include Btau = 0
   !> number of Btau points per line is NSlice_Btau_local (not include the origin but include the end point of every line)
   do ibtau=1, NBTau
      if (NBTau>1) then
         Btau_array(ibtau)= (ibtau-1.0d0)/dble(NBTau-1)*BTauMax 
      else
         Btau_array= 0      
      endif
   enddo ! ibtau

   if (NBTau>1) then
      NSlice_Btau= (NBTau-1)*(Nslice_BTau_Max/(NBTau-1))
      Nslice_Btau_local = Nslice_BTau_Max/(NBTau-1)
   else
      NSlice_Btau= 1
      Nslice_Btau_local = 1
   endif
   if (cpuid == 0) write(stdout, '(a, i6, a, i6)')'>> NSlice_Btau=', NSlice_Btau, ' Nslice_Btau_local=', Nslice_Btau_local

   !> (nband, xyz, Nslice_Btau)
   allocate(vout(2, 3, NSlice_Btau+1))
   allocate(kout(2, 3, NSlice_Btau+1))
   vout= 0d0
   kout= 0d0
   
   knv3= Nk1*Nk2*Nk3
   alpha_tau = 15d0
   Rxyz(1) = Origin_cell%Rua(1)
   Rxyz(2) = Origin_cell%Rub(2)
   Rxyz(3) = Origin_cell%Ruc(3)

   mass(1, 1) = m1x
   mass(1, 2) = m1y
   mass(1, 3) = m1z
   mass(2, 1) = m2x
   mass(2, 2) = m2y
   mass(2, 3) = m2z

   call now(time_start) 
   do ik= 1+ cpuid, knv3, num_cpu
      if (cpuid.eq.0.and. mod(ik/num_cpu, 1000).eq.0) then
         call now(time_end) 
         write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
         ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/1000d0
         time_start= time_end
      endif

      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
       + K3D_vec2_cube*(iky-1)/dble(nk2)  &
       + K3D_vec3_cube*(ikz-1)/dble(nk3)

      kxyz = k(1)*Origin_cell%Kua + k(2)*Origin_cell%Kub + k(3)*Origin_cell%Kuc

      EE(1) = (1-cos(k(1)*twopi))/mass(1, 1)+(1-cos(k(2)*twopi))/mass(1, 2)+(1-cos(k(3)*twopi))/mass(1, 3) 
      EE(2) = (1-cos(k(1)*twopi))/mass(2, 1)+(1-cos(k(2)*twopi))/mass(2, 2)+(1-cos(k(3)*twopi))/mass(2, 3) 
      EE(1) = EE(1)*eV2Hartree + Gap_threshold/2
      EE(2) = EE(2)*eV2Hartree - Gap_threshold/2

      do i =1 ,3
         vk(1, i) = sin(k(i)*twopi)*Rxyz(i)/mass(1, i)*eV2Hartree
         vk(2, i) = sin(k(i)*twopi)*Rxyz(i)/mass(2, i)*eV2Hartree
      enddo 
      ! write(stdout, '(a, f12.6)') 'gap', Gap_threshold/2
      
      ! write(stdout, '(a, 4f12.6)') 'band-1: v_k, EE =', vk(1,:), EE(1)
      ! write(stdout, '(a, 4f12.6)') 'band-2: v_k, EE =', vk(2,:), EE(2)
      
      if (NBTau == 1) then
         do iband = 1, 2
            kout(iband, :, 1) = kxyz
            vout(iband, :, 1) = vk(iband, :)
         enddo
      else 
         do iband = 1, 2
            hstep = alpha_tau*BTauMax/Nslice_Btau
            call ode45(kout(iband, :, :), kxyz, hstep, Nslice_Btau+1, mass(iband, 1), mass(iband, 2), mass(iband, 3))
            do ibtau_local = 1, Nslice_Btau+1
               do i =1,3
                  vout(iband, i, ibtau_local) = sin(kout(iband, i, ibtau_local)*Rxyz(i))*Rxyz(i)/mass(iband, i)*eV2Hartree
               enddo
            enddo
         enddo 
      endif
      

      ! outfileindex = outfileindex+1
      ! open(unit=outfileindex, file='kout.dat')
      ! write(outfileindex, '("#", 3a)')'kx', 'ky', 'kz'
      ! do ibtau_local = 1, Nslice_Btau
      !    write(outfileindex, '(3f12.6)') kout(1, :, ibtau_local)
      ! enddo
      ! write(outfileindex,'(a)') ''
      ! do ibtau_local = 1, Nslice_Btau
      !    write(outfileindex, '(3f12.6)') kout(2, :, ibtau_local)
      ! enddo
      ! close(outfileindex)

      do iband = 1 ,2
         vkbar(iband, :, 1) = vk(iband, :)
      enddo
      if (NBTau >1) then
         do ibtau = 2, NBTau 
            do iband = 1, 2
               vkbar(iband, :, ibtau) = alpha_tau/(Nslice_Btau_local*(ibtau-1))* &
                           (vout(iband, :, 1)+vout(iband, :, Nslice_Btau_local*(ibtau-1)+1)*exp(-alpha_tau))/2
               do ibtau_local = 2, (Nslice_Btau_local*(ibtau-1))
                  vkbar(iband, :, ibtau) = vkbar(iband, :, ibtau) + alpha_tau/(Nslice_Btau_local*(ibtau-1))* &
                     vout(iband, :, ibtau_local)*exp(-alpha_tau*ibtau_local/(Nslice_Btau_local*(ibtau-1)))
               enddo
            enddo
         enddo
      endif


      !> calculate df/de for each k point and each band
      do ibtau = 1, Nbtau
         do ikt=1, NumT
            KBT= KBT_array(ikt)
            do ie=1, OmegaNum
               mu= mu_array(ie)
               do iband = 1, 2
                  call minusdfde_calc_single(EE(iband), KBT, mu, minusdfde)
                  do i = 1, 3
                     do j = 1, 3
                        L11_mpi(iband, i, j, ie, ikt, ibtau)= L11_mpi(iband, i, j, ie, ikt, ibtau)+ &
                                             minusdfde*vk(iband, i)*vkbar(iband, j, ibtau)/dble(knv3)
                        L12_mpi(iband, i, j, ie, ikt, ibtau)= L12_mpi(iband, i, j, ie, ikt, ibtau)+ &
                                             minusdfde*vk(iband, i)*vkbar(iband, j, ibtau)* &
                                             (EE(iband)-mu)/dble(knv3)
                     enddo
                  enddo 
               enddo
            enddo
         enddo 
      enddo 
      

   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(L11_mpi, L11, size(L11), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(L12_mpi, L12, size(L12), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   L11 = L11_mpi
   L12 = L12_mpi
#endif
   L11 = L11*Echarge**2/Bohr_radius*Hartree_SI/hbar**2/Origin_cell%CellVolume & 
         /kCubeVolume*Origin_cell%ReciprocalCellVolume
   L12 = L12*Echarge**2/Bohr_radius*Hartree_SI**2/hbar**2/Origin_cell%CellVolume &
         /kCubeVolume*Origin_cell%ReciprocalCellVolume

   
   if(cpuid .eq. 0) then
      do ie=1, OmegaNum
         write(muname, '(f12.3)')mu_array(ie)/eV2Hartree
         outfileindex = outfileindex+1
         write(seebeckfilename, '(6a)')'seebeck_twoband','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=seebeckfilename)
         write(outfileindex, '("#",a)')'Seebeck coefficient  (V/K)'
         write(outfileindex, "('#column', i7, i13,9i16)")(i, i=1, 11)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T\cdotps)', 'S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_zx', 'S_zy', 'S_zz'
         do iT = 1, NumT
            do ibtau = 1, NBTau
               L11_inv = L11(1, :, :, ie, iT, ibtau) + L11(2, :, :, ie, iT, ibtau)
               if (abs(det3(L11_inv))>eps6) then
                  call inv_r(3, L11_inv)
                  L12_tmp = L12(1, :, :, ie, iT, ibtau) + L12(2, :, :, ie, iT, ibtau)
                  seebeck = matmul(L11_inv, -L12_tmp/Echarge/T_array(iT))
                  write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                            seebeck(1, 1), seebeck(1, 2), seebeck(1, 3), seebeck(2, 1), seebeck(2, 2), seebeck(2, 3), &
                                                            seebeck(3, 1), seebeck(3, 2), seebeck(3, 3)
               else 
                  write(outfileindex, '(a)') 'No solution because the conductivity tensor is singular'
               endif
            enddo
         write(outfileindex,'(a)') ''
         enddo
         close(outfileindex)

         outfileindex = outfileindex+1
         write(L11filename, '(6a)')'L11_twoband','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=L11filename)
         write(outfileindex, '("#",a)')'L11  '
         write(outfileindex, "('#column', i7, i13,9i16)")(i, i=1, 11)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T\cdotps)', 'L11_xx', 'L11_xy', 'L11_xz', 'L11_yx', 'L11_yy', 'L11_yz', 'L11_zx', 'L11_zy', 'L11_zz'
         do iband = 1, 2   
            write(outfileindex, '("#", a, i5)') 'iband = ', iband
            do iT = 1, NumT
               do ibtau = 1, NBTau
                  write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                            L11(iband, 1, 1, ie, iT, ibtau), L11(iband, 1, 2, ie, iT, ibtau), L11(iband, 1, 3, ie, iT, ibtau), &
                                                            L11(iband, 2, 1, ie, iT, ibtau), L11(iband, 2, 2, ie, iT, ibtau), L11(iband, 2, 3, ie, iT, ibtau), &
                                                            L11(iband, 3, 1, ie, iT, ibtau), L11(iband, 3, 2, ie, iT, ibtau), L11(iband, 3, 3, ie, iT, ibtau)
               enddo
               write(outfileindex,'(a)') ''
            enddo 
            write(outfileindex,'(a)') ''
         enddo 
         close(outfileindex)

         outfileindex = outfileindex+1
         write(L12filename, '(6a)')'L12_twoband','_mu_',&
            trim(adjustl(muname)),'eV.dat'
         open(unit=outfileindex, file=L12filename)
         write(outfileindex, '("#",a)')'L12  '
         write(outfileindex, "('#column', i7, i13,9i16)")(i, i=1, 11)
         write(outfileindex,'("#",a6, a20, 11a16)')'T (K)', 'Btau (T\cdotps)', 'L12_xx', 'L12_xy', 'L12_xz', 'L12_yx', 'L12_yy', 'L12_yz', 'L12_zx', 'L12_zy', 'L12_zz'
         do iband = 1, 2   
            write(outfileindex, '("#", a, i5)') 'iband = ', iband
            do iT = 1, NumT
               do ibtau = 1, NBTau
                  write(outfileindex,'(F14.6,F13.3,9E16.6)')T_array(iT), BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                                                            L12(iband, 1, 1, ie, iT, ibtau), L12(iband, 1, 2, ie, iT, ibtau), L12(iband, 1, 3, ie, iT, ibtau), &
                                                            L12(iband, 2, 1, ie, iT, ibtau), L12(iband, 2, 2, ie, iT, ibtau), L12(iband, 2, 3, ie, iT, ibtau), &
                                                            L12(iband, 3, 1, ie, iT, ibtau), L12(iband, 3, 2, ie, iT, ibtau), L12(iband, 3, 3, ie, iT, ibtau)
               enddo
               write(outfileindex,'(a)') ''
            enddo 
            write(outfileindex,'(a)') ''
         enddo 
         close(outfileindex)

      enddo
   endif

end subroutine two_band_seebeck

subroutine ode45(y, y0, h, num_steps, mx, my, mz)
   use para
   implicit none
   !> f(t, y) = y
   integer, intent(in) :: num_steps
   real(dp), intent(inout) :: y(3, num_steps)
   real(dp), intent(in) :: y0(3), h, mx, my, mz
   real(dp) :: x(num_steps)
   integer :: i
   real(dp) :: k1(3), k2(3), k3(3), k4(3)

   ! allocate(y(num_steps))
   ! allocate(x(num_steps))
   y(:, 1) = y0
   x(1) = 0d0

   ! outfileindex = outfileindex+1
   ! open(unit=outfileindex, file='ode45.dat')
   ! write(outfileindex, '("#", 4a)') 'k1', 'k2', 'k3', 'k4'

   do i = 2, num_steps
      x(i) = x(i-1) + h
      call func_ode(x(i-1), y(:, i-1), k1, mx, my, mz)
      call func_ode(x(i-1) + h/2, y(:, i-1) + h/2*k1, k2, mx, my, mz)
      call func_ode(x(i-1) + h/2, y(:, i-1) + h/2*k2, k3, mx, my, mz)
      call func_ode(x(i-1) + h, y(:, i-1) + h*k3, k4, mx, my, mz)
      y(:, i) = y(:, i-1) + h/6*(k1 + 2*k2 + 2*k3 + k4)
      ! write(outfileindex, '(4E12.4)') k1, k2, k3, k4
   enddo

   
   ! write(outfileindex, '("#", 3a)')'x', 'y', 'y_exact'
   ! do i = 1, num_steps
   !    write(outfileindex, '(3E12.4)') x(i), y(i), sqrt(2*x(i)+1)
   ! enddo
end subroutine ode45

subroutine func_ode(x, y, f, mx, my, mz)
   use para
   implicit none
   real(dp), intent(in) :: x, y(3), mx, my, mz
   real(dp), intent(out) :: f(3)
   real(dp) :: Rxyz(3)

   ! f = x**2/y
   ! f = y - 2*x/y
   Rxyz(1) = Origin_cell%Rua(1)
   Rxyz(2) = Origin_cell%Rub(2)
   Rxyz(3) = Origin_cell%Ruc(3) 

   !> dk/dBtau = -ev \times hat(B)
   !> e=1, dBtau is minus number because it represents the past time
   !> v = sin(kR)*R/mass
   f(1) = sin(y(2)*Rxyz(2))*Rxyz(2)/my*Bdirection(3)*eV2Hartree-sin(y(3)*Rxyz(3))*Rxyz(3)/mz*Bdirection(2)*eV2Hartree
   f(2) = sin(y(3)*Rxyz(3))*Rxyz(3)/mz*Bdirection(1)*eV2Hartree-sin(y(1)*Rxyz(1))*Rxyz(1)/mx*Bdirection(3)*eV2Hartree
   f(3) = sin(y(1)*Rxyz(1))*Rxyz(1)/mx*Bdirection(2)*eV2Hartree-sin(y(2)*Rxyz(2))*Rxyz(2)/my*Bdirection(1)*eV2Hartree

end subroutine func_ode

