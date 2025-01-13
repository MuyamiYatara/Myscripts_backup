

program main
   implicit none
   !@ Here set the Number of wannier orbitals without spin 
   integer :: num_of_orbital= 31

   call generate_HmnR(num_of_orbital)
   
end program main


subroutine generate_HmnR(Num_wann)
   !>> Read in the tight-binding model from wannier90_hr.dat
   !>  The format is defined by the wannier90 software
   ! Constructed by quansheng wu 4/2/2010
   !
   ! Yifei Guan added the sparse hr file parsing June/2018
   ! License: GPL V3


   implicit none

   integer, parameter :: dp = kind(1.0d0)
   character*4 :: c_temp
   integer :: i, j, ir, ia, io 
   integer :: i1, i2, i3, i4, i5
   integer :: n, m, ir0, uo1, uo2, do1, do2
   integer :: nwann, Nrpts, file_unit
   integer, allocatable :: irvec(:, :)
   integer, allocatable :: ndegen(:)

   complex(dp), allocatable :: HmnR_up(:, :, :)
   complex(dp), allocatable :: HmnR_dn(:, :, :)
   complex(dp), allocatable :: HmnR_all(:, :, :)
   
   
   integer, intent(in) :: Num_wann
   real(dp) :: static_potential
   real(dp) :: tot, rh, ih

   
   open(12, file="up_hr.dat", status='OLD')
   open(13, file="down_hr.dat", status='OLD')

      !> for Normal HmnR obtained from Wannier90 or sparse HmnR

      !> skip a comment line
      read(12, *)
      read(13, *)

      !> number of Wannier orbitals in the hr file
      nwann=0
      read(12, *)nwann
      read(13, *)
      if (nwann==0) then
         stop "ERROR : num_wann is zero in hr file"
      endif
  
      !> number of lattice vectors taken into account
      read(12, *)Nrpts
      read(13, *)
      if (.not. allocated(HmnR_all)) allocate(HmnR_all(Num_wann*2,Num_wann*2,nrpts))
      allocate(ndegen(Nrpts))
      allocate(irvec(3,Nrpts))
      ndegen=1

      !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
      read(12, *)(ndegen(i), i=1, Nrpts)
      read(13, *)(ndegen(i), i=1, Nrpts)
      ir=0
      
      !@ generate the total hr dat file, up and down sub space remains same as the original file, non-diagonal terms which have different spin index are all set to zero.
      
      !@ we use the spin order follow the convention of wannier 1.2 , which is up up dn dn.

      !@ directly get the value from up and down spin hr dat file
      do ir=1,Nrpts
         do n=1,nwann
            do m=1,nwann
               read(12,*)i1, i2, i3, uo1, uo2, rh, ih
               irvec(1,ir)=i1
               irvec(2,ir)=i2
               irvec(3,ir)=i3
               HmnR_all(uo1, uo2, ir)=dcmplx(rh,ih)
               read(13,*)i1, i2, i3, do1, do2, rh, ih
               HmnR_all(do1+Num_wann, do2+Num_wann, ir)=dcmplx(rh,ih)
               HmnR_all(uo1, do2+Num_wann, ir)=dcmplx(0d0,0d0)
               HmnR_all(do1+Num_wann, uo2, ir)=dcmplx(0d0,0d0)
            end do
         enddo
      enddo

      !@ write the file of total hr dat
      file_unit = 66
      open(file_unit, file="spin_wo_soc_hr.dat", form='formatted', status='unknown')
      write(file_unit,*) "combined by up and down spin hr.dat"
      write(file_unit,*) Num_wann*2
      write(file_unit,*) Nrpts
      write(file_unit,'(15I5)') (ndegen(i),i=1,Nrpts)
      do ir=1,Nrpts
         do i=1,Num_wann*2
            do j=1,Num_wann*2
               write( file_unit,'(5I5,2F12.6)') irvec(:, ir), j, i, HmnR_all(j,i,ir)
            end do
         end do
      end do

end subroutine generate_HmnR

