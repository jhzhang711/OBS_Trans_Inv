!=======================================================================
!   SEIS_FILO: 
!   SEISmological tools for Flat Isotropic Layered structure in the Ocean
!   Copyright (C) 2019 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
module cls_vmodel
  use cls_line_text
  use, intrinsic :: iso_fortran_env
  implicit none 
  
  type vmodel
     private
     integer :: nlay
     double precision, allocatable :: vp(:), vs(:), rho(:), h(:)
     double precision, allocatable :: qp(:), qs(:)
   
   contains
     ! setter
     procedure :: set_nlay => vmodel_set_nlay
     procedure :: set_vp => vmodel_set_vp
     procedure :: set_vs => vmodel_set_vs
     procedure :: set_rho => vmodel_set_rho
     procedure :: set_h => vmodel_set_h
     procedure :: set_qp => vmodel_set_qp
     procedure :: set_qs => vmodel_set_qs
     ! getter
     procedure :: get_nlay => vmodel_get_nlay
     procedure :: get_vp => vmodel_get_vp
     procedure :: get_vs => vmodel_get_vs
     procedure :: get_rho => vmodel_get_rho
     procedure :: get_h => vmodel_get_h
     procedure :: get_qp => vmodel_get_qp
     procedure :: get_qs => vmodel_get_qs
     ! utilities
     procedure :: read_file => vmodel_read_file
     procedure :: set_example_ocean => vmodel_set_example_ocean
     procedure :: set_example_land => vmodel_set_example_land
     procedure :: display => vmodel_display
     procedure :: vp2rho_brocher => vmodel_vp2rho_brocher
     procedure :: vp2vs_brocher => vmodel_vp2vs_brocher
     procedure :: vs2vp_brocher => vmodel_vs2vp_brocher
     procedure :: sphere2flat => vmodel_sphere2flat
     procedure :: setvmodel => vmodel_setvmodel
  end type vmodel
  
  interface vmodel
     module procedure init_vmodel
  end interface vmodel
  
contains
  
  !---------------------------------------------------------------------

  type(vmodel) function init_vmodel()
    
    init_vmodel%nlay = -999
    
    return 
  end function init_vmodel
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_set_nlay(self, nlay)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: nlay

    ! Reset if already set
    if (self%nlay > 0) then
       deallocate(self%vp, self%vs, self%rho, self%h)
    end if
    
    self%nlay = nlay
    allocate(self%vp(nlay), self%vs(nlay), self%rho(nlay), self%h(nlay))
    allocate(self%qp(nlay))
    allocate(self%qs(nlay))
    
    self%vp(1:nlay) = -999.d0
    self%vs(1:nlay) = -999.d0
    self%rho(1:nlay) = -999.d0
    self%h(1:nlay) = -999.d0
    self%qp(1:nlay) = 1.0d20
    self%qs(1:nlay) = 1.0d20


    return 
  end subroutine vmodel_set_nlay

  !---------------------------------------------------------------------

  subroutine vmodel_set_vp(self, i, vp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: vp
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (set_vp)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    !if (vp < 0) then
    !   write(0,*) "ERROR: invalid Vp (set_vp)"
    !   write(0,*) "     : Vp=", vp
    !   stop
    !end if
    self%vp(i) = vp

    return 
  end subroutine vmodel_set_vp
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_vs(self, i, vs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: vs
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (set_vs)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    !if (vs < 0 .and. i /= 1) then
    !   write(0,*) "ERROR: invalid Vs (set_vs)"
    !   write(0,*) "     : Vs=", vs
    !   stop
    !end if
    self%vs(i) = vs

    return 
  end subroutine vmodel_set_vs
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_set_rho(self, i, rho)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: rho
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (set_rho)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    !if (rho < 0) then
    !   write(0,*) "ERROR: invalid density (set_rho)"
    !   write(0,*) "     : rho=", rho
    !   stop
    !end if

    self%rho(i) = rho

    return 
  end subroutine vmodel_set_rho
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_h(self, i, h)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: h
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (set_h)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    !if (h < 0.d0 .and. i /= self%nlay) then
    !   write(0,*) "ERROR: invalid thickness (set_h)"
    !   write(0,*) "     : h=", h, ", i=", i
    !   stop
    !end if
    self%h(i) = h

    return 
  end subroutine vmodel_set_h

  !---------------------------------------------------------------------

  subroutine vmodel_set_qp(self, i, qp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: qp
    
    self%qp(i) = qp

    return 
  end subroutine vmodel_set_qp
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_qs(self, i, qs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: qs
    
    self%qs(i) = qs

    return 
  end subroutine vmodel_set_qs
  
  !---------------------------------------------------------------------

  
  integer function vmodel_get_nlay(self) result(nlay)
    class(vmodel), intent(inout) :: self
    
    nlay = self%nlay
    
    return 
  end function vmodel_get_nlay

  !---------------------------------------------------------------------
  
  double precision function vmodel_get_vp(self, i) result(vp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (get_vp)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    vp = self%vp(i)
    
    return 
  end function vmodel_get_vp
  
  !---------------------------------------------------------------------
  
  double precision function vmodel_get_vs(self, i) result(vs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (get_vs)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    vs = self%vs(i)
    
    return 
  end function vmodel_get_vs
  
  !---------------------------------------------------------------------

  double precision function vmodel_get_rho(self, i) result(rho)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (get_rho)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    rho = self%rho(i)
    
    return 
  end function vmodel_get_rho

  !---------------------------------------------------------------------
  
  double precision function vmodel_get_h(self, i) result(h)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    !if(i < 1 .or. i > self%nlay) then
    !   write(0,*) "ERROR: out of range (get_h)"
    !   write(0,*) "     : i=", i
    !   stop
    !end if
    h = self%h(i)
    
    return 
  end function vmodel_get_h
  
  !---------------------------------------------------------------------

  double precision function vmodel_get_qp(self, i) result(qp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    qp = self%qp(i)
    
    return 
  end function vmodel_get_qp
  
  !---------------------------------------------------------------------

  double precision function vmodel_get_qs(self, i) result(qs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    qs = self%qs(i)
    
    return 
  end function vmodel_get_qs
  
  !---------------------------------------------------------------------

  subroutine vmodel_display(self, io, is_attenuative)
    class(vmodel), intent(inout) :: self
    integer, intent(in), optional :: io
    logical, intent(in), optional :: is_attenuative
    integer :: i,  i_unit
    
    i_unit = output_unit
    if (present(io)) then
       i_unit = io
    end if

    block 
      logical :: flag 
      if (present(is_attenuative)) then
         flag = is_attenuative
      else
         flag = .false.
      end if
      do i = 1, self%nlay
         if (flag) then
            write(i_unit,'(I5, 6F10.3)') &
                 & i, self%vp(i), self%vs(i), &
                 & self%rho(i), self%h(i), &
                 & self%qp(i), self%qs(i)
         else
            write(i_unit,'(I5, 4F10.3)') &
                 & i, self%vp(i), self%vs(i), &
                 & self%rho(i), self%h(i)
         end if
      end do
    end block
    
    return 
  end subroutine vmodel_display
  
  !---------------------------------------------------------------------

  subroutine vmodel_read_file(self, vmod_in, is_attenuative)
    class(vmodel), intent(inout) :: self
    character(*), intent(in) :: vmod_in
    logical, intent(in), optional :: is_attenuative
    integer :: ierr, io, i, nlay
    type(line_text) :: lt
    character(line_max) :: line
    
    
    write(*,*)"Reading velocity model from ", trim(vmod_in)
    
    open(newunit = io, file = vmod_in, iostat = ierr, &
         & status = 'old')
    if (ierr /= 0) then
       write(0, *)"ERROR: cannot open ", trim(vmod_in)
       !call mpi_finalize(ierr)
       ! not necesarry call mpi_finalize here because this routine 
       ! is only used for forward computation
       stop
    end if
    
    
    ! Get # of layers
    do 
       read(io, '(a)')line
       lt = line_text(line, ignore_space = .false.)
       line = lt%get_line()
       if (len_trim(line) == 0) cycle
       read(line, *)nlay
       exit
    end do
    call self%set_nlay(nlay) ! <= allocate Vp, Vs, Rho, H, qp, qs
    
    block
      logical :: flag
      if (present(is_attenuative)) then
         flag = is_attenuative
      else 
         flag = .false.
      end if

      do i = 1, self%nlay
         do 
            read(io, '(a)')line
            lt = line_text(line, ignore_space = .false.)
            line = lt%get_line()
            if (len_trim(line) == 0) cycle
            if (flag) then
               read(line, *) &
                    & self%vp(i), self%vs(i), self%rho(i), &
                    & self%h(i), self%qp(i), self%qs(i)
            else
               read(line, *) &
                    & self%vp(i), self%vs(i), self%rho(i), self%h(i)
            end if
            exit
         end do
      end do
      close(io)
      call self%display(is_attenuative=flag)
    end block


    
    write(*,*)


    return 
  end subroutine vmodel_read_file

  !---------------------------------------------------------------------
  
  subroutine vmodel_set_example_ocean(self)
    class(vmodel), intent(inout) :: self
    
    call self%set_nlay(3)
    
    ! Ocean
    call self%set_vp(1, 1.5d0)
    call self%set_vs(1, -1.d0)
    call self%set_rho(1, 1.d0)
    call self%set_h(1, 3.d0)

    ! Sediment
    call self%set_vp(2, 1.6d0)
    call self%vp2vs_brocher(2)
    call self%vp2rho_brocher(2)
    call self%set_h(2, 1.d0)

    ! Basement
    call self%set_vp(3, 3.d0)
    call self%vp2vs_brocher(3)
    call self%vp2rho_brocher(3)
    call self%set_h(3, 100.d0)

    return 
  end subroutine vmodel_set_example_ocean
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_example_land(self)
    class(vmodel), intent(inout) :: self
    
    call self%set_nlay(3)
    
    ! Ocean
    call self%set_vp(1, 0.255d0)
    call self%set_vs(1, 0.150d0)
    call self%set_rho(1, 1.00d0)
    call self%set_h(1, 0.002d0)

    ! Sediment
    call self%set_vp(2, 0.340d0)
    call self%set_vs(2, 0.200d0)
    call self%set_rho(2, 1.00d0)
    call self%set_h(2, 0.002d0)

    ! Basement
    call self%set_vp(3, 0.510d0)
    call self%set_vs(3, 0.300d0)
    call self%set_rho(3, 1.00d0)
    call self%set_h(3, 100.d0)

    return 
  end subroutine vmodel_set_example_land
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_vp2rho_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4, a5
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vp2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vp(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1
    a5 = a4 * a1

    self%rho(i) = &
         & 1.6612d0 * a1 - 0.4721d0 * a2 + 0.0671d0 * a3 - &
         & 0.0043d0 * a4 + 0.000106d0 * a5
    
    return 
  end subroutine vmodel_vp2rho_brocher
    
  !---------------------------------------------------------------------

  subroutine vmodel_vp2vs_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vs2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vp(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1

    self%vs(i) = &
         & 0.7858d0 - 1.2344d0 * a1 + 0.7949d0 * a2 - &
         & 0.1238d0 * a3 + 0.0064d0 * a4
    
    return 
  end subroutine vmodel_vp2vs_brocher
    
  !---------------------------------------------------------------------

  subroutine vmodel_vs2vp_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vs2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vs(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1

    self%vp(i) = &
         & 0.9409d0 + 2.0947d0 * a1 - 0.8206d0 * a2 + &
         & 0.2683d0 * a3 - 0.0251d0 * a4
    
    return 
  end subroutine vmodel_vs2vp_brocher
    
  !---------------------------------------------------------------------
  
  subroutine vmodel_sphere2flat(self, r_earth, vm_out)
    class(vmodel), intent(in) :: self
    class(vmodel), intent(out) :: vm_out
    double precision, intent(in) :: r_earth
    double precision :: r, ztop, zbot
    integer :: i

    call vm_out%set_nlay(self%nlay)
    r = r_earth
    do i= 1, self%nlay
       ztop = r_earth * log(r_earth / r)
       r = r - 0.5d0 * self%h(i) ! Get radius at layer center
       vm_out%vp(i) = self%vp(i) * r_earth / r
       !write(*,*)i, self%vp(i), vm_out%vp(i)
       vm_out%vs(i) = self%vs(i) * r_earth / r
       vm_out%rho(i) = self%rho(i) * r / r_earth
       r = r - 0.5d0 * self%h(i)
       zbot = r_earth * log(r_earth / r)
       vm_out%h(i) = zbot - ztop
       vm_out%qp(i) = self%qp(i)
       vm_out%qs(i) = self%qs(i)


    end do

    return 
  end subroutine vmodel_sphere2flat

  !---------------------------------------------------------------------
  subroutine vmodel_setvmodel(self,layer,vp,vs,rho,thk,vm_out)
    class(vmodel), intent(in) :: self
    class(vmodel), intent(out) :: vm_out
    integer, intent(in) :: layer
    real, intent(in) :: vp(layer),vs(layer),rho(layer),thk(layer)
    integer :: i
 
	
        call vm_out%set_nlay(layer)
    	do i = 1, layer
    		call vm_out%set_vp(i, real(vp(i),8))
		call vm_out%set_vs(i, real(vs(i),8))
		call vm_out%set_rho(i, real(rho(i),8))
		call vm_out%set_h(i, real(thk(i),8))
		call vm_out%set_qp(i, 1.0d20)
		call vm_out%set_qs(i, 1.0d20)
		!write(*,*) i,vm_out%h(i),vm_out%vp(i),vm_out%vs(i),vm_out%rho(i)
    	end do
    return
  end subroutine vmodel_setvmodel

  !---------------------------------------------------------------------

end module cls_vmodel
