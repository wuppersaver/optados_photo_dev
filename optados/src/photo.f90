!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!
! This file is part of OptaDOS
!
! OptaDOS - For obtaining electronic structure properties based on
!             integrations over the Brillouin zone
! Copyright (C) 2011  Andrew J. Morris,  R. J. Nicholls, C. J. Pickard
!                         and J. R. Yates
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!===============================================================================
module od_photo
  !! This is the module for calculating the photoemission.
  use od_constants, only: dp

  implicit none
  private
  public :: photo_calculate

  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: pdos_weights_atoms
  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: pdos_weights_boxes
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: projected_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: foptical_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :) :: weighted_jdos
  real(kind=dp), allocatable, public, dimension(:, :) :: absorp_layer
  real(kind=dp), allocatable, public, dimension(:, :, :) :: pdos_weights_k_band
  real(kind=dp), allocatable, public, save :: E(:)
  real(kind=dp), allocatable, dimension(:, :, :) :: imfp_val
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: electron_esc
  real(kind=dp), dimension(:, :), allocatable :: I_layer
  real(kind=dp), allocatable, dimension(:, :) :: reflect_photo
  real(kind=dp), allocatable, dimension(:, :) :: absorp_photo
  real(kind=dp), allocatable, dimension(:, :) :: refract
  real(kind=dp), allocatable, dimension(:)  :: reflect
  real(kind=dp), allocatable, dimension(:) :: absorp
  real(kind=dp), dimension(:), allocatable :: thickness_atom
  real(kind=dp), dimension(:), allocatable :: thickness_layer
  real(kind=dp), dimension(:), allocatable :: volume_atom
  real(kind=dp)                            :: box_height
  real(kind=dp)                            :: box_volume
  integer, dimension(:), allocatable       :: box_atom
  integer, dimension(:), allocatable       :: atoms_per_box
  integer                                  :: num_boxes
  real(kind=dp)                            :: photo_slab_volume
  real(kind=dp)                            :: slab_middle_ref
  real(kind=dp)                            :: cell_area
  real(kind=dp), dimension(:), allocatable :: atom_imfp
  real(kind=dp), dimension(:, :, :), allocatable :: band_imfp
  integer :: first_atom_second_l, last_atom_secondlast_l
  real(kind=dp), dimension(:), allocatable :: boxes_top_z_coord
  real(kind=dp), dimension(:, :), allocatable :: new_atom_coordinates
  real(kind=dp), allocatable, dimension(:, :, :, :) :: phi_arpes
  real(kind=dp), allocatable, dimension(:, :, :, :) :: theta_arpes
  real(kind=dp), allocatable, dimension(:, :, :, :) :: theta_arpes_internal
  real(kind=dp), allocatable, dimension(:, :, :, :) :: E_kinetic
  real(kind=dp), allocatable, dimension(:, :, :, :) :: E_transverse
  real(kind=dp), allocatable, dimension(:) :: t_energy
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: weighted_temp
  integer :: max_energy = -1
  real(kind=dp), allocatable, dimension(:, :, :, :)    :: qe_osm
  real(kind=dp), allocatable, dimension(:, :, :, :)    :: te_osm
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: qe_tsm
  real(kind=dp), allocatable, dimension(:, :, :, :) :: te_tsm
  real(kind=dp) :: mean_te
  real(kind=dp) :: total_qe
  real(kind=dp), allocatable, dimension(:) :: layer_qe
  integer, dimension(:), allocatable :: atom_order
  integer, dimension(:), allocatable :: atoms_per_layer
  real(kind=dp) :: work_function_eff
  real(kind=dp) :: evacuum
  real(kind=dp) :: evacuum_eff
  real(kind=dp) :: total_field_emission
  real(kind=dp), allocatable, dimension(:, :, :) :: field_emission
  integer, allocatable, dimension(:) :: layer
  integer :: N_geom
  integer :: max_atoms
  integer :: max_layer
  real(kind=dp) :: q_weight
  ! Added by Felix Mildner, 12/2022 and later
  integer, allocatable, dimension(:)  :: index_energy
  integer                             :: number_energies, current_energy_index, current_photo_energy_index
  real(kind=dp)                       :: temp_photon_energy, time_a, time_b
  integer, allocatable, dimension(:, :):: min_index_unocc
  ! The Free Electron Matrix (FEM) elements are calculated for a specific E_fermi offset, workfct and photon
  ! energies in Castep. Thus we must read it from the file and ensure they are compatible with the parameters
  ! used for the OptaDOS run.
  ! fem_energy_info: energy_count, energy_min, energy_step, energy_fermi, energy_workfct
  integer                             :: energy_count
  real(kind=dp)                       :: energy_min, energy_step, energy_fermi, energy_workfct
  logical                             :: new_geom_choice = .True. ! hard coded choice of geometry definition
  ! Allowing debug output makes the calculation a lot slower since a very hot if statement is not optimised out druing compilation.
  logical                             :: enable_debug_output = .True. ! hard coded extra printing
contains

  subroutine photo_calculate
    !! Main subroutine calling all the other subroutine steps.
    use od_electronic, only: elec_dealloc_optical, elec_pdos_read, elec_read_optical_mat, &
                             efermi, efermi_set, elec_read_foptical_mat, elec_dealloc_pdos
    use od_jdos_utils, only: jdos_utils_calculate, setup_energy_scale
    use od_comms, only: on_root
    use od_parameters, only: photo_work_function, photo_model, photo_elec_field, photo_output, photo_photon_sweep, &
                             photo_photon_min, jdos_spacing, photo_photon_energy, iprint, devel_flag
    use od_dos_utils, only: dos_utils_set_efermi, dos_utils_calculate_at_e, dos_utils_deallocate
    use od_io, only: stdout, io_error, io_time
    use od_pdos, only: pdos_calculate

    implicit none

    integer :: i

    if (on_root) then
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                             Photoemission Calculation                      +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    if (index(devel_flag, 'geom_analysis') > 0) then
      write (stdout, '(1x,a78)') '+           Only performing the analysis of the supplied geometry!           +'
      call analyse_geometry
      call elec_pdos_read
      call make_pdos_weights_atoms
      call elec_dealloc_pdos
      return
    end if

    if (.not. efermi_set) then
      call dos_utils_set_efermi
      call dos_utils_deallocate
    end if

    ! Identify layers
    call analyse_geometry
    call calc_band_info
    call calc_photon_energies

    if (.not. index(photo_model, 'ds_like_pe') > 0) then
      call elec_read_optical_mat
      ! THIS PART COMES FROM THE PDOS MODULE
      ! read in the pdos weights
      call elec_pdos_read
      call make_pdos_weights_atoms
      call elec_dealloc_pdos

      ! Calculate the optical properties of the slab
      call calc_photo_optics
      if (index(devel_flag, 'output_ome_itof') > 0) return

      call calc_absorp_layer
    end if

    ! Electric field and field emission
    if (photo_elec_field .gt. 0.0_dp) then
      call effect_wf
      call calc_field_emission
    else
      evacuum_eff = efermi + photo_work_function
      work_function_eff = photo_work_function
    end if

    if (photo_photon_sweep) then
      do i = 1, number_energies
        time_a = io_time()
        temp_photon_energy = photo_photon_min + (i - 1)*jdos_spacing
        if (on_root) write (stdout, '(1x,a50,f8.4,a20)') '+--------------- Starting Photoemission Sweep with', temp_photon_energy, &
          ' eV ---------------+'
        current_photo_energy_index = i
        current_energy_index = index_energy(i)
        ! Calculate the photoemission angles theta/phi and transverse energy
        ! We will not need that when calculating the simplified model
        if (.not. index(photo_model, 'ds_like_pe') > 0) then
          call calc_angle

          !Calculate the electron escape length
          call calc_electron_esc

          call bulk_emission
        end if

        !Calculate the QE
        !Three-step-model
        if (index(photo_model, '3step') > 0) then
          call calc_three_step_model
          !One-step-model
        elseif (index(photo_model, '1step') > 0) then
          !Read the one-step matrix elements
          if (.not. allocated(foptical_matrix_weights)) call elec_read_foptical_mat
          !Calculate the one-step optical matrix
          call make_foptical_weights
          !Calculate QE
          call calc_one_step_model
          ! Simplified DS like model
        elseif (index(photo_model, 'ds_like_pe') > 0) then
          call calc_ds_like_model
        end if

        !Weight the contribution of each electron
        !to the transverse energy spread according to their QE
        call weighted_mean_te

        call write_qe_data
        ! Only call the binding energy gaussian broadening and file printing if necessary
        if (.not. index(photo_output, 'off') > 0) then
          !Broaden ouputs using a gaussian function
          if (index(photo_output, 'e_bind') > 0) then
            call binding_energy_broadening
          end if
          ! Write either a binding energy output with after Gaussian broadening or the reduced QE tensor
          call write_qe_output_files
        end if
        time_b = io_time()
        if (on_root .and. iprint > 1) then
          write (stdout, '(1x,a44,15x,f11.3,a8)') '+ Time to calculate Photoemission sweep step', time_b - time_a, ' (sec) +'
        end if
      end do
    else
      temp_photon_energy = photo_photon_energy
      current_photo_energy_index = 1
      current_energy_index = index_energy(1)
      ! Calculate the photoemission angles theta/phi and transverse energy
      ! We will not need that when calculating the simplified model
      if (.not. index(photo_model, 'ds_like_pe') > 0) then
        call calc_angle

        !Calculate the electron escape length
        call calc_electron_esc

        call bulk_emission
      end if

      !Calculate the QE
      !Three-step-model
      if (index(photo_model, '3step') > 0) then
        call calc_three_step_model
        !One-step-model
      elseif (index(photo_model, '1step') > 0) then
        !Read the one-step matrix elements
        if (.not. allocated(foptical_matrix_weights)) call elec_read_foptical_mat
        !Calculate the one-step optical matrix
        call make_foptical_weights
        !Calculate QE
        call calc_one_step_model
        ! Simplified DS like model
      elseif (index(photo_model, 'ds_like_pe') > 0) then
        call calc_ds_like_model
      end if

      !Weight the contribution of each electron
      !to the transverse energy spread according to their QE
      call weighted_mean_te
      call write_qe_data

      ! Only call the binding energy gaussian broadening and file printing if necessary
      if (.not. index(photo_output, 'off') > 0) then
        !Broaden ouputs using a gaussian function
        if (index(photo_output, 'e_bind') > 0 .and. temp_photon_energy .ge. (photo_work_function - 0.1)) then
          call binding_energy_broadening
        end if
        !Write either a binding energy output with after Gaussian broadening
        call write_qe_output_files
      end if

    end if
    ! Deallocate the rest that was needed for the photoemission calcs
    call photo_deallocate

    if (on_root) write (stdout, '(1x,a78)') '| End of Photoemission Calculation                                           |'

  end subroutine photo_calculate

  subroutine analyse_geometry
    !* This subroutine identifies and defines a set of layers
    ! and finds the layer that each of the atoms belongs to
    ! or it defines a set of boxes with a height = interlayer distance
    ! at the middle of the slab. The choice between the two is
    ! hardcoded in the variable "new_geom_choice" (True -> Boxes)
    use od_constants, only: dp, periodic_table_name, periodic_table_vdw, deg_to_rad
    use od_cell, only: num_atoms, atoms_pos_cart_photo, atoms_label_tmp, cell_volume, real_lattice
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    use od_parameters, only: photo_imfp_value, photo_slab_max, photo_slab_min, iprint
    implicit none
    integer :: atom_1, atom_2, i, atom_index, temp, first, ierr, atom, ic, counter
    real(kind=dp), allocatable, dimension(:) :: vdw_radii
    real(kind=dp)                            :: z_temp, z_zero = 0.0_dp
    real(kind=dp)                            :: diff_temp, diff_top = 10000.0_dp, diff_bottom = 10000.0_dp
    integer, dimension(2)                    :: indices_top_bottom
    real(kind=dp), dimension(2)              :: mean_heights = 0.0_dp

    allocate (atom_order(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of atom_order failed')

    allocate (vdw_radii(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of atom_order failed')

    allocate (box_atom(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of box_atom failed')
    box_atom = 0

    allocate (layer(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of layer failed')
    do i = 1, num_atoms
      atom_order(i) = i
    end do

    ! Check that we have gamma = 90 deg as that is currently assumed for a lot of calculations!!
    if (real_lattice(3, 1) .gt. 0.000001_dp .and. real_lattice(3, 2) .gt. 0.000001_dp) then
      call io_error('ERROR: analyse_geometry - The c axis is not parallel to the cart. z axis - not currently implemented!')
    end if

    !SORTING ALGORITHM to sort the atoms by z-coordinate
    do atom_1 = 1, num_atoms - 1
      first = atom_order(atom_1)
      do atom_2 = atom_1 + 1, num_atoms
        atom_index = atom_1
        if (atoms_pos_cart_photo(3, atom_order(atom_2)) .gt. atoms_pos_cart_photo(3, first)) then
          first = atom_order(atom_2)
          atom_index = atom_2
        end if

        if (atom_index /= atom_1) then
          temp = atom_order(atom_1)
          atom_order(atom_1) = atom_order(atom_index)
          atom_order(atom_index) = temp
        end if
      end do
    end do

    ! Capitalise the first letter of the atomic label for later
    do atom = 1, num_atoms
      ic = ichar(atoms_label_tmp(atom_order(atom)) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms_label_tmp(atom_order(atom)) (1:1) = char(ic + ichar('Z') - ichar('z'))
    end do

    ! DEFINE THE LAYER FOR EACH ATOM
    ! Assume that a new layer starts if the atom type changes or
    ! the atom is more than 0.5 Angstrom lower than the current layer
    i = 1
    layer(1) = 1
    do atom = 2, num_atoms
      if ((trim(atoms_label_tmp(atom_order(atom))) .ne. trim(atoms_label_tmp(atom_order(atom - 1)))) .or. &
          (abs(atoms_pos_cart_photo(3, atom_order(atom)) - atoms_pos_cart_photo(3, atom_order(atom - 1))) .gt. 0.50)) then
        i = i + 1
      end if
      layer(atom) = i
    end do

    ! --------------------------------------------------------------------------------------------
    ! *    The following code was added in Nov 2023 to test out a new layer assignment scheme    *
    ! *    A set of boxes with the height of the central slab layer distance is created and the  *
    ! *    atoms are sorted into those boxes by their z-coordinate.                              *
    ! --------------------------------------------------------------------------------------------
    ! determine the cell area, and photo_slab_volume for later use
    cell_area = cell_volume/real_lattice(3, 3)
    photo_slab_volume = (photo_slab_max - photo_slab_min)*cell_area
    ! determine the approximate middle of slab as reference
    slab_middle_ref = (photo_slab_max + photo_slab_min)/2
    ! find the nearest two atoms to the middle and determine their layers
    indices_top_bottom = 1
    do atom = 1, num_atoms
      diff_temp = atoms_pos_cart_photo(3, atom_order(atom)) - slab_middle_ref
      ! write(*,*) 'diff_temp(',atom,')',diff_temp
      if (abs(diff_temp) .lt. 0.1) then
        indices_top_bottom(1) = atom_order(atom)
        indices_top_bottom(2) = atom_order(atom + 1)
        slab_middle_ref = atoms_pos_cart_photo(3, atom_order(atom)) - 1
        exit
      end if
      if (diff_temp .gt. 0.0_dp) then
        if (diff_temp .lt. diff_top) then
          indices_top_bottom(1) = atom_order(atom)
          diff_top = diff_temp
        end if
      end if
      if (diff_temp .lt. 0.0_dp) then
        if (abs(diff_temp) .lt. diff_bottom) then
          indices_top_bottom(2) = atom_order(atom)
          diff_bottom = abs(diff_temp)
        end if
      end if
    end do
    ! find potential atoms in the vicinity of the top and bottom atom within 0.5 A
    ! and determing the mean z-coordinate of them (to get mean z-coord of a layer of atoms)
    ! write (*,*) 'indices_top_bottom',indices_top_bottom
    do i = 1, 2
      counter = 0
      diff_top = atoms_pos_cart_photo(3, indices_top_bottom(i)) + 0.5
      diff_bottom = atoms_pos_cart_photo(3, indices_top_bottom(i)) - 0.5
      do atom = 1, num_atoms
        if (atoms_pos_cart_photo(3, atom_order(atom)) .gt. diff_bottom .and. &
            atoms_pos_cart_photo(3, atom_order(atom)) .lt. diff_top) then
          counter = counter + 1
          mean_heights(i) = mean_heights(i) + atoms_pos_cart_photo(3, atom_order(atom))
        end if
      end do
      mean_heights(i) = mean_heights(i)/counter
    end do
    ! determine the box height + box_volume + new slab middle reference
    box_height = mean_heights(1) - mean_heights(2)
    slab_middle_ref = sum(mean_heights)/2
    box_volume = box_height*cell_area
    ! determine the number of boxes we need until we have reached the top of the slab
    num_boxes = ceiling((atoms_pos_cart_photo(3, atom_order(1)) - slab_middle_ref)/box_height)
    if (num_boxes .eq. 0) num_boxes = 1
    ! set up box top points as middle_reference + n(1...)*box_height
    if (.not. allocated(boxes_top_z_coord)) then
      allocate (boxes_top_z_coord(num_boxes))
    end if
    if (.not. allocated(atoms_per_box)) then
      allocate (atoms_per_box(num_boxes))
    end if
    atoms_per_box = 0
    do i = 1, num_boxes
      boxes_top_z_coord(i) = slab_middle_ref + (num_boxes + 1 - i)*box_height
    end do
    ! put each of the atoms into a box
    do i = 1, num_boxes
      counter = 0
      diff_top = boxes_top_z_coord(i)
      diff_bottom = boxes_top_z_coord(i) - box_height
      do atom = 1, num_atoms
        if (atoms_pos_cart_photo(3, atom_order(atom)) .gt. diff_bottom .and. &
            atoms_pos_cart_photo(3, atom_order(atom)) .lt. diff_top) then
          counter = counter + 1
          box_atom(atom) = i
        end if
      end do
      atoms_per_box(i) = counter
    end do
    max_atoms = sum(atoms_per_box)
    box_atom(max_atoms + 1) = num_boxes + 1
    if (on_root) then
      if (iprint .gt. 1) then
        write (stdout, 420) '+', 'box height [Ang] = ', box_height, ',', '# of boxes = ', num_boxes, '+'
420     format(1x, a1, 5x, a19, F13.9, a1, 12x, a13, I4, 9x, a1)
        write (stdout, 421) '+', '# of atoms in each box:', (atoms_per_box(i), i=1, num_boxes)
421     format(1x, a1, 5x, a23, 99(1x, I2))
        write (stdout, 422) '+', 'volume of each box = ', box_volume, ' [Ang**3]', '+'
422     format(1x, a1, 5x, a21, F14.9, a9, 27x, a1)
      end if
      write (stdout, '(1x,a78)') '+------------------------------- Atomic Order  ------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Box     |         Atom Z-Coordinate (Ang)          |'

      do atom = 1, num_atoms
        if (box_atom(atom) .gt. 0 .and. (box_atom(atom) .ne. num_boxes + 1)) then
          write (stdout, '(1x,a3,a2,8x,i3,11x,i3,18x,F12.7,a18)') "|  ", trim(atoms_label_tmp(atom_order(atom))), &
            atom_order(atom), box_atom(atom), atoms_pos_cart_photo(3, atom_order(atom)), "|"
        end if
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, 226) '|  Max number of atoms:', max_atoms, '  Total number of boxes:', num_boxes, '   |'
    end if

    ! This is the older code to define layers and sort atoms into their respective layers
    ! --------------------------------------------------------------------------------------------

    if (on_root) then
      write (stdout, '(1x,a78)') '+------------------------------- Atomic Order  ------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Layer   |         Atom Z-Coordinate (Ang)          |'

      do atom = 1, num_atoms
        write (stdout, '(1x,a3,a2,8x,i3,11x,i3,18x,F12.7,a18)') "|  ", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
          layer(atom), atoms_pos_cart_photo(3, atom_order(atom)), "|"
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    !CALCULATE THE MAX LAYER (HALF SLAB)
    max_layer = ((layer(num_atoms) + 1)/2)

    !CALCULATE THE MAX ATOM (HALF SLAB)
    max_atoms = 0
    do atom = 1, num_atoms
      if (layer(atom) .le. max_layer) then
        max_atoms = max_atoms + 1
      end if
    end do

    if (on_root) then
      write (stdout, 226) '|  Max number of atoms:', max_atoms, '   Max  number of layers:', max_layer, '   |'
226   format(1x, a23, I12, 1x, a25, 1x, I12, a4)

      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    allocate (thickness_atom(max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of thickness_atom failed')
    thickness_atom = 0.0_dp

    allocate (volume_atom(max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of thickness_atom failed')
    thickness_atom = 0.0_dp

    !CALCULATE HOW MANY ATOMS PER LAYER THERE ARE
    allocate (atoms_per_layer(max_layer), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of atoms_per_layer failed')
    atoms_per_layer = 1
    do atom = 2, max_atoms
      if (layer(atom) .eq. layer(atom - 1)) then
        atoms_per_layer(layer(atom)) = atoms_per_layer(layer(atom)) + 1
      end if
    end do

    ! Retreive the van-der-Waals radii from the constants
    do atom_1 = 1, num_atoms
      do i = 1, 109
        if (atoms_label_tmp(atom_order(atom_1)) .eq. periodic_table_name(i)) then
          vdw_radii(atom_1) = periodic_table_vdw(i)
        end if
      end do
    end do

    ! Redefine z=0 (i.e. surface) to the top part of the atom sticking out of the surface the most
    do atom = 1, max_atoms
      if (layer(atom) .gt. 3) exit
      z_temp = atoms_pos_cart_photo(3, atom_order(atom)) + vdw_radii(atom)*sin(45*deg_to_rad)
      if (z_temp .gt. z_zero) then
        z_zero = z_temp
      end if
    end do

    ! Calculate the thickness of each of the layers
    if (max_layer .lt. 2) then
      ! TODO: Check if this is a reasonable estimate for the single layer case
      do atom_1 = 1, max_atoms
        thickness_atom = vdw_radii(atom_1)*2
      end do
    else
      ! Setting thickness_atom for the first atom in layer # 1 and finding the first atom in the second layer
      do atom = 2, max_atoms
        if (layer(atom) .gt. 1) then
          ! write (stdout, *) 'z_zero', z_zero, '1', atoms_pos_cart_photo(3, atom_order(1)), 'second layer', &
          ! atoms_pos_cart_photo(3, atom_order(atom))
          thickness_atom(1) = z_zero - ((atoms_pos_cart_photo(3, atom_order(1)) + atoms_pos_cart_photo(3, atom_order(atom)))/2)
          first_atom_second_l = atom
          exit
        end if
      end do
      ! Setting thickness_atom for the rest of the atoms in the first layer
      do i = 2, first_atom_second_l - 1
        thickness_atom(i) = ((atoms_pos_cart_photo(3, atom_order(i)) - &
                              atoms_pos_cart_photo(3, atom_order(first_atom_second_l)))/2)*2
      end do
      ! Setting thickness_atom for the last atom in the last layer and finding the last atom in the second to last layer
      do i = 1, max_atoms
        if (layer(max_atoms - i) .lt. layer(max_atoms)) then
          thickness_atom(max_atoms) = (abs(atoms_pos_cart_photo(3, atom_order(max_atoms)) - &
                                           atoms_pos_cart_photo(3, atom_order(max_atoms - i)))/2)*2
          last_atom_secondlast_l = max_atoms - i
          exit
        end if
      end do
      ! Setting thickness_atom for the atoms in the last layer, but not for the last atom
      do i = last_atom_secondlast_l + 1, max_atoms - 1
        thickness_atom(i) = (abs(atoms_pos_cart_photo(3, atom_order(i)) - &
                                 atoms_pos_cart_photo(3, atom_order(last_atom_secondlast_l)))/2)*2
      end do
      ! Setting thickness_atom for the atoms in between the first and last layers
      ! Formula -> abs( (z(last atom in n-1th layer) - z(atom) )/2 ) + abs( (z(first atom in n+1th layer) - z(atom) )/2 )
      do atom = first_atom_second_l, last_atom_secondlast_l
        thickness_atom(atom) = abs((atoms_pos_cart_photo(3, atom_order(sum(atoms_per_layer(1:layer(atom) - 1)))) &
                                    - atoms_pos_cart_photo(3, atom_order(atom)))/2) + &
                               abs((atoms_pos_cart_photo(3, atom_order(atom)) - &
                                    atoms_pos_cart_photo(3, atom_order(sum(atoms_per_layer(1:layer(atom))) + 1)))/2)
      end do
    end if

    allocate (thickness_layer(max_layer), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of thickness_layer failed')
    thickness_layer = 0.0_dp

    do atom = 1, max_atoms
      thickness_layer(layer(atom)) = thickness_layer(layer(atom)) + thickness_atom(atom)
    end do
    do i = 1, max_layer
      thickness_layer(i) = thickness_layer(i)/atoms_per_layer(i)
    end do

    do atom = 1, max_atoms
      volume_atom(atom) = (thickness_layer(layer(atom))*cell_area)/atoms_per_layer(layer(atom))
    end do
    if (on_root) then
      write (stdout, '(1x,a78)') '+--------------------- Geometric Analysis of Structure ----------------------+'
      write (stdout, '(1x,a78)') '| Atom | Atom Order | Layer | Layer Thickness | used vdW-rad  | calc. volume |'

      ! Write out the atomic volumes
      do atom = 1, max_atoms
        write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
          layer(atom), thickness_layer(layer(atom)), vdw_radii(atom), volume_atom(atom), "    |"
225     format(1x, a1, a4, 6x, I3, 8x, I3, 6x, E14.6E3, 3x, F11.4, 3x, F11.4, a5)
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    !TEST IF THE SUPPLIED IMFP LIST IS LONG ENOUGH
    if (allocated(photo_imfp_value) .and. size(photo_imfp_value, 1) .gt. 1 .and. size(photo_imfp_value, 1) .lt. max_layer) then
      call io_error('The supplied list of layer dependent imfp values is less than the calculated max_layer. Check input!')
    end if

  end subroutine analyse_geometry

  subroutine calc_band_info
    !===============================================================================
    ! This subroutine determines useful indices of band energies for later use in
    ! the QE and MTE calculation to reduce loop times.
    ! This relies on an IMPORTANT assumption: the bands file is ordered by energy
    ! and not by band number (e.g. after being processed by bands2orbitals)
    ! Felix Mildner, 28th March 2023
    !===============================================================================
    use od_electronic, only: efermi, band_energy, nbands, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: iprint
    use od_io, only: stdout, io_time, io_error
    implicit none
    integer         :: N_k, N_spin, n_eigen, ierr
    real(kind=dp)   :: time0, time1

    time0 = io_time()

    allocate (min_index_unocc(nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_band_info - allocation of min_index_unocc failed')

    do N_k = 1, num_kpoints_on_node(my_node_id)  ! Loop over kpoints
      do N_spin = 1, nspins                           ! Loop over spins
        do n_eigen = 2, nbands                        ! Loop over bands
          ! TODO: Test if this is the behaviour we want and or if we have to change the condition
          if (band_energy(n_eigen - 1, N_spin, N_k) .gt. band_energy(n_eigen, N_spin, N_k)) then
            call io_error('Error: the band energies in the .bands file used are NOT ORDERED CORRECTLY (i.e. by increasing energy) &
            & which will give WRONG RESULTS with the current code!')
          end if
        end do
      end do
    end do

    do N_k = 1, num_kpoints_on_node(my_node_id)  ! Loop over kpoints
      do N_spin = 1, nspins                           ! Loop over spins
        do n_eigen = 1, nbands                        ! Loop over bands
          ! TODO: Test if this is the behaviour we want and or if we have to change the condition
          if (band_energy(n_eigen, N_spin, N_k) .gt. efermi) then
            min_index_unocc(N_spin, N_k) = n_eigen
            exit
          end if
        end do
      end do
    end do

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a36,23x,f11.3,a8)') '+ Time to calculate Band Energy Info', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_band_info

  subroutine calc_photon_energies
    use od_constants, only: dp
    use od_parameters, only: photo_photon_sweep, photo_photon_min, photo_photon_max, jdos_spacing, photo_photon_energy
    use od_io, only: io_error
    implicit none
    real(kind=dp)        ::   num_energies, temp
    integer              ::   ierr, i

    if (photo_photon_sweep) then
      num_energies = (photo_photon_max - photo_photon_min)/jdos_spacing
      number_energies = int(num_energies) + 1
      if (photo_photon_max - photo_photon_min .eq. 0.0_dp) then
        number_energies = 1
      else if (mod(num_energies, 1.0_dp) > 1.0E-10_dp) then
        number_energies = number_energies + 1
        if (abs(mod(num_energies, 1.0_dp) - 1) > 1.0E-10_dp) &
          call io_error('Error: calc_photon_energies - given photon sweep min/max values do not give integer # of photon steps')
      end if
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of index_energy failed')
      do i = 1, number_energies
        temp = (i - 1)*jdos_spacing + photo_photon_min
        ! Account for E = 0.0
        index_energy(i) = int(temp/jdos_spacing) + 1
      end do
    else
      number_energies = 1
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of index_energy failed')
      ! Account for E = 0.0
      index_energy(number_energies) = int(photo_photon_energy/jdos_spacing) + 1
    end if

  end subroutine calc_photon_energies

  subroutine make_pdos_weights_atoms
    !!This subroutine is equivalent to pdos_merge of pdos.F90, but only for atoms
    use od_electronic, only: pdos_orbital, pdos_weights, pdos_mwab, nspins
    use od_cell, only: num_kpoints_on_node, num_atoms, cell_calc_kpoint_r_cart, kpoint_r_cart
    use od_comms, only: my_node_id, on_root
    use od_io, only: io_error, stdout, seedname, io_date
    use od_parameters, only: devel_flag
    implicit none
    character(len=9) :: ctime             ! Temp. time string
    character(len=11):: cdate             ! Temp. date string
    integer :: N_k, N_spin, n_eigen, np, ierr, atom, box, i, i_max, pdos_unit = 32

    allocate (pdos_weights_atoms(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id), num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')

    allocate (pdos_weights_k_band(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_k_band failed')

    pdos_weights_atoms = 0.0_dp
    pdos_weights_k_band = 0.0_dp

    if (new_geom_choice) then
      allocate (pdos_weights_boxes(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id), num_boxes), stat=ierr)
      if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')
      pdos_weights_boxes = 0.0_dp
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, pdos_mwab%nbands
          i = 1
          do np = 1, pdos_mwab%norbitals
            if (np .gt. 1) then
              if (pdos_orbital%rank_in_species(np) .ne. pdos_orbital%rank_in_species(np - 1)) then
                i = i + 1
              end if
            end if
            pdos_weights_atoms(n_eigen, N_spin, N_k, i) = &
              pdos_weights_atoms(n_eigen, N_spin, N_k, i) + &
              pdos_weights(np, n_eigen, N_k, N_spin)
          end do
        end do
      end do
    end do
    i_max = i
    do atom = 1, num_atoms
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, pdos_mwab%nbands
            if (pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) .lt. 0.0_dp) then
              pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) = 0.0_dp
            end if
            pdos_weights_k_band(n_eigen, N_spin, N_k) = pdos_weights_k_band(n_eigen, N_spin, N_k) + &
                                                        pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom))
          end do
        end do
      end do
    end do
    if (new_geom_choice) then
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, pdos_mwab%nbands
              if (pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) .lt. 0.0_dp) then
                pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) = 0.0_dp
              end if
              pdos_weights_boxes(n_eigen, N_spin, N_k, box_atom(atom)) = &
                pdos_weights_boxes(n_eigen, N_spin, N_k, box_atom(atom)) + &
                pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom))
            end do
          end do
        end do
      end do
    end if

    if (index(devel_flag, 'output_pdos_weights') > 0 .and. on_root) then
      call cell_calc_kpoint_r_cart
      write (stdout, '(a78)') "+---------------- Printing K-Points in Cartesian Coordinates ----------------+"
      i = 0
      do N_k = 1, num_kpoints_on_node(my_node_id)
        write (stdout, '(1x,I4,4x,3(1x,E22.15))') i, kpoint_r_cart(:, N_k)
        i = i + 1
      end do
      call io_date(cdate, ctime)
      ! write out atomic/box weights
      if (new_geom_choice) then
        open (unit=pdos_unit, action='write', file=trim(seedname)//'_pdos_boxes.dat')
        write (pdos_unit, '(1x,a28)') '############################'
        write (pdos_unit, *) '# OptaDOS Photoemission: Printing PDOS-Boxes-Weights on ', cdate, ' at ', ctime
        write (pdos_unit, '(1x,a19,1x,a99)') '# PDOS weights for', seedname
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of PDOS Bands :', size(pdos_weights_boxes, 1)
        write (pdos_unit, '(1x,a24,1x,I2)') '# Number of Spins      :', size(pdos_weights_boxes, 2)
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of K-points   :', size(pdos_weights_boxes, 3)
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of Boxes      :', size(pdos_weights_boxes, 4)
        write (pdos_unit, '(1x,a45)') '# F U L L _ P D O S _ B O X _ W E I G H T S'
        write (pdos_unit, '(1x,a28)') '############################'
        do box = 1, num_boxes
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (pdos_unit, '(9999(1x,es24.16))') (pdos_weights_boxes(n_eigen, N_spin, N_k, box), n_eigen=1, pdos_mwab%nbands)
            end do
          end do
        end do
      else
        open (unit=pdos_unit, action='write', file=trim(seedname)//'_pdos_atoms.dat')
        write (pdos_unit, '(1x,a28)') '############################'
        write (pdos_unit, *) '# OptaDOS Photoemission: Printing PDOS-Atoms-Weights on ', cdate, ' at ', ctime
        write (pdos_unit, '(1x,a19,1x,a99)') '# PDOS weights for', seedname
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of PDOS Bands :', size(pdos_weights_atoms, 1)
        write (pdos_unit, '(1x,a24,1x,I2)') '# Number of Spins      :', size(pdos_weights_atoms, 2)
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of K-points   :', size(pdos_weights_atoms, 3)
        write (pdos_unit, '(1x,a24,1x,I4)') '# Number of Atoms      :', size(pdos_weights_atoms, 4)
        write (pdos_unit, '(1x,a45)') '# F U L L _ P D O S _ A T O M _ W E I G H T S'
        write (pdos_unit, '(1x,a28)') '############################'
        do atom = 1, num_atoms
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (pdos_unit, '(9999(1x,es24.16))') (pdos_weights_atoms(n_eigen, N_spin, N_k, atom), n_eigen=1, pdos_mwab%nbands)
            end do
          end do
        end do
      end if
      close (unit=pdos_unit)
      ! Write out the k-band weights
      open (unit=pdos_unit, action='write', file=trim(seedname)//'_pdos_k_band.dat')
      write (pdos_unit, '(1x,a28)') '############################'
      write (pdos_unit, *) '# OptaDOS Photoemission: Printing PDOS-Weights-K-Band on ', cdate, ' at ', ctime
      write (pdos_unit, '(1x,a19,1x,a99)') '# PDOS weights for', seedname
      write (pdos_unit, '(1x,a24,1x,I4)') '# Number of PDOS Bands :', size(pdos_weights_k_band, 1)
      write (pdos_unit, '(1x,a24,1x,I2)') '# Number of Spins      :', size(pdos_weights_k_band, 2)
      write (pdos_unit, '(1x,a24,1x,I4)') '# Number of K-points   :', size(pdos_weights_k_band, 3)
      write (pdos_unit, '(1x,a45)') '# F U L L _ P D O S _ K _ B A N D _ W E I G H T S'
      write (pdos_unit, '(1x,a28)') '############################'
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          write (pdos_unit, '(9999(1x,es24.16))') (pdos_weights_k_band(n_eigen, N_spin, N_k), n_eigen=1, pdos_mwab%nbands)
        end do
      end do
      close (unit=pdos_unit)
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------ Printing pDOS_weights_atoms -----------------------+'
      write (stdout, 125) shape(pdos_weights_atoms)
      write (stdout, 125) i_max, pdos_mwab%nbands, num_kpoints_on_node(my_node_id), nspins
125   format(4(1x, I4))
      write (stdout, '(9999(es15.8))') ((((pdos_weights_atoms(n_eigen, N_spin, N_k, i), N_spin=1, nspins) &
                                          , n_eigen=1, pdos_mwab%nbands), N_k=1, num_kpoints_on_node(my_node_id)), i=1, i_max)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
      write (stdout, '(1x,a78)') '+----------------------- Printing pDOS_weights_k_band -----------------------+'
      write (stdout, 124) shape(pdos_weights_k_band)
      write (stdout, 124) pdos_mwab%nbands, num_kpoints_on_node(my_node_id), nspins
124   format(3(1x, I4))
      write (stdout, '(9999(es15.8))') (((pdos_weights_k_band(n_eigen, N_spin, N_k), &
                                          N_k=1, num_kpoints_on_node(my_node_id)), N_spin=1, nspins), n_eigen=1, pdos_mwab%nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if
  end subroutine make_pdos_weights_atoms

  subroutine calc_photo_optics
    !! This subroutine calculates the projected optical characteristics for each layer.
    use od_optics, only: make_weights, calc_epsilon_2, calc_epsilon_1, calc_refract, calc_absorp, calc_reflect, &
                         epsilon, refract, absorp, reflect, intra, write_absorp, write_epsilon, write_reflect, write_refract
    use od_io, only: stdout, io_error, io_time, seedname, io_date
    use od_electronic, only: elec_read_optical_mat, nbands, nspins, efermi, elec_dealloc_optical, elec_read_band_gradient, &
                             nbands, nspins, band_energy
    use od_cell, only: num_kpoints_on_node, num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r
    use od_jdos_utils, only: jdos_utils_calculate, jdos_nbins, setup_energy_scale, jdos_deallocate, E
    use od_comms, only: comms_bcast, on_root, my_node_id
    use od_parameters, only: optics_intraband, jdos_spacing, devel_flag, iprint, jdos_max_energy, photo_model
    use od_dos_utils, only: dos_utils_calculate_at_e
    use od_constants, only: epsilon_0, e_charge
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: dos_matrix_weights
    real(kind=dp), allocatable, dimension(:, :) :: weighted_dos_at_e
    real(kind=dp), allocatable, dimension(:, :) :: dos_at_e
    integer :: N_k, N2, N_spin, n_eigen, n_eigen_final, atom, ierr, energy, box, initial
    integer :: jdos_bin, i, s, is, idos, wjdos_unit = 23, ome_unit = 32
    real(kind=dp)    :: time0, time1
    ! logical, dimension(3) :: gam
    character(len=9) :: ctime             ! Temp. time string
    character(len=11):: cdate             ! Temp. date string
    character(len=3) :: atom_s
    character(len=4) :: initial_s

    time0 = io_time()

    if (new_geom_choice) then
      allocate (absorp_photo(num_boxes, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')
    else
      allocate (absorp_photo(max_atoms, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')
    end if

    if (new_geom_choice) then
      allocate (reflect_photo(num_boxes, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')
    else
      allocate (reflect_photo(max_atoms, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')
    end if

    call make_weights(matrix_weights)
    N_geom = size(matrix_weights, 5)
    call elec_dealloc_optical

    if (index(devel_flag, 'output_ome_itof') > 0 .and. on_root) then
      call io_date(cdate, ctime)
      do N_k = 1, size(kpoint_r, 2)
        if (all(abs(kpoint_r(:, N_k)) < 1.0E-10_dp)) then
          is = N_k
          exit
        end if
      end do
      write (stdout, *) 'The gamma point was determined to be - ', is
      i = index(devel_flag, 'output_ome_itof')
      read (devel_flag(i + 16:i + 20), *) initial
      write (initial_s, '(I4)') initial
      write (stdout, '(1x,a27,I4,a32)') 'Outputting OMEs for band # ', initial, ' to the bands above it at Gamma.'
      open (unit=ome_unit, action='write', file=trim(seedname)//'_OMEs_from_band_'//trim(adjustl(initial_s))//'.dat')
      write (ome_unit, '(1x,a28)') '############################'
      write (ome_unit, *) '# OptaDOS Photoemission: Printing PDOS-Atoms-Weights on ', cdate, ' at ', ctime
      write (ome_unit, '(1x,a16,1x,a99)') '# OM weights for', seedname
      write (ome_unit, '(1x,a23,1x,I4)') '# Initial Band Choice :', initial
      write (ome_unit, '(1x,a23,1x,F15.7)') '# Band Energy        : ', (band_energy(initial, 1, is) - efermi)
      write (ome_unit, '(1x,a28)') '############################'
      do n_eigen = initial + 1, nbands
        write (ome_unit, '(1x,a6,I4,1x,E20.12E3,1x,F15.7)') 'Band #', n_eigen, matrix_weights(initial, n_eigen, is, 1, 1),&
        & (band_energy(n_eigen, 1, is) - efermi)
      end do
      close (unit=ome_unit)
      return
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+-------------------------- Printing Matrix Weights -------------------------+'
      write (stdout, 126) shape(matrix_weights)
      write (stdout, 126) nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
      write (stdout, '(9999(es15.8))') (((((matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2), N2=1, N_geom), &
                                           N_spin=1, nspins), N_k=1, num_kpoints_on_node(my_node_id)), &
                                         n_eigen_final=1, nbands), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    allocate (projected_matrix_weights(nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of projected_matrix_weights failed')
    if (new_geom_choice) then
      if (.not. index(devel_flag, 'ds_like_pe') > 0) then
        do box = 1, num_boxes                           ! Loop over boxes
          !
          if (iprint > 1 .and. on_root) then
            write (stdout, 145) '+------------------------ Starting BOX  # ', box, ' of ', num_boxes, ' ------------------------+'
          end if
          ! (Re-)Setting the weights for new box
          projected_matrix_weights = 0.0_dp

          do N2 = 1, N_geom
            do N_k = 1, num_kpoints_on_node(my_node_id)    ! Loop over kpoints
              do N_spin = 1, nspins                    ! Loop over spins
                do n_eigen = 1, nbands               ! Loop over state 1
                  do n_eigen_final = n_eigen, nbands    ! Loop over state 2
                    if (band_energy(n_eigen, N_spin, N_k) > efermi .and. n_eigen /= n_eigen_final) cycle
                    if (band_energy(n_eigen_final, N_spin, N_k) < efermi .and. n_eigen /= n_eigen_final) cycle
                    if (pdos_weights_k_band(n_eigen, N_spin, N_k) .eq. 0.0_dp) then
                      cycle
                    end if
                    projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2) = &
                      matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2)* &
                      (pdos_weights_boxes(n_eigen, N_spin, N_k, box)/pdos_weights_k_band(n_eigen, N_spin, N_k))
                  end do                        ! Loop over state 2
                end do                            ! Loop over state 1
              end do                                ! Loop over spins
            end do                                    ! Loop over kpoints
          end do

          if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
            write (stdout, '(1x,a37,I3,a38)') '+-------------------------------Atom-', atom, &
              '-------------------------------------+'
            write (stdout, '(1x,a78)') '+--------------------- Printing Projected Matrix Weights --------------------+'
            write (stdout, 126) shape(projected_matrix_weights)
            write (stdout, 126) nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
            write (stdout, '(9999(es15.8))') (((((projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2), &
                                                  N2=1, N_geom), N_spin=1, nspins), N_k=1, num_kpoints_on_node(my_node_id)), &
                                               n_eigen_final=1, nbands), n_eigen=1, nbands)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          end if

          ! Send matrix element to jDOS routine and get weighted jDOS back
          call jdos_utils_calculate(projected_matrix_weights, weighted_jdos=weighted_jdos)

          if (on_root .and. iprint .gt. 2) then
            N_geom = size(matrix_weights, 5)
            write (atom_s, '(I3)') box + 100
            open (unit=wjdos_unit, action='write', file=trim(seedname)//'_weighted_jdos_'//trim(adjustl(atom_s))//'.dat')
            write (wjdos_unit, '(1x,a28)') '############################'
            write (wjdos_unit, '(1x,a19,1x,a99)') '# Weighted JDOS for', seedname
            write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# maximum JDOS energy :', jdos_max_energy, '[eV]'
            write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# JDOS step size      :', jdos_spacing, '[eV]'
            write (wjdos_unit, '(1x,a28)') '############################'
            do is = 1, nspins
              write (wjdos_unit, *) 'Spin Channel :', is
              do idos = 1, jdos_nbins
                write (wjdos_unit, *) E(idos), ' , ', sum(weighted_jdos(idos, is, 1:N_geom))
              end do
            end do
            close (unit=wjdos_unit)
          end if

          if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
            write (stdout, '(1x,a78)') '+------------------------ Printing Weighted Joint-DOS -----------------------+'
            write (stdout, 124) shape(weighted_jdos)
            write (stdout, 124) jdos_nbins, nspins, N_geom
            write (stdout, '(9999(es15.8))') (((weighted_jdos(jdos_bin, N_spin, N2), N2=1, N_geom), N_spin=1, nspins) &
                                              , jdos_bin=1, jdos_nbins)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          end if

          if (optics_intraband) then
            allocate (dos_matrix_weights(size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of dos_matrix_weights failed')
            allocate (dos_at_e(3, nspins), stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of dos_at_e failed')
            allocate (weighted_dos_at_e(nspins, size(matrix_weights, 5)), stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of weighted_dos_at_e failed')
            dos_at_e = 0.0_dp
            weighted_dos_at_e = 0.0_dp
            do N_k = 1, size(matrix_weights, 5)
              do N2 = 1, nbands
                dos_matrix_weights(N_k, N2, :, :) = matrix_weights(N2, N2, :, :, N_k)
              end do
            end do
            call dos_utils_calculate_at_e(efermi, dos_at_e, dos_matrix_weights, weighted_dos_at_e)
            weighted_dos_at_e = weighted_dos_at_e/atoms_per_box(box)
          end if

          if (on_root) then
            if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0 .and. optics_intraband) then
              write (stdout, '(1x,a36,f8.4,a34)') '+------------------------ E_Fermi = ', efermi, &
                '---------------------------------+'
              write (stdout, '(1x,a78)') '+------------------------ Printing DOS Matrix Weights -----------------------+'
              write (stdout, 125) shape(dos_matrix_weights)
              write (stdout, 125) size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins
              write (stdout, '(9999(es15.8))') ((((dos_matrix_weights(n_eigen, n_eigen_final, N_k, s), s=1, nspins), N_k=1, &
                                                  num_kpoints_on_node(my_node_id)), n_eigen_final=1, nbands), n_eigen=1, &
                                                size(matrix_weights, 5))
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
              write (stdout, '(1x,a78)') '+--------------------------- Printing DOS @ Energy --------------------------+'
              write (stdout, '(9(es15.8))') ((dos_at_e(i, s), i=1, 3), s=1, nspins)
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
              write (stdout, '(1x,a78)') '+----------------------- Printing Weighted DOS @ Energy ---------------------+'
              write (stdout, '(9999(es15.8))') ((weighted_dos_at_e(s, n_eigen), s=1, nspins), n_eigen=1, size(matrix_weights, 5))
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
            end if

            ! Calculate epsilon_2
            call calc_epsilon_2(weighted_jdos, weighted_dos_at_e, box_volume)

            ! Calculate epsilon_1
            call calc_epsilon_1

            ! Calculate other optical properties
            call calc_refract
            call calc_absorp
            call calc_reflect

            if (iprint .gt. 2) then
              call write_epsilon(box + 100, photo_at_e=dos_at_e, photo_volume=box_volume)
              call write_refract(box + 100, photo_volume=box_volume)
              call write_absorp(box + 100, photo_volume=box_volume)
              call write_reflect(box + 100, photo_volume=box_volume)
            end if

            do energy = 1, number_energies
              absorp_photo(box, energy) = absorp(index_energy(energy))
              reflect_photo(box, energy) = reflect(index_energy(energy))
            end do

            if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0) then
              write (stdout, '(1x,a78)') '+-------------------- Printing Material Optical Properties ------------------+'
              write (stdout, '(1x,a78)') '+--------------------------- Printing Epsilon Array -------------------------+'
              write (stdout, 125) shape(epsilon)
              if (.not. optics_intraband) then
                write (stdout, '(9999(E17.8E3))') (((epsilon(jdos_bin, N_k, N2, 1), jdos_bin=1, jdos_nbins), N_k=1, 2), &
                                                   N2=1, N_geom)
              else
                write (stdout, '(9999(E17.8E3))') ((((epsilon(jdos_bin, N_k, N2, i), jdos_bin=1, jdos_nbins), N_k=1, 2), &
                                                    N2=1, N_geom), i=1, 3)
              end if
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

              write (stdout, '(1x,a78)') '+----------------------------- Printing Absorption --------------------------+'
              write (stdout, '(99(E17.8E3))') (absorp_photo(atom, energy), energy=1, number_energies)
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

              write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection --------------------------+'
              write (stdout, '(99(E17.8E3))') (reflect_photo(atom, energy), energy=1, number_energies)
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
            end if
            if (enable_debug_output .and. iprint .gt. 2) then
              write (stdout, '(1x,a78)') '+----------------------------- Printing Absorption - box --------------------+'
              write (stdout, '(99(E17.8E3))') (absorp_photo(box, energy), energy=1, number_energies)
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

              write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection - box --------------------+'
              write (stdout, '(99(E17.8E3))') (reflect_photo(box, energy), energy=1, number_energies)
              write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
            end if
            ! Deallocate extra arrays produced in the case of using optics_intraband
            deallocate (epsilon, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate epsilon')
            deallocate (refract, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate refract')
            deallocate (absorp, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate absorp')
            deallocate (reflect, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate reflect')
            if (optics_intraband) then
              deallocate (intra, stat=ierr)
              if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate intra')
            end if
          end if
          if (optics_intraband) then
            deallocate (dos_matrix_weights, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_matrix_weights')
            deallocate (dos_at_e, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_at_e')
            deallocate (weighted_dos_at_e, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_dos_at_e')
          end if
          call jdos_deallocate
          deallocate (weighted_jdos, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_jdos')
          ! deallocate (E, stat=ierr)
          ! if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate E')
        end do                                        ! Loop over boxes
        call comms_bcast(absorp_photo(1, 1), num_boxes*number_energies)
        call comms_bcast(reflect_photo(1, 1), num_boxes*number_energies)
      end if
    else
      if (.not. index(devel_flag, 'ds_like_pe') > 0) then
      do atom = 1, max_atoms                           ! Loop over atoms
        !
        if (iprint > 1 .and. on_root) then
          write (stdout, 145) '+------------------------ Starting Atom # ', atom, ' of ', max_atoms, ' ------------------------+'
145       format(1x, a42, I3, a4, I3, a26)
        end if
        ! (Re-)Setting the weights for new atom
        projected_matrix_weights = 0.0_dp

        do N2 = 1, N_geom
          do N_k = 1, num_kpoints_on_node(my_node_id)    ! Loop over kpoints
            do N_spin = 1, nspins                    ! Loop over spins
              do n_eigen = 1, nbands               ! Loop over state 1
                do n_eigen_final = n_eigen, nbands    ! Loop over state 2
                  if (band_energy(n_eigen, N_spin, N_k) > efermi .and. n_eigen /= n_eigen_final) cycle
                  if (band_energy(n_eigen_final, N_spin, N_k) < efermi .and. n_eigen /= n_eigen_final) cycle
                  if (pdos_weights_k_band(n_eigen, N_spin, N_k) .eq. 0.0_dp) then
                    ! write (stdout,'(I2,1x,I4,1x,I1,1x,I3)') atom, N_k, N_spin, n_eigen
                    ! write (stdout,'(99(ES19.12))') pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)),&
                    ! pdos_weights_k_band(n_eigen, N_spin, N_k)
                    ! call FLUSH()
                    cycle
                  end if
                  projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2) = &
                    matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2)* &
                    (pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom))/pdos_weights_k_band(n_eigen, N_spin, N_k))
                end do                        ! Loop over state 2
              end do                            ! Loop over state 1
            end do                                ! Loop over spins
          end do                                    ! Loop over kpoints
        end do

        if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
          write (stdout, '(1x,a37,I3,a38)') '+-------------------------------Atom-', atom, '-------------------------------------+'
          write (stdout, '(1x,a78)') '+--------------------- Printing Projected Matrix Weights --------------------+'
          write (stdout, 126) shape(projected_matrix_weights)
          write (stdout, 126) nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
126       format(5(1x, I4))
          write (stdout, '(9999(es15.8))') (((((projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2), N2=1, N_geom), &
                                               N_spin=1, nspins), N_k=1, num_kpoints_on_node(my_node_id)), &
                                             n_eigen_final=1, nbands), n_eigen=1, nbands)
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
        end if

        ! Send matrix element to jDOS routine and get weighted jDOS back
        call jdos_utils_calculate(projected_matrix_weights, weighted_jdos=weighted_jdos)

        if (on_root) then
          N_geom = size(matrix_weights, 5)
          write (atom_s, '(I2)') atom
          open (unit=wjdos_unit, action='write', file=trim(seedname)//'_weighted_jdos_'//trim(adjustl(atom_s))//'.dat')
          write (wjdos_unit, '(1x,a28)') '############################'
          write (wjdos_unit, '(1x,a19,1x,a99)') '# Weighted JDOS for', seedname
          write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# maximum JDOS energy :', jdos_max_energy, '[eV]'
          write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# JDOS step size      :', jdos_spacing, '[eV]'
          write (wjdos_unit, '(1x,a28)') '############################'
          do is = 1, nspins
            write (wjdos_unit, *) 'Spin Channel :', is
            do idos = 1, jdos_nbins
              write (wjdos_unit, *) E(idos), ' , ', sum(weighted_jdos(idos, is, 1:N_geom))
            end do
          end do
          close (unit=wjdos_unit)
        end if

        if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
          write (stdout, '(1x,a78)') '+------------------------ Printing Weighted Joint-DOS -----------------------+'
          write (stdout, 124) shape(weighted_jdos)
          write (stdout, 124) jdos_nbins, nspins, N_geom
124       format(3(1x, I4))
          write (stdout, '(9999(es15.8))') (((weighted_jdos(jdos_bin, N_spin, N2), N2=1, N_geom), N_spin=1, nspins) &
                                            , jdos_bin=1, jdos_nbins)
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
        end if

        if (optics_intraband) then
          allocate (dos_matrix_weights(size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of dos_matrix_weights failed')
          allocate (dos_at_e(3, nspins), stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of dos_at_e failed')
          allocate (weighted_dos_at_e(nspins, size(matrix_weights, 5)), stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of weighted_dos_at_e failed')
          dos_at_e = 0.0_dp
          weighted_dos_at_e = 0.0_dp
          do N_k = 1, size(matrix_weights, 5)
            do N2 = 1, nbands
              dos_matrix_weights(N_k, N2, :, :) = matrix_weights(N2, N2, :, :, N_k)
            end do
          end do
          call dos_utils_calculate_at_e(efermi, dos_at_e, dos_matrix_weights, weighted_dos_at_e)
          weighted_dos_at_e = weighted_dos_at_e/atoms_per_layer(layer(atom))
        end if

        if (on_root) then
          if (index(devel_flag, 'print_qe_constituents') > 0 .and. optics_intraband) then
            write (stdout, '(1x,a36,f8.4,a34)') '+------------------------ E_Fermi = ', efermi, '---------------------------------+'
            write (stdout, '(1x,a78)') '+------------------------ Printing DOS Matrix Weights -----------------------+'
            write (stdout, 125) shape(dos_matrix_weights)
            write (stdout, 125) size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins
125         format(4(1x, I4))
            write (stdout, '(9999(es15.8))') ((((dos_matrix_weights(n_eigen, n_eigen_final, N_k, s), s=1, nspins), N_k=1, &
                                                num_kpoints_on_node(my_node_id)), n_eigen_final=1, nbands), n_eigen=1, &
                                              size(matrix_weights, 5))
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
            write (stdout, '(1x,a78)') '+--------------------------- Printing DOS @ Energy --------------------------+'
            write (stdout, '(9(es15.8))') ((dos_at_e(i, s), i=1, 3), s=1, nspins)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
            write (stdout, '(1x,a78)') '+----------------------- Printing Weighted DOS @ Energy ---------------------+'
            write (stdout, '(9999(es15.8))') ((weighted_dos_at_e(s, n_eigen), s=1, nspins), n_eigen=1, size(matrix_weights, 5))
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          end if

          ! Calculate epsilon_2
          call calc_epsilon_2(weighted_jdos, weighted_dos_at_e, volume_atom(atom))

          ! Calculate epsilon_1
          call calc_epsilon_1

          ! Calculate other optical properties
          call calc_refract
          call calc_absorp
          call calc_reflect

          if (iprint .gt. 2) then
            call write_epsilon(atom, photo_at_e=dos_at_e, photo_volume=volume_atom(atom))
            call write_refract(atom, photo_volume=volume_atom(atom))
            call write_absorp(atom, photo_volume=volume_atom(atom))
            call write_reflect(atom, photo_volume=volume_atom(atom))
          end if

          do energy = 1, number_energies
            absorp_photo(atom, energy) = absorp(index_energy(energy))
            reflect_photo(atom, energy) = reflect(index_energy(energy))
          end do

          if (index(devel_flag, 'print_qe_constituents') > 0) then
            write (stdout, '(1x,a78)') '+-------------------- Printing Material Optical Properties ------------------+'
            write (stdout, '(1x,a78)') '+--------------------------- Printing Epsilon Array -------------------------+'
            write (stdout, 125) shape(epsilon)
            if (.not. optics_intraband) then
              write (stdout, '(9999(E17.8E3))') (((epsilon(jdos_bin, N_k, N2, 1), jdos_bin=1, jdos_nbins), N_k=1, 2), N2=1, N_geom)
            else
              write (stdout, '(9999(E17.8E3))') ((((epsilon(jdos_bin, N_k, N2, i), jdos_bin=1, jdos_nbins), N_k=1, 2), &
                                                  N2=1, N_geom), i=1, 3)
            end if
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

            write (stdout, '(1x,a78)') '+----------------------------- Printing Absorption --------------------------+'
            write (stdout, '(99(E17.8E3))') (absorp_photo(atom, energy), energy=1, number_energies)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

            write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection --------------------------+'
            write (stdout, '(99(E17.8E3))') (reflect_photo(atom, energy), energy=1, number_energies)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          end if
          if (iprint .gt. 2) then
            write (stdout, '(1x,a78)') '+----------------------------- Printing Absorption --------------------------+'
            write (stdout, '(99(E17.8E3))') (absorp_photo(atom, energy), energy=1, number_energies)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

            write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection --------------------------+'
            write (stdout, '(99(E17.8E3))') (reflect_photo(atom, energy), energy=1, number_energies)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          end if
          ! Deallocate extra arrays produced in the case of using optics_intraband
          deallocate (epsilon, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate epsilon')
          deallocate (refract, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate refract')
          deallocate (absorp, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate absorp')
          deallocate (reflect, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate reflect')
          if (optics_intraband) then
            deallocate (intra, stat=ierr)
            if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate intra')
          end if
        end if
        if (optics_intraband) then
          deallocate (dos_matrix_weights, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_matrix_weights')
          deallocate (dos_at_e, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_at_e')
          deallocate (weighted_dos_at_e, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_dos_at_e')
        end if
        call jdos_deallocate
        deallocate (weighted_jdos, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_jdos')
        ! deallocate (E, stat=ierr)
        ! if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate E')
      end do                                        ! Loop over atoms
      call comms_bcast(absorp_photo(1, 1), max_atoms*number_energies)
      call comms_bcast(reflect_photo(1, 1), max_atoms*number_energies)
      end if
    end if

    ! Deallocating this out of the loop to reduce memory operations - could lead to higher memory consumption
    deallocate (projected_matrix_weights, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate projected_matrix_weights')

    if (index(photo_model, '1step') > 0) then
      deallocate (matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate matrix_weights')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a52,7x,f11.3,a8)') '+ Time to calculate Photoemission Optical Properties', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_photo_optics

  subroutine calc_absorp_layer
    !!This subroutine calculates the absorption coefficient for a specific layer
    ! use od_cell, only: atoms_pos_cart_photo
    use od_jdos_utils, only: jdos_nbins
    use od_parameters, only: devel_flag, iprint
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    implicit none
    real(kind=dp) :: I_0
    integer :: atom, box, i, ierr, num_layer

    if (.not. new_geom_choice) then
      allocate (I_layer(max_layer + 1, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of I_layer failed')
      I_layer = 0.0_dp
      I_0 = 1.0_dp
      I_layer = 1.0_dp

      do i = 1, number_energies
        I_layer(1, i) = I_0 - reflect_photo(1, i)
      end do
      if (max_layer .gt. 1) then
        do atom = first_atom_second_l, max_atoms
          do i = 1, number_energies
            ! Thickness_atom is the thickness of each layer representing the path length the light takes through the layer
            I_layer(layer(atom), i) = I_layer(layer(atom) - 1, i)* &
                                      exp(-(absorp_photo(atom, i)*thickness_atom(atom)*1E-10))
            if (I_layer(layer(atom), i) .lt. 0.0_dp) then
              I_layer(layer(atom), i) = 0.0_dp
            end if
          end do
        end do
      end if
      ! Since we later combine the bulk slab emission probability (contains already light intensity) into the
      ! layer by layer emission probability array (does not contain light intensity), we have to set the
      ! intensity value artifically to 1.0 to have it not influence the final value.
      ! We are only ever accessing I_layer to max_atoms, so this has no effect on the rest.
      I_layer(layer(max_atoms + 1), 1:number_energies) = 1.0_dp
    else
      allocate (I_layer(num_boxes + 1, number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of I_layer failed')
      I_layer = 0.0_dp

      I_0 = 1.0_dp
      I_layer = 1.0_dp

      do i = 1, number_energies
        I_layer(1, i) = I_0 - reflect_photo(1, i)
      end do
      ! Calculate the unreflected portion of incoming light
      do i = 1, number_energies
        I_layer(1, i) = I_0 - reflect_photo(1, i)
      end do
      ! If we have more than one box with atoms in it, calculate the incident light intensity for each
      if (num_boxes .gt. 1) then
        do box = 2, num_boxes
          do i = 1, number_energies
            I_layer(box, i) = I_layer(box - 1, i)* &
                              exp(-(absorp_photo(box, i)*box_height*1E-10))
            if (I_layer(box, i) .lt. 0.0_dp) I_layer(box, i) = 0.0_dp
          end do
        end do
      end if
      ! Since we later combine the bulk slab emission probability (contains already light intensity) into the
      ! layer by layer emission probability array (does not contain light intensity), we have to set the
      ! intensity value artifically to 1.0 to have it not influence the final value.
      ! We are only ever accessing I_layer to max_atoms, so this has no effect on the rest.
      I_layer(box_atom(max_atoms + 1), 1:number_energies) = 1.0_dp
    end if

    if (allocated(reflect_photo)) then
      deallocate (reflect_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate reflect_photo')
    end if

    if (allocated(atoms_per_layer)) then
      deallocate (atoms_per_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate atoms_per_layer')
    end if
    if (on_root .and. iprint > 2) then
      write (stdout, '(1x,a78)') '+----------------------- Printing Intensity per Layer -----------------------+'
      write (stdout, '(9999(es15.8))') ((I_layer(num_layer, i), num_layer=1, max_layer), i=1, number_energies)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if
    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------- Printing Intensity per Layer -----------------------+'
      write (stdout, '(1x,I4,1x,I4,1x)') jdos_nbins, max_layer
      write (stdout, '(9999(es15.8))') ((I_layer(num_layer, i), num_layer=1, max_layer), i=1, number_energies)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

  end subroutine calc_absorp_layer

  subroutine effect_wf

    !photo_elec_field given in eV/A

    use od_parameters, only: photo_work_function, photo_elec_field
    use od_electronic, only: efermi
    use od_constants, only: pi, epsilon_zero
    implicit none

    work_function_eff = photo_work_function - sqrt(photo_elec_field/(4*pi*epsilon_zero*1E-4))

    evacuum_eff = work_function_eff + efermi

  end subroutine effect_wf

  subroutine calc_field_emission
    !!*This subroutine calculates the Schottky effect
    ! parameter photo_elec_field given in V/m
    use od_cell, only: num_kpoints_on_node
    use od_parameters, only: photo_work_function, photo_elec_field, photo_temperature
    use od_electronic, only: efermi, band_energy, nbands, nspins
    use od_io, only: io_error
    use od_comms, only: my_node_id, comms_reduce
    use od_constants, only: pi, epsilon_zero, kB, e_charge, b_factor, p1, p2, p3, p4, q1, q2, q3, q4
    implicit none
    integer :: ierr
    real(kind=dp), allocatable, dimension(:, :, :) :: field_energy
    !real(kind=dp), allocatable, dimension(:, :, :) :: tunnel_prob
    real(kind=dp), allocatable, dimension(:, :, :) :: G
    real(kind=dp), allocatable, dimension(:, :, :) :: temp_emission
    real(kind=dp) :: fermi_dirac, barrier_height, argument, exponent
    real(kind=dp) :: l_prime, p_term, q_term, v_function, transmission_prob
    integer :: N_k, N_spin, n_eigen

    allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of field_emission failed')
    field_emission = 0.0_dp

    allocate (field_energy(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of field_energy failed')
    field_energy = 0.0_dp

    allocate (G(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of G failed')
    G = 0.0_dp

    allocate (temp_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of temp_emission failed')
    temp_emission = 0.0_dp

    evacuum = efermi + photo_work_function

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          barrier_height = photo_work_function - (band_energy(n_eigen, N_spin, N_k) - efermi)
          field_energy(n_eigen, N_spin, N_k) = abs(evacuum - band_energy(n_eigen, N_spin, N_k))
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a exponential value of ~1E+/-250
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 575.0_dp) then
            fermi_dirac = 0.0_dp
          elseif (argument .lt. -575.0_dp) then
            fermi_dirac = 1.0_dp
          else
            fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          if ((photo_elec_field/(4.0_dp*pi*epsilon_zero*1E-4_dp)*photo_elec_field) .lt. (field_energy(n_eigen, N_spin, N_k)**2)) &
            then
            if (barrier_height .le. 0.0_dp) then
              field_emission(n_eigen, N_spin, N_k) = 1.0_dp
            else
              l_prime = (e_charge/4*pi*epsilon_zero*1E-4_dp)*photo_elec_field*(1/barrier_height**2)
              p_term = 1.0_dp + (p1*l_prime) + (p2*l_prime**2.0_dp) + (p3*l_prime**3.0_dp) + (p4*l_prime**4.0_dp)
              q_term = q1 + (q2*l_prime) + (q3*l_prime**2.0_dp) + (q4*l_prime**3.0_dp)
              v_function = (1.0_dp - l_prime)*p_term + q_term*l_prime*log(l_prime)

              exponent = -1.0_dp*v_function*b_factor*(barrier_height**(2.0_dp/3.0_dp))*(1.0_dp/photo_elec_field)
              if (exponent .lt. -575.0_dp) then
                transmission_prob = 0.0_dp
              else
                transmission_prob = exp(exponent)
              end if
              field_emission(n_eigen, N_spin, N_k) = transmission_prob
            end if
          end if
          temp_emission(n_eigen, N_spin, N_k) = field_emission(n_eigen, N_spin, N_k)*fermi_dirac
        end do
      end do
    end do

    total_field_emission = sum(temp_emission(1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id)))/cell_area
    call comms_reduce(total_field_emission, 1, "SUM")

    deallocate (field_energy, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate field_energy')

    deallocate (G, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate G')

    deallocate (temp_emission, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate temp_emission')

    !ADD COMMENT
  end subroutine calc_field_emission

  !===============================================================================
  subroutine calc_angle
    !*******=======================================================================
    ! This subroutine calculates the photoemission angles theta and phi
    ! Theta: angle between the photoemitted electron and the surface normal
    ! Phi: angle between the photoemission direction and the x axis
    ! orig. Victor Chang, 7th February 2020
    ! parts rewritten Felix Mildner, after Mar 2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart
    use od_electronic, only: nbands, nspins, band_energy, band_gradient, elec_read_band_gradient, elec_read_band_curvature, &
                             band_curvature, photo_spectral_func, elec_read_spec_function
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_model, photo_momentum, devel_flag, iprint, photo_sf_max_vectors
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout, io_time
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: hbar, ev_to_j, j_to_ev, e_mass, rad_to_deg
    implicit none
    integer :: N_k, N_spin, n_eigen, ierr, sf_maxvec, gdx

    real(kind=dp), allocatable, dimension(:, :, :, :):: E_x
    real(kind=dp), allocatable, dimension(:, :, :, :):: E_y
    real(kind=dp) :: tol = 1.0E-10_dp
    real(kind=dp) :: time0, time1

    time0 = io_time()
    sf_maxvec = photo_sf_max_vectors

    if (.not. allocated(E_transverse)) then
      allocate (E_transverse(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_transverse failed')
    end if
    E_transverse = 0.0_dp

    if (.not. allocated(theta_arpes)) then
      allocate (theta_arpes(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of theta_arpes failed')
    end if
    ! Impossible value as default that is equal to no emission
    theta_arpes = 91.0_dp

    if (.not. allocated(theta_arpes_internal)) then
      allocate (theta_arpes_internal(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of theta_arpes_internal failed')
    end if
    ! Impossible value as default that is equal to no emission
    theta_arpes_internal = 91.0_dp

    if (.not. allocated(phi_arpes)) then
      allocate (phi_arpes(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of phi_arpes failed')
    end if
    phi_arpes = 0.0_dp
    if (.not. allocated(E_kinetic)) then
      allocate (E_kinetic(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_kinetic failed')
    end if
    E_kinetic = 0.0_dp

    allocate (E_x(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_x failed')
    E_x = 0.0_dp

    allocate (E_y(sf_maxvec, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_y failed')
    E_y = 0.0_dp

    if (index(photo_momentum, 'kp') > 0) then
      call elec_read_band_gradient
      call elec_read_band_curvature
    end if
    if (index(photo_momentum, 'operator') > 0) then
      call elec_read_band_gradient
    end if

    if (index(photo_momentum, 'crystal') > 0) call cell_calc_kpoint_r_cart

    if (index(photo_momentum, 'specfn') > 0) then
      call elec_read_spec_function(sf_maxvec)
    else
      if (.not. allocated(photo_spectral_func)) then
        allocate (photo_spectral_func(3, 1, nbands, nspins, num_kpoints_on_node(my_node_id)))
        photo_spectral_func = 1.0_dp
      end if
    end if

    if ((index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) .or. &
        (index(devel_flag, 'print_qe_matrix_full') > 0 .and. on_root) &
        .or. (index(devel_flag, 'print_qe_matrix_reduced') > 0 .and. on_root)) then
      call cell_calc_kpoint_r_cart
      write (stdout, '(a78)') "+---------------- Printing K-Points in Cartesian Coordinates ----------------+"
      do N_k = 1, num_kpoints_on_node(my_node_id)
        write (stdout, '(3(1x,E22.15))') kpoint_r_cart(:, N_k)
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          do gdx = 1, photo_sf_max_vectors
            if (index(photo_momentum, 'kp') > 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = abs &
                                               (0.5_dp*(1/(band_curvature(n_eigen, 1, 1, N_k, N_spin)*ev_to_j*1E-20/(hbar**2)))* &
                                                (band_gradient(n_eigen, 1, N_k, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = abs &
                                               (0.5_dp*(1/(band_curvature(n_eigen, 2, 2, N_k, N_spin)*ev_to_j*1E-20/(hbar**2)))* &
                                                (band_gradient(n_eigen, 2, N_k, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
            end if
            if (index(photo_momentum, 'specfn') > 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))* &
                                                ((photo_spectral_func(1, gdx, n_eigen, N_spin, N_k)*1E+10)**2))*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))* &
                                                ((photo_spectral_func(2, gdx, n_eigen, N_spin, N_k)*1E+10)**2))*j_to_ev
            end if
            if (index(photo_momentum, 'crystal') > 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(1, N_k)*1E+10)**2))*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(2, N_k)*1E+10)**2))*j_to_ev
            end if
            if (index(photo_momentum, 'operator') > 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = abs &
                                               (0.5_dp*e_mass* &
                                                (band_gradient(n_eigen, 1, N_k, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = abs &
                                               (0.5_dp*e_mass* &
                                                (band_gradient(n_eigen, 2, N_k, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
            end if
            E_transverse(gdx, n_eigen, N_spin, N_k) = E_x(gdx, n_eigen, N_spin, N_k) + E_y(gdx, n_eigen, N_spin, N_k)
          end do
        end do
      end do
    end do

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          do gdx = 1, photo_sf_max_vectors
            if ((abs(E_x(gdx, n_eigen, N_spin, N_k)) .lt. tol) .and. (abs(E_y(gdx, n_eigen, N_spin, N_k)) .lt. tol)) then
              phi_arpes(gdx, n_eigen, N_spin, N_k) = 0.0_dp
            elseif ((abs(E_y(gdx, n_eigen, N_spin, N_k)) .lt. tol)) then
              phi_arpes(gdx, n_eigen, N_spin, N_k) = 90.0_dp
            else
              phi_arpes(gdx, n_eigen, N_spin, N_k) = atan(E_x(gdx, n_eigen, N_spin, N_k)/E_y(gdx, n_eigen, N_spin, N_k))*rad_to_deg
            end if
          end do
        end do
      end do
    end do
    ! theta is the angle between emitted electron and the surface normal
    ! 3 Step Model - calculating the final energy of the electrons as the FINAL STATE ENERGY
    if (index(photo_model, '3step') > 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            ! total kinetic energy after emission from the excited band
            E_kinetic(1:photo_sf_max_vectors, n_eigen, N_spin, N_k) = (band_energy(n_eigen, N_spin, N_k) - evacuum_eff)
            do gdx = 1, photo_sf_max_vectors
              if (E_kinetic(gdx, n_eigen, N_spin, N_k) .lt. E_transverse(gdx, n_eigen, N_spin, N_k)) cycle
              ! Angle of electron outside material, after passing the surface and loosing E(work_function)
              ! acos(E_ortho/E_kinetic)
              theta_arpes(gdx, n_eigen, N_spin, N_k) = (acos((E_kinetic(gdx, n_eigen, N_spin, N_k) - &
                                                              E_transverse(gdx, n_eigen, N_spin, N_k))/ &
                                                             E_kinetic(gdx, n_eigen, N_spin, N_k)))*rad_to_deg
              ! angle of electron within material, before passing the surface
              theta_arpes_internal(gdx, n_eigen, N_spin, N_k) = (acos((E_kinetic(gdx, n_eigen, N_spin, N_k) + work_function_eff &
                                                                       - E_transverse(gdx, n_eigen, N_spin, N_k))/ &
                                                                      (E_kinetic(gdx, n_eigen, N_spin, N_k) + work_function_eff))) &
                                                                *rad_to_deg
            end do
          end do
        end do
      end do
      ! 1 Step Model - calculating the final energy as E_band + E_photon
    else if (index(photo_model, '1step') > 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            do gdx = 1, photo_sf_max_vectors
              ! total kinetic energy after emission and passing through work function potential step
              E_kinetic(gdx, n_eigen, N_spin, N_k) = (band_energy(n_eigen, N_spin, N_k) + temp_photon_energy - evacuum_eff)

              ! E_kinetic is the final kinetic energy of the electron after emission
              if (E_kinetic(gdx, n_eigen, N_spin, N_k) .lt. E_transverse(gdx, n_eigen, N_spin, N_k)) cycle
              ! Angle of electron outside material, after passing the surface and loosing E(work_function)
              ! acos(E_ortho/E_kinetic)
              theta_arpes(gdx, n_eigen, N_spin, N_k) = (acos((E_kinetic(gdx, n_eigen, N_spin, N_k) - &
                                                              E_transverse(gdx, n_eigen, N_spin, N_k))/ &
                                                             E_kinetic(gdx, n_eigen, N_spin, N_k)))*rad_to_deg
              ! Angle of electron within material, before passing the surface
              theta_arpes_internal(gdx, n_eigen, N_spin, N_k) = (acos((E_kinetic(gdx, n_eigen, N_spin, N_k) + work_function_eff &
                                                                       - E_transverse(gdx, n_eigen, N_spin, N_k)) &
                                                                      /(E_kinetic(gdx, n_eigen, N_spin, N_k) + &
                                                                        work_function_eff)))*rad_to_deg
            end do
          end do
        end do
      end do
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------ Printing Transverse Energy ------------------------+'
      write (stdout, '(3(1x,I4))') shape(E_transverse)
      write (stdout, '(3(1x,I4))') nbands, num_kpoints_on_node(my_node_id), nspins
      write (stdout, '(9999(es15.8))') ((((E_transverse(gdx, n_eigen, N_spin, N_k), gdx=1, photo_sf_max_vectors), &
                                          N_spin=1, nspins), N_k=1, num_kpoints_on_node(my_node_id)), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    deallocate (E_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_y')

    deallocate (E_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_x')

    if (allocated(band_curvature)) then
      deallocate (band_curvature, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate band_curvature')
    end if

    if (allocated(band_gradient)) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate band_gradient')
    end if

    if (allocated(kpoint_r_cart)) then
      deallocate (kpoint_r_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate kpoint_r_cart')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate Photoemission Angle', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_angle

  subroutine calc_electron_esc
    !! This subroutine calculates the electron escape probability for each of the layers
    use od_constants, only: dp, deg_to_rad, bohr2ang, H2eV, pi
    use od_electronic, only: nbands, nspins, band_energy
    use od_cell, only: num_kpoints_on_node, atoms_pos_cart_photo, atoms_label_tmp
    use od_io, only: io_error, stdout, io_time
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_imfp_value, photo_imfp_choice, devel_flag, iprint, photo_sf_max_vectors
    implicit none
    integer :: atom, N_k, N_spin, n_eigen, ierr, i, gdx
    real(kind=dp) :: tolerance
    real(kind=dp) :: exponent, time0, time1, scale_factor, scaled_x, g1, g2

    tolerance = 0.1E-11_dp
    time0 = io_time()
    allocate (new_atom_coordinates(3, max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of new_atom_coordinates failed')

    !Redefine new z coordinates where the first layer is at z=0
    new_atom_coordinates = atoms_pos_cart_photo
    do atom = 1, max_atoms
      new_atom_coordinates(3, atom_order(atom)) = atoms_pos_cart_photo(3, atom_order(atom)) - &
                                                  (atoms_pos_cart_photo(3, atom_order(1)))
    end do

    if (.not. allocated(electron_esc)) then
      allocate (electron_esc(photo_sf_max_vectors, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    end if
    electron_esc = 0.0_dp

    if (.not. allocated(atom_imfp)) then
      allocate (atom_imfp(max_atoms), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_electron_esc_list - allocation of atom_imfp failed')
    end if
    atom_imfp = 0.0_dp
    ! TODO: Must rename this to agree with box thicknesses
    if (index(photo_imfp_choice, 'layers') > 0) then
      ! Calculate the mean thickness of the atoms in a layer
      do atom = 1, max_atoms
        thickness_layer(layer(atom)) = thickness_layer(layer(atom)) + thickness_atom(atom)
      end do
      do i = 1, max_layer
        thickness_layer(i) = thickness_layer(i)/atoms_per_layer(i)
      end do
      if (on_root) then
        write (stdout, '(1x,a78)') '+--------------- User Supplied and Calculated IMFP Constants ----------------+'
        write (stdout, '(1x,a78)') '| Atom | Atom Order | Layer | Layer Thickness | User Input IMFP | Calc. IMFP |'
      end if

      ! Calculate the layer dependent imfp constant as a list for each layer
      do atom = 1, max_atoms
        do i = 1, layer(atom)
          atom_imfp(atom) = atom_imfp(atom) + thickness_layer(i)*photo_imfp_value(i)
        end do
        atom_imfp(atom) = atom_imfp(atom)/sum(thickness_layer(1:i - 1))
        if (on_root) then
          write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
            layer(atom), thickness_layer(layer(atom)), photo_imfp_value(layer(atom)), atom_imfp(atom), "    |"
225       format(1x, a1, a4, 6x, I3, 8x, I3, 6x, E14.6E3, 3x, F11.4, 3x, F11.4, a5)
        end if
      end do
      if (on_root) write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    else if (index(photo_imfp_choice, 'const') > 0) then
      atom_imfp = photo_imfp_value(1)
    else if (index(photo_imfp_choice, 'curve') > 0) then
      if (.not. allocated(band_imfp)) then
        allocate (band_imfp(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_electron_esc_list - allocation of atom_imfp failed')
      end if
      band_imfp = 0.0_dp
      scale_factor = 6.9_dp
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            scaled_x = (band_energy(n_eigen, N_spin, N_k)/scale_factor) + 1
            if ((1.0_dp - scaled_x) > 1E-10_dp) cycle
            g1 = LOG(scaled_x - 1.0_dp) + ((8.0_dp/3.0_dp) - 2.0_dp*LOG(2.0_dp))
            if (scaled_x < 2.0_dp) then
              g2 = g2 + (2.0_dp/3.0_dp)*(SQRT(2.0_dp - scaled_x)**(3.0_dp))
              g2 = g2 + (2.0_dp*SQRT(2.0_dp - scaled_x))
              g2 = g2 + LOG(ABS((SQRT(2.0_dp - scaled_x) - 1.0_dp)/(SQRT(2.0_dp - scaled_x) + 1.0_dp)))
            else
              g2 = 0.0_dp
            end if
            band_imfp(n_eigen, N_spin, N_k) = bohr2ang*(4.0_dp*pi/3.0_dp)*(scaled_x/(g1 - g2))*(SQRT(2.0_dp*scale_factor/H2eV))
            ! write (stdout, *) "scaled_x", scaled_x, "g1", g1, "g2", g2, band_imfp(n_eigen, N_spin, N_k)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+------------------ IMFP Values from Energy Dependent Curve -----------------+'
      write (stdout, *) 'min', minval(band_imfp), 'max', maxval(band_imfp)
    end if

    if ((index(photo_imfp_choice, 'const') > 0) .or. (index(photo_imfp_choice, 'layers') > 0)) then
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              do gdx = 1, photo_sf_max_vectors
                ! is the emission possible?
                if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. tolerance) then
                  ! The electron's kinetic energy inside the material is higher, than after the emission
                  ! through the surface. Thus follows an angle closer to normal direction and one needs
                  ! the internal theta angle.
                  exponent = (new_atom_coordinates(3, atom_order(atom))/ &
                              cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/atom_imfp(atom)
                  if (exponent .gt. -575.0_dp) then
                    electron_esc(gdx, n_eigen, N_spin, N_k, atom) = exp(exponent)
                  else
                    electron_esc(gdx, n_eigen, N_spin, N_k, atom) = 0.0_dp
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    else if ((index(photo_imfp_choice, 'curve') > 0)) then
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              do gdx = 1, photo_sf_max_vectors
                if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. tolerance) then
                  exponent = (new_atom_coordinates(3, atom_order(atom))/ &
                              cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/band_imfp(n_eigen, N_spin, N_k)
                  if ((exponent .gt. -575.0_dp) .and. (exponent .lt. 575.0_dp)) then
                    electron_esc(gdx, n_eigen, N_spin, N_k, atom) = exp(exponent)
                  else if (exponent .gt. -575.0_dp) then
                    electron_esc(gdx, n_eigen, N_spin, N_k, atom) = 1.0_dp
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------- Printing P(Escape) per Layer -----------------------+'
      write (stdout, 125) shape(electron_esc)
      write (stdout, 125) nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms
125   format(4(1x, I4))
      write (stdout, '(9999(es15.8))') (((((electron_esc(gdx, n_eigen, N_spin, N_k, atom), gdx=1, photo_sf_max_vectors), &
                                           atom=1, max_atoms), N_k=1, num_kpoints_on_node(my_node_id)), &
                                         N_spin=1, nspins), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate Photoemission Escape', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_electron_esc

  subroutine bulk_emission
    !! This subroutine calculates the contribution from the approximated bulk material
    use od_constants, only: dp, deg_to_rad
    use od_electronic, only: nbands, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root, comms_reduce, comms_bcast
    use od_parameters, only: photo_imfp_value, photo_imfp_choice, photo_bulk_cutoff, iprint, photo_sf_max_vectors
    use od_io, only: io_error, io_time, stdout
    implicit none
    real(kind=dp), dimension(:), allocatable :: bulk_light_tmp
    integer :: N_k, N_spin, n_eigen, i, num_layers, ierr, gdx
    real(kind=dp) :: exponent, time0, time1, band_imfp_max

    time0 = io_time()

235 format(1x, a1, 5x, a8, I3, 5x, a10, E13.6E2, 2x, a8, E13.6E2, 9x, a1)
    if (.not. new_geom_choice) then

      if (index(photo_imfp_choice, 'layers') > 0) then
        num_layers = int((atom_imfp(max_atoms)*photo_bulk_cutoff)/thickness_atom(max_atoms))
      else if (index(photo_imfp_choice, 'const') > 0) then
        num_layers = int((photo_imfp_value(1)*photo_bulk_cutoff)/thickness_atom(max_atoms))
      else if (index(photo_imfp_choice, 'curve') > 0) then
        band_imfp_max = maxval(band_imfp)
        call comms_reduce(band_imfp_max, 1, 'MAX')
        call comms_bcast(band_imfp_max, 1)
        num_layers = min(5000, int((band_imfp_max*photo_bulk_cutoff)/thickness_atom(max_atoms)))
      end if

      allocate (bulk_light_tmp(num_layers), stat=ierr)
      if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_light_tmp failed')
      bulk_light_tmp = 0.0_dp

      bulk_light_tmp(1) = I_layer(layer(max_atoms), current_photo_energy_index)* &
                          exp(-(absorp_photo(max_atoms, current_photo_energy_index)*thickness_atom(max_atoms)*1E-10))
      do i = 2, num_layers
        bulk_light_tmp(i) = bulk_light_tmp(i - 1)* &
                            exp(-(absorp_photo(max_atoms, current_photo_energy_index)*i*thickness_atom(max_atoms)*1E-10))
      end do

      if ((index(photo_imfp_choice, 'layers') > 0) .or. (index(photo_imfp_choice, 'const') > 0)) then
        do i = 1, num_layers
          do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
            do N_spin = 1, nspins                    ! Loop over spins
              do n_eigen = 1, nbands
                do gdx = 1, photo_sf_max_vectors
                  if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                    exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms)/ &
                                cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/atom_imfp(max_atoms)
                    ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                    ! As exp(-575) is ~1E-250, this should be more than enough precision.
                    if (exponent .gt. -575.0_dp) then
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                        electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      else if (index(photo_imfp_choice, 'curve') > 0) then
        do i = 1, num_layers
          do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
            do N_spin = 1, nspins                    ! Loop over spins
              do n_eigen = 1, nbands
                do gdx = 1, photo_sf_max_vectors
                  if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                    exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms)/ &
                                cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/band_imfp(n_eigen, N_spin, N_k)
                    ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                    ! As exp(-575) is ~1E-250, this should be more than enough precision.
                    if (exponent .gt. -575.0_dp) then
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                        electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end if ! If statement

      if (iprint .gt. 1 .and. on_root) then
        ! write out the bulk properties
        write (stdout, '(1x,a78)') '+----------------------- Approximated Bulk Slab Info ------------------------+'
        ! write out num_layers
        write (stdout, '(1x,a1,5x,a18,1x,a1,1x,I5,45x,a1)') '+', 'Number Bulk layers', '=', num_layers, '+'
        ! write out the total volume + volume per layer
        write (stdout, '(1x,a1,5x,a14,5x,a1,1x,F10.4,40x,a1)') '+', 'Vol. per layer', '=', &
          thickness_atom(max_atoms)*cell_area, '+'
        write (stdout, '(1x,a1,5x,a12,7x,a1,1x,F10.4,40x,a1)') '+', 'Total Volume', '=', num_layers* &
          thickness_atom(max_atoms)*cell_area, '+'
        write (stdout, '(1x,a78)') '+---- P_esc values for an electron with E = E_fermi and E_transverse = 0 ----+'
        ! write out bulk_light_tmp
        if (num_layers .lt. 6) then
          do i = 1, num_layers
            exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms))/atom_imfp(max_atoms)
            ! This makes sure, that exp(exponent) does not underflow the dp fp value.
            ! As exp(-575) is ~1E-250, this should be more than enough precision.
            if (exponent .gt. -575.0_dp) then
              exponent = exp(exponent)
            else
              exponent = 0.0_dp
            end if
            write (stdout, 235) '+', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '+'
          end do
        else
          do i = 1, num_layers
            exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms))/atom_imfp(max_atoms)
            ! This makes sure, that exp(exponent) does not underflow the dp fp value.
            ! As exp(-575) is ~1E-250, this should be more than enough precision.
            if (exponent .gt. -575.0_dp) then
              exponent = exp(exponent)
            else
              exponent = 0.0_dp
            end if
            if (i .le. 3 .or. i .gt. num_layers - 3) then
              write (stdout, 235) '+', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '+'
            elseif (i .eq. 4) then
              write (stdout, '(1x,a1,35x,a6,35x,a1)') '+', '......', '+'
            end if
          end do
        end if
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if

    else ! choice of the new geometry definition (hardcoded)

      if (index(photo_imfp_choice, 'layers') > 0) then
        num_layers = int((atom_imfp(max_atoms)*photo_bulk_cutoff)/box_height)
      else if (index(photo_imfp_choice, 'const') > 0) then
        num_layers = int((photo_imfp_value(1)*photo_bulk_cutoff)/box_height)
      else if (index(photo_imfp_choice, 'curve') > 0) then
        band_imfp_max = maxval(band_imfp)
        call comms_reduce(band_imfp_max, 1, 'MAX')
        call comms_bcast(band_imfp_max, 1)
        num_layers = min(5000, int((band_imfp_max*photo_bulk_cutoff)/box_height))
      end if

      allocate (bulk_light_tmp(num_layers), stat=ierr)
      if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_light_tmp failed')
      bulk_light_tmp = 0.0_dp

      bulk_light_tmp(1) = I_layer(box_atom(max_atoms), current_photo_energy_index)* &
                          exp(-(absorp_photo(box_atom(max_atoms), current_photo_energy_index)*box_height*1E-10))
      do i = 2, num_layers
        bulk_light_tmp(i) = bulk_light_tmp(i - 1)* &
                            exp(-(absorp_photo(box_atom(max_atoms), current_photo_energy_index)*i*box_height*1E-10))
      end do

      if ((index(photo_imfp_choice, 'layers') > 0) .or. (index(photo_imfp_choice, 'const') > 0)) then
        do i = 1, num_layers
          do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
            do N_spin = 1, nspins                    ! Loop over spins
              do n_eigen = 1, nbands
                do gdx = 1, photo_sf_max_vectors
                  if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                    exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*box_height/ &
                                cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/atom_imfp(max_atoms)
                    ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                    ! As exp(-575) is ~1E-250, this should be more than enough precision.
                    if (exponent .gt. -575.0_dp) then
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                        electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      else if (index(photo_imfp_choice, 'curve') > 0) then
        do i = 1, num_layers
          do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
            do N_spin = 1, nspins                    ! Loop over spins
              do n_eigen = 1, nbands
                do gdx = 1, photo_sf_max_vectors
                  if (cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                    exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*box_height/ &
                                cos(theta_arpes_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/band_imfp(n_eigen, N_spin, N_k)
                    ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                    ! As exp(-575) is ~1E-250, this should be more than enough precision.
                    if (exponent .gt. -575.0_dp) then
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                        electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end if ! If statement

      if (iprint .gt. 1 .and. on_root) then
        ! write out the bulk properties
        write (stdout, '(1x,a78)') '+---------------------- Bulk Approximation Slab Info ------------------------+'
        ! write out num_layers
        write (stdout, '(1x,a1,5x,a18,1x,a1,1x,I5,45x,a1)') '+', 'Number Bulk layers', '=', num_layers, '+'
        ! write out the total volume + volume per layer
        write (stdout, '(1x,a1,5x,a14,5x,a1,1x,F10.4,40x,a1)') '+', 'Vol. per layer', '=', box_volume, '+'
        write (stdout, '(1x,a1,5x,a12,7x,a1,1x,F10.4,40x,a1)') '+', 'Total Volume', '=', num_layers*box_volume, '+'
        write (stdout, '(1x,a78)') '+---- P_esc values for an electron with E = E_fermi and E_transverse = 0 ----+'
        ! write out bulk_light_tmp
        if (num_layers .lt. 6) then
          do i = 1, num_layers
            exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*box_height)/atom_imfp(max_atoms)
            ! This makes sure, that exp(exponent) does not underflow the dp fp value.
            ! As exp(-575) is ~1E-250, this should be more than enough precision.
            if (exponent .gt. -575.0_dp) then
              exponent = exp(exponent)
            else
              exponent = 0.0_dp
            end if
            write (stdout, 235) '+', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '+'
          end do
        else
          do i = 1, num_layers
            exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*box_height)/atom_imfp(max_atoms)
            ! This makes sure, that exp(exponent) does not underflow the dp fp value.
            ! As exp(-575) is ~1E-250, this should be more than enough precision.
            if (exponent .gt. -575.0_dp) then
              exponent = exp(exponent)
            else
              exponent = 0.0_dp
            end if
            if (i .le. 3 .or. i .gt. num_layers - 3) then
              write (stdout, 235) '+', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '+'
            elseif (i .eq. 4) then
              write (stdout, '(1x,a1,35x,a6,35x,a1)') '+', '......', '+'
            end if
          end do
        end if ! If statement printing of slab light intensities formatting
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if ! If statement extra printing of slab data
    end if ! If statement for geometry choice

    deallocate (bulk_light_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate bulk_light_tmp')

    deallocate (new_atom_coordinates, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate new_atom_coordinates')

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a38,21x,f11.3,a8)') '+ Time to calculate Bulk Photoemission', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine bulk_emission

  subroutine calc_ds_like_model
    !*===============================================================================
    ! This subroutine calculates the QE using a simplified model following a Dowell-
    ! Schmerge like Model by Saha et al.
    ! Felix Mildner, May 2024
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient, &
                             elec_read_band_curvature
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, comms_bcast
    use od_parameters, only: photo_temperature, devel_flag, iprint, num_exclude_bands, &
                             exclude_bands, photo_model
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, io_time, seedname, io_date
    use od_jdos_utils, only: jdos_utils_calculate, setup_energy_scale
    use od_constants, only: pi, kB, inv_sqrt_two_pi
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:) :: qe_k_temp
    real(kind=dp) :: width, norm_vac, qe_factor, argument, time0, time1, final_fd, initial_fd, excess_energy
    integer :: N_k, N_spin, n_eigen, n_eigen_final, ierr, i, qe_unit, token, inode
    character(len=10)                           :: char_e
    character(len=99)                           :: filename
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    width = (1.0_dp/11604.45_dp)*photo_temperature
    qe_factor = 1.0_dp/(cell_area)
    norm_vac = inv_sqrt_two_pi/width

    time0 = io_time()

    if (.not. allocated(qe_tsm)) then
      allocate (qe_tsm(nbands, nbands, nspins, num_kpoints_on_node(my_node_id), 3), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of qe_tsm failed')
    end if
    qe_tsm = 0.0_dp

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    call photo_calculate_delta(delta_temp, .false.)

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+-------------------------- Calculating DS Like QE --------------------------+'
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a an exponent value of ~1E+/-250
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 575.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -575.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if
        end do
      end do
    end do

    call setup_energy_scale(E)
    i = 0
    if (on_root) write (stdout, *) '***   Calculating a simplified Dowell Schmerge like model for PE   ***'

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
          if (num_exclude_bands .gt. 1) then
            if (any(exclude_bands == n_eigen_final)) then
              cycle
            end if
          end if
          excess_energy = band_energy(n_eigen_final, N_spin, N_k) - evacuum_eff
          excess_energy = max(excess_energy, 0.0_dp)
          final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
          do n_eigen = 1, n_eigen_final - 1
            initial_fd = fermi_dirac(n_eigen, N_spin, N_k)
            ! Calculating the QE denominator
            qe_tsm(n_eigen, n_eigen_final, N_spin, N_k, 1) = delta_temp(n_eigen, n_eigen_final, N_spin, N_k)* &
                                                             electrons_per_state*kpoint_weight(N_k)* &
                                                             final_fd*initial_fd
            ! Calculating the QE numerator and MTE denominator
            qe_tsm(n_eigen, n_eigen_final, N_spin, N_k, 2) = delta_temp(n_eigen, n_eigen_final, N_spin, N_k)* &
                                                             electrons_per_state*kpoint_weight(N_k)* &
                                                             final_fd*initial_fd*excess_energy
            ! Calculating the MTE numerator
            qe_tsm(n_eigen, n_eigen_final, N_spin, N_k, 3) = delta_temp(n_eigen, n_eigen_final, N_spin, N_k)* &
                                                             electrons_per_state*kpoint_weight(N_k)* &
                                                             final_fd*initial_fd*excess_energy**2
          end do
        end do
      end do
    end do

    if (allocated(delta_temp)) then
      deallocate (delta_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - failed to deallocate delta_temp')
    end if

    if (allocated(fermi_dirac)) then
      deallocate (fermi_dirac, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - failed to deallocate fermi_dirac')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a41,18x,f11.3,a8)') '+ Time to calculate DS like Photoemission', time1 - time0, ' (sec) +'
    end if

    if (index(devel_flag, 'print_kpt_qe_data') > 0) then
      if (on_root) then
        qe_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_k_point_QE.dat'
        write (stdout, *) 'opening file'
        open (unit=qe_unit, action='write', file=filename)
        write (qe_unit, *) '# The k point dependent QE values'
        call io_date(cdate, ctime)
        write (qe_unit, *) '## OptaDOS Photoemission: Printing QE K point Data on ', cdate, ' at ', ctime
      end if

      allocate (qe_k_temp(num_kpoints_on_node(0)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_three_step_model - failed to allocate qe_k_temp on root')
      token = -1

      ! allocate and sum the 3step qe matrix on non-root
      if (.not. on_root) then
        do N_k = 1, num_kpoints_on_node(my_node_id)
          qe_k_temp(N_k) = sum(qe_tsm(:, :, :, N_k, :))
        end do
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that node
        call comms_send(qe_k_temp(1), num_kpoints_on_node(my_node_id), 0)
        ! - send token back to root node
        call comms_send(token, 1, 0)
      end if

      if (on_root) then
        do inode = 1, num_nodes - 1
          ! - send to the token to notes in turn
          call comms_send(token, 1, inode)
          ! - receive the qe_matrix from the other notes and write it to the file
          call comms_recv(qe_k_temp(1), num_kpoints_on_node(inode), inode)
          ! write out the qe_matrix to the file
          do N_k = 1, num_kpoints_on_node(inode)
            write (qe_unit, *) qe_k_temp(N_k)
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (qe_unit, *) sum(qe_tsm(:, :, :, N_k, :))
        end do
        close (unit=qe_unit)
      end if
      deallocate (qe_k_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - failed to deallocate qe_k_temp')
    end if

  end subroutine calc_ds_like_model

  !===============================================================================
  subroutine calc_three_step_model
    !*===============================================================================
    ! This subroutine calculates the QE using the three step model.
    ! Victor Chang, 7th February 2020
    ! edited by Felix Mildner, 03/2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient, &
                             elec_read_band_curvature, transmit_prob, elec_read_transmit_prob, photo_spectral_func
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, comms_bcast
    use od_parameters, only: scissor_op, photo_temperature, devel_flag, photo_photon_sweep, iprint, num_exclude_bands, &
                             exclude_bands, photo_model, photo_sf_max_vectors, photo_output
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, io_time, seedname, io_date
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB, inv_sqrt_two_pi
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:) :: qe_k_temp
    real(kind=dp) :: width, norm_vac, vacuum_gauss, transverse_gauss, qe_factor, argument, &
                     time0, time1, final_fd, initial_fd, temp_contribution
    integer :: N_k, N2, N_spin, n_eigen_init, n_eigen_final, atom, ierr, i, gdx, qe_unit, token, inode
    character(len=10)                           :: char_e
    character(len=99)                           :: filename
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    width = (1.0_dp/11604.45_dp)*photo_temperature
    qe_factor = 1.0_dp/(cell_area)
    norm_vac = inv_sqrt_two_pi/width

    time0 = io_time()

    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of field_emission failed')
      field_emission = 0.0_dp
    end if

    if (.not. allocated(qe_tsm)) then
      allocate (qe_tsm(nbands, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of qe_tsm failed')
    end if
    qe_tsm = 0.0_dp

    if (.not. allocated(te_tsm)) then
      allocate (te_tsm(nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of te_tsm failed')
    end if
    te_tsm = 0.0_dp

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    if (index(devel_flag, 'no_transmit') > 0) then
      if (.not. allocated(transmit_prob)) then
        allocate (transmit_prob(nbands, nspins, num_kpoints_on_node(my_node_id)))
        if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
      end if
      transmit_prob = 1.0_dp
    else
      call elec_read_transmit_prob()
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+----------------- Printing Matrix Weights in 3Step Function ----------------+'
      write (stdout, '(5(1x,I4))') shape(matrix_weights)
      write (stdout, '(5(1x,I4))') nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
      do N2 = 1, N_geom
        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(es15.8))') ((matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, N2), &
                                                n_eigen_final=1, nbands), n_eigen_init=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    call photo_calculate_delta(delta_temp, .false.)

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 3Step QE ---------------------------+'
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+---------------------- Printing Delta Function Values ----------------------+'
      write (stdout, '(5(1x,I4))') shape(delta_temp)
      write (stdout, '(5(1x,I4))') nbands, nbands, num_kpoints_on_node(my_node_id), nspins
      do N_spin = 1, nspins
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (stdout, '(99999(es15.8))') ((delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k), &
                                              n_eigen_final=1, nbands), n_eigen_init=1, nbands)
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) &
      then
      i = 17 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
      write (stdout, '(1x,a78)') '+------------ Printing list of values going into 3step QE Values ------------+'
      write (stdout, '(14(1x,a17))') 'calced_qe_value', 'initial_state_energy', 'final_state_energy', 'spectral_func', &
        'matrix_weights', &
        'delta_temp', 'electron_esc', 'kpoint_weight', 'I_layer', 'transverse_gauss', 'vacuum_gauss', 'fermi_dirac', &
        'pdos_weights_atoms', 'pdos_weights_k_band'
      write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', max_atoms, nbands, nbands, nspins, num_kpoints_on_node(my_node_id), i
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen_init = 1, nbands
          argument = (band_energy(n_eigen_init, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a an exponent value of ~1E+/-250
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 575.0_dp) then
            fermi_dirac(n_eigen_init, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -575.0_dp) then
            fermi_dirac(n_eigen_init, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen_init, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if
        end do
      end do
    end do

    do atom = 1, max_atoms
      if (iprint > 2 .and. on_root) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', &
          "Calculating atom ", atom, " of", max_atoms, "<-- QE-3S |"
      end if
      do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
            if (num_exclude_bands .gt. 1) then
              if (any(exclude_bands == n_eigen_final)) then
                cycle
              end if
            end if
            final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
            do n_eigen_init = 1, n_eigen_final - 1
              initial_fd = fermi_dirac(n_eigen_init, N_spin, N_k)
              do gdx = 1, photo_sf_max_vectors
                ! evacuum_eff = efermi + photo_work_function
                ! Is (photon_energy - transverse energy) > (work_function - E_field_lowering)
                ! Is the final kinetic energy > 0?
                if ((temp_photon_energy - E_transverse(gdx, n_eigen_init, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                  transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen_init, N_spin, N_k)), &
                                              width, (evacuum_eff - efermi))/norm_vac
                else
                  transverse_gauss = 1.0_dp
                end if
                if (band_energy(n_eigen_final, N_spin, N_k) .lt. evacuum_eff) then
                  vacuum_gauss = gaussian(band_energy(n_eigen_final, N_spin, N_k) + &
                                          scissor_op, width, evacuum_eff)/norm_vac
                else
                  vacuum_gauss = 1.0_dp
                end if

                !! this could be checked if it has an impact on the final value
                ! if (band_energy(n_eigen_final, N_spin, N_k) .lt. efermi) cycle
                if (enable_debug_output .and. index(devel_flag, 'reduced_pe') > 0) then
                  if (index(devel_flag, 'projected_pe') > 0) then
                    qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) = &
                      matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1)* &
                      delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k)* &
                      electrons_per_state*kpoint_weight(N_k)* &
                      (pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom))/ &
                       pdos_weights_k_band(n_eigen_init, N_spin, N_k))
                  else
                    qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) = &
                      matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1)* &
                      delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k)* &
                      electrons_per_state*kpoint_weight(N_k)
                  end if
                else
                  if (.not. new_geom_choice) then
                    temp_contribution = &
                      qe_factor*photo_spectral_func(3, gdx, n_eigen_init, N_spin, N_k) &
                      *matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1) &
                      *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                      *electron_esc(gdx, n_eigen_final, N_spin, N_k, atom) &
                      *transmit_prob(n_eigen_final, N_spin, N_k) &
                      *electrons_per_state*kpoint_weight(N_k) &
                      *I_layer(layer(atom), current_photo_energy_index) &
                      *transverse_gauss*vacuum_gauss*initial_fd*final_fd &
                      *(pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                        /pdos_weights_k_band(n_eigen_init, N_spin, N_k)) &
                      *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
                  else
                    temp_contribution = &
                      qe_factor*photo_spectral_func(3, gdx, n_eigen_init, N_spin, N_k) &
                      *matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1) &
                      *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                      *electron_esc(gdx, n_eigen_final, N_spin, N_k, atom) &
                      *transmit_prob(n_eigen_final, N_spin, N_k) &
                      *electrons_per_state*kpoint_weight(N_k) &
                      *(I_layer(box_atom(atom), current_photo_energy_index)) &
                      *transverse_gauss*vacuum_gauss*initial_fd*final_fd &
                      *(pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                        /pdos_weights_k_band(n_eigen_init, N_spin, N_k)) &
                      *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
                  end if
                  qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) = &
                    qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) &
                    + temp_contribution
                  te_tsm(n_eigen_init, N_spin, N_k, atom) = te_tsm(n_eigen_init, N_spin, N_k, atom) &
                                                            + temp_contribution*E_transverse(gdx, n_eigen_init, N_spin, N_k)
                end if
                if (temp_contribution .gt. 0.0_dp) then
                if (enable_debug_output .and. index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
                  write (stdout, '(5(1x,I4))') n_eigen_init, n_eigen_final, N_spin, N_k, atom
                  write (stdout, '(18(1x,E17.9E3))') qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom), &
                    band_energy(n_eigen_init, N_spin, N_k), &
                    band_energy(n_eigen_final, N_spin, N_k), photo_spectral_func(3, gdx, n_eigen_init, N_spin, N_k), &
                    matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1), &
                    delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k), electron_esc(gdx, n_eigen_final, N_spin, N_k, atom), &
                    transmit_prob(n_eigen_final, N_spin, N_k), kpoint_weight(N_k), &
                    I_layer(layer(atom), current_photo_energy_index), transverse_gauss, vacuum_gauss, initial_fd, final_fd, &
                    pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)), pdos_weights_k_band(n_eigen_init, N_spin, N_k)
                end if
                end if
              end do
            end do
          end do
        end do
      end do
    end do

    call photo_calculate_delta(delta_temp, .true.)

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
          final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
          do n_eigen_init = 1, n_eigen_final - 1
            initial_fd = fermi_dirac(n_eigen_init, N_spin, N_k)
            do gdx = 1, photo_sf_max_vectors
              if ((temp_photon_energy - E_transverse(gdx, n_eigen_init, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen_init, N_spin, N_k)), &
                                            width, (evacuum_eff - efermi))/norm_vac
              else
                transverse_gauss = 1.0_dp
              end if
              if ((band_energy(n_eigen_init, N_spin, N_k) + temp_photon_energy) .lt. evacuum_eff) then
                vacuum_gauss = gaussian((band_energy(n_eigen_init, N_spin, N_k) + temp_photon_energy) + &
                                        scissor_op, width, evacuum_eff)/norm_vac
              else
                vacuum_gauss = 1.0_dp
              end if
              temp_contribution = &
                (qe_factor*photo_spectral_func(3, gdx, n_eigen_init, N_spin, N_k) &
                 *matrix_weights(n_eigen_init, n_eigen_final, N_k, N_spin, 1) &
                 *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                 *electron_esc(gdx, n_eigen_final, N_spin, N_k, max_atoms + 1) &
                 *transmit_prob(n_eigen_final, N_spin, N_k) &
                 *electrons_per_state*kpoint_weight(N_k) &
                 *transverse_gauss*vacuum_gauss*initial_fd*final_fd &
                 *(pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(max_atoms)) &
                   /pdos_weights_k_band(n_eigen_init, N_spin, N_k))) &
                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
              qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) = &
                qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) &
                + temp_contribution
              te_tsm(n_eigen_init, N_spin, N_k, atom) = te_tsm(n_eigen_init, N_spin, N_k, atom) &
                                                        + temp_contribution*E_transverse(gdx, n_eigen_init, N_spin, N_k)
            end do
          end do
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (allocated(delta_temp)) then
      deallocate (delta_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate delta_temp')
    end if

    if (allocated(fermi_dirac)) then
      deallocate (fermi_dirac, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate fermi_dirac')
    end if

    if (allocated(transmit_prob) .and. index(photo_output, 'e_bind') .eq. 0) then
      deallocate (transmit_prob, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate transmit_prob')
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_qe_matrix_full') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------- Printing Full 3step QE Matrix ----------------------+'
      write (stdout, '(6(1x,I4))') shape(qe_tsm)
      write (stdout, '(6(1x,I4))') photo_sf_max_vectors, nbands, nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1
      do atom = 1, max_atoms + 1
        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(ES16.8E3))') ((qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom), &
                                                  n_eigen_final=1, nbands), n_eigen_init=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 3step Photoemission', time1 - time0, ' (sec) +'
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_kpt_qe_data') > 0) then
      if (on_root) then
        qe_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_k_point_QE.dat'
        write (stdout, *) 'opening file'
        open (unit=qe_unit, action='write', file=filename)
        write (qe_unit, *) '# The k point dependent QE values'
        call io_date(cdate, ctime)
        write (qe_unit, *) '## OptaDOS Photoemission: Printing QE K point Data on ', cdate, ' at ', ctime
      end if

      allocate (qe_k_temp(num_kpoints_on_node(0)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_three_step_model - failed to allocate qe_k_temp on root')
      token = -1

      ! allocate and sum the 3step qe matrix on non-root
      if (.not. on_root) then
        do N_k = 1, num_kpoints_on_node(my_node_id)
          qe_k_temp(N_k) = sum(qe_tsm(:, :, :, N_k, :))
        end do
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that node
        call comms_send(qe_k_temp(1), num_kpoints_on_node(my_node_id), 0)
        ! - send token back to root node
        call comms_send(token, 1, 0)
      end if

      if (on_root) then
        do inode = 1, num_nodes - 1
          ! - send to the token to notes in turn
          call comms_send(token, 1, inode)
          ! write(stdout, *) 'sent token to node and receiving data from', inode
          call comms_recv(qe_k_temp(1), num_kpoints_on_node(inode), inode)
          ! write out the qe_matrix to the file
          do N_k = 1, num_kpoints_on_node(inode)
            write (qe_unit, *) qe_k_temp(N_k)
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (qe_unit, *) sum(qe_tsm(:, :, :, N_k, :))
        end do
        close (unit=qe_unit)
      end if
      deallocate (qe_k_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate qe_k_temp')
    end if
  end subroutine calc_three_step_model

  !===============================================================================
  subroutine photo_calculate_delta(delta_temp, calculate_bulk)
    !*===============================================================================
    ! Wrapper around the delta function subroutine to pass correct arguments
    ! Victor Chang, 7th February 2020
    !===============================================================================
    use od_parameters, only: linear, fixed, adaptive, quad, iprint
    use od_electronic, only: elec_read_band_gradient, band_gradient, efermi_set
    use od_comms, only: on_root
    use od_io, only: stdout, io_error, io_time
    ! use od_cell, only: cell_volume
    use od_dos_utils, only: dos_utils_set_efermi
    use od_jdos_utils, only: setup_energy_scale, jdos_deallocate

    implicit none

    real(kind=dp) :: time0, time1
    integer       :: ierr

    logical, intent(in)                                  :: calculate_bulk
    real(kind=dp), intent(out), allocatable, optional    :: delta_temp(:, :, :, :)  !I've added this

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+---------------------- Calculate JDOS DELTA FUNCTION -----------------------+'
    end if

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if (quad .or. linear .or. adaptive) then
      if (.not. allocated(band_gradient)) call elec_read_band_gradient
    end if
    !-------------------------------------------------------------------------------
    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
    if (.not. efermi_set) call dos_utils_set_efermi

    time0 = io_time()

    call setup_energy_scale(E)

    if (fixed) then
      if (calculate_bulk) then
        call calculate_delta('f', delta_temp, .true.)
      else
        call calculate_delta('f', delta_temp, .false.)
      end if
    end if
    if (adaptive) then
      if (calculate_bulk) then
        call calculate_delta('a', delta_temp, .true.)
      else
        call calculate_delta('a', delta_temp, .false.)
      end if
    end if
    if (linear) then
      if (calculate_bulk) then
        call calculate_delta('l', delta_temp, .true.)
      else
        call calculate_delta('l', delta_temp, .false.)
      end if
    end if

    if (quad) then
      call io_error("quadratic broadening not implemented")
    end if

    call jdos_deallocate

    if (allocated(E)) then
      deallocate (E, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_calculate_delta - failed to deallocate E')
    end if

    if (allocated(band_gradient) .and. current_photo_energy_index .eq. number_energies) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_calculate_delta - failed to deallocate band_gradient')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a34,25x,f11.3,1x,a7)') &
        '+ Time to calculate Delta Function', time1 - time0, '(sec) +'
    end if
    !-------------------------------------------------------------------------------

  end subroutine photo_calculate_delta

  !===============================================================================
  subroutine calculate_delta(delta_type, delta_temp, calculate_bulk)
    !*===============================================================================
    ! This subroutine evaluates the delta function between the valence band
    ! and the conduction band using the method specified in the input.
    ! orig. Victor Chang, 7 February 2020
    ! edited by Felix Mildner, after March 2022
    !===============================================================================
    use od_comms, only: my_node_id, on_root
    use od_cell, only: num_kpoints_on_node, kpoint_grid_dim, recip_lattice
    use od_parameters, only: adaptive_smearing, fixed_smearing, iprint, finite_bin_correction, &
                             scissor_op, hybrid_linear_grad_tol, hybrid_linear, exclude_bands, &
                             num_exclude_bands, jdos_max_energy, photo_slab_max
    use od_io, only: io_error, stdout
    use od_electronic, only: band_gradient, nbands, band_energy, nspins, efermi
    use od_jdos_utils, only: jdos_nbins
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_constants, only: pi, inv_sqrt_two_pi
    implicit none

    integer :: ik, is, ib, jb, i, ierr
    real(kind=dp) :: cuml, width, adaptive_smearing_temp
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4), sub_cell_length(1:3)
    real(kind=dp), save                   :: delta_bins

    character(len=1), intent(in)                      :: delta_type
    real(kind=dp), intent(inout), allocatable, optional :: delta_temp(:, :, :, :)
    logical, intent(in)                               :: calculate_bulk

    logical :: linear, fixed, adaptive, force_adaptive
    real(kind=dp) :: half_slab_height, norm_width

    linear = .false.
    fixed = .false.
    adaptive = .false.

    select case (delta_type)
    case ("l")
      linear = .true.
    case ("a")
      adaptive = .true.
    case ("f")
      fixed = .true.
    case default
      call io_error(" ERROR : unknown jdos_type in calculate_delta ")
    end select

    width = 0.0_dp
    delta_bins = jdos_max_energy/real(jdos_nbins - 1, dp)
    half_slab_height = photo_slab_max - slab_middle_ref

    if (linear .or. adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:), dp)/2.0_dp
    if (adaptive .or. hybrid_linear) then
      do i = 1, 2
        sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
      end do
      if (calculate_bulk) then
        sub_cell_length(3) = sqrt(recip_lattice(3, 1)**2 + recip_lattice(3, 1)**2 + (pi/box_height)**2)*step(3)
      else
        sub_cell_length(3) = sqrt(recip_lattice(3, 1)**2 + recip_lattice(3, 1)**2 + (pi/half_slab_height)**2)*step(3)
      end if
      adaptive_smearing_temp = adaptive_smearing*sum(sub_cell_length)/3.0_dp
    end if

    if (fixed) width = fixed_smearing

    if (.not. allocated(delta_temp)) then
      allocate (delta_temp(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_delta - allocation of delta_temp failed')
    end if
    delta_temp = 0.0_dp

    do ik = 1, num_kpoints_on_node(my_node_id)
      if (iprint > 2 .and. on_root) then
        if (mod(real(ik, dp), 10.0_dp) == 0.0_dp) write (stdout, '(1x,a1,a38,i4,a3,i4,1x,a14,2x,a11)') '^', &
          "Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id), 'on this node.', "<-- Delta |"
      end if
      do is = 1, nspins
        occ_states: do ib = 1, nbands
          if (num_exclude_bands > 0) then
            if (any(exclude_bands == ib)) cycle
          end if
          if (band_energy(ib, is, ik) .ge. efermi) cycle occ_states
          unocc_states: do jb = 1, nbands
            if (band_energy(jb, is, ik) .lt. efermi) cycle unocc_states
            if (linear .or. adaptive) grad(:) = band_gradient(jb, :, ik, is) - band_gradient(ib, :, ik, is)

            ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
            ! adaptive smearing (and take advantage of FBCS if required).
            force_adaptive = .false.
            if (.not. fixed) then
              if (hybrid_linear .and. (hybrid_linear_grad_tol > sqrt(dot_product(grad, grad)))) force_adaptive = .true.
              if (linear .and. .not. force_adaptive) call doslin_sub_cell_corners(grad, step, band_energy(jb, is, ik) - &
                                                                                  band_energy(ib, is, ik) + scissor_op, EV)
              if (adaptive .or. force_adaptive) width = sqrt(dot_product(grad, grad))*adaptive_smearing_temp
            end if
            ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
            ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
            if (finite_bin_correction .and. (width < delta_bins)) width = delta_bins
            norm_width = inv_sqrt_two_pi/width

            ! The linear method has a special way to calculate the integrated dos
            ! we have to take account for this here.
            if (linear .and. .not. force_adaptive) then
              delta_temp(ib, jb, is, ik) = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), E(current_energy_index), cuml)
            else
              delta_temp(ib, jb, is, ik) = gaussian((band_energy(jb, is, ik) - band_energy(ib, is, ik)) + scissor_op, width, &
                                                    E(current_energy_index))
            end if

          end do unocc_states
        end do occ_states
      end do
    end do

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine calculate_delta

  !===============================================================================
  subroutine make_foptical_weights
    !*===============================================================================
    ! This subroutine calculates the optical matrix elements for the one step
    ! photoemission model.
    ! orig. Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after April 2024
    !===============================================================================
    use od_constants, only: dp, hbar, e_mass
    use od_electronic, only: nbands, nspins, num_electrons, electrons_per_state, foptical_mat, fem_energy_info, efermi
    use od_cell, only: num_kpoints_on_node, cell_get_symmetry, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_parameters, only: optics_geom, optics_qdir, legacy_file_format, devel_flag, photo_photon_sweep, jdos_spacing,&
     & photo_work_function
    use od_io, only: io_error, stdout
    use od_comms, only: my_node_id, on_root

    implicit none

    complex(kind=dp), dimension(3) :: g
    real(kind=dp), dimension(3) :: qdir, qdir1, qdir2
    real(kind=dp), dimension(2) :: num_occ
    real(kind=dp) :: q_weight1, q_weight2, factor, energy_max, tolerance = 0.0001_dp
    integer :: N_k, i, j, N_in, N_spin, N2, N3, n_eigen, num_symm, ierr, energy_index

    if (.not. legacy_file_format .and. index(devel_flag, 'old_filename') > 0) then
      num_symm = 0
      call cell_get_symmetry
    end if
    num_symm = num_crystal_symmetry_operations

    num_occ = 0.0_dp
    do N_spin = 1, nspins
      num_occ(N_spin) = num_electrons(N_spin)
    end do

    if (electrons_per_state == 2) then
      num_occ(1) = num_occ(1)/2.0_dp
    end if

    ! fem_energy_info: energy_count, energy_min, energy_step, energy_fermi, energy_workfct
    energy_count = int(fem_energy_info(1))
    energy_min = fem_energy_info(2)
    energy_step = fem_energy_info(3)
    energy_fermi = fem_energy_info(4)
    energy_workfct = fem_energy_info(5)
    energy_max = energy_min + energy_step*(energy_count - 1)
    ! Check the inputs are compatible
    ! Are the jdos_step and energy_step compatible?
    ! Check this specifically for photon_sweep, as that is quite important, otherwise check if the current energy can be
    ! reached using the input step
    if (jdos_spacing .lt. energy_step) then
      if (on_root) then
        write (stdout, *) 'jdos_spacing = ', jdos_spacing, '1step energy steps for OMEs:', energy_step
        write (stdout, *) 'The jdos_spacing is smaller than the supplied energy_step from the .fem_bin and thus incompatible!'
      end if
      call io_error('The jdos_spacing is smaller than the supplied energy_step from the .fem_bin and thus incompatible!')
    end if
    ! If energy_step is lt jdos_spacing - is the mod==0?
    if (energy_step .lt. jdos_spacing) then
      if (abs(modulo(jdos_spacing, energy_step)) .gt. tolerance) then
        if (on_root) then
          write (stdout, *) 'jdos_spacing = ', jdos_spacing, '1step energy steps for OMEs:', energy_step
          write (stdout, *) 'The jdos_spacing and energy_step for 1step OMEs are not a multiple of each other!'
        end if
        call io_error('The jdos_spacing and energy_step for 1step OMEs are not a multiple of each other!')
      end if
    end if
    ! Is the current photon_energy within the bounds of the energy_min and energy_max values?
    if (temp_photon_energy .gt. energy_max .or. temp_photon_energy .lt. energy_min) then
      if (on_root) then
        write (stdout, *) 'current E_photon = ', temp_photon_energy, 'energy bounds for 1step OMEs:' &
          , energy_min, '->', energy_max
        write (stdout, *) 'The current photon energy is out of the min->max range of the 1step OMEs!'
      end if
      call io_error('The current photon energy is out of the min->max range of the 1step OMEs!')
    end if
    ! Is the fermi_energy within error?
    if (abs(energy_fermi - efermi) .gt. tolerance) then
      if (on_root) then
        write (stdout, *) 'optados E_fermi:', efermi, '1step OME E_fermi:', energy_fermi
        write (stdout, *) 'The Fermi Energy calculated in OptaDOS and supplied from the .fem_bin are incompatible!'
      end if
      call io_error('The Fermi Energy calculated in OptaDOS and supplied from the .fem_bin are incompatible!')
    end if
    ! Is the energy_workfct within error?
    if (abs(energy_workfct - photo_work_function) .gt. tolerance) then
      if (on_root) then
        write (stdout, *) 'optados workfct:', photo_work_function, '1step OME workfct:', energy_workfct
        write (stdout, *) 'The Workfct from OptaDOS input and supplied from the .fem_bin are incompatible!'
      end if
      call io_error('The Workfct from OptaDOS input and supplied from the .fem_bin are incompatible!')
    end if

    ! Calculate the correct energy index in foptical_mat to use for the population of foptical_matrix_weights
    energy_index = nint(((temp_photon_energy - energy_min)/energy_step)) + 1
    if (on_root .and. enable_debug_output) write (stdout, *) 'energy_index:', energy_index

    if (.not. allocated(foptical_matrix_weights)) then
      allocate (foptical_matrix_weights(nbands, num_kpoints_on_node(my_node_id), nspins, N_geom), stat=ierr)
      if (ierr /= 0) call io_error('Error: make_foptical_weights - allocation of foptical_matrix_weights failed')
    end if
    foptical_matrix_weights = 0.0_dp

    if (index(optics_geom, 'polar') > 0) then
      qdir = optics_qdir
      q_weight = ((qdir(1)**2) + (qdir(2)**2) + (qdir(3)**2))**0.5_dp
      if (q_weight < 0.001_dp) call io_error("Error:  please check optics_qdir, norm close to zero")
    end if

    if (index(optics_geom, 'unpolar') > 0) then
      !TO CHANGE WHEN THE light_direction IS CORRECTED
      if (optics_qdir(3) .lt. 1E-06) then
        qdir1(1) = 0.0_dp
        qdir1(2) = 0.0_dp
        qdir1(3) = 1.0_dp
      else
        qdir1(1) = 1.0_dp
        qdir1(2) = 1.0_dp
        qdir1(3) = -(optics_qdir(1) + optics_qdir(2))/optics_qdir(3)
      end if
      qdir2(1) = (optics_qdir(2)*qdir1(3)) - (optics_qdir(3)*qdir1(2))
      qdir2(2) = (optics_qdir(3)*qdir1(1)) - (optics_qdir(1)*qdir1(3))
      qdir2(3) = (optics_qdir(1)*qdir1(2)) - (optics_qdir(2)*qdir1(1))
      q_weight1 = ((qdir1(1)**2) + (qdir1(2)**2) + (qdir1(3)**2))**0.5_dp
      q_weight2 = ((qdir2(1)**2) + (qdir2(2)**2) + (qdir2(3)**2))**0.5_dp
    end if

    N_in = 1  ! 0 = no inversion, 1 = inversion
    g = 0.0_dp

    do N_k = 1, num_kpoints_on_node(my_node_id)                 ! Loop over kpoints
      do N_spin = 1, nspins                                   ! Loop over spins
        do n_eigen = 1, nbands                                ! Loop over state
          factor = 1.0_dp/(temp_photon_energy**2)
          if (index(optics_geom, 'unpolar') > 0) then
            if (num_symm == 0) then
              g(1) = (((qdir1(1)*foptical_mat(n_eigen, 1, energy_index, N_k, N_spin)) + &
                       (qdir1(2)*foptical_mat(n_eigen, 2, energy_index, N_k, N_spin)) + &
                       (qdir1(3)*foptical_mat(n_eigen, 3, energy_index, N_k, N_spin)))/q_weight1)
              g(2) = (((qdir2(1)*foptical_mat(n_eigen, 1, energy_index, N_k, N_spin)) + &
                       (qdir2(2)*foptical_mat(n_eigen, 2, energy_index, N_k, N_spin)) + &
                       (qdir2(3)*foptical_mat(n_eigen, 3, energy_index, N_k, N_spin)))/q_weight2)
              foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) = &
                0.5_dp*factor*(real(g(1)*conjg(g(1)), dp) + real(g(2)*conjg(g(2)), dp))
            else ! begin unpolar symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  ! Calculating foptical_matrix_weights contribution for qdir1
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir1(j))
                    end do
                  end do
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, 1, energy_index, N_k, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, 2, energy_index, N_k, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, 3, energy_index, N_k, N_spin)))/q_weight1)
                  foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*real(g(1)*conjg(g(1)), dp)*factor
                  g(1) = 0.0_dp
                  ! Calculating foptical_matrix_weights contribution for qdir2
                  do i = 1, 3 ! if I include an extra variable I can merge this and the last do loops
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir2(j))
                    end do
                  end do
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, 1, energy_index, N_k, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, 2, energy_index, N_k, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, 3, energy_index, N_k, N_spin)))/q_weight2)
                  foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*real(g(1)*conjg(g(1)), dp)*factor
                end do
              end do
            end if !end unpolar symmetric
          elseif (index(optics_geom, 'polar') > 0) then
            if (num_symm == 0) then
              g(1) = (((qdir(1)*foptical_mat(n_eigen, nbands + 1, 1, N_k, N_spin)) + &
                       (qdir(2)*foptical_mat(n_eigen, nbands + 1, 2, N_k, N_spin)) + &
                       (qdir(3)*foptical_mat(n_eigen, nbands + 1, 3, N_k, N_spin)))/q_weight)
              foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) = factor*real(g(1)*conjg(g(1)), dp)
            else !begin polar symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))* &
                                (crystal_symmetry_operations(j, i, N2)*optics_qdir(j))
                    end do
                  end do
                  g(1) = 0.0_dp
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, 1, energy_index, N_k, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, 2, energy_index, N_k, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, 3, energy_index, N_k, N_spin)))/q_weight)
                  foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, N_k, N_spin, N_geom) + &
                    (1.0_dp/Real((num_symm*(N_in + 1)), dp))*factor*real(g(1)*conjg(g(1)), dp)
                end do
              end do
            end if ! end polar symmetric
          end if ! end photo_geom
        end do ! loop over state 1
      end do ! loop over spins
    end do ! loop over kpoints

    if (allocated(foptical_mat) .and. current_photo_energy_index .eq. number_energies) then
      deallocate (foptical_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error: make_foptical_weights - failed to deallocate foptical_mat')
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+------------------------- Printing Free OM Weights -------------------------+'
      write (stdout, 126) shape(foptical_matrix_weights)
      write (stdout, 126) nbands + 1, nbands + 1, num_kpoints_on_node(my_node_id), nspins, N_geom
126   format(5(1x, I4))
      do N2 = 1, N_geom
        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(es15.8))') (foptical_matrix_weights(n_eigen, N_k, N_spin, N2), n_eigen=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

  end subroutine make_foptical_weights

  !===============================================================================
  subroutine calc_one_step_model
    !===============================================================================
    ! This subroutine calculates the QE using a one step model.
    ! orig. Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after April 2024
    !===============================================================================

    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient,&
    & elec_read_band_curvature, photo_spectral_func
    use od_comms, only: my_node_id, num_nodes
    use od_parameters, only: scissor_op, photo_temperature, devel_flag, photo_photon_sweep, &
                             iprint, photo_model, photo_sf_max_vectors
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_comms, only: on_root, comms_recv, comms_send
    use od_io, only: stdout, io_error, io_file_unit, io_time, seedname, io_date
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB, inv_sqrt_two_pi
    implicit none
    integer :: N_k, N_spin, n_eigen, atom, ierr, i, gdx, kpt_total, inode, token, qe_unit

    real(kind=dp) :: width, norm_vac, vacuum_gauss, transverse_gauss, qe_factor, argument, time0, time1
    real(kind=dp) :: temp_contribution
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:) :: qe_k_temp
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    if (index(devel_flag, 'write_fem_matrix') > 0) then
      kpt_total = sum(num_kpoints_on_node(0:num_nodes - 1))
      call write_distributed_fem_data(kpt_total)
    end if

    qe_factor = 1.0_dp/(cell_area)
    width = (1.0_dp/11604.45_dp)*photo_temperature
    norm_vac = inv_sqrt_two_pi/width

    time0 = io_time()
    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 1Step QE ---------------------------+'
    end if

    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of field_emission failed')
    end if
    field_emission = 0.0_dp

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    if (.not. allocated(qe_osm)) then
      allocate (qe_osm(nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of qe_osm failed')
    end if
    qe_osm = 0.0_dp

    if (.not. allocated(te_osm)) then
      allocate (te_osm(nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of te_osm failed')
    end if
    te_osm = 0.0_dp

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a an exponent value of ~1E+/-100
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      i = 13 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
      write (stdout, '(1x,a78)') '+------------ Printing list of values going into 1step QE Values ------------+'
      write (stdout, '(10(1x,a17))') 'calced_qe_value', 'foptical_matrix_weights', 'electron_esc', 'kpoint_weight',&
      & 'I_layer', 'transverse_gauss', 'vacuum_gauss', 'fermi_dirac', 'pdos_weights_atoms', 'pdos_weights_k_band'
      write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', i, max_atoms, nbands, nspins, num_kpoints_on_node(my_node_id)
    end if
    do atom = 1, max_atoms + 1
      if (iprint > 2 .and. on_root .and. (atom .le. max_atoms)) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', "Calculating atom ", atom, " of", max_atoms, "<-- QE-1S |"
      end if
      do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            do gdx = 1, photo_sf_max_vectors
              if ((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)), &
                                            width, (evacuum_eff - efermi))/norm_vac
              else
                transverse_gauss = 1.0_dp
              end if
              if ((band_energy(n_eigen, N_spin, N_k) + temp_photon_energy) .lt. evacuum_eff) then
                vacuum_gauss = gaussian((band_energy(n_eigen, N_spin, N_k) + temp_photon_energy) + &
                                        scissor_op, width, evacuum_eff)/norm_vac
              else
                vacuum_gauss = 1.0_dp
              end if
              temp_contribution = (qe_factor*photo_spectral_func(3, gdx, n_eigen, N_spin, N_k) &
                                   *foptical_matrix_weights(n_eigen, N_k, N_spin, 1) &
                                   *(electron_esc(gdx, n_eigen, N_spin, N_k, atom)) &
                                   *electrons_per_state*kpoint_weight(N_k) &
                                   *(I_layer(layer(atom), current_photo_energy_index)) &
                                   *transverse_gauss*vacuum_gauss*fermi_dirac(n_eigen, N_spin, N_k) &
                                   *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                                     /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                                  *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
              qe_osm(n_eigen, N_spin, N_k, atom) = qe_osm(n_eigen, N_spin, N_k, atom) &
                                                   + temp_contribution
              te_osm(n_eigen, N_spin, N_k, atom) = te_osm(n_eigen, N_spin, N_k, atom) &
                                                   + temp_contribution*E_transverse(gdx, n_eigen, N_spin, N_k)
            end do
            if (enable_debug_output .and. index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
              write (stdout, '(4(1x,I4))') atom, n_eigen, N_spin, N_k
              write (stdout, '(10(7x,E17.9E3))') qe_osm(n_eigen, N_spin, N_k, atom), &
                foptical_matrix_weights(n_eigen, N_k, N_spin, 1), &
                electron_esc(gdx, n_eigen, N_spin, N_k, atom), kpoint_weight(N_k), &
                I_layer(layer(atom), current_photo_energy_index), transverse_gauss, vacuum_gauss, &
                fermi_dirac(n_eigen, N_spin, N_k), &
                pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)), pdos_weights_k_band(n_eigen, N_spin, N_k)
            end if
          end do
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if ((index(devel_flag, 'print_qe_constituents') > 0 .or. index(devel_flag, 'print_qe_matrix_full') > 0) .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------- Printing 1step QE Matrix -------------------------+'
      write (stdout, 125) shape(qe_osm)
      write (stdout, 125) nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1
125   format(4(1x, I4))
      do atom = 1, max_atoms + 1
        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(9999(ES16.8E3))') (qe_osm(n_eigen, N_spin, N_k, atom), n_eigen=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (enable_debug_output .and. index(devel_flag, 'print_kpt_qe_data') > 0) then
      if (on_root) then
        qe_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_k_point_QE.dat'
        write (stdout, *) 'opening file'
        open (unit=qe_unit, action='write', file=filename)
        write (qe_unit, *) '# The k point dependent QE values'
        call io_date(cdate, ctime)
        write (qe_unit, *) '## OptaDOS Photoemission: Printing QE K point Data on ', cdate, ' at ', ctime
      end if

      allocate (qe_k_temp(num_kpoints_on_node(0)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_one_step_model - failed to allocate qe_k_temp on root')
      token = -1

      ! allocate and sum the 3step qe matrix on non-root
      if (.not. on_root) then
        do N_k = 1, num_kpoints_on_node(my_node_id)
          qe_k_temp(N_k) = sum(qe_osm(:, :, N_k, :))
        end do
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that node
        call comms_send(qe_k_temp(1), num_kpoints_on_node(my_node_id), 0)
        ! - send token back to root node
        call comms_send(token, 1, 0)
      end if

      if (on_root) then
        do inode = 1, num_nodes - 1
          ! - send to the token to notes in turn
          call comms_send(token, 1, inode)
          ! - receive the qe_matrix from the other notes and write it to the file
          call comms_recv(qe_k_temp(1), num_kpoints_on_node(inode), inode)
          ! write out the qe_matrix to the file
          do N_k = 1, num_kpoints_on_node(inode)
            write (qe_unit, *) qe_k_temp(N_k)
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (qe_unit, *) sum(qe_osm(:, :, N_k, :))
        end do
        close (unit=qe_unit)
      end if
      deallocate (qe_k_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - failed to deallocate qe_k_temp')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 1step Photoemission', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_one_step_model

  !===============================================================================
  subroutine weighted_mean_te
    !*===============================================================================
    ! This subroutine calculates the weighted arithmetic mean transverse energy
    ! sum(QE*mte)/(total QE)
    ! orig. Victor Chang, 7 February 2020
    ! edited by Felix Mildner, after June 2023
    !===============================================================================
    use od_cell, only: cell_calc_kpoint_r_cart
    use od_electronic, only: elec_read_band_gradient, elec_read_band_curvature
    use od_comms, only: on_root, comms_reduce, comms_bcast
    use od_parameters, only: photo_model, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: io_error, io_file_unit, io_time, stdout
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: inv_sqrt_two_pi
    implicit none
    real(kind=dp), allocatable, dimension(:) :: layer_e_transverse
    real(kind=dp)                            :: time0, time1, qe_term1, qe_term2, mte_term1, mte_term2
    integer                                  :: atom, ierr

    time0 = io_time()

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Calculating MTE ------------------------------+'
    end if

    if (.not. allocated(layer_qe)) then
      allocate (layer_qe(max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of layer_qe failed')
    end if
    layer_qe = 0.0_dp

    if (.not. allocated(layer_e_transverse)) then
      allocate (layer_e_transverse(max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of layer_e_transverse failed')
    end if
    layer_e_transverse = 0.0_dp

    if (index(photo_model, '3step') > 0) then
      do atom = 1, max_atoms + 1
        ! Calculate the qe contribution of each atom/layer
        layer_qe(atom) = sum(qe_tsm(:, :, :, :, atom))
        layer_e_transverse(atom) = sum(te_tsm(:, :, :, atom))
      end do

      call comms_reduce(layer_e_transverse(1), max_atoms + 1, 'SUM')

      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(layer_qe(1), max_atoms + 1, 'SUM')
      ! Calculate the total QE
      total_qe = sum(layer_qe(1:(max_atoms + 1)))
      mean_te = sum(te_tsm(:, :, :, :))
      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(mean_te, 1, 'SUM')
      call comms_bcast(total_qe, 1)

      if (total_qe .gt. 0.0_dp) then
        mean_te = mean_te/total_qe
      else
        mean_te = 0.0_dp
      end if

      deallocate (te_tsm, stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - failed to deallocate te_tsm')

      deallocate (layer_e_transverse, stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - failed to deallocate layer_e_transverse')

    elseif (index(photo_model, '1step') > 0) then
      do atom = 1, max_atoms + 1
        ! Calculate the qe contribution of each atom/layer
        layer_qe(atom) = sum(qe_osm(:, :, :, atom))
      end do

      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(layer_qe(1), max_atoms + 1, 'SUM')
      ! Calculate the total QE
      total_qe = sum(layer_qe)
      call comms_bcast(total_qe, 1)

      ! Calculate the sum of transverse E from all the bands and k-points on node
      mean_te = sum(te_osm)
      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(mean_te, 1, 'SUM')

      if (total_qe .gt. 0.0_dp) then
        mean_te = mean_te/total_qe
      else
        mean_te = 0.0_dp
      end if

      deallocate (te_osm, stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - failed to deallocate te_osm')

    else if (index(photo_model, 'ds_like_pe') > 0) then

      qe_term1 = sum(qe_tsm(:, :, :, :, 2))
      call comms_reduce(qe_term1, 1, 'SUM')
      qe_term2 = sum(qe_tsm(:, :, :, :, 1))
      call comms_reduce(qe_term2, 1, 'SUM')
      mte_term1 = sum(qe_tsm(:, :, :, :, 3))
      call comms_reduce(mte_term1, 1, 'SUM')
      mte_term2 = sum(qe_tsm(:, :, :, :, 2))
      call comms_reduce(mte_term2, 1, 'SUM')
      total_qe = qe_term1/qe_term2
      mean_te = 0.5*(mte_term1/mte_term2)

    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a23,36x,f11.3,a8)') '+ Time to calculate MTE', time1 - time0, ' (sec) +'
    end if

  end subroutine weighted_mean_te

  subroutine write_qe_data
    !*===============================================================================
    ! This subroutine writes the calculated Photoemission data to the output file.
    ! The contents of this routine used to be part of the subroutine weighted_mean_te, but were moved
    ! here to make the subroutine names more representative of their function.
    !===============================================================================
    use od_cell, only: cell_calc_kpoint_r_cart, atoms_label_tmp
    use od_comms, only: on_root
    use od_parameters, only: photo_work_function, photo_elec_field, photo_model
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout
    use od_jdos_utils, only: jdos_utils_calculate
    integer :: atom

    if (on_root) then
      write (stdout, '(1x,a78)') '+------------------------------ Photoemission -------------------------------+'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, 223) '| Work Function     ', photo_work_function, &
        'eV         Photon Energy', temp_photon_energy, 'eV   |'
      write (stdout, 224) '| Effective Work Function', work_function_eff, &
        'eV         Electric Field', photo_elec_field, 'V/A  |'

      if (index(photo_model, '3step') > 0) then
        write (stdout, '(1x,a78)') '| Final State : Bloch State                                                  |'
      elseif (index(photo_model, '1step') > 0) then
        write (stdout, '(1x,a78)') '| Final State : Free Electron State                                          |'
      end if
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Layer   |             Quantum Efficiency           |'
      if (index(photo_model, 'ds_like_pe') > 0) then
        write (stdout, '(1x,a78)') '|       ********** Calculated a simplified DS like PE model **********       |'
        write (stdout, 227) '| Total Quantum Efficiency (electrons/photon):', total_qe, '   |'

        write (stdout, 228) '| Weighted Mean Transverse Energy (eV):', mean_te, '      |'
      else
        do atom = 1, max_atoms
          write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
            layer(atom), layer_qe(atom), "      |"
        end do
        write (stdout, 226) "| Bulk", layer_qe(max_atoms + 1), &
        &"      |"

        write (stdout, 227) '| Total Quantum Efficiency (electrons/photon):', total_qe, '   |'

        write (stdout, 228) '| Weighted Mean Transverse Energy (eV):', mean_te, '      |'
      end if

      if (photo_elec_field .gt. 0.0_dp) then
        write (stdout, 229) '| Total field emission (electrons/A^2):', total_field_emission, '      |'
      end if

      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
223 format(1x, a20, f15.4, 1x, a24, f11.4, a7)
224 format(1x, a25, f10.4, 1x, a25, f10.4, a7)
225 format(1x, a1, a4, 8x, I3, 10x, I3, 16x, E24.16E3, 2x, a7)
226 format(1x, a6, 38x, E25.16E3, 2x, a7)
227 format(1x, a46, E25.16E3, a7)
228 format(1x, a39, 7x, E25.16E3, a7)
229 format(1x, a39, 7x, E25.16E3, a7)
  end subroutine write_qe_data

  subroutine binding_energy_broadening
    !===============================================================================
    !* This subroutine applies a Gaussian broadening to the binding energy
    ! orig. Victor Chang, 7 February 2020
    ! edited Felix Mildner, after August 2024
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, photo_spectral_func, transmit_prob
    use od_parameters, only: photo_work_function, photo_model, photo_theta_min, photo_theta_max, photo_temperature, &
    & photo_phi_min, photo_phi_max, photo_bindenergy_broadening, photo_sf_max_vectors, scissor_op, iprint, num_exclude_bands, &
    & exclude_bands
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :, :) :: binding_temp
    real(kind=dp) :: time0, time1

    real(kind=dp) :: final_fd, initial_fd
    integer :: N_k, N_spin, n_eigen, n_eigen_final, atom, e_scale, gdx, ierr
    integer :: middle_idx, width_idx

    real(kind=dp) :: temp_contribution, norm_vac, qe_factor, width, argument, transverse_gauss
    real(kind=dp), allocatable, dimension(:, :, :, :) :: fermi_dirac

    time0 = io_time()
    qe_factor = 1.0_dp/(cell_area)
    width = (1.0_dp/11604.45_dp)*photo_temperature
    norm_vac = inv_sqrt_two_pi/width

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(2, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    max_energy = int((temp_photon_energy - photo_work_function)*1000) + 100

    if (max_energy .lt. 0) return

    allocate (t_energy(max_energy), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_broadening - allocation of t_energy failed')
    t_energy = 0.0_dp

    allocate (weighted_temp(max_energy, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_broadening - allocation of weighted_temp failed')
    weighted_temp = 0.0_dp

    allocate (binding_temp(max_energy, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_broadening - allocation of binding_temp failed')
    binding_temp = 0.0_dp

    do e_scale = 1, max_energy
      t_energy(e_scale) = real(e_scale - 1, dp)/1000
    end do

    do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          middle_idx = ceiling((efermi - band_energy(n_eigen, N_spin, N_k))/0.001)
          width_idx = ceiling((photo_bindenergy_broadening*10)/0.001)
          ! do e_scale = max(middle_idx - width_idx, 1), min(middle_idx + width_idx, max_energy)
          do e_scale = 1, max_energy
            binding_temp(e_scale, n_eigen, N_spin, N_k) = &
              gaussian((efermi - band_energy(n_eigen, N_spin, N_k)), photo_bindenergy_broadening, t_energy(e_scale))
          end do
        end do
      end do
    end do
    if (index(photo_model, '3step') > 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
            ! This is a bit of an arbitrary condition, but it turns out
            ! that this corresponds to a an exponent value of ~1E+/-100
            ! and this cutoff condition saves us from running into arithmetic
            ! issues when computing fermi_dirac due to possible underflow.
            if (argument .gt. 230.0_dp) then
              fermi_dirac(1, n_eigen, N_spin, N_k) = 0.0_dp
            elseif (argument .lt. -230.0_dp) then
              fermi_dirac(1, n_eigen, N_spin, N_k) = 1.0_dp
            else
              fermi_dirac(1, n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
            end if
            ! normally called vacuum_gauss
            if ((band_energy(n_eigen, N_spin, N_k)) .lt. evacuum_eff) then
              fermi_dirac(2, n_eigen, N_spin, N_k) = gaussian((band_energy(n_eigen, N_spin, N_k)) + &
                                                              scissor_op, width, evacuum_eff)/norm_vac
            else
              fermi_dirac(2, n_eigen, N_spin, N_k) = 1.0_dp
            end if
          end do
        end do
      end do

      call photo_calculate_delta(delta_temp, .false.)

      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen_final = 2, nbands
              if (num_exclude_bands .gt. 1) then
                if (any(exclude_bands == n_eigen_final)) then
                  cycle
                end if
              end if
              final_fd = 1 - fermi_dirac(1, n_eigen_final, N_spin, N_k)
              do n_eigen = 1, n_eigen_final - 1
                initial_fd = fermi_dirac(1, n_eigen, N_spin, N_k)
                middle_idx = ceiling((efermi - band_energy(n_eigen, N_spin, N_k))/0.001)
                width_idx = ceiling((photo_bindenergy_broadening*10)/0.001)
                do gdx = 1, photo_sf_max_vectors
                  if (theta_arpes(gdx, n_eigen_final, N_spin, N_k) .ge. photo_theta_min .and. &
                      theta_arpes(gdx, n_eigen_final, N_spin, N_k) .le. photo_theta_max) then
                    if (phi_arpes(gdx, n_eigen_final, N_spin, N_k) .ge. photo_phi_min .and. &
                        phi_arpes(gdx, n_eigen_final, N_spin, N_k) .le. photo_phi_max) then
                      ! evacuum_eff = efermi + photo_work_function
                      ! Is (photon_energy - transverse energy) > (work_function - E_field_lowering)
                      ! Is the final kinetic energy > 0?
                      if ((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                        transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)), &
                                                    width, (evacuum_eff - efermi))/norm_vac
                      else
                        transverse_gauss = 1.0_dp
                      end if
                      temp_contribution = &
                        qe_factor*photo_spectral_func(3, gdx, n_eigen, N_spin, N_k) &
                        *matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, 1) &
                        *delta_temp(n_eigen, n_eigen_final, N_spin, N_k) &
                        *electron_esc(gdx, n_eigen_final, N_spin, N_k, atom) &
                        *transmit_prob(n_eigen_final, N_spin, N_k) &
                        *electrons_per_state*kpoint_weight(N_k) &
                        *(I_layer(box_atom(atom), current_photo_energy_index)) &
                        !                             vacuum_gauss
                        *transverse_gauss*fermi_dirac(2, n_eigen_final, N_spin, N_k) &
                        *initial_fd*final_fd &
                        *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                          /pdos_weights_k_band(n_eigen, N_spin, N_k)) &
                        *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
                      ! do e_scale = max(middle_idx - width_idx, 1), min(middle_idx + width_idx, max_energy)
                      do e_scale = 1, max_energy
                        weighted_temp(e_scale, n_eigen, N_spin, N_k, atom) = weighted_temp(e_scale, n_eigen, N_spin, N_k, atom) + &
                                                                             binding_temp(e_scale, n_eigen, N_spin, N_k) &
                                                                             *temp_contribution
                      end do
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

      call photo_calculate_delta(delta_temp, .true.)

      do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen_final = 2, nbands
            if (num_exclude_bands .gt. 1) then
              if (any(exclude_bands == n_eigen_final)) then
                cycle
              end if
            end if
            final_fd = 1 - fermi_dirac(1, n_eigen_final, N_spin, N_k)
            do n_eigen = 1, n_eigen_final - 1
              middle_idx = ceiling((efermi - band_energy(n_eigen, N_spin, N_k))/0.001)
              width_idx = ceiling((photo_bindenergy_broadening*10)/0.001)
              initial_fd = fermi_dirac(1, n_eigen, N_spin, N_k)
              do gdx = 1, photo_sf_max_vectors
                if (theta_arpes(gdx, n_eigen_final, N_spin, N_k) .ge. photo_theta_min .and. &
                    theta_arpes(gdx, n_eigen_final, N_spin, N_k) .le. photo_theta_max) then
                  if (phi_arpes(gdx, n_eigen_final, N_spin, N_k) .ge. photo_phi_min .and. &
                      phi_arpes(gdx, n_eigen_final, N_spin, N_k) .le. photo_phi_max) then
                    ! evacuum_eff = efermi + photo_work_function
                    ! Is (photon_energy - transverse energy) > (work_function - E_field_lowering)
                    ! Is the final kinetic energy > 0?
                    if ((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                      transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)), &
                                                  width, (evacuum_eff - efermi))/norm_vac
                    else
                      transverse_gauss = 1.0_dp
                    end if

                    temp_contribution = &
                      (qe_factor*photo_spectral_func(3, gdx, n_eigen, N_spin, N_k) &
                       *matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, 1) &
                       *delta_temp(n_eigen, n_eigen_final, N_spin, N_k) &
                       *electron_esc(gdx, n_eigen_final, N_spin, N_k, max_atoms + 1) &
                       *transmit_prob(n_eigen_final, N_spin, N_k) &
                       *electrons_per_state*kpoint_weight(N_k) &
                       !                             vacuum_gauss
                       *transverse_gauss*fermi_dirac(2, n_eigen_final, N_spin, N_k) &
                       *initial_fd*final_fd &
                       *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(max_atoms)) &
                         /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                      *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
                    ! do e_scale = max(middle_idx - width_idx, 1), min(middle_idx + width_idx, max_energy)
                    do e_scale = 1, max_energy
                      weighted_temp(e_scale, n_eigen, N_spin, N_k, max_atoms + 1) = &
                        weighted_temp(e_scale, n_eigen, N_spin, N_k, max_atoms + 1) + &
                        binding_temp(e_scale, n_eigen, N_spin, N_k)*temp_contribution
                    end do
                  end if
                end if
              end do
            end do
          end do
        end do
      end do

      if (allocated(transmit_prob)) then
        deallocate (transmit_prob, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate transmit_prob')
      end if

    elseif (index(photo_model, '1step') > 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
            ! This is a bit of an arbitrary condition, but it turns out
            ! that this corresponds to a an exponent value of ~1E+/-100
            ! and this cutoff condition saves us from running into arithmetic
            ! issues when computing fermi_dirac due to possible underflow.
            if (argument .gt. 230.0_dp) then
              fermi_dirac(1, n_eigen, N_spin, N_k) = 0.0_dp
            elseif (argument .lt. -230.0_dp) then
              fermi_dirac(1, n_eigen, N_spin, N_k) = 1.0_dp
            else
              fermi_dirac(1, n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
            end if
            ! normally called vacuum_gauss
            if ((band_energy(n_eigen, N_spin, N_k) + temp_photon_energy) .lt. evacuum_eff) then
              fermi_dirac(2, n_eigen, N_spin, N_k) = gaussian((band_energy(n_eigen, N_spin, N_k) + temp_photon_energy) + &
                                                              scissor_op, width, evacuum_eff)/norm_vac
            else
              fermi_dirac(2, n_eigen, N_spin, N_k) = 1.0_dp
            end if
          end do
        end do
      end do

      do atom = 1, max_atoms + 1
        do N_k = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              middle_idx = ceiling((efermi - band_energy(n_eigen, N_spin, N_k))/0.001)
              width_idx = ceiling((photo_bindenergy_broadening*10)/0.001)
              do gdx = 1, photo_sf_max_vectors
                if (theta_arpes(gdx, n_eigen, N_spin, N_k) .ge. photo_theta_min .and. &
                    theta_arpes(gdx, n_eigen, N_spin, N_k) .le. photo_theta_max) then
                  if (phi_arpes(gdx, n_eigen, N_spin, N_k) .ge. photo_phi_min .and. &
                      phi_arpes(gdx, n_eigen, N_spin, N_k) .le. photo_phi_max) then
                    if ((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)) .le. (evacuum_eff - efermi)) then
                      transverse_gauss = gaussian((temp_photon_energy - E_transverse(gdx, n_eigen, N_spin, N_k)), &
                                                  width, (evacuum_eff - efermi))/norm_vac
                    else
                      transverse_gauss = 1.0_dp
                    end if
                    temp_contribution = (qe_factor*photo_spectral_func(3, gdx, n_eigen, N_spin, N_k) &
                                         *foptical_matrix_weights(n_eigen, N_k, N_spin, 1) &
                                         *electron_esc(gdx, n_eigen, N_spin, N_k, atom) &
                                         *electrons_per_state*kpoint_weight(N_k) &
                                         *I_layer(layer(atom), current_photo_energy_index) &
                                         !                                  vacuum_gauss
                                         *transverse_gauss*fermi_dirac(2, n_eigen, N_spin, N_k) &
                                         *fermi_dirac(1, n_eigen, N_spin, N_k) &
                                         *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                                           /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                                        *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
                    do e_scale = max(middle_idx - width_idx, 1), min(middle_idx + width_idx, max_energy)
                      weighted_temp(e_scale, n_eigen, N_spin, N_k, atom) = weighted_temp(e_scale, n_eigen, N_spin, N_k, atom) + &
                                                                           binding_temp(e_scale, n_eigen, N_spin, N_k) &
                                                                           *temp_contribution
                    end do
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    deallocate (binding_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_broadening - failed to deallocate binding_temp')

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a46,13x,f11.3,a8)') '+ Time to calculate broadened binding energies', time1 - time0, ' (sec) +'
    end if
  end subroutine binding_energy_broadening

  subroutine write_qe_output_files
    !*===============================================================================
    ! This subroutine writes either the transverse energy or the binding energy
    ! after the Gaussian broadening has been applied.
    ! orig. Victor Chang, 7 February 2020
    ! edited Felix Mildner, after April 2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nbands, nspins
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_reduce, comms_bcast
    use od_io, only: io_error, seedname, io_file_unit, io_date, io_time, stdout
    use od_parameters, only: photo_output, photo_model, photo_work_function, iprint, devel_flag, &
                             photo_theta_min, photo_theta_max, photo_phi_min, photo_phi_max
    implicit none
    integer :: atom, ierr, e_scale, binding_unit, matrix_unit
    integer :: N_k, N_spin, n_eigen, kpt_total, band_num

    real(kind=dp), allocatable, dimension(:, :) :: qe_atom
    real(kind=dp) :: time0, time1, total_weighted, qe_norm
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (index(photo_output, 'qe_matrix') > 0) then
      call cell_calc_kpoint_r_cart
      kpt_total = sum(num_kpoints_on_node(0:num_nodes - 1))
      if (num_nodes .gt. 1) then
        call write_distributed_qe_data(kpt_total)
      else
        matrix_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
          filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix_final.dat'
        else
          filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix.dat'
        end if
        open (unit=matrix_unit, action='write', file=filename)
        call io_date(cdate, ctime)
        write (matrix_unit, *) '## OptaDOS Photoemission: Printing QE Matrix on ', cdate, ' at ', ctime
        write (matrix_unit, *) '## Seedname: ', trim(seedname)
        write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
        write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
        if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
          write (matrix_unit, *) '## Writing the contributions of excitations into the !!FINAL!! states'
        end if
        write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'

        ! Printing out the info on root_node

        if (index(photo_model, '3step') > 0) then
          if (index(devel_flag, 'single') > 0) then
            n_eigen = len_trim(devel_flag)
            read (devel_flag(n_eigen - 2:n_eigen), *) band_num
            write (matrix_unit, '(1x,a42,1x,I3)') '## Writing contributions into final band #', band_num
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                      & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N_k = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') &
                    (qe_tsm(n_eigen, band_num, N_spin, N_k, atom), n_eigen=1, nbands)
                end do
              end do
            end do
          else if (index(devel_flag, 'final') > 0) then
            write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
            write (matrix_unit, *) '## at a certain k-point, spin, and atom'
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                      & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N_k = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') &
                    (sum(qe_tsm(1:nbands, n_eigen, N_spin, N_k, atom)), n_eigen=1, nbands)
                end do
              end do
            end do
          else
            write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
            write (matrix_unit, *) '## at a certain k-point, spin, and atom'
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                      & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N_k = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') &
                    (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N_k, atom)), n_eigen=1, nbands)
                end do
              end do
            end do
          end if
        elseif (index(photo_model, '1step') > 0) then
          write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
          write (matrix_unit, *) '## at a certain k-point, spin, and atom'
          write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                    & ')'
          do atom = 1, max_atoms + 1
            if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') (qe_osm(n_eigen, N_spin, N_k, atom), n_eigen=1, nbands)
              end do
            end do
          end do
        end if
        close (unit=matrix_unit)
      end if
    end if

    if (index(photo_output, 'e_bind') > 0 .and. max_energy .gt. 0) then

      allocate (qe_atom(max_atoms + 1, max_energy), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - allocation of qe_atom failed')
      qe_atom = 0.0_dp
      do e_scale = 1, max_energy !loop over binding energy
        do atom = 1, max_atoms + 1
          qe_atom(atom, e_scale) = &
            sum(weighted_temp(e_scale, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), atom))
        end do
      end do

      call comms_reduce(qe_atom(1, 1), max_energy*(max_atoms + 1), "SUM")

      total_weighted = sum(qe_atom(:, :))
      call comms_reduce(total_weighted, 1, "SUM")
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_qe/total_weighted
      else
        qe_norm = 1.0_dp
      end if

      qe_atom = qe_atom*qe_norm

      if (on_root) then
        binding_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))// &
                   '_binding_energy.dat'
        open (unit=binding_unit, action='write', file=filename)
        call io_date(cdate, ctime)
        write (binding_unit, '(1x,a60,a9,a4,a11)') '## OptaDOS Photoemission: Printing Broadened Binding Energy on ',&
        & cdate, ' at ', ctime
        write (binding_unit, '(1x,a13,a80)') '## Seedname: ', trim(seedname)
        write (binding_unit, '(1x,a24,a12)') '## Photoemission Model: ', trim(photo_model)
        write (binding_unit, '(1x,a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
        write (binding_unit, '(1x,a64,2(1x,f7.2))') '## Emission angle theta min, max (w.r.t. surface normal) [deg]: ', &
          photo_theta_min, photo_theta_max
        write (binding_unit, '(1x,a54,2(1x,f7.2))') '## Emission angle phi min, max (w.r.t. x-axis) [deg]: ', &
          photo_phi_min, photo_phi_max
        write (binding_unit, '(1x,a34,f9.5)') '## Fermi Energy Ekin offset [eV]: ', (temp_photon_energy - photo_work_function)
        write (binding_unit, '(1x,a66,1x,a50)') '## Binding Energy (EB) [eV] | Total QE from sum(atoms + bulk) @ EB',&
        &'| Contributions from: atom1 | atom2 | ... | bulk |'

        do e_scale = 1, max_energy
          write (binding_unit, '(1x,ES13.6E2,2x,ES25.12E3,1x,999(1x,ES25.12E3))') t_energy(e_scale), &
            sum(qe_atom(1:max_atoms + 1, e_scale)), qe_atom(1:max_atoms + 1, e_scale)
        end do

        close (unit=binding_unit)
      end if
    end if

    if (allocated(weighted_temp)) then
      deallocate (weighted_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate weighted_temp')
    end if

    if (allocated(qe_atom)) then
      deallocate (qe_atom, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate qe_atom')
    end if

    if (allocated(t_energy)) then
      deallocate (t_energy, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate t_energy')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a32,27x,f11.3,a8)') '+ Time to write the output files', time1 - time0, ' (sec) +'
    end if

  end subroutine write_qe_output_files

  subroutine write_distributed_qe_data(kpt_total)
    !* This subroutine writes the distributed qe tensor to a single file.
    ! To save on required memory the output file is accessed by each MPI process in turn
    ! and writes its values/contents one after the other.
    ! F. Mildner, June 2023
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nspins, nbands
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_bcast
    use od_io, only: io_error, io_file_unit, io_date, io_time, seedname
    use od_parameters, only: photo_model, devel_flag

    implicit none
    real(kind=dp), dimension(:, :, :), allocatable :: qe_mat_temp
    real(kind=dp), dimension(:, :, :, :), allocatable :: tsm_reduced
    real(kind=dp), dimension(:, :, :, :), allocatable :: osm_reduced
    integer, intent(in)                         :: kpt_total
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    integer:: N_k, N_spin, n_eigen, atom, token, matrix_unit, ierr, inode

    ! On root open file and write header

    if (on_root) then
      ! Writing header to output file
      write (char_e, '(F7.3)') temp_photon_energy
      if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix_final.dat'
      else
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix.dat'
      end if
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)
      write (matrix_unit, *) '## OptaDOS Photoemission: Printing QE Matrix on ', cdate, ' at ', ctime
      write (matrix_unit, *) '## Seedname: ', trim(seedname)
      write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
      write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
      write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'
      write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
      write (matrix_unit, *) '## at a certain k-point, spin, and atom'
      write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, &
        max_atoms + 1, ')'
      allocate (qe_mat_temp(nbands, nspins, num_kpoints_on_node(0)), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate qe_mat_temp on root')
      token = -1
    end if

    ! allocate and sum the 3step qe matrix on non-root
    if (.not. on_root) then
      if (index(photo_model, '3step') > 0) then
        allocate (tsm_reduced(nbands, nspins, num_kpoints_on_node(0), max_atoms + 1), stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate tsm_reduced')
        if (index(devel_flag, 'final') > 0) then
          tsm_reduced = sum(qe_tsm, dim=1)
        else
          tsm_reduced = sum(qe_tsm, dim=2)
        end if
      else if (index(photo_model, '1step') > 0) then
        allocate (osm_reduced(nbands, nspins, num_kpoints_on_node(0), max_atoms + 1), stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate tsm_reduced')
        osm_reduced = qe_osm
      end if
    end if
    ! For each atom until max_atoms+1
    do atom = 1, max_atoms + 1
      ! On non root nodes
      if (.not. on_root) then
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that specific atom
        if (index(photo_model, '3step') > 0) then
          call comms_send(tsm_reduced(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
        elseif (index(photo_model, '1step') > 0) then
          call comms_send(osm_reduced(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
        end if
        ! - send token back to root node
        call comms_send(token, 1, 0)
        ! On root node
      elseif (on_root) then
        do inode = 1, num_nodes - 1
          ! - send to the token to notes in turn
          call comms_send(token, 1, inode)
          ! - receive the qe_matrix from the other notes and write it to the file
          call comms_recv(qe_mat_temp(1, 1, 1), nbands*nspins*num_kpoints_on_node(inode), inode)
          ! write out the qe_matrix to the file
          do N_k = 1, num_kpoints_on_node(inode)
            do N_spin = 1, nspins
              write (matrix_unit, '(9999(ES16.8E3))') (qe_mat_temp(n_eigen, N_spin, N_k), n_eigen=1, nbands)
            end do
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        if (index(photo_model, '3step') > 0) then
          if (index(devel_flag, 'final') > 0) then
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') &
                  (sum(qe_tsm(1:nbands, n_eigen, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          else
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') &
                  (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          end if
        elseif (index(photo_model, '1step') > 0) then
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (matrix_unit, '(9999(ES16.8E3))') &
                (qe_osm(n_eigen, N_spin, N_k, atom), n_eigen=1, nbands)
            end do
          end do
        end if
        ! Write header for bulk contrib using root node
        if (atom .eq. max_atoms) write (matrix_unit, '(1x,a21)') '## Bulk Contribution:'
      end if
    end do
    if (on_root) then
      close (unit=matrix_unit)
      deallocate (qe_mat_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to deallocate qe_mat_temp')
    elseif (.not. on_root) then
      if (index(photo_model, '3step') > 0) then
        deallocate (tsm_reduced, stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to deallocate tsm_reduced')
      end if
    end if
  end subroutine write_distributed_qe_data

  subroutine write_distributed_fem_data(kpt_total)
    !***************************************************************
    ! This subroutine writes the distributed free electron matrix element tensor to a single file.
    ! To save on required memory the output file is accessed by each MPI process in turn
    ! and writes its values/contents one after the other.
    ! F. Mildner, June 2023

    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nspins, nbands
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_bcast
    use od_io, only: io_error, io_file_unit, io_date, io_time, seedname
    use od_parameters, only: photo_model

    implicit none
    real(kind=dp), dimension(:, :, :), allocatable :: fem_mat_temp
    ! real(kind=dp), dimension(:, :, :, :), allocatable :: tsm_reduced
    integer, intent(in)                         :: kpt_total
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    integer:: N_k, N_spin, n_eigen, token, matrix_unit, ierr, inode

    ! On root open file and write header

    if (on_root) then
      ! Writing header to output file
      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_fem_matrix.dat'
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)
      write (matrix_unit, *) '## OptaDOS Photoemission: Printing OME Matrix on ', cdate, ' at ', ctime
      write (matrix_unit, *) '## Seedname: ', trim(seedname)
      write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
      write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
      write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'
      write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
      write (matrix_unit, *) '## at a certain k-point, spin, and atom'
      write (matrix_unit, '(1x,a31,3(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, kpt_total, nspins, ')'
      allocate (fem_mat_temp(nbands, num_kpoints_on_node(0), nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate fem_mat_temp on root')
      token = -1
    end if

    ! On non root nodes
    if (.not. on_root) then
      ! - wait for the token
      call comms_recv(token, 1, 0)
      ! - send the respective qe_matrix for that specific atom
      call comms_send(foptical_matrix_weights(1, 1, 1, 1), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
      ! - send token back to root node
      call comms_send(token, 1, 0)
      ! On root node
    elseif (on_root) then
      do inode = 1, num_nodes - 1
        ! - send to the token to notes in turn
        call comms_send(token, 1, inode)
        ! - receive the qe_matrix from the other notes and write it to the file
        call comms_recv(fem_mat_temp(1, 1, 1), nbands*nspins*num_kpoints_on_node(inode), inode)
        ! write out the qe_matrix to the file
        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(inode)
            write (matrix_unit, '(9999(ES16.8E3))') (fem_mat_temp(n_eigen, N_k, N_spin), n_eigen=1, nbands)
          end do
        end do
        ! - receive the token from a node
        call comms_recv(token, 1, inode)
      end do
      ! - write root qe_matrix elements
      do N_spin = 1, nspins
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (matrix_unit, '(9999(ES16.8E3))') (foptical_matrix_weights(n_eigen, N_k, N_spin, 1), n_eigen=1, nbands)
        end do
      end do
      close (unit=matrix_unit)
      deallocate (fem_mat_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to deallocate fem_mat_temp')
    end if
  end subroutine write_distributed_fem_data

  subroutine photo_deallocate
    !***************************************************************
    ! This subroutine deallocates all the quantities which have not
    ! been deallocated yet

    use od_io, only: io_error
    use od_electronic, only: foptical_mat
    implicit none
    integer :: ierr

    if (allocated(phi_arpes)) then
      deallocate (phi_arpes, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate phi_arpes')
    end if

    if (allocated(theta_arpes)) then
      deallocate (theta_arpes, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate theta_arpes')
    end if

    if (allocated(theta_arpes_internal)) then
      deallocate (theta_arpes_internal, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate theta_arpes_internal')
    end if

    if (allocated(refract)) then
      deallocate (refract, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate refract')
    end if

    if (allocated(absorp)) then
      deallocate (absorp, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate absorp')
    end if

    if (allocated(electron_esc)) then
      deallocate (electron_esc, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate electron_esc')
    end if

    if (allocated(layer)) then
      deallocate (layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate layer')
    end if

    if (allocated(layer_qe)) then
      deallocate (layer_qe, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate layer_qe')
    end if

    if (allocated(imfp_val)) then
      deallocate (imfp_val, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate imfp_val')
    end if

    if (allocated(reflect)) then
      deallocate (reflect, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate reflect')
    end if

    if (allocated(matrix_weights)) then
      deallocate (matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate matrix_weights')
    end if

    if (allocated(E_transverse)) then
      deallocate (E_transverse, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate E_transverse')
    end if

    if (allocated(absorp_photo)) then
      deallocate (absorp_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate absorp_photo')
    end if

    if (allocated(atom_order)) then
      deallocate (atom_order, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate atom_order')
    end if

    if (allocated(atoms_per_layer)) then
      deallocate (atoms_per_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate atoms_per_layer')
    end if

    if (allocated(thickness_atom)) then
      deallocate (thickness_atom, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate thickness_atom')
    end if

    if (allocated(pdos_weights_atoms)) then
      deallocate (pdos_weights_atoms, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_atoms')
    end if

    if (allocated(pdos_weights_k_band)) then
      deallocate (pdos_weights_k_band, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_k_band')
    end if

    if (allocated(index_energy)) then
      deallocate (index_energy, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate index_energy')
    end if

    if (allocated(I_layer)) then
      deallocate (I_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate I_layer')
    end if

    if (allocated(thickness_layer)) then
      deallocate (thickness_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate thickness_layer')
    end if

    if (allocated(E_kinetic)) then
      deallocate (E_kinetic, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate E_kinetic')
    end if

    if (allocated(field_emission)) then
      deallocate (field_emission, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate field_emission')
    end if

    ! if (allocated(bulk_prob)) then
    !   deallocate (bulk_prob, stat=ierr)
    !   if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate bulk_prob')
    ! end if

    if (allocated(qe_tsm)) then
      deallocate (qe_tsm, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate qe_tsm')
    end if

    if (allocated(qe_osm)) then
      deallocate (qe_osm, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate qe_osm')
    end if

    if (allocated(foptical_mat)) then
      deallocate (foptical_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate foptical_mat')
    end if

    if (allocated(foptical_matrix_weights)) then
      deallocate (foptical_matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate foptical_matrix_weights')
    end if

  end subroutine photo_deallocate

end module od_photo
