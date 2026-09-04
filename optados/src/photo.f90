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
! along with this program.  If not, see .lt.http://www.gnu.org/licenses/.gt..
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
  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: photo_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: projected_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :) :: foptical_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :) :: weighted_jdos
  real(kind=dp), allocatable, public, dimension(:, :) :: absorp_layer
  real(kind=dp), allocatable, public, dimension(:, :, :) :: pdos_weights_k_band
  real(kind=dp), allocatable, public, save :: E(:)
  real(kind=dp), allocatable, dimension(:, :, :) :: imfp_val
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: electron_esc
  real(kind=dp), dimension(:, :), allocatable :: I_layer
  real(kind=dp), allocatable, dimension(:, :) :: reflect_photo
  ! Per-box dielectric function, (box, energy, 1:2) for eps_1 and eps_2.
  ! Kept because the slab reflectivity has to be formed by averaging epsilon
  ! and not by averaging anything computed from it.
  real(kind=dp), allocatable, dimension(:, :, :) :: epsilon_photo
  real(kind=dp), allocatable, dimension(:, :) :: absorp_photo
  real(kind=dp), allocatable, dimension(:, :) :: refract
  real(kind=dp), allocatable, dimension(:)  :: reflect
  real(kind=dp), allocatable, dimension(:) :: absorp
  real(kind=dp), allocatable, dimension(:) :: box_heights
  real(kind=dp), allocatable, dimension(:) :: box_volumes
  integer, dimension(:), allocatable       :: box_atom
  integer, dimension(:), allocatable       :: atoms_per_box
  integer                                  :: num_boxes
  real(kind=dp)                            :: slab_middle_ref
  !! Centroid spacing from one explicit layer to the next: the optical path the
  !! light travels between them, which is not a box height.
  real(kind=dp), dimension(:), allocatable :: layer_gap
  !! Repeat spacing of the substrate, which bulk_emission steps by. Taken from
  !! the layers at and below the explicit region, not from the deepest box.
  real(kind=dp)                            :: bulk_repeat
  !! Confinement length of the slab along z: the top of the slab down to its
  !! middle. Stands in for the k_z sub-cell length, which a slab does not have.
  !! Set once in analyse_geometry and used by both calculate_delta and the
  !! photoemission JDOS, so the two cannot drift apart.
  real(kind=dp)                            :: slab_half_height
  real(kind=dp)                            :: cell_area
  real(kind=dp), dimension(:), allocatable :: atom_imfp
  real(kind=dp), dimension(:, :, :), allocatable :: band_imfp
  real(kind=dp), dimension(:), allocatable :: boxes_top_z_coord
  logical                                  :: single_layer
  real(kind=dp), dimension(:, :), allocatable :: new_atom_coordinates
  !! Fraction of the symmetry images of the transverse momentum whose azimuth
  !! falls inside the acceptance wedge. Used by the outputs that have no azimuth
  !! axis and therefore no symmetry loop of their own; the map and tensor
  !! routines test each image directly with phi_accepted instead.
  real(kind=dp), allocatable, dimension(:, :, :, :) :: phi_accept_frac
  real(kind=dp), allocatable, dimension(:, :, :, :) :: theta_arpes
  real(kind=dp), allocatable, dimension(:, :, :, :) :: theta_internal
  real(kind=dp), allocatable, dimension(:, :, :, :) :: E_kinetic
  real(kind=dp), allocatable, dimension(:, :, :, :) :: E_transverse
  real(kind=dp), allocatable, dimension(:) :: bind_energy
  real(kind=dp), allocatable, dimension(:, :) :: weighted_be_atom
  real(kind=dp)                               :: total_be_contribs
  real(kind=dp)                               :: total_be_kmat_contribs
  real(kind=dp), allocatable, dimension(:, :) :: ekin_k_matrix
  real(kind=dp), allocatable, dimension(:, :) :: kxky_matrix
  real(kind=dp), allocatable, dimension(:, :, :) :: p_tensor
  integer, dimension(3) :: max_bin_p
  integer :: max_energy = -1
  real(kind=dp), allocatable, dimension(:, :, :, :)    :: qe_osm
  real(kind=dp), allocatable, dimension(:, :, :, :)    :: te_osm
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: qe_tsm
  ! Dowell-Schmerge like model. Kept apart from qe_tsm, whose last index is
  ! the box a transition belongs to: the DS model was reusing that index to
  ! mean "which of five quantities", so the same module array had two
  ! unrelated shapes, and slots 4 and 5 spent an nbands x nbands x nspins x nk
  ! array on two scalars.
  real(kind=dp), allocatable, dimension(:, :, :, :) :: ds_qe_den, ds_qe_num, ds_mte_num
  real(kind=dp) :: ds_dos_mte_num, ds_dos_mte_den
  real(kind=dp), allocatable, dimension(:, :, :, :) :: te_tsm
  real(kind=dp), allocatable, dimension(:, :, :, :) :: gkgrid_weight
  integer :: photo_gkmax
  real(kind=dp) :: mean_te
  real(kind=dp) :: total_qe
  real(kind=dp), allocatable, dimension(:) :: layer_qe
  integer, dimension(:), allocatable :: atom_order
  real(kind=dp) :: work_function_eff
  !! The step an escaping electron climbs, used for the refraction at the
  !! surface. photo_inner_potential when it is given, otherwise the work
  !! function, which is what the code always used and is too small by roughly a
  !! factor of three for a metal.
  real(kind=dp) :: surface_barrier
  logical       :: barrier_warned = .false.
  real(kind=dp) :: evacuum
  real(kind=dp) :: evacuum_eff
  real(kind=dp) :: total_field_emission
  real(kind=dp), allocatable, dimension(:, :, :) :: field_emission
  integer :: N_geom
  integer :: max_atoms
  integer :: max_bin_e, max_bin_k
  real(kind=dp) :: max_e_kinetic, max_k_transverse, plot_extra_upper = 1.0_dp
  ! Margin added to the transverse momentum axis of the E_kin vs p maps, in 1/A.
  real(kind=dp) :: k_extra_padding = 0.2_dp
  real(kind=dp) :: q_weight
  ! Added by Felix Mildner, 12/2022 and later
  integer, allocatable, dimension(:)  :: index_energy
  integer                             :: number_energies, current_energy_index, current_photo_energy_index
  real(kind=dp)                       :: temp_photon_energy, time_a, time_b
  integer, allocatable, dimension(:, :):: min_index_unocc
  ! The coherency tensor is calculated in Castep over a window of final-state
  ! energies. That window is described by fem_energy_info, read from the file by
  ! od_electronic, in the order
  !   n_Ef, Ef_min, Ef_step, Ef_broadening, Ef_origin
  ! all in eV. The lookup is referenced to the vacuum level rather than to a
  ! Fermi energy or work function baked into the file, so nothing here has to
  ! agree with a value Castep chose - only the window has to reach far enough
  ! for the photon energies being swept.
contains

  subroutine photo_calculate
    !! Main subroutine calling all the other subroutine steps.
    use od_electronic, only: elec_dealloc_optical, elec_pdos_read, elec_read_optical_mat, &
      efermi, efermi_set, elec_read_foptical_mat, elec_dealloc_pdos
    use od_jdos_utils, only: jdos_utils_calculate, setup_energy_scale
    use od_comms, only: on_root, comms_bcast
    use od_parameters, only: photo_work_function, photo_model, photo_elec_field, photo_output, photo_energy_sweep, &
      photo_photon_min, jdos_spacing, photo_photon_energy, photo_momentum, iprint
    use od_dos_utils, only: dos_utils_set_efermi, dos_utils_calculate_at_e, dos_utils_deallocate
    use od_io, only: stdout, io_error, io_time
    use od_pdos, only: pdos_calculate

    implicit none

    integer :: idx, token

    if (on_root) then
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                             Photoemission Calculation                      +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    if (.not. efermi_set) then
      call dos_utils_set_efermi
      call dos_utils_deallocate
    end if

    ! Identify layers
    call analyse_geometry
    call calc_band_info
    call calc_photon_energies

    if (index(photo_model, 'dosds') .eq. 0) then
      call elec_read_optical_mat
      call elec_pdos_read
      call make_pdos_weights_atoms
      call elec_dealloc_pdos

      ! Calculate the optical properties of the slab
      call calc_photo_optics
      call calc_absorp_layer
    end if

    ! Electric field and field emission
    if (photo_elec_field .gt. 1.0E-30_dp) then
      call effective_wf
    else
      evacuum_eff = efermi + photo_work_function
      work_function_eff = photo_work_function
    end if

    do idx = 1, number_energies
      time_a = io_time()
      if (photo_energy_sweep) then
        temp_photon_energy = photo_photon_min + (idx - 1)*jdos_spacing
      else
        temp_photon_energy = photo_photon_energy
      end if
      if (on_root) write (stdout, '(1x,a47,f8.4,a23)') '+------------------ Starting Photoemission with', temp_photon_energy, &
        ' eV ------------------+'
      current_photo_energy_index = idx
      current_energy_index = index_energy(idx)

      if (photo_elec_field .gt. 1.0E-30_dp) call calc_field_emission

      ! Three-step-model
      if (index(photo_model, '3step') .gt. 0) then
        ! Calculate the photoemission angles theta/phi and transverse energy
        call calc_angle
        ! Calculate the electron escape length
        call calc_electron_esc
        call bulk_emission
        ! Calculate QE
        call calc_three_step_model

        ! One-step-model
      elseif (index(photo_model, '1step') .gt. 0) then
        ! Calculate the photoemission angles theta/phi and transverse energy
        call calc_angle
        ! Calculate the electron escape length
        call calc_electron_esc
        call bulk_emission
        ! Read the one-step matrix elements
        if (.not. allocated(foptical_matrix_weights)) call elec_read_foptical_mat
        ! Calculate the one-step optical matrix
        call make_foptical_weights
        ! Calculate QE
        call calc_one_step_model

        ! Simplified DS like model
      elseif (index(photo_model, 'dosds') .gt. 0) then
        call calc_ds_like_model
      end if

      !Weight the contribution of each electron
      !to the transverse energy spread according to their QE
      call weighted_mean_te
      if (on_root) then
        call write_qe_results
        token = 1
      end if
      ! quick comms synchronisation
      call comms_bcast(token, 1)

      ! Only call the binding energy gaussian broadening and file printing if necessary
      if (index(photo_output, 'off') == 0) then
        !Broaden ouputs using a gaussian function
        if (index(photo_output, 'bindenergy_curve') .gt. 0) call binding_energy_curve
        if (index(photo_output, 'ekin_ptrans_map') .gt. 0) then
          if (index(photo_momentum, 'gkgrid') .gt. 0) then
            call kinetic_energy_momentum_map_gkgrid
          else
            call kinetic_energy_momentum_map
          end if
        end if
        if (index(photo_output, 'const_bindenergy_p_map') .gt. 0) then
          if (index(photo_momentum, 'gkgrid') .gt. 0) then
            call const_binding_energy_map_gkgrid
          else
            call const_binding_energy_map
          end if
        end if
        if (index(photo_output, 'p_tensor') .gt. 0) then
          if (index(photo_momentum, 'gkgrid') .gt. 0) then
            call full_momentum_tensor_gkgrid
          else
            call full_momentum_tensor
          end if
        end if
        if (index(photo_output, 'qe_tensor') .gt. 0) call write_qe_tensor
      end if
      time_b = io_time()
      if (on_root .and. iprint .gt. 1) then
        write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate Photoemission (step)', time_b - time_a, ' (sec) +'
      end if
    end do
    ! Deallocate the rest that was needed for the photoemission calcs
    call photo_deallocate

    if (on_root) write (stdout, '(1x,a78)') '| End of Photoemission Calculation                                           |'
  end subroutine photo_calculate

  subroutine analyse_geometry
    !* This subroutine identifies and defines a set of boxes,
    ! that represent layers, with a height = interlayer distance
    ! at the middle of the slab. All atoms are then sorted into
    ! these boxes for later use.
    use od_constants, only: dp, periodic_table_name, deg_to_rad
    use od_cell, only: num_atoms, atoms_pos_cart_photo, atoms_label_tmp, cell_volume, real_lattice
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    use od_parameters, only: photo_imfp_value, photo_slab_max, photo_slab_min, photo_slab_middle, photo_layers_tops, iprint, &
      photo_slab_mode, SLAB_MODE_LAYERS, photo_slab_middle_set
    implicit none
    integer :: ierr, atom, counter, i, j, ic, atom_index, first, temp, atom_1, atom_2
    integer :: n_layers, n_species_seen, isp
    integer, allocatable, dimension(:)       :: layer_index
    real(kind=dp)                            :: z_top, z_middle
    real(kind=dp), allocatable, dimension(:) :: box_boundaries
    integer, allocatable, dimension(:, :)    :: layer_species_count
    character(len=10), allocatable, dimension(:) :: species_seen
    logical                                  :: same_species_set, same_species_counts
    logical                                  :: symmetric_stack
    character(len=80)                        :: comp_str, temp_str
    real(kind=dp)                            :: h_min, h_max, h_mean
    real(kind=dp)                            :: diff_temp, current_top, diff_top = 10000.0_dp, diff_bottom = 10000.0_dp
    real(kind=dp)                            :: max_gap, typical_gap, layer_tol
    real(kind=dp)                            :: z_middle_tol, wrap_gap
    character(len=78)                        :: box_msg
    real(kind=dp), allocatable, dimension(:) :: z_gaps, large_gaps, layer_centroid
    integer, allocatable, dimension(:)       :: atoms_in_layer
    logical                                  :: layers_from_input
    ! Below this the atoms are all at the same height and no interlayer
    ! spacing can be inferred from the structure (e.g. a graphene monolayer).
    real(kind=dp), parameter                 :: min_layer_gap = 0.1_dp
    allocate (atom_order(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of atom_order failed')

    do i = 1, num_atoms
      atom_order(i) = i
    end do

    allocate (box_atom(num_atoms + 1), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of box_atom failed')
    box_atom = 1000

    ! real_lattice is stored column-wise, so (3,1) and (3,2) are the z-components
    ! of a and b. Both must vanish for the slab normal to be cartesian z, which
    ! a lot of what follows assumes, and which makes cell_area = cell_volume /
    ! real_lattice(3,3) exact. Either one being non-zero breaks that, and a tilt
    ! in the negative direction breaks it just as badly as a positive one, so the
    ! test has to be .or. on the magnitudes.
    if (abs(real_lattice(3, 1)) .gt. 1.0E-5_dp .or. abs(real_lattice(3, 2)) .gt. 1.0E-5_dp) then
      call io_error('Error: analyse_geometry - the a and b lattice vectors are not in the cart. xy plane &
      &(the slab normal is not along z) - not currently implemented!')
    end if

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

    ! --------------------------------------------------------------------------------------------
    ! *    The following code was added in Nov 2023 to test out a new layer assignment scheme    *
    ! *    A set of boxes with the height of the central slab layer distance is created and the  *
    ! *    atoms are sorted into those boxes by their z-coordinate.                              *
    ! --------------------------------------------------------------------------------------------
    ! determine the cell area for later use
    cell_area = cell_volume/real_lattice(3, 3)

    ! User has given the slab center and top coordinates for each of
    ! the layers. We infer that the upper surface of the slab is the
    ! surface layer's top coordinate (in parameter.f90).
    layers_from_input = (photo_slab_mode .eq. SLAB_MODE_LAYERS)

    ! ------------------------------------------------------------------------
    ! One geometry, however it was supplied.
    !
    ! The slab runs from z_top down to z_middle and is cut by num_boxes+1 planes
    !
    !     z_top = b(0) > b(1) > ... > b(num_boxes) = z_middle
    !
    ! so box i spans [b(i), b(i-1)], the boxes tile the slab exactly and their
    ! volumes add up to it. Everything below z_middle is covered by the bulk
    ! extrapolation. photo_layers_tops, when given, supplies the interior planes;
    ! otherwise they are the midpoints between neighbouring layer centroids,
    ! which tiles for any spacing and reduces to the old boxes when the spacings
    ! are equal.
    !
    ! The atoms are clustered into layers either way, because three different
    ! lengths are needed and only one of them is a box height:
    !
    !   box_heights  the extent of a box, which normalises its epsilon
    !   layer_gap    the centroid spacing, the optical path from the atoms of
    !                one layer to those of the next
    !   bulk_repeat  the substrate interlayer spacing, which bulk_emission steps
    !                by below the explicit region
    !
    ! Conflating those three is what let a relaxed surface, an unequal spacing or
    ! a user-chosen slab middle leak into quantities that have nothing to do with
    ! them.
    ! ------------------------------------------------------------------------

    ! ---- cluster the atoms into layers -------------------------------------
    allocate (z_gaps(max(num_atoms - 1, 1)), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of z_gaps failed')
    z_gaps = 0.0_dp
    do atom = 1, num_atoms - 1
      z_gaps(atom) = atoms_pos_cart_photo(3, atom_order(atom)) - &
                     atoms_pos_cart_photo(3, atom_order(atom + 1))
    end do

    max_gap = 0.0_dp
    if (num_atoms .gt. 1) max_gap = maxval(z_gaps(1:num_atoms - 1))

    ! Does the slab straddle the cell boundary?
    !
    ! Under 3D periodic boundaries a slab is a connected set of atoms separated
    ! from its own image by vacuum, and the machinery below assumes that vacuum
    ! is the gap that wraps round the cell -- in other words that the slab sits
    ! in one piece somewhere inside it. Shift it so that half sits at the bottom
    ! and half at the top and the vacuum becomes an *interior* gap in the sorted
    ! list, at which point it is taken for an interlayer spacing: on a 5-layer
    ! Cu(100) slab that gave one box instead of three, two atoms assigned to no
    ! box at all, a box volume five times too large, and a bulk repeat of
    ! 25.5 Ang for a crystal whose spacing is 1.78 -- with no error and a
    ! plausible looking quantum efficiency.
    !
    ! The test is exact: compare the largest interior gap against the one that
    ! wraps round. Whichever is bigger is the vacuum. Verified over 24 rigid
    ! shifts of an 8-layer slab through a cell, where it flags precisely the
    ! positions at which the layer clustering breaks and no others.
    !
    ! Refusing is the conservative half of the fix. The general remedy is to
    ! re-origin the atoms just above the larger gap, which restores a contiguous
    ! slab for any position; that also has to carry photo_slab_min, _max,
    ! _middle and photo_layers_tops through the same shift, which is a decision
    ! about what frame those keywords are written in, so it is left for now.
    if (num_atoms .gt. 1) then
      wrap_gap = (atoms_pos_cart_photo(3, atom_order(num_atoms)) + real_lattice(3, 3)) &
                 - atoms_pos_cart_photo(3, atom_order(1))
      if (max_gap .gt. wrap_gap) then
        if (on_root) then
          write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
          write (stdout, '(1x,a78)') '! Error: the structure is split across the cell boundary in z. The largest   !'
          write (stdout, '(1x,a78)') '! gap between atoms lies inside the sorted list rather than wrapping round   !'
          write (stdout, '(1x,a78)') '! the cell, so part of the slab sits at the bottom and part at the top.      !'
          write (stdout, '(1x,a78)') '! The vacuum would then be read as an interlayer spacing, giving too few     !'
          write (stdout, '(1x,a78)') '! layers, atoms belonging to no box, and a bulk repeat the size of the       !'
          write (stdout, '(1x,a78)') '! vacuum -- silently, with a plausible number at the end of it.              !'
          write (stdout, '(1x,a78)') '! Translate the structure so the slab sits in one piece inside the cell.     !'
          write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
          write (stdout, '(1x,a46,1x,f10.4,20x,a1)') '|  Largest gap between atoms  (Ang)          :', max_gap, '|'
          write (stdout, '(1x,a46,1x,f10.4,20x,a1)') '|  Gap wrapping round the cell (Ang)         :', wrap_gap, '|'
        end if
        call io_error('Error: analyse_geometry - the slab is split across the cell boundary in z. '// &
                      'Translate it so it sits in one piece.')
      end if
    end if

    if (num_atoms .eq. 1 .or. max_gap .lt. min_layer_gap) then
      ! All atoms at essentially the same height, so the thickness cannot be
      ! inferred and has to come from the user's slab bounds.
      single_layer = .true.
      layer_tol = huge(1.0_dp)
    else
      single_layer = .false.
      ! Typical interlayer spacing = median of the "large" gaps. Gaps below 10%
      ! of the largest are intra-layer rumpling and are excluded, so a layer
      ! holding many atoms cannot drag the median towards zero.
      allocate (large_gaps(num_atoms - 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of large_gaps failed')
      counter = 0
      do atom = 1, num_atoms - 1
        if (z_gaps(atom) .gt. 0.1_dp*max_gap) then
          counter = counter + 1
          large_gaps(counter) = z_gaps(atom)
        end if
      end do
      do i = 2, counter
        diff_temp = large_gaps(i)
        j = i - 1
        do while (j .ge. 1)
          if (large_gaps(j) .le. diff_temp) exit
          large_gaps(j + 1) = large_gaps(j)
          j = j - 1
        end do
        large_gaps(j + 1) = diff_temp
      end do
      typical_gap = large_gaps((counter + 1)/2)
      layer_tol = 0.5_dp*typical_gap
      deallocate (large_gaps, stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - deallocation of large_gaps failed')
    end if

    allocate (layer_index(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of layer_index failed')
    n_layers = 1
    layer_index(1) = 1
    do atom = 2, num_atoms
      if (z_gaps(atom - 1) .gt. layer_tol) n_layers = n_layers + 1
      layer_index(atom) = n_layers
    end do

    allocate (layer_centroid(n_layers), atoms_in_layer(n_layers), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of layer_centroid failed')
    layer_centroid = 0.0_dp
    atoms_in_layer = 0
    do atom = 1, num_atoms
      i = layer_index(atom)
      layer_centroid(i) = layer_centroid(i) + atoms_pos_cart_photo(3, atom_order(atom))
      atoms_in_layer(i) = atoms_in_layer(i) + 1
    end do
    do i = 1, n_layers
      layer_centroid(i) = layer_centroid(i)/real(atoms_in_layer(i), dp)
    end do

    ! ---- the three planes that define the slab ------------------------------
    z_top = photo_slab_max
    if (single_layer) then
      ! Nothing beneath a single layer to extrapolate, so the whole slab is the
      ! one box and bulk_emission returns immediately.
      z_middle = photo_slab_min
    else if (layers_from_input .or. photo_slab_middle_set) then
      z_middle = photo_slab_middle
    else
      ! Half the slab, stated geometrically rather than as a layer count. Close
      ! to the old (n_layers+1)/2 for an evenly spaced slab, and no longer
      ! dependent on how many layers happen to be present.
      z_middle = 0.5_dp*(photo_slab_max + photo_slab_min)
    end if
    if (z_middle .ge. z_top) &
      call io_error('Error: analyse_geometry - the bottom of the explicit region is not below '// &
                    'its top. Check photo_slab_middle against photo_slab_max.')
    if (z_top .lt. layer_centroid(1)) &
      call io_error('Error: analyse_geometry - the slab surface lies below the topmost layer. '// &
                    'Check photo_slab_max, or photo_layers_tops(1), against the structure.')

    ! ---- the box boundaries -------------------------------------------------
    if (layers_from_input) then
      num_boxes = size(photo_layers_tops, 1)
    else if (single_layer) then
      num_boxes = 1
    else
      ! A layer lying *on* the dividing plane is treated explicitly. The old
      ! layer-count form of this, num_boxes = (n_layers + 1)/2, could not get
      ! this wrong because it never compared a coordinate; stating the split
      ! geometrically instead made the middle layer of a slab with an odd number
      ! of layers a coin toss, because that layer sits exactly on the default
      ! z_middle = 0.5*(photo_slab_max + photo_slab_min). The same 13 layer slab
      ! gave 7 explicit boxes at one vacuum spacing and 6 at another, decided by
      ! rounding in the fourth decimal of the bounds on a margin of 1E-7 Ang.
      !
      ! The tolerance is a hundredth of the typical layer spacing, so it can only
      ! ever catch a layer that is coincident with the plane: the nearest
      ! genuinely distinct layer is half a spacing away, fifty times further out.
      ! For an evenly spaced slab this reproduces (n_layers + 1)/2 exactly, while
      ! an unevenly spaced one still splits on geometry rather than on how many
      ! layers the inference happened to find.
      if (photo_slab_middle_set) then
        z_middle_tol = 0.01_dp*typical_gap
        num_boxes = 0
        do i = 1, n_layers
          if (layer_centroid(i) .gt. z_middle - z_middle_tol) num_boxes = num_boxes + 1
        end do
      else
        ! Half the inferred layers, counted rather than measured. Stating the
        ! split geometrically instead -- counting the centroids above
        ! 0.5*(photo_slab_max + photo_slab_min) -- makes it inherit the
        ! arbitrariness of both bounds, which are the user's choice of where the
        ! electron density has fallen off rather than a property of the
        ! structure. Moving photo_slab_max by 0.5 Ang moves that plane by
        ! 0.25 Ang and reassigns a whole layer between the explicit region and
        ! the bulk extrapolation. The layer count cannot do that, and for an odd
        ! number of layers the +1 puts the middle layer on the explicit side
        ! instead of leaving it exactly on the plane to be decided by rounding.
        num_boxes = (n_layers + 1)/2
      end if
      if (num_boxes .gt. 0 .and. on_root .and. iprint .gt. 1) then
        if (layer_centroid(num_boxes) .le. z_middle) then
          write (box_msg, '(a1,5x,a,i0,a)') '|', 'Layer ', num_boxes, &
            ' sits on the explicit/bulk plane; kept explicit'
          box_msg(78:78) = '|'
          write (stdout, '(1x,a78)') box_msg
        end if
      end if
      if (num_boxes .lt. 1) &
        call io_error('Error: analyse_geometry - the bottom of the explicit region lies above '// &
                      'every layer, so no layer would be treated explicitly.')
      if (num_boxes .ge. n_layers) &
        call io_error('Error: analyse_geometry - the bottom of the explicit region lies below '// &
                      'every layer, so nothing is left for the bulk extrapolation to stand '// &
                      'for. It must leave at least one layer beneath it.')
      ! Snap the last boundary to the same midpoint rule as the others. Without
      ! this the deepest box is truncated wherever photo_slab_middle happens to
      ! fall, so its height -- and therefore its volume, its epsilon and its
      ! optical constants -- depends on the exact value rather than on the
      ! material. photo_slab_middle then does what it says: it selects which
      ! layer is the last explicit one.
      z_middle = 0.5_dp*(layer_centroid(num_boxes) + layer_centroid(num_boxes + 1))
    end if

    allocate (box_boundaries(0:num_boxes), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of box_boundaries failed')
    box_boundaries(0) = z_top
    box_boundaries(num_boxes) = z_middle
    if (layers_from_input) then
      do i = 1, num_boxes - 1
        box_boundaries(i) = photo_layers_tops(i + 1)
      end do
    else
      do i = 1, num_boxes - 1
        box_boundaries(i) = 0.5_dp*(layer_centroid(i) + layer_centroid(i + 1))
      end do
    end if

    if (.not. allocated(box_heights)) then
      allocate (box_heights(num_boxes), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of box_heights failed')
    end if
    if (.not. allocated(box_volumes)) then
      allocate (box_volumes(num_boxes), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of box_volumes failed')
    end if
    if (.not. allocated(boxes_top_z_coord)) then
      allocate (boxes_top_z_coord(num_boxes), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of boxes_top_z_coord failed')
    end if
    if (.not. allocated(atoms_per_box)) then
      allocate (atoms_per_box(num_boxes), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of atoms_per_box failed')
    end if

    do i = 1, num_boxes
      box_heights(i) = box_boundaries(i - 1) - box_boundaries(i)
      boxes_top_z_coord(i) = box_boundaries(i - 1)
      if (box_heights(i) .le. 0.0_dp) &
        call io_error('Error: analyse_geometry - the box boundaries do not decrease. Check '// &
                      'photo_layers_tops, and that photo_slab_middle lies below the last of them.')
    end do
    box_volumes = box_heights*cell_area
    slab_middle_ref = z_middle
    slab_half_height = z_top - z_middle

    ! ---- the two lengths that are not box heights ---------------------------
    allocate (layer_gap(max(num_boxes - 1, 1)), stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of layer_gap failed')
    layer_gap = 0.0_dp
    do i = 1, min(num_boxes - 1, n_layers - 1)
      layer_gap(i) = layer_centroid(i) - layer_centroid(i + 1)
    end do

    if (single_layer .or. n_layers .le. num_boxes) then
      bulk_repeat = box_heights(num_boxes)
    else
      ! Median of the spacings at and below the explicit region: that is the
      ! substrate's repeat, and it is what the bulk extrapolation walks in.
      ! z_gaps is reused as scratch here, its clustering job being done.
      counter = 0
      do i = max(num_boxes, 1), n_layers - 1
        counter = counter + 1
        z_gaps(counter) = layer_centroid(i) - layer_centroid(i + 1)
      end do
      do i = 2, counter
        diff_temp = z_gaps(i)
        j = i - 1
        do while (j .ge. 1)
          if (z_gaps(j) .le. diff_temp) exit
          z_gaps(j + 1) = z_gaps(j)
          j = j - 1
        end do
        z_gaps(j + 1) = diff_temp
      end do
      bulk_repeat = z_gaps((counter + 1)/2)
    end if

    ! ---- put the atoms in their boxes ---------------------------------------
    ! Scanning from the top means an atom sitting exactly on a shared plane goes
    ! to the upper box rather than being claimed twice or dropped. Anything below
    ! z_middle keeps num_boxes+1 and belongs to the bulk.
    box_atom = num_boxes + 1
    atoms_per_box = 0
    do atom = 1, num_atoms
      current_top = atoms_pos_cart_photo(3, atom_order(atom))
      if (current_top .gt. z_top) then
        if (on_root) write (stdout, '(1x,a,i0,a,f12.7,a,f12.7)') &
          'Error: atom ', atom_order(atom), ' at z = ', current_top, &
          ' lies above the slab surface at z = ', z_top
        call io_error('Error: analyse_geometry - an atom lies above the slab surface. Check '// &
                      'photo_slab_max, or photo_layers_tops(1), against the structure.')
      end if
      do i = 1, num_boxes
        if (current_top .le. box_boundaries(i - 1) .and. current_top .ge. box_boundaries(i)) then
          box_atom(atom) = i
          atoms_per_box(i) = atoms_per_box(i) + 1
          exit
        end if
      end do
    end do

    if (.not. layers_from_input) then
      ! ----------------------------------------------------------------------
      ! Can the inferred-layer model describe this structure at all?
      !
      ! Everything downstream assumes the layers are equivalent repeats of one
      ! material: num_boxes keeps the top half of them whatever they contain,
      ! box_volumes gives each the same kind of extent, one epsilon per box sets
      ! the optics of that box - the topmost one alone fixes the reflectivity of
      ! the whole surface - and bulk_emission repeats the deepest explicit box
      ! downwards as though everything below it were more of the same.
      !
      ! None of that survives an interface, an adsorbate, or a compound whose
      ! planes alternate in species. An upright CO on Cu, for instance, clusters
      ! into a one-atom O layer above a one-atom C layer, and the oxygen's
      ! dielectric function then becomes the reflectivity of the copper surface.
      ! Say so and stop, rather than return a confident wrong answer.
      ! ----------------------------------------------------------------------
      allocate (species_seen(num_atoms), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of species_seen failed')
      n_species_seen = 0
      do atom = 1, num_atoms
        isp = 0
        do i = 1, n_species_seen
          if (trim(species_seen(i)) .eq. trim(atoms_label_tmp(atom_order(atom)))) isp = i
        end do
        if (isp .eq. 0) then
          n_species_seen = n_species_seen + 1
          species_seen(n_species_seen) = atoms_label_tmp(atom_order(atom))
        end if
      end do

      allocate (layer_species_count(n_species_seen, n_layers), stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - allocation of layer_species_count failed')
      layer_species_count = 0
      do atom = 1, num_atoms
        do i = 1, n_species_seen
          if (trim(species_seen(i)) .eq. trim(atoms_label_tmp(atom_order(atom)))) then
            layer_species_count(i, layer_index(atom)) = layer_species_count(i, layer_index(atom)) + 1
          end if
        end do
      end do

      ! Which species are present is the fatal test; how many of each there are
      ! is only a warning, since a vacancy or a reconstruction changes the counts
      ! without making the layers different materials.
      same_species_set = .true.
      same_species_counts = .true.
      do i = 2, n_layers
        do isp = 1, n_species_seen
          if ((layer_species_count(isp, i) .gt. 0) .neqv. (layer_species_count(isp, 1) .gt. 0)) &
            same_species_set = .false.
          if (layer_species_count(isp, i) .ne. layer_species_count(isp, 1)) same_species_counts = .false.
        end do
      end do

      ! A symmetric sandwich is describable even though its layers are not all
      ! the same material. If the composition sequence reads the same from either
      ! face -- film, substrate, film -- then the middle layer of the stack is
      ! the middle of the substrate, and num_boxes = (n_layers + 1)/2 puts the
      ! deepest explicit box exactly there. bulk_emission then replicates the
      ! centre of the sandwich downwards, which is the same approximation the
      ! code already makes for a homogeneous symmetric slab.
      !
      ! That is the whole content of the test: symmetry is precisely the
      ! condition under which halving the stack lands on the right material. A
      ! film grown on one face only fails it, and there the deepest explicit box
      ! is film whenever the film is half the stack or more, so the bulk term
      ! stands in for a substrate it has never seen.
      !
      ! The test is on composition, not position, so relaxation of the two faces
      ! does not break it -- only a genuine asymmetry does, such as an adsorbate
      ! or a vacancy on one side.
      symmetric_stack = .false.
      if (.not. same_species_set .and. n_layers .gt. 2) then
        symmetric_stack = .true.
        do i = 1, n_layers/2
          do isp = 1, n_species_seen
            if (layer_species_count(isp, i) .ne. layer_species_count(isp, n_layers + 1 - i)) &
              symmetric_stack = .false.
          end do
        end do
      end if

      if (symmetric_stack .and. on_root .and. iprint .gt. 1) then
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        write (stdout, '(1x,a78)') '| The layers are not all one material, but the stack reads the same from     |'
        write (stdout, '(1x,a78)') '| either face, so it is treated as a symmetric sandwich. The explicit        |'
        write (stdout, '(1x,a78)') '| region is the top half and the bulk term stands for the middle layer,      |'
        write (stdout, '(1x,a78)') '| which is the centre of the slab. Set photo_layers_tops and                 |'
        write (stdout, '(1x,a78)') '| photo_slab_middle instead if that is not what is wanted.                   |'
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if

      if (.not. same_species_set .and. .not. symmetric_stack) then
        if (on_root) then
          write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
          write (stdout, '(1x,a78)') '| The inferred layers do not all contain the same species, so this           |'
          write (stdout, '(1x,a78)') '| structure is an interface, an adsorbate system, or a compound with         |'
          write (stdout, '(1x,a78)') '| alternating planes. The layer model cannot describe it: it takes the       |'
          write (stdout, '(1x,a78)') '| top half of the stack as the explicit region whatever is in it, gives      |'
          write (stdout, '(1x,a78)') '| the topmost layer the optical constants of the whole surface, and          |'
          write (stdout, '(1x,a78)') '| extrapolates the deepest explicit layer downwards as though the            |'
          write (stdout, '(1x,a78)') '| substrate were made of it.                                                 |'
          write (stdout, '(1x,a78)') '|                                                                            |'
          write (stdout, '(1x,a78)') '| Layer | Composition, top layer first                                       |'
          do i = 1, n_layers
            comp_str = ' '
            do isp = 1, n_species_seen
              if (layer_species_count(isp, i) .gt. 0) then
                write (temp_str, '(1x,a,i0)') trim(species_seen(isp)), layer_species_count(isp, i)
                comp_str = trim(comp_str)//trim(temp_str)
              end if
            end do
            write (stdout, '(1x,a1,i6,1x,a1,1x,a67,a1)') '|', i, '|', adjustl(comp_str), '|'
          end do
          write (stdout, '(1x,a78)') '|                                                                            |'
          write (stdout, '(1x,a78)') '| Set the layers yourself with photo_layers_tops and photo_slab_middle.      |'
          write (stdout, '(1x,a78)') '| Those tile the slab exactly and let you group a full repeat unit of        |'
          write (stdout, '(1x,a78)') '| material into one layer rather than one plane of one species.              |'
          write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        end if
        call io_error('Error: analyse_geometry - the inferred layers are not all the same material. '// &
                      'Use photo_layers_tops and photo_slab_middle to define them explicitly.')
      end if

      if (.not. same_species_counts .and. on_root) then
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        write (stdout, '(1x,a78)') '! Warning: the inferred layers hold the same species but not the same        !'
        write (stdout, '(1x,a78)') '! numbers of them, as a vacancy or a reconstruction would give. Each layer   !'
        write (stdout, '(1x,a78)') '! still gets the same kind of box volume, so their optical constants are     !'
        write (stdout, '(1x,a78)') '! not on an equal footing.                                                   !'
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
      end if

      ! The layer spacings feed box_volumes and the optical path lengths, and the
      ! boxes only tile the slab exactly when consecutive spacings are equal.
      h_min = minval(layer_centroid(1:n_layers - 1) - layer_centroid(2:n_layers))
      h_max = maxval(layer_centroid(1:n_layers - 1) - layer_centroid(2:n_layers))
      h_mean = (layer_centroid(1) - layer_centroid(n_layers))/real(n_layers - 1, dp)
      if ((h_max - h_min) .gt. 0.01_dp*h_mean .and. on_root) then
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        write (temp_str, '(a,f7.4,a,f7.4,a)') 'Warning: the inferred layer spacings vary, from ', &
          h_min, ' to ', h_max, ' Ang.'
        write (stdout, '(1x,a1,1x,a74,1x,a1)') '!', adjustl(temp_str), '!'
        write (stdout, '(1x,a78)') '! Consecutive boxes then overlap or leave a gap by half that difference,     !'
        write (stdout, '(1x,a78)') '! and their volumes no longer add up to the slab. Give photo_layers_tops     !'
        write (stdout, '(1x,a78)') '! instead if that matters for this structure.                                !'
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
      end if

      deallocate (layer_species_count, species_seen, stat=ierr)
      if (ierr /= 0) call io_error('Error: analyse_geometry - deallocation of the composition check arrays failed')
    end if

    deallocate (box_boundaries, layer_index, z_gaps, layer_centroid, atoms_in_layer, stat=ierr)
    if (ierr /= 0) call io_error('Error: analyse_geometry - deallocation of the geometry work arrays failed')

    max_atoms = sum(atoms_per_box)

    ! An empty box divides by zero in the optics and leaves a hole in the
    ! layer-by-layer light attenuation, and it means the supplied boundaries do
    ! not match the structure. Stop rather than silently produce nonsense.
    do i = 1, num_boxes
      if (atoms_per_box(i) .eq. 0) then
        if (on_root) write (stdout, '(1x,a,i0,a,i0,a)') &
          'Error: layer/box ', i, ' of ', num_boxes, ' contains no atoms.'
        call io_error('Error: analyse_geometry - empty layer/box. Check photo_layers_tops '// &
                      'or photo_slab_min/photo_slab_max against the structure.')
      end if
    end do

    ! The bulk slab is given the artificial box index num_boxes + 1, which the QE
    ! calculation uses to index I_layer.
    box_atom(max_atoms + 1) = num_boxes + 1

    if (on_root) then
      ! Only the upper and lower surface have been given (SLAB_MODE_BOUNDS)
      ! so with debug printing on, the user gets the inferred box height, number of
      ! boxes/layers and the # of atoms in each layer/box
      if (.not. layers_from_input) then
        if (photo_slab_middle_set) then
          write (stdout, 229) '|  Explicit region ends at        (Ang)   :         ', photo_slab_middle, '|'
        else
          write (stdout, '(1x,a78)') &
            '|  Explicit region: top half of the inferred layers (no photo_slab_middle)   |'
        end if
      end if
      if (iprint .gt. 2 .and. .not. layers_from_input) then
        write (stdout, 420) '+', 'box height (Ang) = ', box_heights(1), ',', '# of boxes = ', num_boxes, '+'
420     format(1x, a1, 5x, a19, F13.9, a1, 12x, a13, I4, 9x, a1)
        write (stdout, 421) '+', '# of atoms in each box:', (atoms_per_box(i), i=1, num_boxes)
421     format(1x, a1, 5x, a23, 99(1x, I2))
      end if
      write (stdout, '(1x,a78)') '+------------------------------- Atomic Order  ------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  | Box/Layer |         Atom Z-Coordinate (Ang)          |'

      do atom = 1, num_atoms
        if ((box_atom(atom) .le. num_boxes)) then
          write (stdout, '(1x,a3,a2,8x,i3,11x,i3,18x,F12.7,a18)') "|  ", trim(atoms_label_tmp(atom_order(atom))), &
            atom_order(atom), box_atom(atom), atoms_pos_cart_photo(3, atom_order(atom)), "|"
        else
          write (stdout, '(1x,a3,a2,8x,i3,14x,18x,F12.7,a18)') "|  ", trim(atoms_label_tmp(atom_order(atom))), &
            atom_order(atom), atoms_pos_cart_photo(3, atom_order(atom)), "|"
        end if
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, 226) '|  Max number of atoms:', max_atoms, '  Total number of boxes:', num_boxes, '   |'
      if (layers_from_input) then
        write (stdout, '(1x,a1,76x,a1)') '|', '|'
        do i = 1, num_boxes
          write (stdout, 228) '|  Volume of box for layer ', i, ' (Ang^3) :         ', box_volumes(i), '|'
        end do
        write (stdout, '(1x,a1,76x,a1)') '|', '|'
        do i = 1, num_boxes
          write (stdout, 228) '|  Top z-coord. for layer #', i, ' (Ang)   :         ', boxes_top_z_coord(i), '|'
        end do
        write (stdout, 229) '|  Z-coord. for slab middle       (Ang)   :         ', photo_slab_middle, '|'
      else
        write (stdout, 227) '|  Volume of box for layer selection (Ang^3) :           ', box_volumes(1), '      |'
      end if
      write (stdout, 229) '|  Slab confinement length        (Ang)   :         ', slab_half_height, '|'
      write (stdout, 229) '|  Bulk repeat spacing            (Ang)   :         ', bulk_repeat, '|'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
226 format(1x, a23, I12, 1x, a25, 1x, I12, a4)
227 format(1x, a57, f14.5, a7)
228 format(1x, a27, i6, a19, f14.5, 11x, a1)
229 format(1x, a52, f14.5, 11x, a1)

    ! Test if the supplied IMFP list has same length as # of layers
    ! Otherwise, we run out of imfp values for layers
    if ((size(photo_imfp_value, 1) .gt. 1) .and. &
        (size(photo_imfp_value, 1) .ne. num_boxes)) then
      call io_error('Error : the # supplied IMFP values does not match the # layers. Check input!')
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
    ! If every band at a (spin, k) lies below E_F the search below never assigns,
    ! so start from nbands + 1: that makes every "is this an unoccupied band"
    ! test false and every "loop over the unoccupied bands" zero-trip, rather
    ! than leaving the entry undefined.
    min_index_unocc = nbands + 1

    do N_k = 1, num_kpoints_on_node(my_node_id)  ! Loop over kpoints
      do N_spin = 1, nspins
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
      do N_spin = 1, nspins
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
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a36,23x,f11.3,a8)') '+ Time to calculate Band Energy Info', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_band_info

  subroutine calc_photon_energies
    use od_constants, only: dp
    use od_parameters, only: photo_energy_sweep, photo_photon_min, photo_photon_max, jdos_spacing, photo_photon_energy
    use od_io, only: io_error
    implicit none
    real(kind=dp)        ::   num_energies, temp
    integer              ::   ierr, i

    if (photo_energy_sweep) then
      num_energies = (photo_photon_max - photo_photon_min)/jdos_spacing
      if (photo_photon_max - photo_photon_min .lt. 1.0e-12_dp) then
        number_energies = 1
      else if (mod(num_energies, 1.0_dp) .gt. 1.0E-10_dp) then
        ! The bounds must span a whole number of jdos_spacing steps, because
        ! each photon energy is mapped onto a JDOS bin index below.  A spacing
        ! such as 0.05 has no exact binary representation, so the ratio may sit
        ! just below an integer rather than on it; only that case is accepted.
        if (abs(mod(num_energies, 1.0_dp) - 1) .gt. 1.0E-10_dp) &
          call io_error('Error: calc_photon_energies - given photon sweep min/max values do not give integer # of photon steps')
      end if
      number_energies = int(num_energies) + 1
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photon_energies - allocation of index_energy failed')
      do i = 1, number_energies
        temp = (i - 1)*jdos_spacing + photo_photon_min
        ! Account for E = 0.0
        index_energy(i) = int(temp/jdos_spacing) + 1
      end do
      ! We only have one photon energy to do the calculation for.
    else
      number_energies = 1
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photon_energies - allocation of index_energy failed')
      ! Account for E = 0.0
      index_energy(number_energies) = int(photo_photon_energy/jdos_spacing) + 1
    end if
  end subroutine calc_photon_energies

  subroutine make_pdos_weights_atoms
    !!This subroutine is equivalent to pdos_merge of pdos.F90, but only for atoms
    use od_electronic, only: pdos_orbital, pdos_weights, pdos_mwab, nspins
    use od_cell, only: num_kpoints_on_node, num_atoms, cell_calc_kpoint_r_cart, kpoint_r_cart
    use od_comms, only: my_node_id, on_root
    use od_io, only: io_error, stdout, seedname, io_date, io_file_unit
    use od_parameters, only: devel_flag
    implicit none
    character(len=9) :: ctime             ! Temp. time string
    character(len=11):: cdate             ! Temp. date string
    integer :: N_k, N_spin, n_eigen, np, ierr, atom, box, i, i_max, pdos_unit
    integer, allocatable, dimension(:) :: orbital_atom

    allocate (pdos_weights_atoms(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id), num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')

    allocate (pdos_weights_k_band(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_k_band failed')

    pdos_weights_atoms = 0.0_dp
    pdos_weights_k_band = 0.0_dp

    allocate (pdos_weights_boxes(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id), num_boxes), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')
    pdos_weights_boxes = 0.0_dp

    ! The orbital -> atom map depends only on the orbital table, so build it once
    ! here rather than re-deriving it inside the k-point, spin and band loops. A
    ! new atom starts wherever either the ion number or the species changes:
    ! testing the rank alone is not enough, because the rank restarts at 1 for
    ! every species, so a species holding exactly one ion leaves the rank
    ! unchanged across the boundary and its successor would be folded into it,
    ! shifting every atom index after that against atoms_pos_cart_photo. The
    ! (species, rank) pair is what projection_utils and core already key on.
    allocate (orbital_atom(pdos_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of orbital_atom failed')
    i = 1
    orbital_atom(1) = 1
    do np = 2, pdos_mwab%norbitals
      if ((pdos_orbital%rank_in_species(np) .ne. pdos_orbital%rank_in_species(np - 1)) .or. &
          (pdos_orbital%species_no(np) .ne. pdos_orbital%species_no(np - 1))) i = i + 1
      orbital_atom(np) = i
    end do
    i_max = i

    ! Walking the orbitals has to recover exactly the atoms in the cell. Anything
    ! else means the .pdos_bin orbital ordering and the -out.cell atom ordering
    ! have parted company, so every atom-resolved quantity below would be
    ! attributed to the wrong site. Checked before the accumulation rather than
    ! after it, because i_max > num_atoms would otherwise be written past the end
    ! of pdos_weights_atoms first.
    if (i_max .ne. num_atoms) then
      if (on_root) write (stdout, '(1x,a,i0,a,i0,a)') &
        'Error: the pdos orbitals map onto ', i_max, ' atoms but the cell has ', num_atoms, '.'
      call io_error('Error: make_pdos_weights_atoms - the .pdos_bin orbital ordering does not '// &
                    'match the -out.cell atom ordering.')
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, pdos_mwab%nbands
          do np = 1, pdos_mwab%norbitals
            pdos_weights_atoms(n_eigen, N_spin, N_k, orbital_atom(np)) = &
              pdos_weights_atoms(n_eigen, N_spin, N_k, orbital_atom(np)) + &
              pdos_weights(np, n_eigen, N_k, N_spin)
          end do
        end do
      end do
    end do

    deallocate (orbital_atom, stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - failed to deallocate orbital_atom')
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
    ! We need the pdos contributions for each box to calculate the optical properties for
    ! each box representing a layer. The values are summed up for all the atoms in that
    ! specific box.
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

    if (index(devel_flag, 'output_pdos_weights') .gt. 0 .and. on_root) then
      call cell_calc_kpoint_r_cart
      write (stdout, '(a78)') "+---------------- Printing K-Points in Cartesian Coordinates ----------------+"
      i = 0
      do N_k = 1, num_kpoints_on_node(my_node_id)
        write (stdout, '(1x,I4,4x,3(1x,E22.15))') i, kpoint_r_cart(:, N_k)
        i = i + 1
      end do
      call io_date(cdate, ctime)
      ! write out atomic/box weights
      pdos_unit = io_file_unit()
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
      close (unit=pdos_unit)

      pdos_unit = io_file_unit()
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
      close (unit=pdos_unit)

      ! Write out the k-band weights
      pdos_unit = io_file_unit()
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

    if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root) then
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
    use od_io, only: stdout, io_error, io_time, seedname, io_date, io_file_unit
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
    integer :: N_k, N2, N_spin, n_eigen, n_eigen_final, ierr, energy, box
    integer :: jdos_bin, i, s, is, idos, wjdos_unit, initial, ome_unit
    real(kind=dp)    :: time0, time1
    character(len=3) :: atom_s
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    character(len=4) :: initial_s

    time0 = io_time()

    allocate (absorp_photo(num_boxes, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')

    allocate (reflect_photo(num_boxes, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')
    allocate (epsilon_photo(num_boxes, number_energies, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of epsilon_photo failed')
    epsilon_photo = 0.0_dp

    ! Advanced and tricky user option to read the optical properties from a previous run
    ! Must be used with caution, since currently no checking of parameters is performed!
    if (index(devel_flag, 'optics_restart') .gt. 0) then
      call setup_energy_scale(E)
      if (on_root) then
        if (.not. allocated(absorp)) then
          allocate (absorp(jdos_nbins), stat=ierr)
          if (ierr /= 0) call io_error("Error: calc_photo_optics cannot allocate absorp")
        end if
        if (.not. allocated(reflect)) then
          allocate (reflect(jdos_nbins), stat=ierr)
          if (ierr /= 0) call io_error("Error: calc_photo_optics cannot allocate reflect")
        end if
        call read_absorp_file
        call read_reflect_file
      end if
      call comms_bcast(absorp_photo(1, 1), num_boxes*number_energies)
      call comms_bcast(reflect_photo(1, 1), num_boxes*number_energies)

      time1 = io_time()
      if (on_root .and. iprint .gt. 1) then
        write (stdout, '(1x,a47,12x,f11.3,a8)') '+ Time to read Photoemission Optical Properties', time1 - time0, ' (sec) +'
      end if

      call make_weights(matrix_weights)
      call elec_dealloc_optical

      if (index(photo_model, '3step') .gt. 0 .or. index(photo_model, 'dosds') .gt. 0) then
        ! Flip the kpt and spin indices in the matrix_weights array for contiguous memory access later
        allocate (photo_matrix_weights(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)))
        if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of photo_matrix_weights failed')

        do N_spin = 1, nspins
          do N_k = 1, num_kpoints_on_node(my_node_id)
            photo_matrix_weights(:, :, N_spin, N_k) = matrix_weights(:, :, N_k, N_spin, 1)
          end do
        end do
      end if
      ! get rid of the old, now unnecessary array - either because we have the 1step model,
      ! or we have transferred the relevant data to photo_matrix_weights
      deallocate (matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate photo_matrix_weights')
      N_geom = 1
      return
    end if

    call make_weights(matrix_weights)
    N_geom = size(matrix_weights, 5)
    call elec_dealloc_optical

    if (index(devel_flag, 'output_ome_itof') .gt. 0 .and. on_root) then
      is = -1
      call io_date(cdate, ctime)
      do N_k = 1, size(kpoint_r, 2)
        if (all(abs(kpoint_r(:, N_k)) .lt. 1.0E-10_dp)) then
          is = N_k
          exit
        end if
      end do
      if (is .lt. 0) call io_error('Error: this devel_flag should be run in serial. No gamma point found on root.')
      write (stdout, *) 'The gamma point was determined to be - ', is
      i = index(devel_flag, 'output_ome_itof')
      read (devel_flag(i + 16:i + 20), *) initial
      write (initial_s, '(I4)') initial
      write (stdout, '(1x,a27,I4,a32)') 'Outputting OMEs for band # ', initial, ' to the bands above it at Gamma.'
      ome_unit = io_file_unit()
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
    end if

    if (index(photo_model, 'dosds') .eq. 0) then
      allocate (projected_matrix_weights(nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of projected_matrix_weights failed')
      do box = 1, num_boxes                           ! Loop over boxes
        !
        if (iprint .gt. 1 .and. on_root) then
          write (stdout, 145) '+--------------------- Starting BOX/Layer  # ', box, ' of ', num_boxes, ' ---------------------+'
        end if
        ! (Re-)Setting the weights for new box
        projected_matrix_weights = 0.0_dp

        do N2 = 1, N_geom
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              do n_eigen = 1, nbands               ! Loop over state 1
                do n_eigen_final = n_eigen, nbands    ! Loop over state 2
                  if (band_energy(n_eigen, N_spin, N_k) .gt. efermi .and. n_eigen /= n_eigen_final) cycle
                  if (band_energy(n_eigen_final, N_spin, N_k) .lt. efermi .and. n_eigen /= n_eigen_final) cycle
                  if (pdos_weights_k_band(n_eigen, N_spin, N_k) .eq. 0.0_dp) then
                    cycle
                  end if
                  projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2) = &
                    matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2)* &
                    (pdos_weights_boxes(n_eigen, N_spin, N_k, box)/pdos_weights_k_band(n_eigen, N_spin, N_k))
                end do                        ! Loop over state 2
              end do                            ! Loop over state 1
            end do
          end do
        end do

        if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root) then
          write (stdout, '(1x,a37,I3,a38)') '+-------------------------------Atom-', box, &
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
        call jdos_utils_calculate(projected_matrix_weights, weighted_jdos=weighted_jdos, &
                                  slab_half_height=slab_half_height)

        if (on_root .and. iprint .gt. 2) then
          write (atom_s, '(I3)') box + 100
          wjdos_unit = io_file_unit()
          open (unit=wjdos_unit, action='write', file=trim(seedname)//'_weighted_jdos_'//trim(adjustl(atom_s))//'.dat')
          write (wjdos_unit, '(1x,a28)') '############################'
          write (wjdos_unit, '(1x,a19,1x,a99)') '# Weighted JDOS for', seedname
          write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# maximum JDOS energy :', jdos_max_energy, '[eV]'
          write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# JDOS step size      :', jdos_spacing, '[eV]'
          write (wjdos_unit, '(1x,a28)') '############################'
          do is = 1, nspins
            write (wjdos_unit, *) 'Spin Channel :', is
            do idos = 1, jdos_nbins
              write (wjdos_unit, *) E(idos), ' , ', sum(weighted_jdos(idos, is, 1:size(matrix_weights, 5)))
            end do
          end do
          close (unit=wjdos_unit)
        end if

        if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root) then
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
          do N_geom = 1, size(matrix_weights, 5)
            do n_eigen = 1, nbands
              dos_matrix_weights(N_geom, n_eigen, :, :) = matrix_weights(n_eigen, n_eigen, :, :, N_geom)
            end do
          end do
          N_geom = size(matrix_weights, 5)
          call dos_utils_calculate_at_e(efermi, dos_at_e, dos_matrix_weights, weighted_dos_at_e)
          ! No per-atom normalisation here: weighted_dos_at_e is the sum over the
          ! atoms in the box, exactly like weighted_jdos, and calc_epsilon_2
          ! divides both by the same box_volumes(box). Dividing only the
          ! intraband term by atoms_per_box would give the two terms different
          ! normalisations and would break lateral supercells, where doubling
          ! the in-plane cell correctly doubles box_volumes but must not also
          ! introduce a spurious 1/N_atoms.
        end if

        if (on_root) then
          if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. optics_intraband) then
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
          call calc_epsilon_2(weighted_jdos, weighted_dos_at_e, box_volumes(box))

          ! Calculate epsilon_1
          call calc_epsilon_1

          ! Calculate other optical properties
          call calc_refract
          call calc_absorp
          call calc_reflect

          if (iprint .gt. 2) then
            call write_epsilon(box, photo_at_e=dos_at_e, photo_volume=box_volumes(box))
            call write_refract(box, photo_volume=box_volumes(box))
            call write_absorp(box, photo_volume=box_volumes(box))
            call write_reflect(box, photo_volume=box_volumes(box))
          end if

          do energy = 1, number_energies
            absorp_photo(box, energy) = absorp(index_energy(energy))
            reflect_photo(box, energy) = reflect(index_energy(energy))
            epsilon_photo(box, energy, 1) = epsilon(index_energy(energy), 1, 1, 1)
            epsilon_photo(box, energy, 2) = epsilon(index_energy(energy), 2, 1, 1)
          end do

          if (index(devel_flag, 'print_qe_constituents') .gt. 0) then
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
            write (stdout, '(99(E17.8E3))') (absorp_photo(box, energy), energy=1, number_energies)
            write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

            write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection --------------------------+'
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
      end do                                        ! Loop over boxes
      call comms_bcast(absorp_photo(1, 1), num_boxes*number_energies)
      call comms_bcast(reflect_photo(1, 1), num_boxes*number_energies)
    end if
145 format(1x, a45, I3, a4, I3, a23)
124 format(3(1x, I4))
125 format(4(1x, I4))
126 format(5(1x, I4))

    ! Deallocating this out of the loop to reduce memory operations - could lead to higher memory consumption
    deallocate (projected_matrix_weights, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate projected_matrix_weights')
    if (index(photo_model, '3step') .gt. 0 .or. index(photo_model, 'dosds') .gt. 0) then
      ! Flip the kpt and spin indices in the matrix_weights array for contiguous memory access later
      allocate (photo_matrix_weights(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)))
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of photo_matrix_weights failed')

      do N_spin = 1, nspins
        do N_k = 1, num_kpoints_on_node(my_node_id)
          photo_matrix_weights(:, :, N_spin, N_k) = matrix_weights(:, :, N_k, N_spin, 1)
        end do
      end do
    end if
    ! get rid of the old, now unnecessary array - either because we have the 1step model,
    ! or we have transferred the relevant data to photo_matrix_weights
    deallocate (matrix_weights, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate photo_matrix_weights')

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a52,7x,f11.3,a8)') '+ Time to calculate Photoemission Optical Properties', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_photo_optics

  subroutine read_absorp_file
    !*******=======================================================================
    ! ***Experimental subroutine not intended for use by the normal user.***
    ! The absorption files are read without any checks, so inherently
    ! unsafe for now!
    ! This subroutine reads in a series of absorption coefficient curves
    ! from a number of appropriately named files. This way the relevant
    ! optical data for photoemission can be read as a checkpoint. This can
    ! be used to for example calculate the photoemission for a set of
    ! k-points along a bandstructure path with the optical properties of
    ! a MP grid like k-point distribution, as that is expected to have better
    ! convergence.
    ! Written by F Mildner, Mar 2025
    !===============================================================================
    use od_optics, only: absorp
    use od_jdos_utils, only: jdos_nbins
    use od_io, only: seedname, io_file_unit, io_error

    integer :: absorp_unit, box, i, N, ierr, energy
    character(len=3) :: box_char
    character(len=100) :: dummya
    absorp_unit = io_file_unit()

    do box = 1, num_boxes

      write (box_char, '(I0.3)') box
      open (unit=absorp_unit, file=trim(seedname)//'_absorption_photo_box_'//trim(adjustl(box_char))//'.dat', iostat=ierr)
      if (ierr /= 0) call io_error('Error: Could not open absorption curve .dat file for box #'//trim(adjustl(box_char)))
      ! skip header
      do i = 1, 50
        read (absorp_unit, *) dummya
        if (index(dummya, '#') .eq. 0) exit
      end do
      do N = 2, jdos_nbins
        read (absorp_unit, '(1x,a37,1x,es37.30)') dummya, absorp(N)
      end do
      close (unit=absorp_unit)

      do energy = 1, number_energies
        absorp_photo(box, energy) = absorp(index_energy(energy))
      end do

    end do
  end subroutine read_absorp_file

  subroutine read_reflect_file
    !*******=======================================================================
    ! ***Experimental subroutine not intended for use by the normal user.***
    ! The reflectivity files are read without any checks, so inherently
    ! unsafe for now!
    ! This subroutine reads in a series of reflection coefficient curves
    ! from a number of appropriately named files. This way the relevant
    ! optical data for photoemission can be read as a checkpoint. This can
    ! be used to for example calculate the photoemission for a set of
    ! k-points along a bandstructure path with the optical properties of
    ! a MP grid like k-point distribution, as that is expected to have better
    ! convergence.
    ! Written by F Mildner, Mar 2025
    !===============================================================================
    use od_optics, only: reflect
    use od_jdos_utils, only: jdos_nbins
    use od_io, only: seedname, io_file_unit, io_error

    integer :: reflect_unit, box, i, N, ierr, energy
    character(len=3) :: box_char
    character(len=100) :: dummya

    reflect_unit = io_file_unit()

    do box = 1, num_boxes
      write (box_char, '(I0.3)') box
      open (unit=reflect_unit, file=trim(seedname)//'_reflection_photo_box_'//trim(adjustl(box_char))//'.dat', iostat=ierr)
      if (ierr /= 0) call io_error('Error: Could not open absorption curve .dat file for box #'//trim(adjustl(box_char)))
      ! skip header
      do i = 1, 50
        read (reflect_unit, *) dummya
        if (index(dummya, '#') .eq. 0) exit
      end do
      do N = 2, jdos_nbins
        read (reflect_unit, '(1x,a37,1x,es37.30)') dummya, reflect(N)
      end do
      close (unit=reflect_unit)

      do energy = 1, number_energies
        reflect_photo(box, energy) = reflect(index_energy(energy))
      end do

    end do
  end subroutine read_reflect_file

  subroutine calc_absorp_layer
    !*******=======================================================================
    ! This subroutine calculates the absorption coefficient for all defined/
    ! inferred layers in the slab structure.
    !=======================================================================
    use od_io, only: io_error
    implicit none
    real(kind=dp) :: I_0
    integer :: box, i, ierr

    allocate (I_layer(num_boxes + 1, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of I_layer failed')
    I_layer = 1.0_dp

    I_0 = 1.0_dp

    ! Calculate the unreflected portion of incoming light. The reflectivity is
    ! that of the slab as a stratified medium, not of the surface box on its own
    ! -- see slab_reflectivity. Layer 1 then receives all of the transmitted
    ! light: there is no material above it, so nothing can attenuate the beam
    ! before it reaches the first plane of atoms, and every deeper layer is
    ! reached from there by the centroid-to-centroid recursion below.
    do i = 1, number_energies
      if (allocated(epsilon_photo)) then
        I_layer(1, i) = I_0 - slab_reflectivity(i)
      else
        ! read_reflect_file supplies R per box and no epsilon, so there is
        ! nothing to average and the old surface-box value is all there is.
        I_layer(1, i) = I_0 - reflect_photo(1, i)
      end if
    end do
    ! If we have more than one box with atoms in it, calculate the incident light intensity for each
    !
    ! The step from the atoms of layer box-1 to those of layer box crosses the
    ! material lying between them, which is not the box the light arrives in. The
    ! step is layer_gap(box-1), the spacing between the two layers' centroids,
    ! and the boundary between the boxes lies inside it, so the light crosses part
    ! of box-1 and part of box. Hence the mean of the two absorption coefficients
    ! over that one spacing.
    !
    ! This is not a no-op for a geometrically uniform slab: the boxes carry
    ! different projected optical responses, so the surface layer absorbs
    ! differently from the interior even when the spacings are identical.
    ! bulk_emission keeps a single coefficient below the explicit region, which is
    ! right there - every step it takes is through repeats of the same deepest box.
    if (num_boxes .gt. 1) then
      do box = 2, num_boxes
        do i = 1, number_energies
          I_layer(box, i) = I_layer(box - 1, i)* &
                            exp(-(0.5_dp*(absorp_photo(box - 1, i) + absorp_photo(box, i)) &
                                  *layer_gap(box - 1)*1E-10))
          if (I_layer(box, i) .lt. 0.0_dp) I_layer(box, i) = 0.0_dp
        end do
      end do
    end if
    ! Since we later combine the bulk slab emission probability (contains already light intensity) into the
    ! layer by layer emission probability array (does not contain light intensity), we have to set the
    ! intensity value artifically to 1.0 to have it not influence the final value.
    ! We are only ever accessing I_layer to max_atoms, so this has no effect on the rest.
    I_layer(box_atom(max_atoms + 1), 1:number_energies) = 1.0_dp

    if (allocated(reflect_photo)) then
      deallocate (reflect_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate reflect_photo')
    end if
    if (allocated(epsilon_photo)) then
      deallocate (epsilon_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate epsilon_photo')
    end if
  end subroutine calc_absorp_layer

  function slab_reflectivity(energy) result(reflectivity)
    !*******=======================================================================
    ! Reflectivity of the slab, seen as a semi-infinite stratified medium.
    !
    ! A Fresnel reflectivity samples the material over the absorption depth
    ! 1/alpha, which for Cu at 6 eV is around 100 Ang against an explicit region
    ! of 13 Ang. Taking it from the surface box alone therefore sets the reflected
    ! fraction from the optical response of the layer that accounts for a few per
    ! cent of the light's interaction, and it is the box whose height is the most
    ! arbitrary, being fixed by photo_slab_max. Measured: moving photo_slab_max by
    ! 0.5 Ang changed R by 12 % and the quantum efficiency by 8 %, because
    ! eps_2 goes as 1/box_volume.
    !
    ! So average the dielectric function over depth, weighted by how much of the
    ! light each box actually sees, and take the Fresnel expression of the
    ! average. Three things this has to get right:
    !
    !   * Average epsilon, not R and not n. Kramers-Kronig is linear in eps_2, so
    !     averaging eps_2 and transforming gives the average of the eps_1 values;
    !     n, kappa, alpha and R are all non-linear, so <R>, R(<n>) and R(<eps>)
    !     are three different numbers and only the last one means anything.
    !   * Arithmetic, not harmonic. At normal incidence the field lies in the
    !     layer plane, which is the arithmetic Wiener bound; the harmonic bound is
    !     for a field along the surface normal.
    !   * The tail. Everything below the explicit region is the deepest box
    !     repeated, and for a clean metal it carries most of the weight.
    !
    ! The probe depth depends on the average that is being computed with it, so
    ! iterate. alpha goes as kappa at fixed photon energy, which lets the code's
    ! own absorption coefficient be rescaled rather than rebuilt from constants.
    !===============================================================================
    use od_constants, only: dp
    implicit none
    integer, intent(in) :: energy
    real(kind=dp) :: reflectivity
    integer       :: box, iter
    real(kind=dp) :: probe_depth, depth, weight, weight_sum, sum_1, sum_2
    real(kind=dp) :: eps_1, eps_2, modulus, n_index, kappa, kappa_bulk, alpha

    eps_1 = epsilon_photo(num_boxes, energy, 1)
    eps_2 = epsilon_photo(num_boxes, energy, 2)
    kappa_bulk = kappa_of(eps_1, eps_2)
    alpha = absorp_photo(num_boxes, energy)

    ! Without an absorption coefficient there is no depth to weight over, so
    ! fall back to the deepest box, which is what the weighting tends to.
    if (alpha .le. 0.0_dp .or. kappa_bulk .le. 0.0_dp) then
      reflectivity = fresnel(eps_1, eps_2)
      return
    end if

    do iter = 1, 8
      probe_depth = 1.0E10_dp/alpha        ! 1/alpha in Angstrom
      depth = 0.0_dp
      sum_1 = 0.0_dp; sum_2 = 0.0_dp; weight_sum = 0.0_dp
      do box = 1, num_boxes
        weight = exp(-depth/probe_depth) - exp(-(depth + box_heights(box))/probe_depth)
        sum_1 = sum_1 + epsilon_photo(box, energy, 1)*weight
        sum_2 = sum_2 + epsilon_photo(box, energy, 2)*weight
        weight_sum = weight_sum + weight
        depth = depth + box_heights(box)
      end do
      weight = exp(-depth/probe_depth)     ! the bulk repeat below the explicit region
      sum_1 = sum_1 + epsilon_photo(num_boxes, energy, 1)*weight
      sum_2 = sum_2 + epsilon_photo(num_boxes, energy, 2)*weight
      weight_sum = weight_sum + weight
      if (weight_sum .le. 0.0_dp) exit
      eps_1 = sum_1/weight_sum
      eps_2 = sum_2/weight_sum
      kappa = kappa_of(eps_1, eps_2)
      if (kappa .le. 0.0_dp) exit
      alpha = absorp_photo(num_boxes, energy)*kappa/kappa_bulk
    end do

    reflectivity = fresnel(eps_1, eps_2)
  end function slab_reflectivity

  function kappa_of(eps_1, eps_2) result(kappa)
    !! Extinction coefficient from the dielectric function.
    use od_constants, only: dp
    implicit none
    real(kind=dp), intent(in) :: eps_1, eps_2
    real(kind=dp) :: kappa, modulus
    modulus = sqrt(eps_1**2 + eps_2**2)
    kappa = sqrt(max(0.5_dp*(modulus - eps_1), 0.0_dp))
  end function kappa_of

  function fresnel(eps_1, eps_2) result(reflectivity)
    !! Normal-incidence reflectivity of a semi-infinite medium.
    use od_constants, only: dp
    implicit none
    real(kind=dp), intent(in) :: eps_1, eps_2
    real(kind=dp) :: reflectivity, modulus, n_index, kappa
    modulus = sqrt(eps_1**2 + eps_2**2)
    n_index = sqrt(max(0.5_dp*(modulus + eps_1), 0.0_dp))
    kappa = sqrt(max(0.5_dp*(modulus - eps_1), 0.0_dp))
    reflectivity = (((n_index - 1.0_dp)**2) + kappa**2)/(((n_index + 1.0_dp)**2) + kappa**2)
  end function fresnel

  subroutine effective_wf
    use od_parameters, only: photo_work_function, photo_elec_field
    use od_electronic, only: efermi
    use od_constants, only: pi, epsilon_0, e_charge, j_to_ev, ev_to_j
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    implicit none
    real(kind=dp) :: schottky_lowering

    !photo_elec_field given in V/m
    schottky_lowering = sqrt((e_charge**3*photo_elec_field)/(4*pi*epsilon_0*ev_to_j**2))
    work_function_eff = photo_work_function - schottky_lowering
    ! Nothing downstream survives a barrier that is not positive. evacuum_eff
    ! would drop below the Fermi level, which breaks the assumption that a state
    ! too deep to emit is also too deep to have an escape angle, and in the
    ! fallback branch of the surface barrier it makes theta_internal evaluate
    ! acos(sqrt(negative)). The NaN that follows compares false against every
    ! guard, so the states are dropped silently and the run still prints a
    ! number. Stop here instead.
    if (work_function_eff .le. 0.0_dp) then
      if (on_root) then
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        write (stdout, '(1x,a78)') '! Error: the Schottky lowering from photo_elec_field is at least as large as !'
        write (stdout, '(1x,a78)') '! the work function, so the effective barrier is zero or negative.  That is  !'
        write (stdout, '(1x,a78)') '! the field emission regime, which the three step model does not describe,   !'
        write (stdout, '(1x,a78)') '! and it makes the surface refraction take the square root of a negative     !'
        write (stdout, '(1x,a78)') '! barrier.  Lower photo_elec_field, or set the barrier explicitly with       !'
        write (stdout, '(1x,a78)') '! photo_inner_potential and treat the emission with a field emission model.  !'
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        write (stdout, '(1x,a46,1x,f10.4,20x,a1)') '|  Work function (eV)                        :', &
          photo_work_function, '|'
        write (stdout, '(1x,a46,1x,f10.4,20x,a1)') '|  Schottky lowering (eV)                    :', &
          schottky_lowering, '|'
        write (stdout, '(1x,a46,1x,f10.4,20x,a1)') '|  Effective work function (eV)              :', &
          work_function_eff, '|'
      end if
      call io_error('Error: effective_wf - photo_elec_field lowers the work function to zero or below')
    end if
    evacuum_eff = work_function_eff + efermi
  end subroutine effective_wf

  subroutine calc_field_emission
    !*******=======================================================================
    ! This subroutine calculates the Schottky effect and emission
    ! probabilities of an electron through the surface barrier.
    ! parameter photo_elec_field given in V/A
    ! orig. Victor Chang
    ! updated by Felix Mildner, after Mar 2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node
    use od_parameters, only: photo_work_function, photo_elec_field, photo_temperature, iprint
    use od_electronic, only: efermi, band_energy, nbands, nspins
    use od_io, only: io_error, stdout, io_time
    use od_comms, only: my_node_id, comms_reduce, on_root
    use od_constants, only: pi, epsilon_0, kB, ev_to_j, e_charge
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :) :: temp_emission
    real(kind=dp) :: field_energy_squared, fermi_dirac, barrier_height, argument, exponent
    real(kind=dp) :: transmission_prob
    real(kind=dp)    :: time0, time1
    integer :: N_k, N_spin, n_eigen, ierr

    time0 = io_time()
    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of field_emission failed')
    end if
    field_emission = 0.0_dp

    field_energy_squared = 0.0_dp

    allocate (temp_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of temp_emission failed')
    temp_emission = 0.0_dp

    evacuum = efermi + photo_work_function

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          ! Calculate how much the final electron energy is below
          ! the unmodified work function.
          ! initial band energy + photon energy
          barrier_height = evacuum - band_energy(n_eigen, N_spin, N_k) + temp_photon_energy
          field_energy_squared = (barrier_height)**2
          ! Calculate the fermi dirac occupations
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but exp(+-230) ~ 1E(+-100)
          ! so this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible under/over-flow.
          if (argument .gt. 230.0_dp) then
            fermi_dirac = 0.0_dp
          elseif (argument .lt. -230.0_dp) then
            fermi_dirac = 1.0_dp
          else
            fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          ! Calculating if the "scaled barrier field" - f - is 0 < f < 1 and not f > 1
          ! otherwise the integral borders are not real and the approximation is not defined.
          ! This can happen if the final electron energy is above the barrier, or if the field
          ! is very strong.
          if ((e_charge**3*photo_elec_field)/(4*pi*epsilon_0*ev_to_j**2) .lt. field_energy_squared) then
            if (barrier_height .le. 0.0_dp) then
              field_emission(n_eigen, N_spin, N_k) = 1.0_dp
            else
              call compute_G(barrier_height, photo_elec_field, exponent)
              if ((exponent .lt. -230.0_dp)) then
                transmission_prob = 1.0_dp
              else if (exponent .gt. 230.0_dp) then
                transmission_prob = 0.0_dp
              else
                transmission_prob = exp(-1.0_dp*exponent)
              end if
              field_emission(n_eigen, N_spin, N_k) = transmission_prob
            end if
            ! setting this to 1, so if the scaled barrier field (see long comment above) is outside the
            ! range, which is only achieved if the rounded barrier field is lowered so much that emission
            ! probability is ~1 (high final energy or high electric field (i.e. strongly lowered barrier)).
          else
            field_emission(n_eigen, N_spin, N_k) = 1.0_dp
          end if
          temp_emission(n_eigen, N_spin, N_k) = field_emission(n_eigen, N_spin, N_k)*fermi_dirac
        end do
      end do
    end do

    total_field_emission = sum(temp_emission(1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id)))/cell_area
    call comms_reduce(total_field_emission, 1, "SUM")
    deallocate (temp_emission, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate temp_emission')

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a48,11x,f11.3,a8)') '+ Time to calculate Field Emission Probabilities', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_field_emission

  subroutine compute_G(barrier_eV, F, G)
    use od_io, only: stdout
    use od_constants, only: dp, pi, e_mass, h_planck, ev_to_j, e_charge, epsilon_0
    implicit none
    real(kind=dp), intent(in)  :: barrier_eV, F
    real(kind=dp), intent(out) :: G
    ! The tolerance for adaptive simpson integrator
    real(kind=dp), parameter :: eps = 1.0e-14_dp
    ! Max recursion depth for adaptive simpson integrator
    integer, parameter       :: max_depth = 30
    real(kind=dp) :: phi, g_e_si, tmp, z1, z2
    real(kind=dp) :: fa, fm, fb, whole, integral
    real(kind=dp) :: mval, Mmid, disc, root

    phi = barrier_eV*ev_to_j
    g_e_si = 4.0_dp*pi*sqrt(2.0_dp*e_mass)/h_planck

    ! analytic turning points
    disc = phi*phi - (e_charge**3*F)/(4.0_dp*pi*epsilon_0)

    if (disc < 0.0_dp) then
      write (stdout, *) 'No real turning points: discriminant < 0'
      write (stdout, *) 'discriminant = ', disc
      stop
    end if

    root = sqrt(disc)
    z1 = (phi - root)/(2.0_dp*e_charge*F)
    z2 = (phi + root)/(2.0_dp*e_charge*F)

    ! Make sure, that z2 > z1
    if (z1 > z2) then
      tmp = z1
      z1 = z2
      z2 = tmp
    end if

    ! Calculate barrier height M at z1 and calculate fa = sqrt(M) if M > 0 or fa = 0 if M < 0
    call barrier_M(z1, phi, F, mval)
    call sqrt_clamped(mval, fa)

    ! Calculate barrier height M at z2 and calculate fb = sqrt(M) if M > 0 or fb = 0 if M < 0
    call barrier_M(z2, phi, F, mval)
    call sqrt_clamped(mval, fb)

    ! Calculate barrier height M at middlepoint z and calculate fm = sqrt(M) if M > 0 or fm = 0 if M < 0
    mval = 0.5_dp*(z1 + z2)
    call barrier_M(mval, phi, F, Mmid)
    call sqrt_clamped(Mmid, fm)

    whole = (z2 - z1)*(fa + 4.0_dp*fm + fb)/6.0_dp

    call adaptive_simpson(z1, z2, fa, fm, fb, whole, phi, F, eps, max_depth, integral)
    G = g_e_si*integral
  end subroutine compute_G

  subroutine barrier_M(z, phi, F, M)
    ! Helper function, which computes the value of the Schottky-Nordheim barrier
    ! function
    use od_constants, only: dp, pi, e_charge, epsilon_0
    implicit none
    real(kind=dp), intent(in)  :: z, phi, F
    real(kind=dp), intent(out) :: M
    real(kind=dp) :: image_term

    image_term = (e_charge**2)/(16.0_dp*pi*epsilon_0*z)
    M = phi - e_charge*F*z - image_term
  end subroutine barrier_M

  subroutine sqrt_clamped(x, y)
    use od_constants, only: dp
    implicit none
    real(kind=dp), intent(in)  :: x
    real(kind=dp), intent(out) :: y

    if (x < 0.0_dp) then
      y = 0.0_dp
    else
      y = sqrt(x)
    end if
  end subroutine sqrt_clamped

  recursive subroutine adaptive_simpson(a, b, fa, fm, fb, whole, phi, F, eps, depth, integral)
    ! adaptive Simpson rule based integrator, which calculates the integral until the
    ! integral is within eps
    use od_constants, only: dp
    implicit none
    real(kind=dp), intent(in)   :: a, b, fa, fm, fb, whole, phi, F, eps
    integer, intent(in)   :: depth
    real(kind=dp), intent(out)  :: integral
    real(kind=dp) :: m, lm, rm
    real(kind=dp) :: flm, frm
    real(kind=dp) :: left, right
    real(kind=dp) :: Mval
    real(kind=dp) :: left_int, right_int

    if (depth <= 0) then
      integral = whole
      return
    end if

    m = 0.5_dp*(a + b)
    lm = 0.5_dp*(a + m)
    rm = 0.5_dp*(m + b)

    call barrier_M(lm, phi, F, Mval)
    call sqrt_clamped(Mval, flm)

    call barrier_M(rm, phi, F, Mval)
    call sqrt_clamped(Mval, frm)

    left = (m - a)*(fa + 4.0_dp*flm + fm)/6.0_dp
    right = (b - m)*(fm + 4.0_dp*frm + fb)/6.0_dp

    if (abs(left + right - whole) < 15.0_dp*eps) then
      integral = left + right + (left + right - whole)/15.0_dp
    else
      call adaptive_simpson(a, m, fa, flm, fm, left, phi, F, eps/2.0_dp, depth - 1, left_int)
      call adaptive_simpson(m, b, fm, frm, fb, right, phi, F, eps/2.0_dp, depth - 1, right_int)
      integral = left_int + right_int
    end if
  end subroutine adaptive_simpson

  !===============================================================================
  function phi_accepted(kx, ky) result(inside)
    !*===============================================================================
    ! Is this transverse momentum inside the azimuthal acceptance wedge?
    !
    ! The wedge is a direction, photo_phi_centre, plus a half width. Testing it as
    ! |phi - centre| <= halfwidth would need wrap-around arithmetic; the angle
    ! between the momentum and the acceptance direction is the same condition
    ! written as a dot product, which wraps for free and needs no atan2:
    !
    !   cos(angle) = (k . u)/|k| >= cos(halfwidth),   u = (cos centre, sin centre)
    !
    ! A half width of 180 deg gives cos = -1 and accepts everything, as it must.
    ! Normal emission has no azimuth at all and is accepted by any wedge.
    !===============================================================================
    use od_constants, only: deg_to_rad
    use od_parameters, only: photo_phi_centre, photo_phi_halfwidth
    implicit none
    real(kind=dp), intent(in) :: kx, ky
    logical                   :: inside
    real(kind=dp)             :: k_norm
    real(kind=dp), parameter  :: k_tol = 1.0E-12_dp

    k_norm = sqrt(kx*kx + ky*ky)
    if (k_norm .lt. k_tol) then
      inside = .true.
      return
    end if
    inside = (kx*cos(photo_phi_centre*deg_to_rad) + ky*sin(photo_phi_centre*deg_to_rad))/k_norm &
             .ge. cos(photo_phi_halfwidth*deg_to_rad)
  end function phi_accepted

  !===============================================================================
  function phi_star_fraction(kx, ky) result(frac)
    !*===============================================================================
    ! The fraction of the symmetry star of this transverse momentum that lies
    ! inside the acceptance wedge.
    !
    ! An irreducible k-point carries the weight of its whole star, and the crystal
    ! symmetry operations rotate the transverse momentum around the surface
    ! normal, so a single irreducible point emits into as many azimuths as it has
    ! images. The outputs that resolve the azimuth place each image separately and
    ! test it with phi_accepted; the ones that do not - the EDC and the maps whose
    ! momentum axis is |k| - have nowhere to put them, and the right weight for
    ! them is the fraction of the star that survives the wedge.
    !
    ! With one operation this is 0 or 1, identical to phi_accepted, so a structure
    ! with no symmetry needs no special case.
    !===============================================================================
    use od_cell, only: num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_parameters, only: devel_flag
    implicit none
    real(kind=dp), intent(in) :: kx, ky
    real(kind=dp)             :: frac
    real(kind=dp)             :: image(2)
    integer                   :: nsymm_op, n_accept, n_symm

    if (index(devel_flag, 'no_symmetry') .gt. 0) then
      frac = 0.0_dp
      if (phi_accepted(kx, ky)) frac = 1.0_dp
      return
    end if

    n_symm = 0
    if (allocated(crystal_symmetry_operations)) n_symm = num_crystal_symmetry_operations
    if (n_symm .lt. 1) then
      ! No symmetry block in the -out.cell, so the point stands only for itself.
      frac = 0.0_dp
      if (phi_accepted(kx, ky)) frac = 1.0_dp
      return
    end if

    n_accept = 0
    do nsymm_op = 1, n_symm
      image = matmul(crystal_symmetry_operations(1:2, 1:2, nsymm_op), (/kx, ky/))
      if (phi_accepted(image(1), image(2))) n_accept = n_accept + 1
    end do
    frac = real(n_accept, dp)/real(n_symm, dp)
  end function phi_star_fraction

  subroutine calc_angle
    !*******=======================================================================
    ! This subroutine calculates the photoemission angles theta and phi
    ! Theta: angle between the photoemitted electron and the surface normal
    ! Phi: angle between the photoemission direction and the x axis
    ! orig. Victor Chang, 7th February 2020
    ! parts rewritten Felix Mildner, after Mar 2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart
    use od_electronic, only: nbands, nspins, band_energy, band_gradient, elec_read_band_gradient, &
      photo_gkgrid, elec_read_gk_grid
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_momentum, devel_flag, iprint, &
      photo_inner_potential, photo_inner_potential_set
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout, io_time
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: hbar, ev_to_j, j_to_ev, e_mass, rad_to_deg
    implicit none
    integer :: N_k, N_spin, n_eigen, ierr, gdx

    real(kind=dp), allocatable, dimension(:, :, :, :):: E_x
    real(kind=dp), allocatable, dimension(:, :, :, :):: E_y
    real(kind=dp) :: tol = 1.0E-10_dp
    real(kind=dp) :: time0, time1

    time0 = io_time()

    if (index(photo_momentum, 'gkgrid') .gt. 0) then
      call elec_read_gk_grid()
      photo_gkmax = size(photo_gkgrid, 2)

      if (.not. allocated(gkgrid_weight)) then
        allocate (gkgrid_weight(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_angle - allocation of gkgrid_weight failed')
      end if
      ! move the important spectral weight into the smaller array for later use
      gkgrid_weight(1:photo_gkmax, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id)) = &
        photo_gkgrid(3, 1:photo_gkmax, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id))
    else
      ! populate the gkgrid_weight array with 1 to not alter the final values
      ! so we do not have to write multiple functions for gkgrid or crystal options
      ! later and can just do the multiplication for either.
      photo_gkmax = 1
      if (.not. allocated(gkgrid_weight)) then
        allocate (gkgrid_weight(1, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_angle - allocation of gkgrid_weight failed')
      end if
      gkgrid_weight = 1.0_dp
    end if

    if (.not. allocated(E_transverse)) then
      allocate (E_transverse(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_transverse failed')
    end if
    E_transverse = 0.0_dp

    if (.not. allocated(theta_arpes)) then
      allocate (theta_arpes(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of theta_arpes failed')
    end if
    ! Impossible value as default that is equal to no emission
    theta_arpes = 91.0_dp

    if (.not. allocated(theta_internal)) then
      allocate (theta_internal(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of theta_internal failed')
    end if
    ! Impossible value as default that is equal to no emission
    theta_internal = 91.0_dp

    if (.not. allocated(phi_accept_frac)) then
      allocate (phi_accept_frac(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of phi_accept_frac failed')
    end if
    ! Nothing accepted until the transverse momentum says otherwise. States that
    ! cannot emit are dropped by theta, which keeps its 91 deg sentinel; do not
    ! rely on this default to exclude them.
    phi_accept_frac = 0.0_dp

    if (.not. allocated(E_kinetic)) then
      allocate (E_kinetic(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_kinetic failed')
    end if
    E_kinetic = 0.0_dp

    allocate (E_x(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_x failed')
    E_x = 0.0_dp

    allocate (E_y(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_y failed')
    E_y = 0.0_dp

    ! The refraction at the surface is set by the depth of the well the electron
    ! climbs out of, which is the inner potential, measured from the bottom of the
    ! free-electron-like final state band. The work function is measured from the
    ! Fermi level and is far smaller - about 4.3 against 13.5 eV for Cu - so using
    ! it bends the electron too little, leaves theta_internal too large, and makes
    ! the escape path 1/cos(theta) too long.
    if (photo_inner_potential_set) then
      surface_barrier = photo_inner_potential
    else
      surface_barrier = work_function_eff
      if (on_root .and. .not. barrier_warned) then
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        write (stdout, '(1x,a78)') '! Warning: photo_inner_potential is not set, so the refraction at the surface !'
        write (stdout, '(1x,a78)') '! falls back to the work function. That is the barrier measured from the      !'
        write (stdout, '(1x,a78)') '! Fermi level rather than from the bottom of the final state band, so it is   !'
        write (stdout, '(1x,a78)') '! far too small - roughly 4.3 against 13.5 eV for Cu - and escape depths come !'
        write (stdout, '(1x,a78)') '! out too short. workfct.py prints an estimate to feed the keyword.           !'
        write (stdout, '(1x,a78)') '!----------------------------------------------------------------------------!'
        barrier_warned = .true.
      end if
    end if

    if (index(photo_momentum, 'crystal') .gt. 0) call cell_calc_kpoint_r_cart

    if ((index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root) .or. &
        (index(devel_flag, 'print_qe_matrix_full') .gt. 0 .and. on_root) &
        .or. (index(devel_flag, 'print_qe_matrix_reduced') .gt. 0 .and. on_root)) then
      call cell_calc_kpoint_r_cart
      write (stdout, '(a78)') "+---------------- Printing K-Points in Cartesian Coordinates ----------------+"
      do N_k = 1, num_kpoints_on_node(my_node_id)
        write (stdout, '(3(1x,E22.15))') kpoint_r_cart(:, N_k)
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          do gdx = 1, photo_gkmax
            if (index(photo_momentum, 'crystal') .gt. 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(1, N_k)*1E+10)**2))*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(2, N_k)*1E+10)**2))*j_to_ev
            end if
            if (index(photo_momentum, 'gkgrid') .gt. 0) then
              E_x(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))* &
                                                ((photo_gkgrid(1, gdx, n_eigen, N_spin, N_k)*1E+10)**2))*j_to_ev
              E_y(gdx, n_eigen, N_spin, N_k) = (((hbar**2)/(2*e_mass))* &
                                                ((photo_gkgrid(2, gdx, n_eigen, N_spin, N_k)*1E+10)**2))*j_to_ev
            end if
            E_transverse(gdx, n_eigen, N_spin, N_k) = E_x(gdx, n_eigen, N_spin, N_k) + E_y(gdx, n_eigen, N_spin, N_k)

            ! The azimuth has to come from the signed transverse momentum, not
            ! from E_x and E_y: those are proportional to kx^2 and ky^2, so they
            ! have lost the signs that atan2 needs for the quadrant, and the
            ! angle they give is distorted as well. One irreducible k-point
            ! stands for its whole symmetry star, and the outputs that use this
            ! array have no azimuth axis to place the images on, so what is
            ! stored is the fraction of the star that lands in the wedge.
            if (index(photo_momentum, 'gkgrid') .gt. 0) then
              phi_accept_frac(gdx, n_eigen, N_spin, N_k) = &
                phi_star_fraction(photo_gkgrid(1, gdx, n_eigen, N_spin, N_k), &
                                  photo_gkgrid(2, gdx, n_eigen, N_spin, N_k))
            else
              phi_accept_frac(gdx, n_eigen, N_spin, N_k) = &
                phi_star_fraction(kpoint_r_cart(1, N_k), kpoint_r_cart(2, N_k))
            end if

            ! Emission angle theta is the angle between emitted
            ! electron vector and the surface normal.
            ! total kinetic energy after emission and passing through work
            ! function potential step
            E_kinetic(gdx, n_eigen, N_spin, N_k) = (band_energy(n_eigen, N_spin, N_k) &
                                                  & + temp_photon_energy - evacuum_eff)
            if (E_kinetic(gdx, n_eigen, N_spin, N_k) .lt. E_transverse(gdx, n_eigen, N_spin, N_k)) cycle
            ! Angle of electron outside material, after passing the surface and loosing E(work_function)
            ! acos(E_normal/E_kinetic)
            ! cos(theta) = k_z/|k|, a ratio of momenta. E_normal/E_kinetic is
            ! k_z^2/k^2, so it is cos^2(theta) and needs the square root before
            ! the acos. Both limits are exact either way, which is why leaving it
            ! out was invisible at normal and at grazing emission.
            theta_arpes(gdx, n_eigen, N_spin, N_k) = (acos(sqrt((E_kinetic(gdx, n_eigen, N_spin, N_k) &
                                                                 - E_transverse(gdx, n_eigen, N_spin, N_k)) &
                                                                /E_kinetic(gdx, n_eigen, N_spin, N_k))))*rad_to_deg
            ! Angle of electron within material, before passing the surface
            theta_internal(gdx, n_eigen, N_spin, N_k) = (acos(sqrt((E_kinetic(gdx, n_eigen, N_spin, N_k) + surface_barrier &
                                                                    - E_transverse(gdx, n_eigen, N_spin, N_k)) &
                                                                   /(E_kinetic(gdx, n_eigen, N_spin, N_k) + &
                                                                     surface_barrier))))*rad_to_deg
          end do ! Gkgrid
        end do ! bands
      end do ! spins
    end do ! k-points

    if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------ Printing Transverse Energy ------------------------+'
      write (stdout, '(3(1x,I4))') shape(E_transverse)
      write (stdout, '(3(1x,I4))') nbands, num_kpoints_on_node(my_node_id), nspins
      write (stdout, '(9999(es15.8))') ((((E_transverse(gdx, n_eigen, N_spin, N_k), gdx=1, photo_gkmax), &
                                          N_spin=1, nspins), N_k=1, num_kpoints_on_node(my_node_id)), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    deallocate (E_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_y')

    deallocate (E_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_x')

    if (allocated(band_gradient)) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate band_gradient')
    end if

    if (allocated(kpoint_r_cart)) then
      deallocate (kpoint_r_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate kpoint_r_cart')
    end if

    if (allocated(photo_gkgrid)) then
      deallocate (photo_gkgrid, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate photo_gkgrid')
    end if

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate Photoemission Angle', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_angle

  subroutine calc_electron_esc
    !*******=======================================================================
    ! This subroutine calculates the electron escape probability for each of the
    ! determined or supplied layers. The emission probability is calculated for the
    ! angle dependent propagation length versus the IMFP dependent scattering prob.
    ! orig. Victor Chang and Bruno Camino
    ! parts rewritten Felix Mildner, after Mar 2023
    !===============================================================================
    use od_constants, only: dp, deg_to_rad, bohr2ang, H2eV, pi
    use od_electronic, only: nbands, nspins, band_energy, efermi
    use od_cell, only: num_kpoints_on_node, atoms_pos_cart_photo, atoms_label_tmp, num_atoms
    use od_io, only: io_error, stdout, io_time
    use od_comms, only: my_node_id, on_root, comms_reduce
    use od_parameters, only: photo_imfp_value, photo_imfp_model, photo_model, iprint
    implicit none
    integer :: atom, N_k, N_spin, n_eigen, ierr, i, gdx
    real(kind=dp) :: tolerance, total_depth
    real(kind=dp) :: exponent, time0, time1, scale_factor, scaled_x, g1, g2
    real(kind=dp) :: band_imfp_min, band_imfp_max

    tolerance = 1.0E-12_dp
    time0 = io_time()
    ! One entry per atom in the cell, not per explicitly treated atom: the copy
    ! below is of the whole array, and the indexing that follows is by the
    ! original atom index through atom_order, which is bounded by num_atoms.
    allocate (new_atom_coordinates(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of new_atom_coordinates failed')

    !Redefine new z coordinates where the first layer is at z=0
    new_atom_coordinates = atoms_pos_cart_photo
    do atom = 1, max_atoms
      new_atom_coordinates(3, atom_order(atom)) = atoms_pos_cart_photo(3, atom_order(atom)) - &
                                                  (atoms_pos_cart_photo(3, atom_order(1)))
    end do

    if (.not. allocated(electron_esc)) then
      allocate (electron_esc(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    end if
    electron_esc = 0.0_dp

    if (.not. allocated(atom_imfp)) then
      allocate (atom_imfp(max_atoms), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_electron_esc_list - allocation of atom_imfp failed')
    end if
    atom_imfp = 0.0_dp
    if (index(photo_imfp_model, 'layers') .gt. 0) then
      if (on_root) then
        write (stdout, '(1x,a78)') '+--------------- User Supplied and Calculated IMFP Constants ----------------+'
        write (stdout, '(1x,a78)') '| Atom | Atom Order | Layer | Layer Thickness | User Input IMFP | Calc. IMFP |'
      end if

      ! Calculate the layer dependent imfp constant as a list for each layer
      do atom = 1, max_atoms
        total_depth = 0.0_dp
        do i = 1, box_atom(atom)
          atom_imfp(atom) = atom_imfp(atom) + box_heights(i)*photo_imfp_value(i)
          total_depth = total_depth + box_heights(i)
        end do
        atom_imfp(atom) = atom_imfp(atom)/total_depth
        if (on_root) then
          write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
            box_atom(atom), box_heights(box_atom(atom)), photo_imfp_value(box_atom(atom)), atom_imfp(atom), "    |"
225       format(1x, a1, a4, 6x, I3, 8x, I3, 6x, E14.6E3, 3x, F11.4, 3x, F11.4, a5)
        end if
      end do
      if (on_root) write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'

    else if (index(photo_imfp_model, 'const') .gt. 0) then
      atom_imfp = photo_imfp_value(1)

      ! This is a Cu specific IMFP curve by Nagy,Echenique - https://www.doi.org/10.1103/PhysRevB.85.115131
    else if (index(photo_imfp_model, 'cu_curve') .gt. 0) then
      if (.not. allocated(band_imfp)) then
        allocate (band_imfp(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_electron_esc_list - allocation of atom_imfp failed')
      end if
      band_imfp = 0.0_dp
      scale_factor = 6.9_dp
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            scaled_x = ((band_energy(n_eigen, N_spin, N_k) + temp_photon_energy - efermi) &
                        /scale_factor) + 1
            if ((1.0_dp - scaled_x) .gt. 1E-10_dp) cycle
            g1 = LOG(scaled_x - 1.0_dp) + ((8.0_dp/3.0_dp) - 2.0_dp*LOG(2.0_dp))
            ! g2 is a property of this state alone. It has to start from zero every
            ! time: accumulating onto whatever the previous state left behind made
            ! the IMFP depend on how many states preceded it, and hence on the band
            ! structure of the particular slab and on the MPI decomposition.
            g2 = 0.0_dp
            if (scaled_x .lt. 2.0_dp) then
              g2 = g2 + (2.0_dp/3.0_dp)*(SQRT(2.0_dp - scaled_x)**(3.0_dp))
              g2 = g2 + (2.0_dp*SQRT(2.0_dp - scaled_x))
              g2 = g2 + LOG(ABS((SQRT(2.0_dp - scaled_x) - 1.0_dp)/(SQRT(2.0_dp - scaled_x) + 1.0_dp)))
            end if
            band_imfp(n_eigen, N_spin, N_k) = bohr2ang*(4.0_dp*pi/3.0_dp)*(scaled_x/(g1 - g2))*(SQRT(2.0_dp*scale_factor/H2eV))
            ! write (stdout, *) "scaled_x", scaled_x, "g1", g1, "g2", g2, band_imfp(n_eigen, N_spin, N_k)
          end do
        end do
      end do
      ! Report the range over the states that actually have an IMFP. The zeros
      ! left behind by the cycle above are structural, not a short mean free
      ! path, and taking them into the minimum only ever printed 0.
      band_imfp_min = minval(band_imfp, mask=band_imfp .gt. 0.0_dp)
      band_imfp_max = maxval(band_imfp)
      call comms_reduce(band_imfp_min, 1, 'MIN')
      call comms_reduce(band_imfp_max, 1, 'MAX')
      if (on_root) then
        write (stdout, '(1x,a78)') '+------------------ IMFP Values from Energy Dependent Curve -----------------+'
        write (stdout, '(1x,a1,5x,a24,1x,a1,1x,E14.6E3,30x,a1)') '|', 'Min. IMFP over emitting ', '=', band_imfp_min, '|'
        write (stdout, '(1x,a1,5x,a24,1x,a1,1x,E14.6E3,30x,a1)') '|', 'Max. IMFP over emitting ', '=', band_imfp_max, '|'
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if
      ! set atom imfp to curve minimum (~2.5 Angstrom) to get
      ! "estimate value" during bulk slab printing
      atom_imfp = 2.50_dp
    end if

    if ((index(photo_imfp_model, 'const') .gt. 0) .or. (index(photo_imfp_model, 'layers') .gt. 0)) then
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              do gdx = 1, photo_gkmax
                ! is the emission possible?
                if (cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. tolerance) then
                  ! The electron's kinetic energy inside the material is higher, than after the emission
                  ! through the surface. It follows an angle closer to normal direction and one has to
                  ! use the internal theta angle here.
                  exponent = (new_atom_coordinates(3, atom_order(atom))/ &
                              cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/atom_imfp(atom)
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
    else if ((index(photo_imfp_model, 'cu_curve') .gt. 0)) then
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              ! The curve leaves band_imfp at zero for every state whose final
              ! energy is below E_F, because those states have no phase space to
              ! propagate through and cannot emit. Their theta_internal is also
              ! left at the 91 degree sentinel, so the cosine test below already
              ! rejects them, but that couples this loop to calc_angle several
              ! hundred lines away and leaves the division relying on the sign
              ! of an infinity. Skip them here instead.
              if (band_imfp(n_eigen, N_spin, N_k) .le. 0.0_dp) cycle
              do gdx = 1, photo_gkmax
                ! is the emission possible?
                if (cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. tolerance) then
                  ! The electron's kinetic energy inside the material is higher, than after the emission
                  ! through the surface. It follows an angle closer to normal direction and one has to
                  ! use the internal theta angle here.
                  ! For the curve a band dependent IMFP value is calculated, not a layer
                  ! dependent one.
                  exponent = (new_atom_coordinates(3, atom_order(atom))/ &
                              cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/band_imfp(n_eigen, N_spin, N_k)
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

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate Photoemission Escape', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_electron_esc

  subroutine bulk_emission
    !*******=======================================================================
    ! This subroutine calculates the product of the electron escape probability and
    ! the light intensity for the bulk approximation slab. The emission probability
    ! is calculated for the angle dependent propagation length versus the IMFP
    ! dependent scattering probability.
    ! orig. Victor Chang and Bruno Camino
    ! parts rewritten Felix Mildner, after Mar 2023
    !===============================================================================
    use od_constants, only: dp, deg_to_rad
    use od_electronic, only: nbands, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root, comms_reduce, comms_bcast
    use od_parameters, only: photo_imfp_value, photo_imfp_model, photo_bulk_cutoff, iprint
    use od_io, only: io_error, io_time, stdout
    implicit none
    real(kind=dp), dimension(:), allocatable :: bulk_light_tmp
    integer :: N_k, N_spin, n_eigen, i, num_layers, ierr, gdx
    real(kind=dp) :: exponent, time0, time1, band_imfp_max

    if (single_layer) then
      deallocate (new_atom_coordinates, stat=ierr)
      if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate new_atom_coordinates')
      return
    end if

    time0 = io_time()
    if (index(photo_imfp_model, 'layers') .gt. 0) then
      num_layers = ceiling((atom_imfp(max_atoms)*photo_bulk_cutoff)/bulk_repeat)
    else if (index(photo_imfp_model, 'const') .gt. 0) then
      num_layers = ceiling((photo_imfp_value(1)*photo_bulk_cutoff)/bulk_repeat)
    else if (index(photo_imfp_model, 'cu_curve') .gt. 0) then
      ! Calculate the emission probability for at most 1000 layers,
      ! since the propagation IMFPs for low energy electrons can be
      ! quite large in the curve case.
      band_imfp_max = maxval(band_imfp)
      call comms_reduce(band_imfp_max, 1, 'MAX')
      call comms_bcast(band_imfp_max, 1)
      num_layers = min(1000, ceiling((band_imfp_max*photo_bulk_cutoff)/bulk_repeat))
    end if

    allocate (bulk_light_tmp(num_layers), stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_light_tmp failed')
    bulk_light_tmp = 0.0_dp

    ! Every step below the explicit region is through another repeat of the same
    ! substrate, so one coefficient and one spacing -- bulk_repeat, not the height
    ! of the deepest box, which is only the distance left between the last
    ! boundary and the middle of the slab.
    bulk_light_tmp(1) = I_layer(box_atom(max_atoms), current_photo_energy_index)* &
                        exp(-(absorp_photo(box_atom(max_atoms), current_photo_energy_index)*bulk_repeat*1E-10))
    do i = 2, num_layers
      bulk_light_tmp(i) = bulk_light_tmp(i - 1)* &
                          exp(-(absorp_photo(box_atom(max_atoms), current_photo_energy_index)*bulk_repeat*1E-10))
    end do

    if ((index(photo_imfp_model, 'layers') .gt. 0) .or. (index(photo_imfp_model, 'const') .gt. 0)) then
      do i = 1, num_layers
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              do gdx = 1, photo_gkmax
                if (cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                  exponent = ((new_atom_coordinates(3, atom_order(max_atoms)) - i*bulk_repeat)/ &
                              cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/atom_imfp(max_atoms)
                  ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                  ! As exp(-230) is ~1E-100, this should be more than enough precision.
                  if (exponent .gt. -230.0_dp) then
                    electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    else if (index(photo_imfp_model, 'cu_curve') .gt. 0) then
      do i = 1, num_layers
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              ! As in calc_electron_esc: no IMFP means the state cannot emit.
              if (band_imfp(n_eigen, N_spin, N_k) .le. 0.0_dp) cycle
              do gdx = 1, photo_gkmax
                if (cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad) .gt. 0.0_dp) then
                  exponent = ((new_atom_coordinates(3, atom_order(max_atoms)) - i*bulk_repeat)/ &
                              cos(theta_internal(gdx, n_eigen, N_spin, N_k)*deg_to_rad))/band_imfp(n_eigen, N_spin, N_k)
                  ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                  ! As exp(-230) is ~1E-100, this should be more than enough precision.
                  if (exponent .gt. -230.0_dp) then
                    electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) = &
                      electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) + exp(exponent)*bulk_light_tmp(i)
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    if (on_root) then
      ! write out the bulk properties
      write (stdout, '(1x,a78)') '+---------------------- Bulk Approximation Slab Info ------------------------+'
      ! write out num_layers
      write (stdout, '(1x,a1,5x,a18,1x,a1,1x,I5,45x,a1)') '|', 'Number Bulk layers', '=', num_layers, '|'
      ! write out the total volume + volume per layer
      write (stdout, '(1x,a1,5x,a14,5x,a1,1x,F10.4,40x,a1)') '|', 'Vol. per layer', '=', box_volumes(num_boxes), '|'
      write (stdout, '(1x,a1,5x,a12,7x,a1,1x,F10.4,40x,a1)') '|', 'Total Volume', '=', num_layers*box_volumes(num_boxes), '|'
      write (stdout, '(1x,a78)') '+---- P_esc values for an electron with E = E_fermi and E_transverse = 0 ----+'
      ! write out bulk_light_tmp
      if (num_layers .lt. 6) then
        do i = 1, num_layers
          exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*bulk_repeat)/atom_imfp(max_atoms)
          ! This makes sure, that exp(exponent) does not underflow the dp fp value.
          ! As exp(-230) is ~1E-100, this should be more than enough precision.
          if (exponent .gt. -575.0_dp) then
            exponent = exp(exponent)
          else
            exponent = 0.0_dp
          end if
          write (stdout, 235) '|', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '|'
        end do
      else
        do i = 1, num_layers
          exponent = (new_atom_coordinates(3, atom_order(max_atoms)) - i*bulk_repeat)/atom_imfp(max_atoms)
          ! This makes sure, that exp(exponent) does not underflow the dp fp value.
          ! As exp(-230) is ~1E-100, this should be more than enough precision.
          if (exponent .gt. -230.0_dp) then
            exponent = exp(exponent)
          else
            exponent = 0.0_dp
          end if
          if (i .le. 3 .or. i .gt. num_layers - 3) then
            write (stdout, 235) '|', 'Layer # ', i, 'I_light = ', bulk_light_tmp(i), 'P_esc = ', exponent, '|'
          elseif (i .eq. 4) then
            write (stdout, '(1x,a1,35x,a6,35x,a1)') '|', '......', '|'
          end if
        end do
      end if ! If statement printing of slab light intensities formatting
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if ! If statement extra printing of slab data

    deallocate (bulk_light_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate bulk_light_tmp')

    deallocate (new_atom_coordinates, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate new_atom_coordinates')

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a38,21x,f11.3,a8)') '+ Time to calculate Bulk Photoemission', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

235 format(1x, a1, 5x, a8, I3, 5x, a10, E13.6E2, 2x, a8, E13.6E2, 9x, a1)
  end subroutine bulk_emission

  subroutine calc_ds_like_model
    !*===============================================================================
    ! This subroutine calculates the QE and MTE using a both the DOS dependent model
    ! and a simplified model following the Dowell-Schmerge like Model by Saha et al.
    ! The simplified model is an approximation, as the contributions are taken from
    ! individual bands, rather than the total DOS of the cell.
    ! Felix Mildner, after May 2024
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_temperature, iprint, num_exclude_bands, &
      exclude_bands, fixed, adaptive, linear
    use od_dos_utils, only: dos_adaptive, dos_fixed, dos_linear, dos_utils_calculate, dos_E => E
    use od_io, only: stdout, io_error, io_time
    use od_jdos_utils, only: jdos_energy_scale => setup_energy_scale
    use od_constants, only: kB
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :)    :: fermi_dirac
    real(kind=dp), allocatable, dimension(:)          :: fd
    real(kind=dp), allocatable, dimension(:)          :: dos_temp
    real(kind=dp) :: argument, time0, time1, final_fd, initial_fd, excess_energy, delta_e, diff, &
                     initial_dos, final_dos, temp_value
    integer :: N_k, N_spin, n_eigen, n_eigen_final, ierr, N_E, delta_index_photon, index_e

    time0 = io_time()

    if (.not. allocated(ds_qe_den)) then
      allocate (ds_qe_den(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of ds_qe_den failed')
      allocate (ds_qe_num(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of ds_qe_num failed')
      allocate (ds_mte_num(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of ds_mte_num failed')
    end if
    ds_qe_den = 0.0_dp
    ds_qe_num = 0.0_dp
    ds_mte_num = 0.0_dp
    ds_dos_mte_num = 0.0_dp
    ds_dos_mte_den = 0.0_dp

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    call photo_calculate_delta(delta_temp, .false.)

    if (iprint .gt. 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+-------------------------- Calculating DS Like QE --------------------------+'
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but exp(+-230) ~ 1E(+-100)
          ! so this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible under/over-flow.
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

    call jdos_energy_scale(E)
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a78)') '|       Simplified Dowell-Schmerge like model, after Saha et al.             |'
    end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
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
            ! The three sums share a common factor; form it once. The QE is
            ! ds_qe_num/ds_qe_den and the MTE is half of ds_mte_num/ds_qe_num, so
            ! the numerator of the first is the denominator of the second.
            temp_value = delta_temp(n_eigen, n_eigen_final, N_spin, N_k)* &
                         electrons_per_state*kpoint_weight(N_k)*final_fd*initial_fd
            ds_qe_den(n_eigen, n_eigen_final, N_spin, N_k) = temp_value
            ds_qe_num(n_eigen, n_eigen_final, N_spin, N_k) = temp_value*excess_energy
            ds_mte_num(n_eigen, n_eigen_final, N_spin, N_k) = temp_value*excess_energy**2
          end do
        end do
      end do
    end do
    call dos_utils_calculate()
    allocate (dos_temp(size(dos_E)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of dos_temp failed')
    ! One of these is always set by the broadening keyword, but say so rather
    ! than leaving dos_temp undefined if that ever stops being true.
    if (fixed) then
      dos_temp = sum(dos_fixed, dim=2)
    else if (adaptive) then
      dos_temp = sum(dos_adaptive, dim=2)
    else if (linear) then
      dos_temp = sum(dos_linear, dim=2)
    else
      call io_error('Error: calc_ds_like_model - no broadening scheme is active, so there is '// &
                    'no density of states to build the DOS based estimate from')
    end if

    ! Fermi occupation on the DOS energy scale, and the bin holding the initial
    ! state that a photon takes to the vacuum level.
    allocate (fd(size(dos_E)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_ds_like_model - allocation of fd failed')
    delta_e = dos_E(2) - dos_E(1)
    diff = 1.0E6_dp
    do N_e = 1, size(dos_E)
      ! evacuum_eff is work_function_eff + efermi, which is the same reference
      ! the band loop above uses.
      if (abs(dos_E(N_e) - evacuum_eff + temp_photon_energy) .lt. diff) then
        diff = abs(dos_E(N_e) - evacuum_eff + temp_photon_energy)
        index_e = N_e
      end if
      argument = (dos_E(N_e) - efermi)/(kB*photo_temperature)
      if (argument .gt. 230.0_dp) then
        fd(N_e) = 0.0_dp
      elseif (argument .lt. -230.0_dp) then
        fd(N_e) = 1.0_dp
      else
        fd(N_e) = 1.0_dp/(exp(argument) + 1.0_dp)
      end if
    end do
    delta_index_photon = int(temp_photon_energy/delta_e)
    initial_fd = fd(index_e)
    ! Walk up the initial-state energy until either the states are empty or the
    ! final state runs off the end of the grid.
    do while ((initial_fd .gt. 1.0E-50_dp) .and. ((index_e + delta_index_photon) .lt. (size(dos_E) - delta_index_photon - 2)))
      initial_fd = fd(index_e)
      final_fd = 1 - fd(index_e + delta_index_photon)
      initial_dos = dos_temp(index_e)
      final_dos = dos_temp(index_e + delta_index_photon)
      excess_energy = dos_E(index_e + delta_index_photon) - evacuum_eff
      temp_value = initial_dos*initial_fd*final_dos*final_fd*excess_energy
      ds_dos_mte_num = ds_dos_mte_num + temp_value*excess_energy
      ds_dos_mte_den = ds_dos_mte_den + temp_value
      index_e = index_e + 1
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
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a41,18x,f11.3,a8)') '+ Time to calculate DS like Photoemission', time1 - time0, ' (sec) +'
    end if

    deallocate (dos_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_ds_like_model - failed to deallocate dos_temp')
    deallocate (fd, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_ds_like_model - failed to deallocate fd')

  end subroutine calc_ds_like_model

  !===============================================================================
  subroutine calc_three_step_model
    !*===============================================================================
    ! This subroutine calculates the QE using the three step model.
    ! Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after 03/2023
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient, &
      elec_read_band_curvature, transmit_prob, elec_read_transmit_prob
    use od_comms, only: my_node_id, on_root, comms_send, comms_recv, comms_bcast
    use od_parameters, only: photo_temperature, devel_flag, photo_energy_sweep, iprint, &
      photo_output, photo_use_tmprob
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, io_time, io_date
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB, inv_sqrt_two_pi
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp) :: width, norm_vac, qe_factor, argument, efinal_temp, e_normal, &
                     time0, time1, final_fd, temp_contribution, gk_factor, te_gk_factor
    integer :: N_k, N_spin, n_eigen, n_eigen_init, n_eigen_final, atom, ierr, gdx

    width = kB*photo_temperature
    qe_factor = 1.0_dp/(cell_area)
    norm_vac = inv_sqrt_two_pi/width

    time0 = io_time()
    ! If no electric field is applied for the run, this has not yet been allocated,
    ! so we do that now and set it to 0 as not to affect the values.
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

    if (.not. allocated(emission_gauss)) then
      allocate (emission_gauss(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of emission_gauss failed')
    end if
    emission_gauss = 0.0_dp

    ! Do we use the classical approximation of the emission probability out of the
    ! surface using the transmission coefficient calculated based on PWs propagating
    ! out of the surface? See S. Huefner, Photoelectron Spectroscopy,
    ! Springer Berlin, Heidelberg, Third, 2003, Equation (6.11)
    if (photo_use_tmprob) then
      call elec_read_transmit_prob()
    else
      if (.not. allocated(transmit_prob)) then
        allocate (transmit_prob(nbands, num_kpoints_on_node(my_node_id), nspins))
        if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of transmit_prob failed')
      end if
      transmit_prob = 1.0_dp
    end if

    if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root .and. .not. photo_energy_sweep) then
      write (stdout, '(1x,a78)') '+----------------- Printing Matrix Weights in 3Step Function ----------------+'
      write (stdout, '(5(1x,I4))') shape(photo_matrix_weights)
      write (stdout, '(5(1x,I4))') nbands, nbands, nspins, num_kpoints_on_node(my_node_id), N_geom
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          write (stdout, '(99999(es15.8))') ((photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k), &
                                              n_eigen_final=1, nbands), n_eigen_init=1, nbands)
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (iprint .gt. 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 3Step QE ---------------------------+'
    end if

    ! TODO : Write an extra function, which writes all the values to an external file
    ! if (index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root .and. .not. photo_energy_sweep) &
    !   then
    !   i = 16 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
    !   write (stdout, '(1x,a78)') '+------------ Printing list of values going into 3step QE Values ------------+'
    !   write (stdout, '(16(1x,a17))') 'calced_qe_value', 'initial_state_energy', 'final_state_energy', 'spectral_func', &
    !     'photo_matrix_weights', &
    !     'delta_temp', 'electron_esc', 'kpoint_weight', 'I_layer', 'emisison_gauss', 'transverse_gauss', 'vacuum_gauss', &
    !     'fermi_dirac','final_fd', 'pdos_weights_atoms', 'pdos_weights_k_band'
    !   write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', max_atoms, nbands, nbands, nspins, num_kpoints_on_node(my_node_id), i
    ! end if

    ! Preparing the fermi dirac occupations and gaussian broadened emission
    ! Heaviside step function for the emission probability
    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but exp(+-230) ~ 1E(+-100)
          ! so this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible under/over-flow.
          if (argument .gt. 230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          ! Calculate the final state energy
          efinal_temp = band_energy(n_eigen, N_spin, N_k) + temp_photon_energy
          do gdx = 1, photo_gkmax
            ! is the energy along the normal .gt. 0?
            ! Include now the vacuum level and transverse energy to get the final energy along normal
            ! evacuum_eff = efermi + work_function_eff
            e_normal = efinal_temp - evacuum_eff - E_transverse(gdx, n_eigen, N_spin, N_k)
            if (e_normal .gt. 0.0_dp) then
              emission_gauss(gdx, n_eigen, N_spin, N_k) = 1.0_dp
            else
              emission_gauss(gdx, n_eigen, N_spin, N_k) = gaussian(e_normal, width, 0.0_dp)/norm_vac
            end if
          end do
        end do
      end do
    end do

    ! calculating the QE and transverse energy for all explicit layers
    call photo_calculate_delta(delta_temp, .false.)
    do atom = 1, max_atoms
      if (iprint .gt. 2 .and. on_root) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', &
          "Calculating atom ", atom, " of", max_atoms, ".lt.-- QE-3S |"
        call flush(stdout)
      end if
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
            final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
            do n_eigen_init = 1, n_eigen_final - 1
              ! do most of the calculation
              temp_contribution = (qe_factor*photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                                   *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k)*transmit_prob(n_eigen_final, N_k, N_spin) &
                                   *electrons_per_state*kpoint_weight(N_k)*(I_layer(box_atom(atom), current_photo_energy_index)) &
                                   *fermi_dirac(n_eigen_init, N_spin, N_k)*final_fd &
                                   *(pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                                     /pdos_weights_k_band(n_eigen_init, N_spin, N_k))) &
                                  *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
              do gdx = 1, photo_gkmax
                ! do the gkgrid_dependent part
                gk_factor = gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, atom) &
                            *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
                te_gk_factor = gk_factor*E_transverse(gdx, n_eigen_init, N_spin, N_k)
                qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) = qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, atom) &
                                                                         + temp_contribution*gk_factor
                te_tsm(n_eigen_init, N_spin, N_k, atom) = te_tsm(n_eigen_init, N_spin, N_k, atom) &
                                                          + temp_contribution*te_gk_factor

                ! if (index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root) then
                !   write (stdout, '(6(1x,I4))') gdx,n_eigen, n_eigen_final, N_spin, N_k, atom
                !   write (stdout, '(14(1x,E17.9E3))') qe_tsm(n_eigen, n_eigen_final, N_spin, N_k, atom),&
                !     band_energy(n_eigen, N_spin, N_k), band_energy(n_eigen_final, N_spin, N_k),&
                !     gkgrid_weight(gdx, n_eigen_init, N_spin, N_k), matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, 1), &
                !     delta_temp(n_eigen, n_eigen_final, N_spin, N_k), electron_esc(gdx, n_eigen_final, N_spin, N_k, atom), &
                !     kpoint_weight(N_k), I_layer(box_atom(atom), current_photo_energy_index), &
                !     emission_gauss(gdx, n_eigen_init, N_spin, N_k), fermi_dirac(n_eigen_init, N_spin, N_k), final_fd,&
                !     pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)), pdos_weights_k_band(n_eigen, N_spin, N_k)
                ! end if
              end do
            end do
          end do
        end do
      end do
    end do

    if (.not. single_layer) then
      ! Calculate the QE and Transverse Energy contributions from the bulk slab approximation
      call photo_calculate_delta(delta_temp, .true.)
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
            final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
            do n_eigen_init = 1, n_eigen_final - 1
              temp_contribution = &
                (qe_factor*photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                 *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                 *transmit_prob(n_eigen_final, N_k, N_spin) &
                 *electrons_per_state*kpoint_weight(N_k) &
                 *fermi_dirac(n_eigen_init, N_spin, N_k)*final_fd &
                 *(pdos_weights_boxes(n_eigen_init, N_spin, N_k, num_boxes) &
                   /pdos_weights_k_band(n_eigen_init, N_spin, N_k))) &
                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
              do gdx = 1, photo_gkmax
                gk_factor = gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                            *emission_gauss(gdx, n_eigen_init, N_spin, N_k) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, max_atoms + 1)
                te_gk_factor = gk_factor*E_transverse(gdx, n_eigen_init, N_spin, N_k)
                qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, max_atoms + 1) = &
                  qe_tsm(n_eigen_init, n_eigen_final, N_spin, N_k, max_atoms + 1) + temp_contribution*gk_factor
                te_tsm(n_eigen_init, N_spin, N_k, max_atoms + 1) = te_tsm(n_eigen_init, N_spin, N_k, max_atoms + 1) &
                                                                   + temp_contribution*te_gk_factor
              end do
            end do
          end do
        end do
      end do
    end if

    ! if (index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root) then
    !   write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    ! end if

    if (allocated(delta_temp)) then
      deallocate (delta_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate delta_temp')
    end if

    if (allocated(fermi_dirac)) then
      deallocate (fermi_dirac, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate fermi_dirac')
    end if

    if (allocated(emission_gauss)) then
      deallocate (emission_gauss, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate emission_gauss')
    end if

    if (allocated(transmit_prob) .and. index(photo_output, 'off') .gt. 0) then
      deallocate (transmit_prob, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate transmit_prob')
    end if

    if (allocated(gkgrid_weight) .and. index(photo_output, 'off') .gt. 0) then
      deallocate (gkgrid_weight, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate gkgrid_weight')
    end if

    if (index(devel_flag, 'print_kpt_qe_data') .gt. 0) call print_3step_kpt_qe

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 3step Photoemission', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_three_step_model

  subroutine print_3step_kpt_qe
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, comms_bcast
    use od_io, only: stdout, io_error, io_file_unit, seedname, io_date
    use od_parameters, only: photo_model
    implicit none
    real(kind=dp), allocatable, dimension(:) :: qe_k_temp
    integer :: N_k, ierr, qe_unit, token, inode
    character(len=10)                           :: char_e
    character(len=99)                           :: filename
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    if (on_root) then
      qe_unit = io_file_unit()
      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_k_point_QE.dat'
      write (stdout, *) 'opening file'
      open (unit=qe_unit, action='write', file=filename)
      write (qe_unit, *) '## The k point dependent QE values'
      call io_date(cdate, ctime)
      write (qe_unit, *) '## OptaDOS Photoemission: Printing QE K point Data on ', cdate, ' at ', ctime
    end if

    allocate (qe_k_temp(num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: print_3step_kpt_qe - failed to allocate qe_k_temp on root')
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
    if (ierr /= 0) call io_error('Error: print_3step_kpt_qe - failed to deallocate qe_k_temp')
  end subroutine print_3step_kpt_qe

  !===============================================================================
  subroutine photo_calculate_delta(delta_temp, calculate_bulk)
    !*===============================================================================
    ! Wrapper around the delta function subroutine to pass correct arguments
    ! Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after March 2023
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

    if (iprint .gt. 1 .and. on_root) then
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
    if (iprint .gt. 1 .and. on_root) then
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
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a34,25x,f11.3,1x,a7)') '+ Time to calculate Delta Function', time1 - time0, '(sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
  end subroutine photo_calculate_delta

  subroutine calculate_delta(delta_type, delta_temp, calculate_bulk)
    !*===============================================================================
    ! This subroutine evaluates the delta function between the initial band
    ! and the final band using the method specified in the input.
    ! The calculate_bulk paramater controls for the correct smearing step
    ! width calculation for either an explicit layer or a set of extrapolated
    ! bulk like layers.
    ! This is an adapted version, where all the bands are taken into account,
    ! which is more in line with how the 1-step model is calculated and also
    ! considers bands above the fermi energy with reduced occupation.
    ! orig. Victor Chang, 7 February 2020
    ! edited by Felix Mildner, after March 2022
    !===============================================================================
    use od_comms, only: my_node_id, on_root
    use od_cell, only: num_kpoints_on_node, kpoint_grid_dim, recip_lattice
    use od_parameters, only: adaptive_smearing, fixed_smearing, iprint, finite_bin_correction, &
      hybrid_linear_grad_tol, hybrid_linear, exclude_bands, &
      num_exclude_bands, jdos_max_energy, photo_slab_max, photo_slab_middle, &
      photo_slab_mode, SLAB_MODE_LAYERS
    use od_io, only: io_error, stdout
    use od_electronic, only: band_gradient, nbands, band_energy, nspins
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
    real(kind=dp) :: norm_width

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

    if (linear .or. adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:), dp)/2.0_dp
    if (adaptive .or. hybrid_linear) then
      do i = 1, 2
        sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
      end do
      if (calculate_bulk) then
        sub_cell_length(3) = sqrt(recip_lattice(3, 1)**2 + recip_lattice(3, 2)**2 + (pi/bulk_repeat)**2)*step(3)
      else
        sub_cell_length(3) = sqrt(recip_lattice(3, 1)**2 + recip_lattice(3, 2)**2 + (pi/slab_half_height)**2)*step(3)
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
      do is = 1, nspins
        do jb = 2, nbands
          if (num_exclude_bands .gt. 0) then
            if (any(exclude_bands == jb)) cycle
          end if
          do ib = 1, jb - 1
            if (linear .or. adaptive) grad(:) = band_gradient(jb, :, ik, is) - band_gradient(ib, :, ik, is)

            ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
            ! adaptive smearing (and take advantage of FBCS if required).
            force_adaptive = .false.
            if (.not. fixed) then
              if (hybrid_linear .and. (hybrid_linear_grad_tol .gt. sqrt(dot_product(grad, grad)))) force_adaptive = .true.
              if (linear .and. .not. force_adaptive) call doslin_sub_cell_corners(grad, step, band_energy(jb, is, ik) - &
                                                                                  band_energy(ib, is, ik), EV)
              if (adaptive .or. force_adaptive) width = sqrt(dot_product(grad, grad))*adaptive_smearing_temp
            end if
            ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
            ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
            if (finite_bin_correction) width = max(width, delta_bins)
            norm_width = inv_sqrt_two_pi/width

            ! The linear method has a special way to calculate the integrated dos
            ! we have to take account for this here.
            if (linear .and. .not. force_adaptive) then
              delta_temp(ib, jb, is, ik) = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), E(current_energy_index), cuml)
            else
              delta_temp(ib, jb, is, ik) = gaussian(band_energy(jb, is, ik) - band_energy(ib, is, ik), width, &
                                                    E(current_energy_index))
            end if
          end do
        end do
      end do
    end do

    if (iprint .gt. 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine calculate_delta

  !===============================================================================
  subroutine make_foptical_weights
    !*===============================================================================
    ! Contract the free-electron coherency tensor with the light polarisation to
    ! give the one-step matrix element weight for each initial state.
    !
    ! The tensor is stored against final-state energy and carries no photon
    ! energy, so the weight for photon energy hw comes from evaluating it at
    ! E_f = E_n + hw and contracting
    !
    !     |q.M|^2 summed over G  =  q_i conjg(q_j) A_ij
    !
    ! The polarisation directions here are real, so only the three diagonal
    ! elements and the real parts of the off-diagonals contribute; the imaginary
    ! parts are stored for a circular polarisation to use later.
    ! orig. Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after April 2024
    !===============================================================================
    use od_constants, only: dp, kB
    use od_electronic, only: nbands, nspins, num_electrons, electrons_per_state, fem_tensor, &
      & fem_energy_info, band_energy, efermi
    use od_cell, only: num_kpoints_on_node, cell_get_symmetry, num_crystal_symmetry_operations, &
      & crystal_symmetry_operations
    use od_parameters, only: optics_geom, optics_qdir, legacy_file_format, devel_flag, photo_energy_sweep, &
      & iprint, photo_temperature
    use od_io, only: io_error, stdout
    use od_comms, only: my_node_id, on_root, comms_reduce

    implicit none

    real(kind=dp), dimension(9) :: tens
    real(kind=dp), dimension(3) :: qdir, qdir1, qdir2
    real(kind=dp), dimension(2) :: num_occ
    real(kind=dp) :: q_weight, q_weight1, q_weight2, factor, e_final, ef_top, e_occ_max
    integer :: N_k, i, j, N_in, N_spin, N2, N3, n_eigen, num_symm, ierr
    integer :: n_ef, n_above

    if (.not. legacy_file_format .and. index(devel_flag, 'old_filename') .gt. 0) then
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

    ! fem_energy_info: n_Ef, Ef_min, Ef_step, Ef_broadening, Ef_origin
    n_ef = nint(fem_energy_info(1))
    ef_top = fem_energy_info(2) + fem_energy_info(3)*real(n_ef, dp)

    ! No compatibility checks are needed against the photon grid, the Fermi
    ! energy or the work function: the tensor holds none of them.  The only way
    ! to get this wrong is to ask for a final-state energy the file does not
    ! cover, which is counted below and raised as an error rather than silently
    ! dropping the states that would have been emitted.
    !
    ! Only occupied initial states are worth checking.  A band well above the
    ! Fermi level carries no electrons to emit, so its final-state energy falling
    ! outside the window costs nothing -- and with the extra bands a spectral run
    ! carries, most of them do.
    e_occ_max = efermi + 10.0_dp*kB*photo_temperature

    if (.not. allocated(foptical_matrix_weights)) then
      allocate (foptical_matrix_weights(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: make_foptical_weights - allocation of foptical_matrix_weights failed')
    end if
    foptical_matrix_weights = 0.0_dp

    if (index(optics_geom, 'polar') .gt. 0) then
      qdir = optics_qdir
      q_weight = ((qdir(1)**2) + (qdir(2)**2) + (qdir(3)**2))**0.5_dp
      if (q_weight .lt. 0.001_dp) call io_error("Error:  please check optics_qdir, norm close to zero")
    end if

    if (index(optics_geom, 'unpolar') .gt. 0) then
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
    factor = 1.0_dp/(temp_photon_energy**2)
    n_above = 0

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands                                ! Loop over state
          ! * The tensor is binned by |k+G|^2/2, the plane-wave energy, whose zero
          !   is the cell-averaged potential -- and that average moves when vacuum
          !   is added, by 2.7 eV over a 10-28 A scan on an otherwise identical
          !   slab.  A photoelectron travels in the vacuum, not in the average
          !   potential, so the final state to look up is the plane wave whose
          !   kinetic energy above the vacuum level matches the excess energy.
          !   Referencing the lookup to evacuum_eff makes it cell independent: the
          !   spread over the converged cells falls from 1.264 eV to 0.070 eV, and
          !   what is left is the work-function spread, which is real.
          e_final = band_energy(n_eigen, N_spin, N_k) + temp_photon_energy - evacuum_eff
          if (e_final .gt. ef_top .and. band_energy(n_eigen, N_spin, N_k) .le. e_occ_max) &
            n_above = n_above + 1
          call fem_tensor_at(n_eigen, N_spin, N_k, e_final, tens)

          if (index(optics_geom, 'unpolar') .gt. 0) then
            if (num_symm == 0) then
              foptical_matrix_weights(n_eigen, N_spin, N_k) = 0.5_dp*factor* &
                                                              (fem_contract(tens, qdir1, q_weight1) &
                                                               + fem_contract(tens, qdir2, q_weight2))
            else ! begin unpolar symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir1(j))
                    end do
                  end do
                  foptical_matrix_weights(n_eigen, N_spin, N_k) = &
                    foptical_matrix_weights(n_eigen, N_spin, N_k) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*factor*fem_contract(tens, qdir, q_weight1)
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir2(j))
                    end do
                  end do
                  foptical_matrix_weights(n_eigen, N_spin, N_k) = &
                    foptical_matrix_weights(n_eigen, N_spin, N_k) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*factor*fem_contract(tens, qdir, q_weight2)
                end do
              end do
            end if !end unpolar symmetric

          elseif (index(optics_geom, 'polar') .gt. 0) then
            if (num_symm == 0) then
              foptical_matrix_weights(n_eigen, N_spin, N_k) = factor*fem_contract(tens, qdir, q_weight)
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
                  foptical_matrix_weights(n_eigen, N_spin, N_k) = &
                    foptical_matrix_weights(n_eigen, N_spin, N_k) + &
                    (1.0_dp/Real((num_symm*(N_in + 1)), dp))*factor*fem_contract(tens, qdir, q_weight)
                end do
              end do
            end if ! end polar symmetric

          elseif (index(optics_geom, 'poly') .gt. 0) then
            ! Polycrystalline: no preferred light direction.  Averaging q_i q_j A_ij
            ! over the three cartesian directions is just the trace over three.
            if (num_symm == 0) then
              foptical_matrix_weights(n_eigen, N_spin, N_k) = (factor/3.0_dp)* &
                                                              (tens(1) + tens(2) + tens(3))
            else ! begin poly symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  do i = 1, 3
                    qdir(i) = ((-1.0_dp)**(N3 + 1))*crystal_symmetry_operations(1, i, N2)
                    qdir1(i) = ((-1.0_dp)**(N3 + 1))*crystal_symmetry_operations(2, i, N2)
                    qdir2(i) = ((-1.0_dp)**(N3 + 1))*crystal_symmetry_operations(3, i, N2)
                  end do
                  foptical_matrix_weights(n_eigen, N_spin, N_k) = &
                    foptical_matrix_weights(n_eigen, N_spin, N_k) + &
                    (1.0_dp/Real((num_symm*(N_in + 1)), dp))*factor* &
                    ((fem_contract(tens, qdir, 1.0_dp) + fem_contract(tens, qdir1, 1.0_dp) &
                      + fem_contract(tens, qdir2, 1.0_dp))/3.0_dp)
                end do
              end do
            end if ! end poly symmetric
          end if ! end photo_geom
        end do ! loop over state 1
      end do ! loop over spins
    end do ! loop over kpoints

    call comms_reduce(n_above, 1, 'SUM')
    if (on_root .and. n_above .gt. 0) then
      write (stdout, *) 'photon energy:', temp_photon_energy, ' top of the stored Ef window:', ef_top
      write (stdout, *) n_above, ' occupied (band, k, spin) states need a final-state energy above'
      write (stdout, *) 'the window stored in the .fem_bin.  Those states would be silently dropped.'
      write (stdout, *) 'Raise SPECTRAL_FEM_EF_MAX and regenerate, or lower the photon energy.'
      call flush(stdout)
      call io_error('Error: photon energy takes states above the stored final-state window')
    end if

    if (index(devel_flag, 'print_qe_constituents') .gt. 0 .and. on_root .and. .not. photo_energy_sweep) then
      write (stdout, '(1x,a78)') '+------------------------- Printing Free OM Weights -------------------------+'
      write (stdout, 126) shape(foptical_matrix_weights)
126   format(5(1x, I4))
      do N_spin = 1, nspins
        do N_k = 1, num_kpoints_on_node(my_node_id)
          write (stdout, '(99999(es15.8))') (foptical_matrix_weights(n_eigen, N_spin, N_k), n_eigen=1, nbands)
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

  end subroutine make_foptical_weights

  !===============================================================================
  subroutine fem_tensor_at(n_eigen, N_spin, N_k, e_final, tens)
    !===============================================================================
    ! The nine stored reals of A_ij at final-state energy e_final, linearly
    ! interpolated between the two neighbouring bins.  Bin ie is centred on
    ! Ef_min + (ie - 1/2)*Ef_step.  Energies outside the stored window return
    ! zero: below it there are no free-electron final states at all, and above
    ! it the caller counts the occurrence and raises an error.
    !===============================================================================
    use od_constants, only: dp
    use od_electronic, only: fem_tensor, fem_energy_info
    implicit none
    integer, intent(in) :: n_eigen, N_spin, N_k
    real(kind=dp), intent(in) :: e_final
    real(kind=dp), dimension(9), intent(out) :: tens
    real(kind=dp) :: x, f
    integer :: ie, n_ef

    n_ef = nint(fem_energy_info(1))
    x = (e_final - fem_energy_info(2))/fem_energy_info(3) + 0.5_dp
    ie = floor(x)
    if (ie .lt. 1 .or. ie .ge. n_ef) then
      tens = 0.0_dp
      return
    end if
    f = x - real(ie, dp)
    tens = (1.0_dp - f)*fem_tensor(:, ie, n_eigen, N_k, N_spin) &
           + f*fem_tensor(:, ie + 1, n_eigen, N_k, N_spin)

  end subroutine fem_tensor_at

  !===============================================================================
  pure function fem_contract(tens, q, qnorm) result(c)
    !===============================================================================
    ! q_i q_j A_ij for a real polarisation direction q, from the nine stored
    ! reals.  The imaginary parts of the off-diagonals cancel in pairs when q is
    ! real, so only six of the nine are used here.
    !===============================================================================
    use od_constants, only: dp
    implicit none
    real(kind=dp), dimension(9), intent(in) :: tens
    real(kind=dp), dimension(3), intent(in) :: q
    real(kind=dp), intent(in) :: qnorm
    real(kind=dp) :: c

    c = (q(1)*q(1)*tens(1) + q(2)*q(2)*tens(2) + q(3)*q(3)*tens(3) &
         + 2.0_dp*(q(1)*q(2)*tens(4) + q(1)*q(3)*tens(6) + q(2)*q(3)*tens(8)))/(qnorm*qnorm)

  end function fem_contract

  !===============================================================================
  subroutine calc_one_step_model
    !===============================================================================
    ! This subroutine calculates the QE using a one step model.
    ! orig. Victor Chang, 7th February 2020
    ! edited by Felix Mildner, after April 2024
    !===============================================================================

    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient,&
    & elec_read_band_curvature
    use od_comms, only: my_node_id, num_nodes
    use od_parameters, only: photo_temperature, devel_flag, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_comms, only: on_root, comms_recv, comms_send, comms_reduce
    use od_io, only: stdout, io_error, io_file_unit, io_time, io_date
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB, inv_sqrt_two_pi
    implicit none
    integer :: N_k, N_spin, n_eigen, atom, ierr, gdx, kpt_total

    real(kind=dp) :: width, norm_vac, qe_factor, argument, time0, time1
    real(kind=dp) :: temp_contribution, efinal_temp, e_normal
    real(kind=dp) :: gk_factor, te_gk_factor
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    ! per-band breakdown of the QE integrand, for devel_flag 'print_1step_terms'
    real(kind=dp), allocatable, dimension(:) :: dbg_m, dbg_fd, dbg_pdos, dbg_gk, dbg_qe

    if (index(devel_flag, 'write_fem_matrix') .gt. 0) then
      kpt_total = sum(num_kpoints_on_node(0:num_nodes - 1))
      call write_distributed_fem_data(kpt_total)
    end if

    qe_factor = 1.0_dp/(cell_area)
    width = (1.0_dp/11604.45_dp)*photo_temperature
    norm_vac = inv_sqrt_two_pi/width

    time0 = io_time()
    if (iprint .gt. 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 1Step QE ---------------------------+'
    end if

    ! If no electric field is applied for the run, this has not yet been allocated,
    ! so we do that now and set it to 0 as not to affect the values.
    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of field_emission failed')
      field_emission = 0.0_dp
    end if

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    if (.not. allocated(emission_gauss)) then
      allocate (emission_gauss(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of fermi_dirac failed')
    end if
    emission_gauss = 0.0_dp

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

    ! TODO: Create extra printing function for both photo models
    ! if (index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root .and. .not. photo_energy_sweep) then
    !   i = 14 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
    !   write (stdout, '(1x,a78)') '+------------ Printing list of values going into 1step QE Values ------------+'
    !   write (stdout, '(14(7x,a17))') 'calced_qe_value', 'contribution', 'band_energy', 'gkgrid_weight', &
    !    'foptical_matrix_weights', &
    !    'electron_esc', 'kpoint_weight', 'I_layer', 'emission_gauss', 'transverse_gauss', 'vacuum_gauss', 'fermi_dirac', &
    !    'pdos_weights_atoms', 'pdos_weights_k_band'
    !   write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', i, max_atoms, nbands, nspins, num_kpoints_on_node(my_node_id)
    ! end if

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but exp(+-230) ~ 1E(+-100)
          ! so this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible under/over-flow.
          if (argument .gt. 230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          ! the final total energy of the electron (E_initial + hw) - e_vacuum
          efinal_temp = band_energy(n_eigen, N_spin, N_k) + temp_photon_energy - evacuum_eff

          ! is the photon energy large enough to allow an emission at this kpoint/k+G
          do gdx = 1, photo_gkmax
            ! Include now the vacuum level and transverse energy to get the final energy along normal
            e_normal = efinal_temp - E_transverse(gdx, n_eigen, N_spin, N_k)
            if (e_normal .gt. 0.0_dp) then
              emission_gauss(gdx, n_eigen, N_spin, N_k) = 1.0_dp
            else
              emission_gauss(gdx, n_eigen, N_spin, N_k) = gaussian(e_normal, width, 0.0_dp)/norm_vac
            end if
          end do
        end do
      end do
    end do

    ! Optional per-band breakdown of the surface layer's QE integrand: one row
    ! per band, so it stays small, unlike a dump of every term at every gk/k.
    if (index(devel_flag, 'print_1step_terms') .gt. 0) then
      allocate (dbg_m(nbands), dbg_fd(nbands), dbg_pdos(nbands), dbg_gk(nbands), dbg_qe(nbands), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of the 1step term breakdown failed')
      dbg_m = 0.0_dp; dbg_fd = 0.0_dp; dbg_pdos = 0.0_dp; dbg_gk = 0.0_dp; dbg_qe = 0.0_dp
    end if

    do atom = 1, max_atoms
      if (iprint .gt. 2 .and. on_root .and. (atom .le. max_atoms)) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', "Calculating atom ", atom, " of", max_atoms, ".lt.-- QE-1S |"
      end if
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            temp_contribution = (qe_factor &
                                 *foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                 *electrons_per_state*kpoint_weight(N_k) &
                                 *(I_layer(box_atom(atom), current_photo_energy_index)) &
                                 *fermi_dirac(n_eigen, N_spin, N_k) &
                                 *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                                   /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                                *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
            if (allocated(dbg_m) .and. atom .eq. 1) then
              dbg_m(n_eigen) = dbg_m(n_eigen) + kpoint_weight(N_k)*foptical_matrix_weights(n_eigen, N_spin, N_k)
              dbg_fd(n_eigen) = dbg_fd(n_eigen) + kpoint_weight(N_k)*fermi_dirac(n_eigen, N_spin, N_k)
              dbg_pdos(n_eigen) = dbg_pdos(n_eigen) + kpoint_weight(N_k) &
                                  *pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                                  /pdos_weights_k_band(n_eigen, N_spin, N_k)
            end if
            do gdx = 1, photo_gkmax
              gk_factor = gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                          *electron_esc(gdx, n_eigen, N_spin, N_k, atom) &
                          *emission_gauss(gdx, n_eigen, N_spin, N_k)
              te_gk_factor = gk_factor*E_transverse(gdx, n_eigen, N_spin, N_k)
              if (allocated(dbg_gk) .and. atom .eq. 1) then
                dbg_gk(n_eigen) = dbg_gk(n_eigen) + kpoint_weight(N_k)*gk_factor
                dbg_qe(n_eigen) = dbg_qe(n_eigen) + temp_contribution*gk_factor
              end if
              qe_osm(n_eigen, N_spin, N_k, atom) = qe_osm(n_eigen, N_spin, N_k, atom) &
                                                   + temp_contribution*gk_factor
              te_osm(n_eigen, N_spin, N_k, atom) = te_osm(n_eigen, N_spin, N_k, atom) &
                                                   + temp_contribution*te_gk_factor
              ! if ((temp_contribution*gk_factor) .gt. 0.0_dp .and. index(devel_flag, 'print_qe_formula_values') .gt. 0 &
              !     .and. on_root) then
              !   write (stdout, '(5(1x,I4))') gdx, n_eigen, N_spin, N_k, atom
              !   write (stdout, '(14(7x,E17.9E3))') qe_osm(n_eigen, N_spin, N_k, atom), temp_contribution*gk_factor, &
              !     band_energy(n_eigen, N_spin, N_k), &
              !     gkgrid_weight(gdx, n_eigen, N_spin, N_k), foptical_matrix_weights(n_eigen, N_spin, N_k), &
              !     electron_esc(gdx, n_eigen, N_spin, N_k, atom), kpoint_weight(N_k), &
              !     I_layer(box_atom(atom), current_photo_energy_index), emission_gauss(gdx, n_eigen, N_spin, N_k), &
              !     transverse_gauss(gdx, n_eigen, N_spin, N_k), vacuum_gauss(n_eigen, N_spin, N_k) &
              !     fermi_dirac(n_eigen, N_spin, N_k), pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)), &
              !     pdos_weights_k_band(n_eigen, N_spin, N_k)
              ! end if
            end do
          end do
        end do
      end do
    end do

    if (.not. single_layer) then
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            temp_contribution = (qe_factor &
                                 *foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                 *electrons_per_state*kpoint_weight(N_k) &
                                 *(I_layer(box_atom(max_atoms + 1), current_photo_energy_index)) &
                                 *fermi_dirac(n_eigen, N_spin, N_k) &
                                 *(pdos_weights_boxes(n_eigen, N_spin, N_k, num_boxes) &
                                   /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                                *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
            do gdx = 1, photo_gkmax
              gk_factor = gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                          *electron_esc(gdx, n_eigen, N_spin, N_k, max_atoms + 1) &
                          *emission_gauss(gdx, n_eigen, N_spin, N_k)
              te_gk_factor = gk_factor*E_transverse(gdx, n_eigen, N_spin, N_k)
              qe_osm(n_eigen, N_spin, N_k, max_atoms + 1) = qe_osm(n_eigen, N_spin, N_k, max_atoms + 1) &
                                                            + temp_contribution*gk_factor
              te_osm(n_eigen, N_spin, N_k, max_atoms + 1) = te_osm(n_eigen, N_spin, N_k, max_atoms + 1) &
                                                            + temp_contribution*te_gk_factor
              ! if ((temp_contribution*gk_factor) .gt. 0.0_dp .and. index(devel_flag, 'print_qe_formula_values') .gt. 0 &
              !     .and. on_root) then
              !   write (stdout, '(5(1x,I4))') gdx, n_eigen, N_spin, N_k, atom
              !   write (stdout, '(14(7x,E17.9E3))') qe_osm(n_eigen, N_spin, N_k, atom), temp_contribution*gk_factor, &
              !     band_energy(n_eigen, N_spin, N_k), &
              !     gkgrid_weight(gdx, n_eigen, N_spin, N_k), foptical_matrix_weights(n_eigen, N_spin, N_k), &
              !     electron_esc(gdx, n_eigen, N_spin, N_k, atom), kpoint_weight(N_k), &
              !     I_layer(box_atom(atom), current_photo_energy_index), emission_gauss(gdx, n_eigen, N_spin, N_k), &
              !     transverse_gauss(gdx, n_eigen, N_spin, N_k), vacuum_gauss(n_eigen, N_spin, N_k) &
              !     fermi_dirac(n_eigen, N_spin, N_k), pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)), &
              !     pdos_weights_k_band(n_eigen, N_spin, N_k)
              ! end if
            end do
          end do
        end do
      end do
    end if
    ! if (index(devel_flag, 'print_qe_formula_values') .gt. 0 .and. on_root) then
    !   write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    ! end if

    if (allocated(dbg_m)) then
      call comms_reduce(dbg_m(1), nbands, 'SUM')
      call comms_reduce(dbg_fd(1), nbands, 'SUM')
      call comms_reduce(dbg_pdos(1), nbands, 'SUM')
      call comms_reduce(dbg_gk(1), nbands, 'SUM')
      call comms_reduce(dbg_qe(1), nbands, 'SUM')
      if (on_root) then
        write (stdout, '(1x,a78)') '+------------- 1step QE integrand, surface layer, per band ------------------+'
        write (stdout, '(1x,a6,5(1x,a17))') 'band', 'foptical', 'fermi_dirac', 'pdos_ratio', 'sum_gk', 'contribution'
        do n_eigen = 1, nbands
          write (stdout, '(1x,i6,5(1x,E17.9E3))') n_eigen, dbg_m(n_eigen), dbg_fd(n_eigen), &
            dbg_pdos(n_eigen), dbg_gk(n_eigen), dbg_qe(n_eigen)
        end do
        write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
      end if
      deallocate (dbg_m, dbg_fd, dbg_pdos, dbg_gk, dbg_qe, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - failed to deallocate the 1step term breakdown')
    end if

    if (index(devel_flag, 'print_kpt_qe_data') .gt. 0) call print_1step_kpt_qe

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 1step Photoemission', time1 - time0, ' (sec) +'
    end if
  end subroutine calc_one_step_model

  subroutine print_1step_kpt_qe
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, comms_bcast
    use od_io, only: stdout, io_error, io_file_unit, seedname, io_date
    use od_parameters, only: photo_model
    implicit none
    real(kind=dp), allocatable, dimension(:) :: qe_k_temp
    integer :: N_k, ierr, qe_unit, token, inode
    character(len=10)                           :: char_e
    character(len=99)                           :: filename
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    if (on_root) then
      qe_unit = io_file_unit()
      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_k_point_QE.dat'
      write (stdout, *) 'opening file'
      open (unit=qe_unit, action='write', file=filename)
      write (qe_unit, *) '## The k point dependent QE values'
      call io_date(cdate, ctime)
      write (qe_unit, *) '## OptaDOS Photoemission: Printing QE K point Data on ', cdate, ' at ', ctime
    end if

    allocate (qe_k_temp(num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: print_1step_kpt_qe - failed to allocate qe_k_temp on root')
    token = -1

    ! allocate and the 1step qe matrix on non-root
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
      ! When all the files are read, the order is nodes#1 -> nodes#max -> root
      ! The root gets all the extra kpoints, that do not fit neatly into the
      ! integer division determined step size.
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
    if (ierr /= 0) call io_error('Error: print_1step_kpt_qe - failed to deallocate qe_k_temp')
  end subroutine print_1step_kpt_qe

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
    real(kind=dp)                            :: time0, time1, qe_term1, qe_term2, mte_term1, mte_term2
    integer                                  :: atom, ierr

    time0 = io_time()

    if (iprint .gt. 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Calculating MTE ------------------------------+'
    end if

    if (.not. allocated(layer_qe)) then
      allocate (layer_qe(max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of layer_qe failed')
    end if
    layer_qe = 0.0_dp

    if (index(photo_model, '3step') .gt. 0) then
      do atom = 1, max_atoms + 1
        ! Calculate the qe contribution of each atom/layer
        layer_qe(atom) = sum(qe_tsm(:, :, :, :, atom))
      end do

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

    elseif (index(photo_model, '1step') .gt. 0) then
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

    else if (index(photo_model, 'dosds') .gt. 0) then

      qe_term1 = sum(ds_qe_num)
      call comms_reduce(qe_term1, 1, 'SUM')
      qe_term2 = sum(ds_qe_den)
      call comms_reduce(qe_term2, 1, 'SUM')
      mte_term1 = sum(ds_mte_num)
      call comms_reduce(mte_term1, 1, 'SUM')
      total_qe = qe_term1/qe_term2
      ! The numerator of the QE is the denominator of the MTE.
      mean_te = 0.5_dp*(mte_term1/qe_term1)

    end if

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a23,36x,f11.3,a8)') '+ Time to calculate MTE', time1 - time0, ' (sec) +'
    end if

  end subroutine weighted_mean_te

  subroutine write_qe_results
    !*===============================================================================
    ! This subroutine writes the calculated Photoemission data to the output file.
    ! The contents of this routine used to be part of the subroutine weighted_mean_te, but were moved
    ! here to make the subroutine names more representative of their function.
    !===============================================================================
    use od_cell, only: cell_calc_kpoint_r_cart, atoms_label_tmp
    use od_parameters, only: photo_work_function, photo_elec_field, photo_model, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout
    use od_jdos_utils, only: jdos_utils_calculate
    integer :: atom

    write (stdout, '(1x,a78)') '+------------------------------ Photoemission -------------------------------+'
    write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    write (stdout, 223) '| Work Function     ', photo_work_function, &
      'eV      Photon Energy   ', temp_photon_energy, 'eV   |'
    if (((photo_elec_field .lt. 1.0e3_dp) .and. (photo_elec_field .gt. 1.0e-3_dp)) .or. (photo_elec_field .eq. 0.0_dp)) then
      write (stdout, 224) '| Effective Work Function', work_function_eff, &
        'eV      Electric Field   ', photo_elec_field, 'V/m  |'
    else
      write (stdout, 235) '| Effective Work Function', work_function_eff, &
        'eV      Electric Field   ', photo_elec_field, 'V/m  |'
    end if

    if (index(photo_model, '3step') .gt. 0) then
      write (stdout, '(1x,a78)') '| Final State : Bloch State                                                  |'
    elseif (index(photo_model, '1step') .gt. 0) then
      write (stdout, '(1x,a78)') '| Final State : Free Electron State                                          |'
    end if
    write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    if (index(photo_model, 'dosds') .gt. 0) then
      write (stdout, '(1x,a78)') '|       **********  Results from DOS/Band dep. DS PE model  **********       |'
      write (stdout, '(1x,a78)') '|       **********       Band based estimate values         **********       |'
      write (stdout, 236) '|       QE from single band contributions   :', total_qe, '   |'

      write (stdout, 236) '|       MTE from single band contrib.  (eV) :', mean_te, '      |'
      write (stdout, '(1x,a78)') '|       **********        DOS based MTE estimate            **********       |'
      write (stdout, 236) '|       MTE estimate from DOS          (eV) :', &
        0.5_dp*ds_dos_mte_num/ds_dos_mte_den, '   |'
    else
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Layer   |             Quantum Efficiency           |'
      ! Larger number of digits for debugging purposes
      if (iprint .gt. 2) then
        do atom = 1, max_atoms
          write (stdout, 231) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
            box_atom(atom), layer_qe(atom), "      |"
        end do
        write (stdout, 232) "| Bulk", layer_qe(max_atoms + 1), &
        &"      |"

        write (stdout, 233) '| Total Quantum Efficiency (electrons/photon):', total_qe, '   |'

        write (stdout, 234) '| Weighted Mean Transverse Energy (eV):', mean_te, '      |'
      else
        do atom = 1, max_atoms
          write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
            box_atom(atom), layer_qe(atom), "      |"
        end do
        write (stdout, 226) "| Bulk", layer_qe(max_atoms + 1), &
        &"      |"

        write (stdout, 227) '| Total Quantum Efficiency (electrons/photon):', total_qe, '   |'

        write (stdout, 228) '| Weighted Mean Transverse Energy (eV):', mean_te, '      |'
      end if
    end if

    if (photo_elec_field .gt. 0.0_dp) then
      ! Larger number of digits for debugging purposes
      if (iprint .gt. 2) then
        write (stdout, 234) '| Total field emission (electrons/A^2):', total_field_emission, '      |'
      else
        write (stdout, 228) '| Total field emission (electrons/A^2):', total_field_emission, '      |'
      end if
    end if

    write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    call flush(stdout)
223 format(1x, a20, f15.4, 1x, a24, f11.4, a7)
224 format(1x, a25, f10.4, 1x, a25, f10.4, a7)
235 format(1x, a25, f10.4, 1x, a25, E10.4, a7)
225 format(1x, a1, a4, 8x, I3, 10x, I3, 16x, E17.4E3, 9x, a7)
226 format(1x, a6, 38x, E18.4E3, 9x, a7)
227 format(1x, a46, E20.4E3, 5x, a7)
236 format(1x, a45, E20.4E3, 6x, a7)
228 format(1x, a39, 7x, E20.4E3, 5x, a7)
231 format(1x, a1, a4, 8x, I3, 10x, I3, 16x, E24.16E3, 2x, a7)
232 format(1x, a6, 38x, E25.16E3, 2x, a7)
233 format(1x, a46, E25.16E3, a7)
234 format(1x, a39, 7x, E25.16E3, a7)
  end subroutine write_qe_results

  subroutine prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    !===============================================================================
    !* This subroutine allocates and fills arrays commonly used in other subroutines
    ! to calculate the broadened emissions, maps and momentum tensors.
    ! These arrays include:
    ! 1. Fermi dirac smeared occupations
    ! 2. ARPES mask arrays for angle constrainted emissions
    ! 3. Heaviside step function with 0 replaced by gaussian smearing
    ! written by Felix Mildner, January 2026
    !===============================================================================
    use od_cell, only: num_kpoints_on_node
    use od_electronic, only: nbands, nspins, band_energy, efermi
    use od_parameters, only: photo_theta_centre, photo_theta_halfwidth, photo_temperature
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast
    use od_io, only: io_error, io_file_unit, io_time, io_date
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi
    implicit none

    real(kind=dp), intent(inout), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), intent(inout), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), intent(inout), allocatable, dimension(:, :, :, :) :: emission_gauss

    real(kind=dp) :: norm_vac, width, argument, efinal_temp, e_normal
    real(kind=dp) :: theta_lo, theta_hi
    integer :: N_k, N_spin, n_eigen, gdx, ierr

    if (.not. allocated(fermi_dirac)) then
      allocate (fermi_dirac(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of fermi_dirac failed')
    end if
    fermi_dirac = 0.0_dp

    if (.not. allocated(arpes_mask)) then
      allocate (arpes_mask(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of arpes_mask failed')
    end if
    arpes_mask = 0.0_dp

    if (.not. allocated(emission_gauss)) then
      allocate (emission_gauss(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of emission_gauss failed')
    end if
    emission_gauss = 0.0_dp

    width = kB*photo_temperature
    norm_vac = inv_sqrt_two_pi/width

    ! Theta is measured from the surface normal and cannot physically exceed
    ! 90 deg, so the acceptance window is clipped to 0-90 rather than allowed to
    ! run past it. That clip is load-bearing: theta_arpes keeps its 91 deg
    ! sentinel for states that cannot emit, and an unclipped window could reach
    ! past 90 and admit them.
    theta_lo = max(photo_theta_centre - photo_theta_halfwidth, 0.0_dp)
    theta_hi = min(photo_theta_centre + photo_theta_halfwidth, 90.0_dp)

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N_k) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but exp(+-230) ~ 1E(+-100)
          ! so this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible under/over-flow.
          if (argument .gt. 230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 0.0_dp
          elseif (argument .lt. -230.0_dp) then
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp
          else
            fermi_dirac(n_eigen, N_spin, N_k) = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          ! Calculate the final state energy - e_vacuum
          efinal_temp = band_energy(n_eigen, N_spin, N_k) + temp_photon_energy - evacuum_eff
          ! Is there enough total energy for this kpt/band for E_normal .gt. 0 after passing through surface potential step
          do gdx = 1, photo_gkmax
            ! Unified condition of emission: is the energy along the normal .gt. 0?
            ! Include now the vacuum level and transverse energy to get the final energy along normal
            ! evacuum_eff = efermi + work_function_eff
            e_normal = efinal_temp - E_transverse(gdx, n_eigen, N_spin, N_k)

            if (e_normal .gt. 0.0_dp) then
              emission_gauss(gdx, n_eigen, N_spin, N_k) = 1.0_dp
            else
              emission_gauss(gdx, n_eigen, N_spin, N_k) = gaussian(e_normal, width, 0.0_dp)/norm_vac
            end if

            ! Polar acceptance only. Theta is invariant under the in-plane
            ! symmetry operations, so testing it at the irreducible k-point is
            ! exact. The azimuth is not, and is handled where the emission
            ! direction is actually known - phi_accepted inside the symmetry
            ! loops, phi_accept_frac everywhere else.
            if (theta_arpes(gdx, n_eigen, N_spin, N_k) .ge. theta_lo .and. &
                theta_arpes(gdx, n_eigen, N_spin, N_k) .le. theta_hi) then
              arpes_mask(gdx, n_eigen, N_spin, N_k) = 1.0_dp
            end if
          end do
        end do
      end do
    end do
  end subroutine prepare_emission_arrays

  subroutine deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    !===============================================================================
    !* This subroutine deallocates the commonly used arrays for emission
    ! calculations in the broadened curves and maps
    ! written by Felix Mildner, January 2026
    !===============================================================================
    use od_io, only: io_error
    implicit none

    real(kind=dp), intent(inout), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), intent(inout), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), intent(inout), allocatable, dimension(:, :, :, :) :: emission_gauss

    integer   :: ierr
    if (allocated(fermi_dirac)) then
      deallocate (fermi_dirac, stat=ierr)
      if (ierr /= 0) call io_error('Error: deallocate_emission_arrays - failed to deallocate fermi_dirac')
    end if

    if (allocated(arpes_mask)) then
      deallocate (arpes_mask, stat=ierr)
      if (ierr /= 0) call io_error('Error: deallocate_emission_arrays - failed to deallocate arpes_mask')
    end if

    if (allocated(emission_gauss)) then
      deallocate (emission_gauss, stat=ierr)
      if (ierr /= 0) call io_error('Error: deallocate_emission_arrays - failed to deallocate emission_gauss')
    end if
  end subroutine

  subroutine binding_energy_curve
    !===============================================================================
    !* This subroutine calculates a binding energy vs contributed QE curve and writes
    ! it to a file. Can be thought of as an energy distribution curve (EDC) in an
    ! ARPES experiment
    ! orig. Victor Chang, 7 February 2020
    ! edited Felix Mildner, after August 2024
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, transmit_prob
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_momentum, photo_phi_centre, &
      photo_phi_halfwidth, photo_bindenergy_broadening, iprint, optics_geom, optics_qdir
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :, :) :: binding_temp
    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :) :: qe_atom
    real(kind=dp) :: time0, time1
    real(kind=dp) :: temp_contribution, gk_factor, qe_factor
    real(kind=dp) :: final_fd, be_temp, qe_contrib
    real(kind=dp) :: total_weighted, qe_norm
    integer :: N_k, N_spin, n_eigen_init, n_eigen, n_eigen_final, atom, e_scale, gdx, ierr
    integer :: idx_center, idx_window, window_width, e_min, e_max
    integer :: binding_unit
    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+---------------- Starting E_binding Curve (EDC) Calculation ----------------+'
      call flush(stdout)
    end if
    ! What is the number of expected energy bins with 0.5 eV margin and one extra bin for 0? bin_width = 0.001 eV
    max_energy = int((temp_photon_energy - work_function_eff)*1000) + 501
    ! If we are too far below workfunction, we can expect
    ! to not have any emission so we can skip doing the work
    if (max_energy .lt. -250) then
      if (on_root) write (stdout, '(1x,a78)') '+---------------- No E_binding Curve Calculated - returning -----------------+'
      return
    end if
    ! We are redoing most of the QE calculation due
    ! to memory constraints, so we need the qe_factor again
    qe_factor = 1.0_dp/(cell_area)
    ! How many std. deviations out from the center should
    ! the Gaussian broadening be summed up?
    window_width = 6
    idx_window = ceiling(photo_bindenergy_broadening*window_width*1000)

    allocate (bind_energy(max_energy), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of bind_energy failed')
    bind_energy = 0.0_dp

    allocate (weighted_be_atom(max_energy, max_atoms + 1), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of weighted_be_atom failed')
    weighted_be_atom = 0.0_dp

    allocate (binding_temp(max_energy, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - allocation of binding_temp failed')
    binding_temp = 0.0_dp

    do e_scale = 1, max_energy
      ! Calculate the bins' energy from -0.5 eV up to the max bin
      bind_energy(e_scale) = (e_scale - 1)*0.001_dp - 0.5_dp
    end do

    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          be_temp = efermi - band_energy(n_eigen, N_spin, N_k)
          idx_center = ceiling(be_temp*1000) + 501
          do e_scale = max(idx_center - idx_window, 1), min(idx_center + idx_window, max_energy)
            binding_temp(e_scale, n_eigen, N_spin, N_k) = &
              gaussian(be_temp, photo_bindenergy_broadening, bind_energy(e_scale))
          end do
        end do
      end do
    end do

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    total_be_contribs = 0.0_dp

    if (index(photo_model, '3step') .gt. 0) then
      ! Emission probability for all explicitly considered atoms
      call photo_calculate_delta(delta_temp, .false.)
      do atom = 1, max_atoms
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
              final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
              do n_eigen_init = 1, n_eigen_final - 1
                idx_center = ceiling((efermi - band_energy(n_eigen_init, N_spin, N_k))*1000) + 501
                e_min = max(idx_center - idx_window, 1)
                e_max = min(idx_center + idx_window, max_energy)
                temp_contribution = &
                  qe_factor*photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                  *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k)*transmit_prob(n_eigen_final, N_k, N_spin) &
                  *electrons_per_state*kpoint_weight(N_k)*(I_layer(box_atom(atom), current_photo_energy_index)) &
                  *fermi_dirac(n_eigen_init, N_spin, N_k)*final_fd &
                  *(pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                    /pdos_weights_k_band(n_eigen_init, N_spin, N_k)) &
                  *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
                do gdx = 1, photo_gkmax
                  gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k)*phi_accept_frac(gdx, n_eigen_init, N_spin, N_k) &
                              *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                              *electron_esc(gdx, n_eigen_init, N_spin, N_k, atom) &
                              *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
                  qe_contrib = temp_contribution*gk_factor
                  total_be_contribs = total_be_contribs + qe_contrib
                  weighted_be_atom(e_min:e_max, atom) = &
                    weighted_be_atom(e_min:e_max, atom) + qe_contrib*binding_temp(e_min:e_max, n_eigen_init, N_spin, N_k)
                end do ! gk
              end do ! band_initial
            end do ! bands_final
          end do ! spins
        end do ! kpts
      end do ! atoms

      ! Emission probability for the extended bulk slab
      call photo_calculate_delta(delta_temp, .true.)
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_final = min_index_unocc(N_spin, N_k), nbands
            final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
            do n_eigen_init = 1, n_eigen_final - 1
              idx_center = ceiling((efermi - band_energy(n_eigen_init, N_spin, N_k))*1000) + 500
              e_min = max(idx_center - idx_window, 1)
              e_max = min(idx_center + idx_window, max_energy)
              temp_contribution = &
                (qe_factor*photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                 *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                 *transmit_prob(n_eigen_final, N_k, N_spin) &
                 *electrons_per_state*kpoint_weight(N_k) &
                 *fermi_dirac(n_eigen_init, N_spin, N_k)*final_fd &
                 *(pdos_weights_boxes(n_eigen_init, N_spin, N_k, num_boxes) &
                   /pdos_weights_k_band(n_eigen_init, N_spin, N_k))) &
                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k))
              do gdx = 1, photo_gkmax
                gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k)*phi_accept_frac(gdx, n_eigen_init, N_spin, N_k) &
                            *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, max_atoms + 1) &
                            *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
                qe_contrib = temp_contribution*gk_factor
                total_be_contribs = total_be_contribs + qe_contrib
                weighted_be_atom(e_min:e_max, max_atoms + 1) = &
                  weighted_be_atom(e_min:e_max, max_atoms + 1) + qe_contrib*binding_temp(e_min:e_max, n_eigen_init, N_spin, N_k)
              end do ! gk
            end do ! band_initial
          end do ! bands_final
        end do ! spins
      end do ! kpts

    elseif (index(photo_model, '1step') .gt. 0) then
      do atom = 1, max_atoms + 1
        do N_k = 1, num_kpoints_on_node(my_node_id)
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              idx_center = ceiling((efermi - band_energy(n_eigen, N_spin, N_k))*1000) + 500
              idx_window = ceiling(photo_bindenergy_broadening*window_width*1000)
              e_min = max(idx_center - idx_window, 1)
              e_max = min(idx_center + idx_window, max_energy)
              temp_contribution = (qe_factor*foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                   *electrons_per_state*kpoint_weight(N_k) &
                                   *I_layer(box_atom(atom), current_photo_energy_index) &
                                   *fermi_dirac(n_eigen, N_spin, N_k) &
                                   *(pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                                     /pdos_weights_k_band(n_eigen, N_spin, N_k))) &
                                  *(1.0_dp + field_emission(n_eigen, N_spin, N_k))
              do gdx = 1, photo_gkmax
                gk_factor = arpes_mask(gdx, n_eigen, N_spin, N_k)*phi_accept_frac(gdx, n_eigen, N_spin, N_k) &
                            *gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                            *electron_esc(gdx, n_eigen, N_spin, N_k, atom) &
                            *emission_gauss(gdx, n_eigen, N_spin, N_k)
                qe_contrib = temp_contribution*gk_factor
                total_be_contribs = total_be_contribs + qe_contrib
                weighted_be_atom(e_min:e_max, atom) = &
                  weighted_be_atom(e_min:e_max, atom) + binding_temp(e_min:e_max, n_eigen, N_spin, N_k)*qe_contrib
              end do ! gk
            end do ! bands
          end do ! spins
        end do ! kpts
      end do ! atoms
    end if

    ! Reverse memory mapping order for smoother printing to file
    allocate (qe_atom(max_atoms + 1, max_energy), stat=ierr)
    if (ierr /= 0) call io_error('Error: write_qe_tensor - allocation of qe_atom failed')
    qe_atom = 0.0_dp
    do e_scale = 1, max_energy
      do atom = 1, max_atoms + 1
        qe_atom(atom, e_scale) = weighted_be_atom(e_scale, atom)
      end do
    end do

    ! Gather data from nodes
    call comms_reduce(qe_atom(1, 1), max_energy*(max_atoms + 1), "SUM")
    call comms_reduce(total_be_contribs, 1, "SUM")

    if (on_root) then
      total_weighted = sum(qe_atom(:, :))
      ! Rescale the broadened contributions array
      ! to the sum of unbroadened contributions
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      qe_atom = qe_atom*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))// &
                 '_bindenergy_curve.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)//'_' &
        //trim(adjustl(char_e))//'_bindenergy_curve.dat'
      binding_unit = io_file_unit()
      open (unit=binding_unit, action='write', file=filename)

      call io_date(cdate, ctime)
      write (binding_unit, '(1x,a60,a11,a4,a9)') '## OptaDOS Photoemission: Printing Broadened Binding Energy on ',&
      & cdate, ' at ', ctime
      write (binding_unit, '(1x,a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (binding_unit, '(1x,a24,a12)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (binding_unit, '(1x,a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (binding_unit, '(1x,a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (binding_unit, '(1x,a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (binding_unit, '(1x,a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (binding_unit, '(1x,a35,f9.5)') '## Binding Energy Broadening [eV]: ', photo_bindenergy_broadening
      write (binding_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg]: ', &
        photo_theta_centre, photo_theta_halfwidth
      write (binding_unit, '(1x,a63,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg]: ', &
        photo_phi_centre, photo_phi_halfwidth
      write (binding_unit, '(1x,a34,f9.5)') '## Fermi Energy Ekin offset [eV]: ', (temp_photon_energy - work_function_eff)
      write (binding_unit, '(1x,a66,1x,a50)') '## Binding Energy (EB) [eV] | Total QE from sum(atoms + bulk) @ EB',&
      &'| Contributions from: atom1 | atom2 | ... | bulk |'

      write (out_string, '(a,I0,"(1x,",a,")")') "1x,ES25.6E2,", max_atoms + 2, "ES25.12E3"
      do e_scale = 1, max_energy
        write (binding_unit, '('//trim(out_string)//')') bind_energy(e_scale), &
          sum(qe_atom(1:max_atoms + 1, e_scale)), qe_atom(1:max_atoms + 1, e_scale)
      end do

      close (unit=binding_unit)
    end if

    deallocate (qe_atom, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - failed to deallocate qe_atom')

    deallocate (binding_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - failed to deallocate binding_temp')

    deallocate (bind_energy, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - failed to deallocate bind_energy')

    deallocate (weighted_be_atom, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_curve - failed to deallocate weighted_be_atom')

    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate binding energy curve', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
  end subroutine binding_energy_curve

  subroutine kinetic_energy_momentum_map
    !*===============================================================================
    ! This subroutine calculates a binding energy vs reciprocal transverse momentum
    ! map of the gaussian broadened band contributions and writes it to a file.
    ! Can be thought of the bandstructure projection along the transverse diagonal
    ! showing the contributions of emitting bands projected onto that diagonal.
    ! written by Felix Mildner, after May 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_grid_dim, recip_lattice
    use od_electronic, only: nbands, nspins
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_momentum, photo_phi_centre, &
      photo_phi_halfwidth, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, optics_geom, optics_qdir
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, e_mass, hbar, ev_to_j
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:)          :: gauss_k

    integer :: i, N_k, N_spin, n_eigen_init, n_eigen, first_final, kdx, edx, ierr, window_width, &
               matrix_unit, k_window, e_window, center_bin_e, center_bin_k, kdx_min, kdx_max, edx_min, edx_max
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp
    real(kind=dp) :: step(1:2), sub_cell_length(1:2)
    real(kind=dp) :: k_broadening, temp_k, min_e, gauss_e, e_temp, qe_contrib, &
                     total_weighted, qe_norm, qe_factor, gk, time0, time1

    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+--------- Starting Binding Energy vs Transverse P Map Calculation ----------+'
      call flush(stdout)
    end if
    ! We are redoing parts of the QE calculation, so we need these factors
    qe_factor = 1.0_dp/(cell_area)
    ! How many standard deviations out from the center should the Gaussian broadening be summed up?
    max_energy = int((temp_photon_energy - work_function_eff)*1000) + 500
    if (max_energy .lt. -250) then
      write (stdout, '(1x,a78)') '+-------------- No E_kin vs p_trans map calculated - returning --------------+'
      return
    end if

    call cell_calc_kpoint_r_cart
    step(:) = 0.5_dp/real(kpoint_grid_dim(1:2), dp)
    do i = 1, 2
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do
    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2))
    ! so that FWHM = the step distance between the kpoints
    k_broadening = sqrt(sub_cell_length(1)**2 + sub_cell_length(2)**2)/(4.70964009_dp)

    ! calculate the number of bins to go left and right of center
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    k_window = window_width*ceiling(k_broadening/photo_pmat_bin_width)
    e_window = window_width*ceiling(photo_bindenergy_broadening/photo_pmat_bin_width)
    ! get the maximum transverse k value
    max_k_transverse = 0.0_dp
    do N_k = 1, num_kpoints_on_node(my_node_id)
      temp_k = sqrt(kpoint_r_cart(1, N_k)**2 + kpoint_r_cart(2, N_k)**2)
      max_k_transverse = max(max_k_transverse, temp_k)
    end do
    call comms_reduce(max_k_transverse, 1, "MAX")
    call comms_bcast(max_k_transverse, 1)
    ! Pad the extent so the gaussian broadening of the outermost contributions
    ! is shown rather than cut off. The padded value is kept, because it is what
    ! the header reports as the extent of the momentum axis.
    max_k_transverse = max_k_transverse + k_extra_padding
    max_bin_k = ceiling(max_k_transverse/photo_pmat_bin_width)

    ! calculating upper bound of energy range with some extra for plotting
    max_e_kinetic = temp_photon_energy - work_function_eff + plot_extra_upper
    ! calculating lower bound of energy range
    ! Restrict lower E_kinetic bound to either -0.25 eV or (minimal E_kinetic - 0.25 eV)
    ! This makes sure the program does not print huge matrices at higher photon energies
    min_e = max(minval(E_kinetic) - 0.25_dp, -0.25_dp)
    call comms_reduce(min_e, 1, 'MIN')
    call comms_bcast(min_e, 1)
    max_bin_e = ceiling((max_e_kinetic - min_e)/photo_pmat_bin_width)

    if (max_bin_e .lt. 0 .or. max_bin_k .lt. 0) then
      write (stdout, '(1x,a78)') '+-------------------- map array size negative - returning -------------------+'
      return
    end if

    ! set up the matrix of energy vs transverse k
    allocate (ekin_k_matrix(max_bin_k, max_bin_e), stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map - allocation of ekin_k_matrix failed')
    ekin_k_matrix = 0.0_dp

    allocate (gauss_k(max_bin_k), stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map - allocation of gauss_k failed')
    gauss_k = 0.0_dp

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    ! Reset the running total of unbroadened contributions. It is a module
    ! variable, so without this it keeps its value from the previous photon
    ! energy of a sweep and qe_norm below normalises the matrix to the QE
    ! summed over every step so far instead of this one.
    total_be_kmat_contribs = 0.0_dp

    if (index(photo_model, '3step') .gt. 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)
        temp_k = sqrt(kpoint_r_cart(1, N_k)**2 + kpoint_r_cart(2, N_k)**2)
        ! calculate the bin position in k and e
        center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
        kdx_min = max(center_bin_k - k_window, 1)
        kdx_max = min(center_bin_k + k_window, max_bin_k)
        gk = (kdx_min - 1)*photo_pmat_bin_width
        do kdx = kdx_min, kdx_max
          gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
          gk = gk + photo_pmat_bin_width
        end do
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            ! The energy and k patches depend on the initial band and the k-point
            ! only, so the sums over final bands and over atoms factor out of the
            ! accumulation. The pair set is unchanged: a pair contributes when
            ! final >= min_index_unocc and init < final.
            first_final = max(min_index_unocc(N_spin, N_k), n_eigen_init + 1)
            qe_contrib = sum(qe_tsm(n_eigen_init, first_final:nbands, N_spin, N_k, 1:max_atoms)) &
                         *arpes_mask(1, n_eigen_init, N_spin, N_k)*phi_accept_frac(1, n_eigen_init, N_spin, N_k)
            total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
            if (abs(qe_contrib) .lt. tiny_contribution) cycle
            center_bin_e = ceiling((E_kinetic(1, n_eigen_init, N_spin, N_k) - min_e)/photo_pmat_bin_width)
            edx_min = max(center_bin_e - e_window, 1)
            edx_max = min(center_bin_e + e_window, max_bin_e)
            e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
            do edx = edx_min, edx_max
              gauss_e = gaussian(E_kinetic(1, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, e_temp)
              e_temp = e_temp + photo_pmat_bin_width
              ekin_k_matrix(kdx_min:kdx_max, edx) = ekin_k_matrix(kdx_min:kdx_max, edx) &
                                                    + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
            end do ! e_kinetic
          end do ! bands_init
        end do ! spins
      end do ! kpts

      do N_k = 1, num_kpoints_on_node(my_node_id)
        temp_k = sqrt(kpoint_r_cart(1, N_k)**2 + kpoint_r_cart(2, N_k)**2)
        ! calculate the bin position in k and e
        center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
        kdx_min = max(center_bin_k - k_window, 1)
        kdx_max = min(center_bin_k + k_window, max_bin_k)
        gk = (kdx_min - 1)*photo_pmat_bin_width
        do kdx = kdx_min, kdx_max
          gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
          gk = gk + photo_pmat_bin_width
        end do
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            first_final = max(min_index_unocc(N_spin, N_k), n_eigen_init + 1)
            qe_contrib = sum(qe_tsm(n_eigen_init, first_final:nbands, N_spin, N_k, max_atoms + 1)) &
                         *arpes_mask(1, n_eigen_init, N_spin, N_k)*phi_accept_frac(1, n_eigen_init, N_spin, N_k)
            total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
            if (abs(qe_contrib) .lt. tiny_contribution) cycle
            center_bin_e = ceiling((E_kinetic(1, n_eigen_init, N_spin, N_k) - min_e)/photo_pmat_bin_width)
            edx_min = max(center_bin_e - e_window, 1)
            edx_max = min(center_bin_e + e_window, max_bin_e)
            e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
            do edx = edx_min, edx_max
              gauss_e = gaussian(E_kinetic(1, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, e_temp)
              e_temp = e_temp + photo_pmat_bin_width
              ekin_k_matrix(kdx_min:kdx_max, edx) = ekin_k_matrix(kdx_min:kdx_max, edx) &
                                                    + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
            end do ! e_kinetic
          end do ! bands_init
        end do ! spins
      end do ! kpts
    end if

    if (index(photo_model, '1step') .gt. 0) then
      do N_k = 1, num_kpoints_on_node(my_node_id)
        temp_k = sqrt(kpoint_r_cart(1, N_k)**2 + kpoint_r_cart(2, N_k)**2)
        ! calculate the bin position in k and e
        center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
        kdx_min = max(center_bin_k - k_window, 1)
        kdx_max = min(center_bin_k + k_window, max_bin_k)
        gk = (kdx_min - 1)*photo_pmat_bin_width
        do kdx = kdx_min, kdx_max
          gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
          gk = gk + photo_pmat_bin_width
        end do
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            ! Neither patch carries an atom index, so sum the atoms here.
            qe_contrib = sum(qe_osm(n_eigen, N_spin, N_k, 1:max_atoms + 1)) &
                         *arpes_mask(1, n_eigen, N_spin, N_k)*phi_accept_frac(1, n_eigen, N_spin, N_k)
            total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
            if (abs(qe_contrib) .lt. tiny_contribution) cycle
            center_bin_e = ceiling((E_kinetic(1, n_eigen, N_spin, N_k) - min_e)/photo_pmat_bin_width)
            edx_min = max(center_bin_e - e_window, 1)
            edx_max = min(center_bin_e + e_window, max_bin_e)
            e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
            do edx = edx_min, edx_max
              gauss_e = gaussian(E_kinetic(1, n_eigen, N_spin, N_k), photo_bindenergy_broadening, e_temp)
              e_temp = e_temp + photo_pmat_bin_width
              ekin_k_matrix(kdx_min:kdx_max, edx) = ekin_k_matrix(kdx_min:kdx_max, edx) &
                                                    + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
            end do ! e_kinetic
          end do ! bands
        end do ! spins
      end do ! kpts
    end if

    call comms_reduce(ekin_k_matrix(1, 1), max_bin_e*max_bin_k, 'SUM')
    call comms_reduce(total_be_kmat_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(ekin_k_matrix(:, :))
      ! Rescale the broadened contributions array
      ! to the sum of unbroadened contributions
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_kmat_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      ekin_k_matrix = ekin_k_matrix*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))// &
                 '_Ekin_ptrans_map.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)//'_' &
        //trim(adjustl(char_e))//'_Ekin_ptrans_map.dat'
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)

      call io_date(cdate, ctime)
      write (matrix_unit, '(a64,a11,a4,a9)') '## OptaDOS Photoemission: Kinetic Energy vs P_transverse matrix ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a14,a)') '## Seedname : ', trim(adjustl(seedname))
      write (matrix_unit, '(a25,a12)') '## Photoemission Model : ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a24,f7.3)') '## Photon Energy [eV] : ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a36,f9.5)') '## Binding Energy Broadening [eV] : ', photo_bindenergy_broadening
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg] : ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg] : ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a35,f9.5)') '## Fermi Energy Ekin offset [eV] : ', max_e_kinetic - plot_extra_upper
      write (matrix_unit, '(a34,f9.5)') '## Max k_transverse value [1/A] : ', max_k_transverse
      write (matrix_unit, '(a20,f9.5)') '## Bin width [eV] : ', photo_pmat_bin_width
      write (matrix_unit, '(a19,2(1x,I10),a2)') '## Matrix Shape : (', max_bin_e, max_bin_k, ' )'

      write (out_string, '(I0,"(1x,",a,")")') max_bin_k, 'ES25.12E3'
      do edx = 1, max_bin_e
        write (matrix_unit, '('//trim(out_string)//')') (ekin_k_matrix(kdx, edx), kdx=1, max_bin_k)
      end do

      close (unit=matrix_unit)
    end if
    ! Safety comms sync
    call comms_bcast(qe_norm, 1)

    deallocate (gauss_k, stat=ierr)
    if (ierr /= 0) call io_error('Error : kinetic_energy_momentum_map - failed to deallocate gauss_k')
    deallocate (ekin_k_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map - failed to deallocate ekin_k_matrix')
    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a46,13x,f11.3,a8)') '+ Time to calculate kinetic energy p_trans map', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if
  end subroutine kinetic_energy_momentum_map

  subroutine kinetic_energy_momentum_map_gkgrid
    !*===============================================================================
    ! This subroutine calculates a kinetic energy vs reciprocal transverse momentum
    ! map of the gaussian broadened band contributions and writes it to a file.
    ! This is the optimised version for the photo_momentum option to allow supercell
    ! calculations. Can be thought of the bandstructure projection along the
    ! transverse diagonal showing the contributions of emitting bands.
    ! written by Felix Mildner, after May 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_weight, &
      kpoint_grid_dim, recip_lattice
    use od_electronic, only: nbands, nspins, electrons_per_state, transmit_prob, &
      photo_gkgrid, elec_read_gk_grid
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_momentum, photo_phi_centre, &
      photo_phi_halfwidth, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, optics_geom, optics_qdir
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, e_mass, hbar, ev_to_j
    implicit none

    integer :: i, N_k, N_spin, n_eigen_init, n_eigen, n_eigen_final, atom, kdx, edx, gdx, ierr
    integer :: window_width, matrix_unit
    integer :: k_window, e_window, center_bin_e, center_bin_k, kdx_min, kdx_max, edx_min, edx_max
    real(kind=dp) :: temp_contribution, gk_factor, qe_factor, gk
    real(kind=dp) :: step(1:2), sub_cell_length(1:2), k_broadening, temp_k, min_e, gauss_e, e_temp
    real(kind=dp) :: final_fd, qe_contrib, total_weighted, qe_norm
    real(kind=dp) :: sum_final, sum_atoms
    ! Contributions below this are dropped rather than broadened into the map:
    ! the patch accumulation is the expensive part and a zero adds nothing.
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp
    real(kind=dp) :: time0, time1

    real(kind=dp), allocatable, dimension(:, :, :) :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:)          :: gauss_k

    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+--------- Starting Binding Energy vs Transverse P Map Calculation ----------+'
      call flush(stdout)
    end if
    ! We are redoing parts of the QE calculation, so we need these factors
    qe_factor = 1.0_dp/(cell_area)
    ! How many standard deviations out from the center should the Gaussian broadening be summed up?
    max_energy = int((temp_photon_energy - work_function_eff)*1000) + 500
    if (max_energy .lt. -250) then
      write (stdout, '(1x,a78)') '+-------------- No E_kin vs p_trans map calculated - returning --------------+'
      return
    end if

    call cell_calc_kpoint_r_cart
    step(:) = 0.5_dp/real(kpoint_grid_dim(1:2), dp)
    do i = 1, 2
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do
    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2))=4.70964...
    ! so that FWHM = the step distance between the kpoints
    k_broadening = sqrt(sub_cell_length(1)**2 + sub_cell_length(2)**2)/(4.70964009_dp)

    ! calculate the number of bins to go left and right of center
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    k_window = window_width*ceiling(k_broadening/photo_pmat_bin_width)
    e_window = window_width*ceiling(photo_bindenergy_broadening/photo_pmat_bin_width)
    ! Get the maximum k as the k of a PW with E_excess+0.5 eV. Unlike the
    ! crystal-momentum variant this cannot be taken from the k-point grid --
    ! the G+k sphere would make the matrices enormous -- so it is bounded here
    ! and then padded by the same allowance, again so that the broadening of
    ! the outermost contributions is not cut off.
    max_k_transverse = sqrt((2*e_mass*((temp_photon_energy - work_function_eff + 0.5)*ev_to_j))/(hbar*hbar))*1E-10
    max_k_transverse = max_k_transverse + k_extra_padding
    max_bin_k = ceiling(max_k_transverse/photo_pmat_bin_width)

    call elec_read_gk_grid()

    ! calculating upper bound of energy range with some extra for plotting
    max_e_kinetic = temp_photon_energy - work_function_eff + plot_extra_upper
    ! calculating lower bound of energy range
    ! Restrict lower E_kinetic bound to either -0.25 eV or minimal E_kinetic
    ! This makes sure the program does not print huge matrices at higher photon energies
    min_e = max(minval(E_kinetic) - 0.25_dp, -0.25_dp)
    call comms_reduce(min_e, 1, 'MIN')
    call comms_bcast(min_e, 1)
    max_bin_e = ceiling((max_e_kinetic - min_e)/photo_pmat_bin_width)

    if (max_bin_e .lt. 0 .or. max_bin_k .lt. 0) then
      write (stdout, '(1x,a78)') '+-------------------- map array size negative - returning -------------------+'
      return
    end if

    ! set up the matrix of energy vs transverse k
    allocate (ekin_k_matrix(max_bin_k, max_bin_e), stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map_gkgrid - allocation of ekin_k_matrix failed')
    ekin_k_matrix = 0.0_dp

    allocate (gauss_k(max_bin_k), stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map_gkgrid - allocation of gauss_k failed')
    gauss_k = 0.0_dp

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    ! Reset the running total of unbroadened contributions. It is a module
    ! variable, so without this it keeps its value from the previous photon
    ! energy of a sweep and qe_norm below normalises the matrix to the QE
    ! summed over every step so far instead of this one.
    total_be_kmat_contribs = 0.0_dp

    if (index(photo_model, '3step') .gt. 0) then

      call photo_calculate_delta(delta_temp, .false.)
      ! The broadening on both axes is binned from E_kinetic and photo_gkgrid,
      ! which are indexed (gdx, band, spin, kpt) only -- neither depends on the
      ! atom nor on the final band. The Gaussian patch is therefore identical
      ! for every atom and every final band, and factors out of both sums:
      !
      !   sum_atom sum_final (contribution * patch) = (sum_atom sum_final contribution) * patch
      !
      ! so the patch is evaluated and accumulated once per (kpt, spin, initial
      ! band, gdx) rather than once per (atom, final band, ...) as well. The
      ! two sums separate further, since the final-band factors do not depend
      ! on the atom or on gdx and the atom factors do not depend on the final
      ! band. Nothing is stored: the saving is in the loop order alone.
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            ! Factors carrying only the initial band.
            temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            ! Sum over final bands: independent of atom and of gdx.
            sum_final = 0.0_dp
            do n_eigen_final = n_eigen_init + 1, nbands
              final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
              sum_final = sum_final &
                          + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
            end do
            if (abs(sum_final) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              ! Everything here is invariant under the atom and final-band sums.
              gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k)*phi_accept_frac(gdx, n_eigen_init, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                          *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
              ! Below threshold emission_gauss is exponentially small and the
              ! mask is hard zero, so most of the grid is skipped outright.
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              ! Sum over atoms: the only atom dependence left is the escape
              ! length, the light intensity and the projection onto that atom.
              sum_atoms = 0.0_dp
              do atom = 1, max_atoms
                sum_atoms = sum_atoms &
                            + I_layer(box_atom(atom), current_photo_energy_index) &
                            *pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, atom)
              end do

              qe_contrib = temp_contribution*sum_final*gk_factor*sum_atoms
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib

              temp_k = sqrt(photo_gkgrid(1, gdx, n_eigen_init, N_spin, N_k)**2 + &
                            photo_gkgrid(2, gdx, n_eigen_init, N_spin, N_k)**2)
              center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
              kdx_min = max(center_bin_k - k_window, 1)
              kdx_max = min(center_bin_k + k_window, max_bin_k)
              gk = (kdx_min - 1)*photo_pmat_bin_width
              do kdx = kdx_min, kdx_max
                gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
                gk = gk + photo_pmat_bin_width
              end do

              center_bin_e = ceiling((E_kinetic(gdx, n_eigen_init, N_spin, N_k) - min_e)/photo_pmat_bin_width)
              edx_min = max(center_bin_e - e_window, 1)
              edx_max = min(center_bin_e + e_window, max_bin_e)
              e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
              do edx = edx_min, edx_max
                gauss_e = gaussian(E_kinetic(gdx, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, &
                                   e_temp)
                e_temp = e_temp + photo_pmat_bin_width
                ekin_k_matrix(kdx_min:kdx_max, edx) = &
                  ekin_k_matrix(kdx_min:kdx_max, edx) + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
              end do ! e_kinetic
            end do ! gk
          end do ! band_init
        end do ! spins
      end do ! kpts

      call photo_calculate_delta(delta_temp, .true.)

      ! Same factorisation as the explicit layers above. There is no atom loop
      ! here -- the extrapolated bulk is one region, indexed max_atoms + 1 --
      ! so only the final-band sum and the patch come out.
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                *pdos_weights_boxes(n_eigen_init, N_spin, N_k, num_boxes) &
                                /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            sum_final = 0.0_dp
            do n_eigen_final = n_eigen_init + 1, nbands
              final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
              sum_final = sum_final &
                          + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
            end do
            if (abs(sum_final) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k)*phi_accept_frac(gdx, n_eigen_init, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                          *electron_esc(gdx, n_eigen_init, N_spin, N_k, max_atoms + 1) &
                          *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              qe_contrib = temp_contribution*sum_final*gk_factor
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib

              temp_k = sqrt(photo_gkgrid(1, gdx, n_eigen_init, N_spin, N_k)**2 &
                            + photo_gkgrid(2, gdx, n_eigen_init, N_spin, N_k)**2)
              center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
              kdx_min = max(center_bin_k - k_window, 1)
              kdx_max = min(center_bin_k + k_window, max_bin_k)
              gk = (kdx_min - 1)*photo_pmat_bin_width
              do kdx = kdx_min, kdx_max
                gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
                gk = gk + photo_pmat_bin_width
              end do

              center_bin_e = ceiling((E_kinetic(gdx, n_eigen_init, N_spin, N_k) - min_e)/photo_pmat_bin_width)
              edx_min = max(center_bin_e - e_window, 1)
              edx_max = min(center_bin_e + e_window, max_bin_e)
              e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
              do edx = edx_min, edx_max
                gauss_e = gaussian(E_kinetic(gdx, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, &
                                   e_temp)
                e_temp = e_temp + photo_pmat_bin_width
                ekin_k_matrix(kdx_min:kdx_max, edx) = &
                  ekin_k_matrix(kdx_min:kdx_max, edx) + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
              end do ! e_kinetic
            end do ! gk
          end do ! band_init
        end do ! spins
      end do ! kpts
    end if

    if (index(photo_model, '1step') .gt. 0) then
      ! One band index here, so only the atom loop factors out of the patch --
      ! the sum over atoms and the extrapolated bulk becomes a scalar.
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            temp_contribution = qe_factor*foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                *electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen, N_spin, N_k)) &
                                /pdos_weights_k_band(n_eigen, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              gk_factor = arpes_mask(gdx, n_eigen, N_spin, N_k)*phi_accept_frac(gdx, n_eigen, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                          *emission_gauss(gdx, n_eigen, N_spin, N_k)
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              sum_atoms = 0.0_dp
              do atom = 1, max_atoms + 1
                sum_atoms = sum_atoms &
                            + I_layer(box_atom(atom), current_photo_energy_index) &
                            *pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                            *electron_esc(gdx, n_eigen, N_spin, N_k, atom)
              end do

              qe_contrib = temp_contribution*gk_factor*sum_atoms
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib

              temp_k = sqrt(photo_gkgrid(1, gdx, n_eigen, N_spin, N_k)**2 + photo_gkgrid(2, gdx, n_eigen, N_spin, N_k)**2)
              center_bin_k = ceiling(temp_k/photo_pmat_bin_width)
              kdx_min = max(center_bin_k - k_window, 1)
              kdx_max = min(center_bin_k + k_window, max_bin_k)
              gk = (kdx_min - 1)*photo_pmat_bin_width
              do kdx = kdx_min, kdx_max
                gauss_k(kdx) = gaussian(temp_k, k_broadening, gk)
                gk = gk + photo_pmat_bin_width
              end do

              center_bin_e = ceiling((E_kinetic(gdx, n_eigen, N_spin, N_k) - min_e)/photo_pmat_bin_width)
              edx_min = max(center_bin_e - e_window, 1)
              edx_max = min(center_bin_e + e_window, max_bin_e)
              e_temp = min_e + (edx_min - 1)*photo_pmat_bin_width
              do edx = edx_min, edx_max
                gauss_e = gaussian(E_kinetic(gdx, n_eigen, N_spin, N_k), photo_bindenergy_broadening, e_temp)
                e_temp = e_temp + photo_pmat_bin_width
                ekin_k_matrix(kdx_min:kdx_max, edx) = &
                  ekin_k_matrix(kdx_min:kdx_max, edx) + gauss_e*gauss_k(kdx_min:kdx_max)*qe_contrib
              end do ! e_kinetic
            end do ! gk
          end do ! bands
        end do ! spins
      end do ! kpts
    end if

    call comms_reduce(ekin_k_matrix(1, 1), max_bin_e*max_bin_k, 'SUM')
    call comms_reduce(total_be_kmat_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(ekin_k_matrix(:, :))
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_kmat_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      ekin_k_matrix = ekin_k_matrix*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))// &
                 '_Ekin_ptrans_map.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)//'_' &
        //trim(adjustl(char_e))//'_Ebind_ptrans_map.dat'
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)

      write (matrix_unit, '(a56,a11,a4,a9)') '## OptaDOS Photoemission: Energy vs P_transverse matrix ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a14,a)') '## Seedname : ', trim(adjustl(seedname))
      write (matrix_unit, '(a25,a12)') '## Photoemission Model : ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a24,f7.3)') '## Photon Energy [eV] : ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a36,f9.5)') '## Binding Energy Broadening [eV] : ', photo_bindenergy_broadening
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg] : ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg] : ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a35,f9.5)') '## Fermi Energy Ekin offset [eV] : ', max_e_kinetic - plot_extra_upper
      write (matrix_unit, '(a34,f9.5)') '## Max k_transverse value [1/A] : ', max_k_transverse
      write (matrix_unit, '(a20,f9.5)') '## Bin width [eV] : ', photo_pmat_bin_width
      write (matrix_unit, '(a19,2(1x,I10),a2)') '## Matrix Shape : (', max_bin_e, max_bin_k, ' )'

      write (out_string, '(I0,"(1x,",a,")")') max_bin_k, 'ES25.12E3'
      do edx = 1, max_bin_e
        write (matrix_unit, '('//trim(out_string)//')') (ekin_k_matrix(kdx, edx), kdx=1, max_bin_k)
      end do

      close (unit=matrix_unit)
    end if
    ! Safety comms sync
    call comms_bcast(qe_norm, 1)

    deallocate (gauss_k, stat=ierr)
    if (ierr /= 0) call io_error('Error : kinetic_energy_momentum_map_gkgrid - failed to deallocate gauss_k')

    deallocate (photo_gkgrid, stat=ierr)
    if (ierr /= 0) call io_error('Error : kinetic_energy_momentum_map_gkgrid - failed to deallocate photo_gkgrid')

    deallocate (ekin_k_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error: kinetic_energy_momentum_map_gkgrid - failed to deallocate ekin_k_matrix')

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate E_kin vs p_trans map', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if
  end subroutine kinetic_energy_momentum_map_gkgrid

  subroutine full_momentum_tensor
    !*===============================================================================
    ! This subroutine calculates the px,py,pz momentum tensor of emitted electrons,
    ! applies a gaussian broadening to each contribution and writes it to a file.
    ! written by Felix Mildner, after May 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight, &
      kpoint_grid_dim, recip_lattice, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_electronic, only: nbands, nspins
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_momentum, photo_phi_centre, &
      photo_phi_halfwidth, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, devel_flag, optics_geom, &
      optics_qdir, photo_momentum
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, e_mass, ev_to_j, hbar
    implicit none

    integer    ::  i, N_k, N_spin, n_eigen_init, n_eigen, ierr
    integer    ::  matrix_unit, total_ks, nsymm_op, x_center, y_center, z_center, xdx, ydx, zdx
    integer    ::  xdx_offset, ydx_offset, zdx_offset, xdx_window, ydx_window, zdx_window
    integer    ::  xdx_min, xdx_max, ydx_min, ydx_max, zdx_min, zdx_max, window_width
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :)    :: fermi_dirac
    real(kind=dp), allocatable, dimension(:, :, :)    :: p_z
    real(kind=dp), allocatable, dimension(:, :)        :: gauss_xy
    real(kind=dp), allocatable, dimension(:)          :: gauss_y, gauss_x
    real(kind=dp) :: step(1:2), sub_cell_length(1:2), temp_mat(2, 2), current_k(2)
    real(kind=dp) :: qe_contrib, gauss_z, total_weighted, qe_norm
    real(kind=dp) :: kx_broadening, ky_broadening, kz_broadening, k_prefactor
    real(kind=dp) :: z_max, z_min, xy_max, wave_prefactor, etemp, min_e
    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    real(kind=dp) :: time0, time1, time2
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+---------------- Starting Full Momentum Tensor Calculation -----------------+'
      call flush(stdout)
    end if
    if (photo_momentum == 'gkgrid') then
      if (on_root) write (stdout, '(1x,a78)') '+----------- Tensor Calculation with Gkgrid scheme not implemented ----------+'
      return
    end if
    ! get kinetic energy at efermi for reference
    max_e_kinetic = temp_photon_energy - work_function_eff
    if (max_e_kinetic .lt. -0.25_dp) then
      if (on_root) write (stdout, '(1x,a78)') '+-------------- Max. E_kinetic .lt. -0.25eV;  no tensor calculated -------------+'
      return
    end if

    wave_prefactor = 2*e_mass/(hbar**2)

    ! calculating lower bound of energy range
    total_ks = kpoint_grid_dim(1)*kpoint_grid_dim(2)
    do i = 1, 2
      step(i) = 0.5_dp/real(kpoint_grid_dim(i), dp)
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do

    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2)) so that
    ! the FWHM = 1/2 the step distance between the kpoints
    kx_broadening = sub_cell_length(1)/(4.70964009_dp)
    ky_broadening = sub_cell_length(2)/(4.70964009_dp)
    ! kz_broadening = photo_bindenergy_broadening[1/A]/(4.70964009_dp)
    kz_broadening = sqrt((2*e_mass*(photo_bindenergy_broadening*ev_to_j))/(hbar*hbar))*1E-10/(4.70964009_dp)

    ! calculate the number of bins to go left and right
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    xdx_window = ceiling(window_width*kx_broadening/photo_pmat_bin_width)
    ydx_window = ceiling(window_width*ky_broadening/photo_pmat_bin_width)
    zdx_window = ceiling(window_width*kz_broadening/photo_pmat_bin_width)

    call cell_calc_kpoint_r_cart
    z_max = sqrt((2*e_mass*((max_e_kinetic + 0.5)*ev_to_j))/(hbar*hbar))*1E-10
    min_e = max(minval(E_kinetic) - 0.25_dp, 0.0_dp)
    z_min = sqrt((2*e_mass*((min_e)*ev_to_j))/(hbar*hbar))*1E-10
    xy_max = min((abs(maxval(kpoint_r_cart(1:2, :))) + 0.5), z_max)

    xdx_offset = ceiling(xy_max/photo_pmat_bin_width)
    ydx_offset = ceiling(xy_max/photo_pmat_bin_width)
    zdx_offset = ceiling((z_max - z_min)/photo_pmat_bin_width) + zdx_window

    call comms_reduce(xdx_offset, 1, "MAX")
    call comms_reduce(ydx_offset, 1, "MAX")
    call comms_bcast(xdx_offset, 1)
    call comms_bcast(ydx_offset, 1)
    max_bin_p(1) = 2*xdx_offset + 1
    max_bin_p(2) = 2*ydx_offset + 1
    max_bin_p(3) = zdx_offset + 1

    allocate (p_z(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - allocation of p_z failed')
    p_z = 10000.0_dp

    allocate (p_tensor(max_bin_p(1), max_bin_p(2), max_bin_p(3)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - allocation of p_tensor failed')
    p_tensor = 0.0_dp

    allocate (gauss_x(max_bin_p(1)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - allocation of gauss_x failed')
    gauss_x = 0.0_dp

    allocate (gauss_y(max_bin_p(2)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - allocation of gauss_y failed')
    gauss_y = 0.0_dp

    allocate (gauss_xy(max_bin_p(1), max_bin_p(2)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - allocation of gauss_xy failed')
    gauss_y = 0.0_dp

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    ! Reset the running total of unbroadened contributions. It is a module
    ! variable, so without this it keeps its value from the previous photon
    ! energy of a sweep and qe_norm below normalises the matrix to the QE
    ! summed over every step so far instead of this one.
    total_be_kmat_contribs = 0.0_dp
    ! array assignment for turning E_kin into electron momentum p
    ! along surface normal, offset by the minimal z_value, included
    ! in the printout
    do N_k = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, nbands
          etemp = E_kinetic(1, n_eigen, N_spin, N_k) - E_transverse(1, n_eigen, N_spin, N_k)
          if (etemp .gt. 0.0_dp) then
            p_z(n_eigen, N_spin, N_k) = sqrt(wave_prefactor*etemp*ev_to_j)*1E-10_dp - z_min
          end if
        end do
      end do
    end do

    if (index(photo_model, '3step') .gt. 0) then
      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, max_bin_p(1))
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, max_bin_p(2))
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands
              ! The x/y patch and p_z carry no atom index, so the sum over atoms
              ! factors out of the accumulation below and is done here instead.
              qe_contrib = sum(qe_tsm(n_eigen_init, 1:nbands, N_spin, N_k, 1:max_atoms))* &
                           arpes_mask(1, n_eigen_init, N_spin, N_k)*k_prefactor
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              z_center = nint((p_z(n_eigen_init, N_spin, N_k))/photo_pmat_bin_width) + 1
              zdx_min = max(z_center - zdx_window, 1)
              zdx_max = min(z_center + zdx_window, max_bin_p(3))
              do zdx = zdx_min, zdx_max
                gauss_z = gaussian(p_z(n_eigen_init, N_spin, N_k), kz_broadening, (zdx - 1)*photo_pmat_bin_width)
                p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) = &
                  p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) &
                  + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_z*qe_contrib
              end do ! pz
            end do ! bands_initial
          end do ! spins
        end do ! kpts
      end do ! symm_ops

      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, max_bin_p(1))
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, max_bin_p(2))
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands
              qe_contrib = sum(qe_tsm(n_eigen_init, 1:nbands, N_spin, N_k, max_atoms + 1))* &
                           arpes_mask(1, n_eigen_init, N_spin, N_k)*k_prefactor
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
              z_center = nint(p_z(n_eigen_init, N_spin, N_k)/photo_pmat_bin_width) + 1
              zdx_min = max(z_center - zdx_window, 1)
              zdx_max = min(z_center + zdx_window, max_bin_p(3))
              do zdx = zdx_min, zdx_max
                gauss_z = gaussian(p_z(n_eigen_init, N_spin, N_k), kz_broadening, (zdx - 1)*photo_pmat_bin_width)
                p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) = &
                  p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) &
                  + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_z*qe_contrib
              end do ! pz
            end do ! bands_initial
          end do ! spins
        end do ! kpts
      end do ! symm_ops
    end if

    if (index(photo_model, '1step') .gt. 0) then

      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, max_bin_p(1))
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, max_bin_p(2))
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              ! No atom index on the patch or p_z, so sum the atoms here.
              qe_contrib = sum(qe_osm(n_eigen, N_spin, N_k, 1:max_atoms + 1))* &
                           arpes_mask(1, n_eigen, N_spin, N_k)*k_prefactor
              total_be_kmat_contribs = total_be_kmat_contribs + qe_contrib
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              z_center = nint(p_z(n_eigen, N_spin, N_k)/photo_pmat_bin_width) + 1
              zdx_min = max(z_center - zdx_window, 1)
              zdx_max = min(z_center + zdx_window, max_bin_p(3))
              do zdx = zdx_min, zdx_max
                gauss_z = gaussian(p_z(n_eigen, N_spin, N_k), kz_broadening, (zdx - 1)*photo_pmat_bin_width)
                p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) = &
                  p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) &
                  + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_z*qe_contrib
              end do ! pz
            end do ! bands
          end do ! spins
        end do ! kpts
      end do ! symm_ops
    end if

    if (on_root .and. iprint .gt. 1) then
      time2 = io_time()
      write (stdout, '(1x,a46,13x,f11.3,a8)') '+ Time to accumulate the momentum tensor', time2 - time0, ' (sec) +'
      call flush(stdout)
    end if
    call comms_reduce(p_tensor(1, 1, 1), max_bin_p(1)*max_bin_p(2)*max_bin_p(3), 'SUM')
    call comms_reduce(total_be_kmat_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(p_tensor(:, :, :))
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_kmat_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      p_tensor = p_tensor*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_ptensor.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)//'_' &
        //trim(adjustl(char_e))//'_ptensor.dat'
      call io_date(cdate, ctime)

      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      write (matrix_unit, '(a59,a11,a4,a9)') '## OptaDOS Photoemission: Printing Full Momentum Tensor on ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (matrix_unit, '(a24,a12)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg]: ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg]: ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a14,f9.5)') '## Bin width: ', photo_pmat_bin_width
      write (matrix_unit, '(a47)') '## Note: x and y are from -k to +k including 0!'
      write (matrix_unit, '(a29,f9.5)') '## p_z value of first z_bin: ', z_min
      write (matrix_unit, '(a19,3(i7,a3))') '## Matrix Shape: ( ', max_bin_p(1), ' , ', max_bin_p(2), ' , ', max_bin_p(3), ' )'

      write (out_string, '(I0,"(",a,")")') max_bin_p(1), 'E9.1E3'
      do zdx = 1, max_bin_p(3)
        do ydx = 1, max_bin_p(2)
          write (matrix_unit, '('//trim(out_string)//')') (p_tensor(xdx, ydx, zdx), xdx=1, max_bin_p(1))
        end do
      end do
      close (unit=matrix_unit)
    end if

    deallocate (p_z, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - failed to deallocate p_z')

    deallocate (p_tensor, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - failed to deallocate p_tensor')

    deallocate (gauss_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - failed to deallocate gauss_x')

    deallocate (gauss_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - failed to deallocate gauss_y')

    deallocate (gauss_xy, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor - failed to deallocate gauss_xy')

    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate full momentum tensor', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if
  end subroutine full_momentum_tensor

  subroutine accumulate_gk_tensor(gdx, n_eigen, N_spin, N_k, qe_contrib, total_ks, &
                                  xdx_offset, ydx_offset, xdx_window, ydx_window, zdx_window, &
                                  kx_broadening, ky_broadening, kz_broadening, &
                                  p_z, gauss_x, gauss_y, gauss_xy, gauss_z_v)
    !*===============================================================================
    ! Lays one already-summed contribution onto the momentum tensor, once for each
    ! crystal symmetry operation.
    !
    ! qe_contrib carries the sums over atoms and final bands and does not depend on
    ! the symmetry operation -- only the transverse momentum rotates -- so the
    ! caller evaluates it once and this routine reuses it. The z Gaussian likewise
    ! depends only on (gdx, band, spin, kpt) and is built before the symmetry loop.
    !===============================================================================
    use od_cell, only: kpoint_weight, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_electronic, only: photo_gkgrid
    use od_parameters, only: photo_pmat_bin_width, devel_flag
    use od_algorithms, only: gaussian
    implicit none

    integer, intent(in)       :: gdx, n_eigen, N_spin, N_k, total_ks
    integer, intent(in)       :: xdx_offset, ydx_offset, xdx_window, ydx_window, zdx_window
    real(kind=dp), intent(in) :: qe_contrib, kx_broadening, ky_broadening, kz_broadening
    real(kind=dp), intent(in) :: p_z(:, :, :, :)
    real(kind=dp), intent(inout) :: gauss_x(:), gauss_y(:), gauss_xy(:, :), gauss_z_v(:)

    integer       :: nsymm_op, xdx, ydx, zdx, x_center, y_center, z_center
    integer       :: xdx_min, xdx_max, ydx_min, ydx_max, zdx_min, zdx_max
    real(kind=dp) :: temp_mat(2, 2), current_k(2), k_prefactor, weighted

    z_center = nint(p_z(gdx, n_eigen, N_spin, N_k)/photo_pmat_bin_width) + 1
    zdx_min = max(z_center - zdx_window, 1)
    zdx_max = min(z_center + zdx_window, size(gauss_z_v))
    if (zdx_min .gt. zdx_max) return
    do zdx = zdx_min, zdx_max
      gauss_z_v(zdx) = gaussian(p_z(gdx, n_eigen, N_spin, N_k), kz_broadening, &
                                (zdx - 1)*photo_pmat_bin_width)
    end do

    do nsymm_op = 1, num_crystal_symmetry_operations
      temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
      if (index(devel_flag, 'no_symmetry') .gt. 0) then
        current_k = photo_gkgrid(1:2, gdx, n_eigen, N_spin, N_k)
        k_prefactor = 1.0_dp/num_crystal_symmetry_operations
      else
        current_k = matmul(temp_mat, photo_gkgrid(1:2, gdx, n_eigen, N_spin, N_k))
        k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
      end if
      ! As above: this is the emission direction of this image, so the wedge
      ! can be tested here and nowhere earlier.
      if (.not. phi_accepted(current_k(1), current_k(2))) cycle
      weighted = qe_contrib*k_prefactor
      total_be_kmat_contribs = total_be_kmat_contribs + weighted

      x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
      xdx_min = max(x_center - xdx_window, 1)
      xdx_max = min(x_center + xdx_window, size(gauss_x))
      if (xdx_min .gt. xdx_max) cycle
      do xdx = xdx_min, xdx_max
        gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
      end do

      y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
      ydx_min = max(y_center - ydx_window, 1)
      ydx_max = min(y_center + ydx_window, size(gauss_y))
      if (ydx_min .gt. ydx_max) cycle
      do ydx = ydx_min, ydx_max
        gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
      end do

      do ydx = ydx_min, ydx_max
        do xdx = xdx_min, xdx_max
          gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
        end do
      end do

      do zdx = zdx_min, zdx_max
        p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) = &
          p_tensor(xdx_min:xdx_max, ydx_min:ydx_max, zdx) &
          + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_z_v(zdx)*weighted
      end do
    end do ! symm_ops
  end subroutine accumulate_gk_tensor

  subroutine full_momentum_tensor_gkgrid
    !*===============================================================================
    ! This subroutine calculates the px,py,pz momentum tensor of emitted electrons,
    ! applies a gaussian broadening to each contribution and writes it to a file.
    ! written by Felix Mildner, after May 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight, &
      kpoint_grid_dim, recip_lattice, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_electronic, only: nbands, nspins, electrons_per_state, transmit_prob, &
      photo_gkgrid, elec_read_gk_grid
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_momentum, photo_phi_centre, &
      photo_phi_halfwidth, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, devel_flag, optics_geom, &
      optics_qdir
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, e_mass, ev_to_j, hbar
    implicit none

    integer    ::  i, N_k, N_spin, n_eigen_init, n_eigen, n_eigen_final, atom, gdx, ierr
    integer    ::  matrix_unit, total_ks, xdx, ydx, zdx
    integer    ::  xdx_offset, ydx_offset, zdx_offset, xdx_window, ydx_window, zdx_window
    integer    ::  window_width
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :, :) :: p_z
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :)    :: fermi_dirac
    real(kind=dp), allocatable, dimension(:)          :: gauss_y, gauss_x, gauss_z_v
    real(kind=dp), allocatable, dimension(:, :)          :: gauss_xy
    real(kind=dp) :: step(1:2), sub_cell_length(1:2)
    real(kind=dp) :: qe_contrib, total_weighted, qe_norm, sum_final, sum_atoms
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp
    real(kind=dp) :: kx_broadening, ky_broadening, kz_broadening
    real(kind=dp) :: final_fd, z_max, xy_max, wave_prefactor, temp_contribution, gk_factor, qe_factor
    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    real(kind=dp) :: time0, time1

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+---------------- Starting Full Momentum Tensor Calculation -----------------+'
      call flush(stdout)
    end if
    ! get kinetic energy at efermi for reference
    max_e_kinetic = temp_photon_energy - work_function_eff
    if (max_e_kinetic .lt. -0.25_dp) then
      if (on_root) write (stdout, '(1x,a78)') '+-------------- Max. E_kinetic .lt. -0.25eV;  no tensor calculated -------------+'
      return
    end if

    qe_factor = 1.0_dp/(cell_area)
    wave_prefactor = 2*e_mass/(hbar**2)

    ! calculating lower bound of energy range
    total_ks = kpoint_grid_dim(1)*kpoint_grid_dim(2)
    do i = 1, 2
      step(i) = 0.5_dp/real(kpoint_grid_dim(i), dp)
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do

    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2)) so that
    ! the FWHM = 1/2 the step distance between the kpoints
    kx_broadening = sub_cell_length(1)/(4.70964009_dp)
    ky_broadening = sub_cell_length(2)/(4.70964009_dp)
    ! kz_broadening = photo_bindenergy_broadening[1/A]/(4.70964009_dp)
    kz_broadening = sqrt((2*e_mass*(photo_bindenergy_broadening*ev_to_j))/(hbar*hbar))*1E-10/(4.70964009_dp)

    ! calculate the number of bins to go left and right
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    xdx_window = ceiling(window_width*kx_broadening/photo_pmat_bin_width)
    ydx_window = ceiling(window_width*ky_broadening/photo_pmat_bin_width)
    zdx_window = ceiling(window_width*kz_broadening/photo_pmat_bin_width)

    call cell_calc_kpoint_r_cart
    max_e_kinetic = temp_photon_energy - work_function_eff
    z_max = sqrt((2*e_mass*((max_e_kinetic + 0.5)*ev_to_j))/(hbar*hbar))*1E-10
    xy_max = min((abs(maxval(kpoint_r_cart(1:2, :))) + 0.5), z_max)
    xdx_offset = ceiling(xy_max/photo_pmat_bin_width)
    ydx_offset = ceiling(xy_max/photo_pmat_bin_width)
    zdx_offset = ceiling(z_max/photo_pmat_bin_width) + zdx_window

    call comms_reduce(xdx_offset, 1, "MAX")
    call comms_reduce(ydx_offset, 1, "MAX")
    call comms_bcast(xdx_offset, 1)
    call comms_bcast(ydx_offset, 1)
    max_bin_p(1) = 2*xdx_offset + 1
    max_bin_p(2) = 2*ydx_offset + 1
    max_bin_p(3) = zdx_offset + 1

    allocate (p_z(photo_gkmax, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of p_z failed')
    p_z = 10000.0_dp

    allocate (p_tensor(max_bin_p(1), max_bin_p(2), max_bin_p(3)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of p_tensor failed')
    p_tensor = 0.0_dp

    allocate (gauss_x(max_bin_p(1)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of gauss_x failed')
    gauss_x = 0.0_dp

    allocate (gauss_y(max_bin_p(2)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of gauss_y failed')
    gauss_y = 0.0_dp

    allocate (gauss_xy(max_bin_p(1), max_bin_p(2)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of gauss_xy failed')
    gauss_xy = 0.0_dp

    allocate (gauss_z_v(max_bin_p(3)), stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - allocation of gauss_z_v failed')
    gauss_z_v = 0.0_dp

    ! Each gkgrid routine frees photo_gkgrid on the way out, so whichever runs
    ! first deallocates it; read it back rather than assuming it is present.
    call elec_read_gk_grid()

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    ! Reset the running total of unbroadened contributions. It is a module
    ! variable, so without this it keeps its value from the previous photon
    ! energy of a sweep and qe_norm below normalises the matrix to the QE
    ! summed over every step so far instead of this one.
    total_be_kmat_contribs = 0.0_dp
    ! conditional array assignment for turning E_kin into electron momentum p
    ! along surface normal
    where (E_kinetic - E_transverse .gt. 0.0_dp)
    p_z = sqrt(wave_prefactor*((E_kinetic - E_transverse)*ev_to_j))*1E-10_dp
    end where

    if (index(photo_model, '3step') .gt. 0) then

      call photo_calculate_delta(delta_temp, .false.)
      ! The transverse momentum comes from photo_gkgrid, so unlike the crystal
      ! scheme the x/y patch depends on gdx and on the initial band as well as
      ! on the k-point. It still carries no atom index and no final band index,
      ! so both sums factor out of the accumulation:
      !
      !   sum_atom sum_final (contribution * patch) = (sum_atom sum_final contribution) * patch
      !
      ! and the patch is built once per (kpt, spin, initial band, gdx, symmetry
      ! operation). The symmetry loop is innermost because the contribution is
      ! invariant under it -- only the momentum rotates -- so the arithmetic is
      ! done once and reused for every operation.
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            sum_final = 0.0_dp
            do n_eigen_final = n_eigen_init + 1, nbands
              final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
              sum_final = sum_final &
                          + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
            end do
            if (abs(sum_final) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                          *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              sum_atoms = 0.0_dp
              do atom = 1, max_atoms
                sum_atoms = sum_atoms &
                            + I_layer(box_atom(atom), current_photo_energy_index) &
                            *pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, atom)
              end do

              qe_contrib = temp_contribution*sum_final*gk_factor*sum_atoms
              if (abs(qe_contrib) .lt. tiny_contribution) cycle

              call accumulate_gk_tensor(gdx, n_eigen_init, N_spin, N_k, qe_contrib, &
                                        total_ks, xdx_offset, ydx_offset, &
                                        xdx_window, ydx_window, zdx_window, &
                                        kx_broadening, ky_broadening, kz_broadening, &
                                        p_z, gauss_x, gauss_y, gauss_xy, gauss_z_v)
            end do ! gk
          end do ! bands_initial
        end do ! spins
      end do ! kpts

      call photo_calculate_delta(delta_temp, .true.)

      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen_init = 1, nbands - 1
            temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                *pdos_weights_boxes(n_eigen_init, N_spin, N_k, num_boxes) &
                                /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            sum_final = 0.0_dp
            do n_eigen_final = n_eigen_init + 1, nbands
              final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
              sum_final = sum_final &
                          + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                          *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
            end do
            if (abs(sum_final) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                          *electron_esc(gdx, n_eigen_init, N_spin, N_k, max_atoms + 1) &
                          *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              qe_contrib = temp_contribution*sum_final*gk_factor
              if (abs(qe_contrib) .lt. tiny_contribution) cycle

              call accumulate_gk_tensor(gdx, n_eigen_init, N_spin, N_k, qe_contrib, &
                                        total_ks, xdx_offset, ydx_offset, &
                                        xdx_window, ydx_window, zdx_window, &
                                        kx_broadening, ky_broadening, kz_broadening, &
                                        p_z, gauss_x, gauss_y, gauss_xy, gauss_z_v)
            end do ! gk
          end do ! bands_initial
        end do ! spins
      end do ! kpts
    end if

    if (index(photo_model, '1step') .gt. 0) then

      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, nbands
            temp_contribution = qe_factor*foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                *electrons_per_state*kpoint_weight(N_k) &
                                *fermi_dirac(n_eigen, N_spin, N_k) &
                                *(1.0_dp + field_emission(n_eigen, N_spin, N_k)) &
                                /pdos_weights_k_band(n_eigen, N_spin, N_k)
            if (abs(temp_contribution) .lt. tiny_contribution) cycle

            do gdx = 1, photo_gkmax
              gk_factor = arpes_mask(gdx, n_eigen, N_spin, N_k) &
                          *gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                          *emission_gauss(gdx, n_eigen, N_spin, N_k)
              if (abs(gk_factor) .lt. tiny_contribution) cycle

              sum_atoms = 0.0_dp
              do atom = 1, max_atoms + 1
                sum_atoms = sum_atoms &
                            + I_layer(box_atom(atom), current_photo_energy_index) &
                            *pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                            *electron_esc(gdx, n_eigen, N_spin, N_k, atom)
              end do

              qe_contrib = temp_contribution*gk_factor*sum_atoms
              if (abs(qe_contrib) .lt. tiny_contribution) cycle

              call accumulate_gk_tensor(gdx, n_eigen, N_spin, N_k, qe_contrib, &
                                        total_ks, xdx_offset, ydx_offset, &
                                        xdx_window, ydx_window, zdx_window, &
                                        kx_broadening, ky_broadening, kz_broadening, &
                                        p_z, gauss_x, gauss_y, gauss_xy, gauss_z_v)
            end do ! gk
          end do ! bands
        end do ! spins
      end do ! kpts
    end if

    call comms_reduce(p_tensor(1, 1, 1), max_bin_p(1)*max_bin_p(2)*max_bin_p(3), 'SUM')
    call comms_reduce(total_be_kmat_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(p_tensor(:, :, :))
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_kmat_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      p_tensor = p_tensor*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_ptensor.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)//'_' &
        //trim(adjustl(char_e))//'_ptensor.dat'
      call io_date(cdate, ctime)

      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      write (matrix_unit, '(a59,a11,a4,a9)') '## OptaDOS Photoemission: Printing Full Momentum Tensor on ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (matrix_unit, '(a24,a12)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg]: ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg]: ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a14,f9.5)') '## Bin width: ', photo_pmat_bin_width
      write (matrix_unit, '(a61)') '## Note: x and y are from -k to +k including 0, z is 0 to kz!'
      write (matrix_unit, '(a19,3(i7,a3))') '## Matrix Shape: ( ', max_bin_p(1), ' , ', max_bin_p(2), ' , ', max_bin_p(3), ' )'

      write (out_string, '(I0,"(",a,")")') max_bin_p(1), 'E9.1E3'
      do zdx = 1, max_bin_p(3)
        do ydx = 1, max_bin_p(2)
          write (matrix_unit, '('//trim(out_string)//')') (p_tensor(xdx, ydx, zdx), xdx=1, max_bin_p(1))
        end do
      end do
      close (unit=matrix_unit)
    end if

    deallocate (p_z, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - failed to deallocate p_z')

    deallocate (p_tensor, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - failed to deallocate p_tensor')

    deallocate (gauss_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - failed to deallocate gauss_x')

    deallocate (gauss_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - failed to deallocate gauss_y')

    deallocate (gauss_xy, stat=ierr)
    if (ierr /= 0) call io_error('Error: full_momentum_tensor_gkgrid - failed to deallocate gauss_xy')

    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate full momentum tensor', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if
  end subroutine full_momentum_tensor_gkgrid

  subroutine const_binding_energy_map
    !*===============================================================================
    ! This subroutine calculates a map of reciprocal space at a specified binding
    ! energy and writes it out to a file.
    ! written by Felix Mildner, after Jan 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight, &
      kpoint_grid_dim, recip_lattice, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_electronic, only: nbands, nspins, band_energy, efermi
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_phi_centre, photo_phi_halfwidth, &
      photo_momentum, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, &
      devel_flag, optics_geom, optics_qdir, photo_const_bindenergy_value
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, ev_to_j, e_mass, hbar
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :)    :: fermi_dirac
    real(kind=dp), allocatable, dimension(:)          :: gauss_y, gauss_x
    real(kind=dp), allocatable, dimension(:, :)        :: gauss_xy
    real(kind=dp) :: step(1:2), sub_cell_length(1:2), gauss_e, temp_mat(2, 2), current_k(2), z_max, xy_max
    real(kind=dp) :: k_prefactor, ref_level, kx_broadening, ky_broadening, qe_contrib, time0, time1
    real(kind=dp) :: total_weighted, qe_norm
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp
    integer    :: i, N_k, N_spin, n_eigen_init, n_eigen, ierr, window_width
    integer    :: matrix_unit, nsymm_op, x_center, y_center, xdx, ydx, xdx_min, xdx_max, ydx_min, ydx_max, px_max, py_max
    integer    :: total_ks, xdx_window, ydx_window, ydx_offset, xdx_offset
    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e, char_ref
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+------------ Starting Constant Binding Energy Map Calculation --------------+'
      call flush(stdout)
    end if

    ! get kinetic energy at efermi for reference
    max_e_kinetic = temp_photon_energy - work_function_eff
    ! calculate total number of k-points in the unreduced grid
    total_ks = kpoint_grid_dim(1)*kpoint_grid_dim(2)
    ! reference binding energy level, at which we want to create the map
    ! gauss_e below compares a band energy, so the level is E_F - E_b, which
    ref_level = efermi - photo_const_bindenergy_value
    do i = 1, 2
      step(i) = 0.5_dp/real(kpoint_grid_dim(i), dp)
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do
    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2)) so that
    ! the FWHM = 1/2 the step distance between the kpoints
    kx_broadening = sub_cell_length(1)/(4.70964009_dp)
    ky_broadening = sub_cell_length(2)/(4.70964009_dp)
    ! k_broadening =  sqrt((2*e_mass*(photo_bindenergy_broadening*0.01_dp*ev_to_j))/(hbar*hbar))*1E-10

    ! calculate the number of bins to go left and right
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    xdx_window = ceiling(window_width*kx_broadening/photo_pmat_bin_width)
    ydx_window = ceiling(window_width*ky_broadening/photo_pmat_bin_width)

    call cell_calc_kpoint_r_cart
    z_max = sqrt((2*e_mass*((max_e_kinetic + 0.5)*ev_to_j))/(hbar*hbar))*1E-10
    xy_max = min((abs(maxval(kpoint_r_cart(1:2, :))) + 0.5), z_max)
    xdx_offset = ceiling(xy_max/photo_pmat_bin_width)
    ydx_offset = ceiling(xy_max/photo_pmat_bin_width)

    call comms_reduce(xdx_offset, 1, "MAX")
    call comms_reduce(ydx_offset, 1, "MAX")
    call comms_bcast(xdx_offset, 1)
    call comms_bcast(ydx_offset, 1)

    px_max = 2*xdx_offset + 1
    py_max = 2*ydx_offset + 1

    ! set up the kx x ky matrix
    allocate (kxky_matrix(px_max, py_max), stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map - allocation of kxky_matrix failed')
    kxky_matrix = 0.0_dp

    if (.not. allocated(gauss_x)) then
      allocate (gauss_x(px_max), stat=ierr)
      if (ierr /= 0) call io_error('Error: const_binding_energy_map - allocation of gauss_x failed')
    end if
    gauss_x = 0.0_dp
    if (.not. allocated(gauss_y)) then
      allocate (gauss_y(py_max), stat=ierr)
      if (ierr /= 0) call io_error('Error: const_binding_energy_map - allocation of gauss_y failed')
    end if
    gauss_y = 0.0_dp

    if (.not. allocated(gauss_xy)) then
      allocate (gauss_xy(px_max, py_max), stat=ierr)
      if (ierr /= 0) call io_error('Error: const_binding_energy_map - allocation of gauss_xy failed')
    end if
    gauss_xy = 0.0_dp

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    total_be_contribs = 0.0_dp

    if (index(photo_model, '3step') .gt. 0) then
      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, px_max)
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, py_max)
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          gauss_xy = 0.0_dp
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands
              gauss_e = gaussian(band_energy(n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, ref_level)
              ! The x/y patch and gauss_e carry no atom index, so the sum over
              ! atoms factors out of the accumulation below.
              qe_contrib = sum(qe_tsm(n_eigen_init, 1:nbands, N_spin, N_k, 1:max_atoms))*k_prefactor &
                           *arpes_mask(1, n_eigen_init, N_spin, N_k)
              total_be_contribs = total_be_contribs + qe_contrib
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) = kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) &
                                                              + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_e*qe_contrib
            end do ! bands_init
          end do ! spins
        end do ! kpts
      end do ! symm_ops

      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, px_max)
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, py_max)
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          gauss_xy = 0.0_dp
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands
              gauss_e = gaussian(band_energy(n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, ref_level)
              qe_contrib = sum(qe_tsm(n_eigen_init, 1:nbands, N_spin, N_k, max_atoms + 1))*k_prefactor &
                           *arpes_mask(1, n_eigen_init, N_spin, N_k)
              total_be_contribs = total_be_contribs + qe_contrib
              kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) = kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) &
                                                              + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_e*qe_contrib
            end do ! bands_init
          end do ! spins
        end do ! kpts
      end do ! symm_ops
    end if

    if (index(photo_model, '1step') .gt. 0) then
      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        do N_k = 1, num_kpoints_on_node(my_node_id)
          if (index(devel_flag, 'no_symmetry') .gt. 0) then
            current_k = kpoint_r_cart(1:2, N_k)
            k_prefactor = 1.0_dp/num_crystal_symmetry_operations
          else
            current_k = matmul(temp_mat, kpoint_r_cart(1:2, N_k))
            k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          end if
          ! current_k is the emission direction of this symmetry image, so the
          ! azimuthal acceptance can finally be tested against something real.
          ! At the irreducible k-point it could not: one point stands for its
          ! whole star, and the star spans many azimuths.
          if (.not. phi_accepted(current_k(1), current_k(2))) cycle
          x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
          xdx_min = max(x_center - xdx_window, 1)
          xdx_max = min(x_center + xdx_window, px_max)
          do xdx = xdx_min, xdx_max
            gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
          end do
          y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
          ydx_min = max(y_center - ydx_window, 1)
          ydx_max = min(y_center + ydx_window, py_max)
          do ydx = ydx_min, ydx_max
            gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
          end do
          gauss_xy = 0.0_dp
          do ydx = ydx_min, ydx_max
            do xdx = xdx_min, xdx_max
              gauss_xy(xdx, ydx) = gauss_x(xdx)*gauss_y(ydx)
            end do
          end do
          do N_spin = 1, nspins
            do n_eigen = 1, nbands
              gauss_e = gaussian(band_energy(n_eigen, N_spin, N_k), photo_bindenergy_broadening, ref_level)
              ! No atom index on the x/y patch, so sum the atoms here.
              qe_contrib = sum(qe_osm(n_eigen, N_spin, N_k, 1:max_atoms + 1)) &
                           *arpes_mask(1, n_eigen, N_spin, N_k)*k_prefactor
              total_be_contribs = total_be_contribs + qe_contrib
              if (abs(qe_contrib) .lt. tiny_contribution) cycle
              kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) = kxky_matrix(xdx_min:xdx_max, ydx_min:ydx_max) &
                                                              + gauss_xy(xdx_min:xdx_max, ydx_min:ydx_max)*gauss_e*qe_contrib
            end do ! bands
          end do ! spins
        end do ! kpts
      end do ! symm_ops
    end if

    call comms_reduce(kxky_matrix(1, 1), px_max*py_max, 'SUM')
    call comms_reduce(total_be_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(kxky_matrix(:, :))
      write (*, *) total_weighted
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      kxky_matrix = kxky_matrix*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      write (char_ref, '(F7.2)') photo_const_bindenergy_value
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_ref_'// &
                 trim(adjustl(char_ref))//'_const_map.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)// &
        '_'//trim(adjustl(char_e))//'_ref_'//trim(adjustl(char_ref))//'_const_map.dat'
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)

      call io_date(cdate, ctime)
      write (matrix_unit, '(a66,a11,a4,a9)') '## OptaDOS Photoemission: Printing Constant Binding Energy Map on ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (matrix_unit, '(a24,a8)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg]: ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg]: ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a38,f9.5)') '## Band Energy of States shown [eV] : ', ref_level
      write (matrix_unit, '(a44,f9.5)') '## Reference Energy of Map (E-E_F)   [eV] : ', photo_const_bindenergy_value
      write (matrix_unit, '(a44,f9.5)') '## Momentum bin width               [1/A] : ', photo_pmat_bin_width
      write (matrix_unit, '(a44,f9.5)') '## Binding energy broadening width   [eV] : ', photo_bindenergy_broadening
      write (matrix_unit, '(a46,i10,a3,i10,a2)') '## Matrix Shape                           : ( ', px_max, ' , ', py_max, ' )'

      write (out_string, '(I0,"(1x,",a,")")') px_max, 'ES25.12E3'
      do ydx = 1, py_max
        write (matrix_unit, '('//trim(out_string)//')') (kxky_matrix(xdx, ydx), xdx=1, px_max)
      end do
      close (unit=matrix_unit)
    end if

    deallocate (kxky_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map - deallocation of kxky_matrix failed')
    deallocate (gauss_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map - failed to deallocate gauss_x')
    deallocate (gauss_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map - failed to deallocate gauss_y')
    deallocate (gauss_xy, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map - failed to deallocate gauss_xy')

    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a45,14x,f11.3,a8)') '+ Time to calculate const. binding energy map', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if

  end subroutine const_binding_energy_map

  subroutine const_binding_energy_map_gkgrid
    !*===============================================================================
    ! This subroutine calculates a map of reciprocal space at a specified binding
    ! energy and writes it out to a file. This is the optimised version for the
    ! photo_momentum option to allow supercell calculations.
    ! written by Felix Mildner, after May 2025
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, kpoint_weight, cell_calc_kpoint_r_cart, &
      kpoint_grid_dim, recip_lattice, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_electronic, only: nbands, nspins, electrons_per_state, transmit_prob, photo_gkgrid, elec_read_gk_grid
    use od_parameters, only: photo_model, photo_theta_centre, photo_theta_halfwidth, photo_temperature, &
      photo_phi_centre, photo_phi_halfwidth, &
      photo_momentum, photo_bindenergy_broadening, iprint, photo_pmat_bin_width, optics_geom, optics_qdir, &
      photo_const_bindenergy_value
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast, on_root
    use od_io, only: io_error, io_file_unit, stdout, io_time, io_date, seedname
    use od_constants, only: inv_sqrt_two_pi, kB, rad_to_deg, twopi, ev_to_j, e_mass, hbar
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp), allocatable, dimension(:, :, :, :) :: arpes_mask
    real(kind=dp), allocatable, dimension(:, :, :, :) :: emission_gauss
    real(kind=dp), allocatable, dimension(:, :, :)    :: fermi_dirac
    real(kind=dp), allocatable, dimension(:)          :: gauss_y, gauss_x
    real(kind=dp) :: step(1:2), sub_cell_length(1:2), gauss_e, temp_mat(2, 2), current_k(2), final_fd, xy_max
    real(kind=dp) :: k_prefactor, ref_level, kx_broadening, ky_broadening, qe_contrib, time0, time1
    real(kind=dp) :: sum_final, sum_atoms
    ! Contributions below this never reach the map, so the patch that would
    ! spread them is not built at all.
    real(kind=dp), parameter :: tiny_contribution = 1.0e-30_dp
    real(kind=dp) :: temp_contribution, gk_factor, norm_vac, qe_factor, width, total_weighted, qe_norm
    integer    :: i, N_k, N_spin, n_eigen_init, n_eigen, n_eigen_final, atom, gdx, ierr, window_width
    integer    :: matrix_unit, nsymm_op, x_center, y_center, xdx, ydx, xdx_min, xdx_max, ydx_min, ydx_max, px_max, py_max
    integer    :: total_ks, xdx_window, ydx_window, ydx_offset, xdx_offset
    character(len=100)                          :: out_string
    character(len=99)                           :: filename
    character(len=10)                           :: char_e, char_ref
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+------------ Starting Constant Binding Energy Map Calculation --------------+'
      call flush(stdout)
    end if

    qe_factor = 1.0_dp/(cell_area)
    width = kB*photo_temperature
    norm_vac = inv_sqrt_two_pi/width

    ! get kinetic energy at efermi for reference
    max_e_kinetic = temp_photon_energy - work_function_eff
    total_ks = kpoint_grid_dim(1)*kpoint_grid_dim(2)
    ref_level = temp_photon_energy - work_function_eff - photo_const_bindenergy_value
    do i = 1, 2
      step(i) = 0.5_dp/real(kpoint_grid_dim(i), dp)
      sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
    end do
    ! diagonal distance between MP points in reciprocal space divided by 2*2*sqrt(2*ln(2)) so that
    ! the FWHM = 1/2 the step distance between the kpoints
    kx_broadening = sub_cell_length(1)/(4.70964009_dp)
    ky_broadening = sub_cell_length(2)/(4.70964009_dp)
    ! k_broadening =  sqrt((2*e_mass*(photo_bindenergy_broadening*0.01_dp*ev_to_j))/(hbar*hbar))*1E-10

    ! calculate the number of bins to go left and right
    ! set to 8 standard deviations (width) of a gaussian function
    window_width = 6
    xdx_window = ceiling(window_width*kx_broadening/photo_pmat_bin_width)
    ydx_window = ceiling(window_width*ky_broadening/photo_pmat_bin_width)

    call cell_calc_kpoint_r_cart
    max_e_kinetic = temp_photon_energy - work_function_eff
    xy_max = sqrt((2*e_mass*((max_e_kinetic + 0.5)*ev_to_j))/(hbar*hbar))*1E-10
    xdx_offset = ceiling(xy_max/photo_pmat_bin_width)
    ydx_offset = ceiling(xy_max/photo_pmat_bin_width)

    call elec_read_gk_grid()

    px_max = 2*xdx_offset + 1
    py_max = 2*ydx_offset + 1

    ! set up the kx x ky matrix
    allocate (kxky_matrix(px_max, py_max), stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - allocation of kxky_matrix failed')
    kxky_matrix = 0.0_dp

    allocate (gauss_x(px_max), stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - allocation of gauss_x failed')
    gauss_x = 0.0_dp

    allocate (gauss_y(py_max), stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - allocation of gauss_y failed')
    gauss_y = 0.0_dp

    call prepare_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)
    total_be_contribs = 0.0_dp
    if (index(photo_model, '3step') .gt. 0) then
      call photo_calculate_delta(delta_temp, .false.)
      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        ! current_k = matmul(temp_mat, photo_gkgrid(1:2, gdx, n_eigen_init, N_spin, N_k))

        ! Everything that sets the k_x/k_y patch is indexed by the initial band
        ! and gdx, and the same is true of gauss_e now that E_kinetic takes the
        ! initial band. Neither depends on the atom nor on the final band, so
        ! both sums factor out of the patch and separate from one another:
        !
        !   sum_atom sum_final (contribution * patch) = (sum_atom sum_final contribution) * patch
        !
        ! gauss_e picks a single binding-energy slice a few hundredths of an eV
        ! wide, so it is negligible for most states; testing it before the patch
        ! is built skips the majority of the work outright.
        do N_k = 1, num_kpoints_on_node(my_node_id)
          k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands - 1
              temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                  *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                  *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                  /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
              if (abs(temp_contribution) .lt. tiny_contribution) cycle

              sum_final = 0.0_dp
              do n_eigen_final = n_eigen_init + 1, nbands
                final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
                sum_final = sum_final &
                            + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                            *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                            *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
              end do
              if (abs(sum_final) .lt. tiny_contribution) cycle

              do gdx = 1, photo_gkmax
                ! Hoisted above the running total: this is the emission direction
                ! of the current symmetry image, and an image outside the azimuthal
                ! wedge must contribute to neither the map nor the normalisation.
                current_k = matmul(temp_mat, photo_gkgrid(1:2, gdx, n_eigen_init, N_spin, N_k))
                if (.not. phi_accepted(current_k(1), current_k(2))) cycle
                gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k) &
                            *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                            *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
                if (abs(gk_factor) .lt. tiny_contribution) cycle

                sum_atoms = 0.0_dp
                do atom = 1, max_atoms
                  sum_atoms = sum_atoms &
                              + I_layer(box_atom(atom), current_photo_energy_index) &
                              *pdos_weights_atoms(n_eigen_init, N_spin, N_k, atom_order(atom)) &
                              *electron_esc(gdx, n_eigen_init, N_spin, N_k, atom)
                end do

                qe_contrib = k_prefactor*temp_contribution*sum_final*gk_factor*sum_atoms
                if (abs(qe_contrib) .lt. tiny_contribution) cycle
                ! The running total excludes gauss_e, as it did before: the map
                ! is normalised against the unsliced quantum efficiency.
                total_be_contribs = total_be_contribs + qe_contrib

                gauss_e = gaussian(E_kinetic(gdx, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, ref_level)
                if (abs(gauss_e*qe_contrib) .lt. tiny_contribution) cycle

                current_k = matmul(temp_mat, photo_gkgrid(1:2, gdx, n_eigen_init, N_spin, N_k))
                x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
                xdx_min = max(x_center - xdx_window, 1)
                xdx_max = min(x_center + xdx_window, px_max)
                do xdx = xdx_min, xdx_max
                  gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
                end do
                y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
                ydx_min = max(y_center - ydx_window, 1)
                ydx_max = min(y_center + ydx_window, py_max)
                do ydx = ydx_min, ydx_max
                  gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
                end do
                do ydx = ydx_min, ydx_max
                  do xdx = xdx_min, xdx_max
                    kxky_matrix(xdx, ydx) = kxky_matrix(xdx, ydx) + gauss_x(xdx)*gauss_y(ydx)*gauss_e*qe_contrib
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      call photo_calculate_delta(delta_temp, .true.)

      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        ! Same factorisation for the extrapolated bulk; no atom loop to hoist,
        ! since the bulk is one region indexed max_atoms + 1.
        do N_k = 1, num_kpoints_on_node(my_node_id)
          k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          do N_spin = 1, nspins
            do n_eigen_init = 1, nbands - 1
              temp_contribution = qe_factor*electrons_per_state*kpoint_weight(N_k) &
                                  *fermi_dirac(n_eigen_init, N_spin, N_k) &
                                  *(1.0_dp + field_emission(n_eigen_init, N_spin, N_k)) &
                                  *pdos_weights_boxes(n_eigen_init, N_spin, N_k, num_boxes) &
                                  /pdos_weights_k_band(n_eigen_init, N_spin, N_k)
              if (abs(temp_contribution) .lt. tiny_contribution) cycle

              sum_final = 0.0_dp
              do n_eigen_final = n_eigen_init + 1, nbands
                ! if (num_exclude_bands .gt. 1) then
                !   if (any(exclude_bands == n_eigen_final)) then
                !     cycle
                !   end if
                ! end if
                final_fd = 1 - fermi_dirac(n_eigen_final, N_spin, N_k)
                sum_final = sum_final &
                            + photo_matrix_weights(n_eigen_init, n_eigen_final, N_spin, N_k) &
                            *delta_temp(n_eigen_init, n_eigen_final, N_spin, N_k) &
                            *transmit_prob(n_eigen_final, N_k, N_spin)*final_fd
              end do
              if (abs(sum_final) .lt. tiny_contribution) cycle

              do gdx = 1, photo_gkmax
                ! Hoisted above the running total: this is the emission direction
                ! of the current symmetry image, and an image outside the azimuthal
                ! wedge must contribute to neither the map nor the normalisation.
                if (.not. phi_accepted(current_k(1), current_k(2))) cycle
                gk_factor = arpes_mask(gdx, n_eigen_init, N_spin, N_k) &
                            *gkgrid_weight(gdx, n_eigen_init, N_spin, N_k) &
                            *electron_esc(gdx, n_eigen_init, N_spin, N_k, max_atoms + 1) &
                            *emission_gauss(gdx, n_eigen_init, N_spin, N_k)
                if (abs(gk_factor) .lt. tiny_contribution) cycle

                qe_contrib = temp_contribution*sum_final*gk_factor*k_prefactor
                if (abs(qe_contrib) .lt. tiny_contribution) cycle
                total_be_contribs = total_be_contribs + qe_contrib

                gauss_e = gaussian(E_kinetic(gdx, n_eigen_init, N_spin, N_k), photo_bindenergy_broadening, ref_level)
                if (abs(gauss_e*qe_contrib) .lt. tiny_contribution) cycle

                x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
                xdx_min = max(x_center - xdx_window, 1)
                xdx_max = min(x_center + xdx_window, px_max)
                do xdx = xdx_min, xdx_max
                  gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
                end do
                y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
                ydx_min = max(y_center - ydx_window, 1)
                ydx_max = min(y_center + ydx_window, py_max)
                do ydx = ydx_min, ydx_max
                  gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
                end do
                do ydx = ydx_min, ydx_max
                  do xdx = xdx_min, xdx_max
                    kxky_matrix(xdx, ydx) = kxky_matrix(xdx, ydx) + gauss_x(xdx)*gauss_y(ydx)*gauss_e*qe_contrib
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    if (index(photo_model, '1step') .gt. 0) then
      do nsymm_op = 1, num_crystal_symmetry_operations
        temp_mat = crystal_symmetry_operations(1:2, 1:2, nsymm_op)
        ! One band index here, so only the atom sum comes out of the patch.
        do N_k = 1, num_kpoints_on_node(my_node_id)
          k_prefactor = kpoint_weight(N_k)*total_ks/num_crystal_symmetry_operations
          do N_spin = 1, nspins
            kxkybands: do n_eigen = 1, nbands
              temp_contribution = qe_factor*foptical_matrix_weights(n_eigen, N_spin, N_k) &
                                  *electrons_per_state*kpoint_weight(N_k) &
                                  *fermi_dirac(n_eigen, N_spin, N_k) &
                                  *(1.0_dp + field_emission(n_eigen, N_spin, N_k)) &
                                  /pdos_weights_k_band(n_eigen, N_spin, N_k)
              if (abs(temp_contribution) .lt. tiny_contribution) cycle

              do gdx = 1, photo_gkmax
                ! Hoisted above the running total: this is the emission direction
                ! of the current symmetry image, and an image outside the azimuthal
                ! wedge must contribute to neither the map nor the normalisation.
                current_k = matmul(temp_mat, photo_gkgrid(1:2, gdx, n_eigen, N_spin, N_k))
                if (.not. phi_accepted(current_k(1), current_k(2))) cycle
                gk_factor = arpes_mask(gdx, n_eigen, N_spin, N_k) &
                            *gkgrid_weight(gdx, n_eigen, N_spin, N_k) &
                            *emission_gauss(gdx, n_eigen, N_spin, N_k)
                if (abs(gk_factor) .lt. tiny_contribution) cycle

                sum_atoms = 0.0_dp
                do atom = 1, max_atoms + 1
                  sum_atoms = sum_atoms &
                              + I_layer(box_atom(atom), current_photo_energy_index) &
                              *pdos_weights_atoms(n_eigen, N_spin, N_k, atom_order(atom)) &
                              *electron_esc(gdx, n_eigen, N_spin, N_k, atom)
                end do

                qe_contrib = temp_contribution*gk_factor*sum_atoms*k_prefactor
                if (abs(qe_contrib) .lt. tiny_contribution) cycle
                total_be_contribs = total_be_contribs + qe_contrib

                gauss_e = gaussian(E_kinetic(gdx, n_eigen, N_spin, N_k), photo_bindenergy_broadening, ref_level)
                if (abs(gauss_e*qe_contrib) .lt. tiny_contribution) cycle

                x_center = nint(current_k(1)/photo_pmat_bin_width) + xdx_offset + 1
                xdx_min = max(x_center - xdx_window, 1)
                xdx_max = min(x_center + xdx_window, px_max)
                do xdx = xdx_min, xdx_max
                  gauss_x(xdx) = gaussian(current_k(1), kx_broadening, (xdx - xdx_offset - 1)*photo_pmat_bin_width)
                end do
                y_center = nint(current_k(2)/photo_pmat_bin_width) + ydx_offset + 1
                ydx_min = max(y_center - ydx_window, 1)
                ydx_max = min(y_center + ydx_window, py_max)
                do ydx = ydx_min, ydx_max
                  gauss_y(ydx) = gaussian(current_k(2), ky_broadening, (ydx - ydx_offset - 1)*photo_pmat_bin_width)
                end do
                do ydx = ydx_min, ydx_max
                  do xdx = xdx_min, xdx_max
                    kxky_matrix(xdx, ydx) = kxky_matrix(xdx, ydx) + gauss_x(xdx)*gauss_y(ydx)*gauss_e*qe_contrib
                  end do
                end do
              end do
            end do kxkybands
          end do
        end do
      end do
    end if

    call comms_reduce(kxky_matrix(1, 1), px_max*py_max, 'SUM')
    call comms_reduce(total_be_contribs, 1, 'SUM')

    if (on_root) then
      total_weighted = sum(kxky_matrix(:, :))
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_be_contribs/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      kxky_matrix = kxky_matrix*qe_norm

      write (char_e, '(F7.3)') temp_photon_energy
      write (char_ref, '(F7.2)') photo_const_bindenergy_value
      filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_ref_'// &
                 trim(adjustl(char_ref))//'_const_map.dat'
      write (stdout, '(1x,a)') '| Writing out to *SEEDNAME*_'//trim(photo_model)// &
        '_'//trim(adjustl(char_e))//'_ref_'//trim(adjustl(char_ref))//'_const_map.dat'
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)

      call io_date(cdate, ctime)
      write (matrix_unit, '(a66,a11,a4,a9)') '## OptaDOS Photoemission: Printing Constant Binding Energy Map on ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (matrix_unit, '(a24,a8)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      write (matrix_unit, '(a72,2(1x,f7.2))') '## Emission angle theta centre, half width (w.r.t. surface normal) [deg]: ', &
        photo_theta_centre, photo_theta_halfwidth
      write (matrix_unit, '(a64,2(1x,f7.2))') '## Emission angle phi centre, half width (w.r.t. x-axis) [deg]: ', &
        photo_phi_centre, photo_phi_halfwidth
      write (matrix_unit, '(a44,f9.5)') '## Kinetic Energy of Electrons shown [eV] : ', ref_level
      write (matrix_unit, '(a44,f9.5)') '## Reference Energy of Map (E-E_F)   [eV] : ', photo_const_bindenergy_value
      write (matrix_unit, '(a44,f9.5)') '## Momentum bin width               [1/A] : ', photo_pmat_bin_width
      write (matrix_unit, '(a44,f9.5)') '## Binding energy broadening width   [eV] : ', photo_bindenergy_broadening
      write (matrix_unit, '(a46,i10,a3,i10,a2)') '## Matrix Shape                           : ( ', px_max, ' , ', py_max, ' )'

      write (out_string, '(I0,"(1x,",a,")")') px_max, 'ES25.12E3'
      do ydx = 1, py_max
        write (matrix_unit, '('//trim(out_string)//')') (kxky_matrix(xdx, ydx), xdx=1, px_max)
      end do
      close (unit=matrix_unit)
    end if

    deallocate (kxky_matrix, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - deallocation of kxky_matrix failed')
    deallocate (gauss_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - failed to deallocate gauss_x')
    deallocate (gauss_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - failed to deallocate gauss_y')
    deallocate (photo_gkgrid, stat=ierr)
    if (ierr /= 0) call io_error('Error: const_binding_energy_map_gkgrid - failed to deallocate photo_gkgrid')
    call deallocate_emission_arrays(fermi_dirac, arpes_mask, emission_gauss)

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a45,14x,f11.3,a8)') '+ Time to calculate const. binding energy map', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      call flush(stdout)
    end if

  end subroutine const_binding_energy_map_gkgrid

  subroutine write_qe_tensor
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
    use od_parameters, only: photo_model, photo_momentum, iprint, devel_flag, optics_geom, optics_qdir
    implicit none

    integer :: atom, matrix_unit
    integer :: N_k, N_spin, n_eigen, kpt_total, band_num
    character(len=99)                           :: filename
    character(len=100)                          :: out_string
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    real(kind=dp) :: time0, time1

    time0 = io_time()

    call cell_calc_kpoint_r_cart
    kpt_total = sum(num_kpoints_on_node(0:num_nodes - 1))
    if (num_nodes .gt. 1) then
      call write_distributed_qe_data(kpt_total)
    else
      matrix_unit = io_file_unit()
      write (char_e, '(F7.3)') temp_photon_energy
      if (index(devel_flag, 'final') .gt. 0 .and. index(photo_model, '3step') .gt. 0) then
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_tensor_final.dat'
      else
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_tensor.dat'
      end if
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)
      write (matrix_unit, '(a53,a11,a4,a9)') '## OptaDOS Photoemission: Printing Full QE tensor on ',&
      & cdate, ' at ', ctime
      write (matrix_unit, '(a13,a)') '## Seedname: ', trim(adjustl(seedname))
      write (matrix_unit, '(a24,a12)') '## Photoemission Model: ', trim(adjustl(photo_model))
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, '(a23,f7.3)') '## Photon Energy [eV]: ', temp_photon_energy
      write (matrix_unit, '(a21,a15)') '## Optics Geometry : ', trim(adjustl(optics_geom))
      write (matrix_unit, '(a39,3(1x,f10.5))') '## Optics q-dir vector [unnormalised] :', optics_qdir(1:3)
      if (index(devel_flag, 'final') .gt. 0 .and. index(photo_model, '3step') .gt. 0) then
        write (matrix_unit, '(a69)') '## Writing the contributions of excitations into the !!FINAL!! states'
      end if
      write (matrix_unit, '(a61,a,a6)') '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'
      ! Printing out the info on root_node
      write (out_string, '(I0,"(1x,",a,")")') nbands, 'ES16.8E3'

      if (index(photo_model, '3step') .gt. 0) then
        if (index(devel_flag, 'single') .gt. 0) then
          n_eigen = len_trim(devel_flag)
          read (devel_flag(n_eigen - 2:n_eigen), *) band_num
          write (matrix_unit, '(a42,1x,I3)') '## Writing contributions into final band #', band_num
          write (matrix_unit, '(a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                       & ')'
          do atom = 1, max_atoms + 1
            if (atom .eq. max_atoms + 1) write (matrix_unit, '(a21)') '## Bulk Contribution:'
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '('//trim(out_string)//')') &
                  (qe_tsm(n_eigen, band_num, N_spin, N_k, atom), n_eigen=1, nbands)
              end do
            end do
          end do
        else if (index(devel_flag, 'final') .gt. 0) then
          write (matrix_unit, '(a79)') '## (Reduced) QE Matrix where each row contains the contributions from each band'
          write (matrix_unit, '(a39)') '## at a certain k-point, spin, and atom'
          write (matrix_unit, '(a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                       & ')'
          do atom = 1, max_atoms + 1
            if (atom .eq. max_atoms + 1) write (matrix_unit, '(a21)') '## Bulk Contribution:'
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '('//trim(out_string)//')') &
                  (sum(qe_tsm(1:nbands, n_eigen, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          end do
        else
          write (matrix_unit, '(a79)') '## (Reduced) QE Matrix where each row contains the contributions from each band'
          write (matrix_unit, '(a39)') '## at a certain k-point, spin, and atom'
          write (matrix_unit, '(a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                                       & ')'
          do atom = 1, max_atoms + 1
            if (atom .eq. max_atoms + 1) write (matrix_unit, '(a21)') '## Bulk Contribution:'
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '('//trim(out_string)//')') &
                  (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          end do
        end if
      elseif (index(photo_model, '1step') .gt. 0) then
        write (matrix_unit, '(a79)') '## (Reduced) QE Matrix where each row contains the contributions from each band'
        write (matrix_unit, '(a39)') '## at a certain k-point, spin, and atom'
        write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                                  & ')'
        do atom = 1, max_atoms + 1
          if (atom .eq. max_atoms + 1) write (matrix_unit, '(a21)') '## Bulk Contribution:'
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (matrix_unit, '('//trim(out_string)//')') (qe_osm(n_eigen, N_spin, N_k, atom), n_eigen=1, nbands)
            end do
          end do
        end do
      end if
      close (unit=matrix_unit)
    end if

    time1 = io_time()
    if (on_root .and. iprint .gt. 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a37,22x,f11.3,a8)') '+ Time to write the qe tensor to file', time1 - time0, ' (sec) +'
    end if

  end subroutine write_qe_tensor

  subroutine write_distributed_qe_data(kpt_total)
    !* This subroutine writes the distributed qe tensor to a single file.
    ! To save on required memory the output file is accessed by each MPI process in turn
    ! and writes its values/contents one after the other.
    ! F. Mildner, June 2023
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nspins, nbands
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_bcast
    use od_io, only: io_error, io_file_unit, io_date, io_time, seedname
    use od_parameters, only: photo_model, photo_momentum, devel_flag

    implicit none
    real(kind=dp), dimension(:, :, :), allocatable :: qe_mat_temp
    real(kind=dp), dimension(:, :, :, :), allocatable :: tsm_reduced
    integer, intent(in)                         :: kpt_total
    character(len=99)                           :: filename
    character(len=100)                          :: out_string
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    integer:: N_k, N_spin, n_eigen, atom, token, matrix_unit, ierr, inode

    ! On root open file and write header
    if (on_root) then
      ! Writing header to output file
      write (char_e, '(F7.3)') temp_photon_energy
      if (index(devel_flag, 'final') .gt. 0 .and. index(photo_model, '3step') .gt. 0) then
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_tensor_final.dat'
      else
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_tensor.dat'
      end if
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)
      write (matrix_unit, *) '## OptaDOS Photoemission: Printing QE Matrix on ', cdate, ' at ', ctime
      write (matrix_unit, *) '## Seedname: ', trim(seedname)
      write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
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
    write (out_string, '(I0,"(1x,",a,")")') nbands, 'ES16.8E3'

    ! allocate and sum the 3step qe matrix on non-root
    if (.not. on_root) then
      if (index(photo_model, '3step') .gt. 0) then
        allocate (tsm_reduced(nbands, nspins, num_kpoints_on_node(0), max_atoms + 1), stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate tsm_reduced')
        if (index(devel_flag, 'final') .gt. 0) then
          tsm_reduced = sum(qe_tsm, dim=1)
        else
          tsm_reduced = sum(qe_tsm, dim=2)
        end if
      end if
    end if
    ! For each atom until max_atoms+1
    do atom = 1, max_atoms + 1
      ! On non root nodes
      if (.not. on_root) then
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that specific atom
        if (index(photo_model, '3step') .gt. 0) then
          call comms_send(tsm_reduced(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
        elseif (index(photo_model, '1step') .gt. 0) then
          call comms_send(qe_osm(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
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
              write (matrix_unit, '('//trim(out_string)//')') (qe_mat_temp(n_eigen, N_spin, N_k), n_eigen=1, nbands)
            end do
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        if (index(photo_model, '3step') .gt. 0) then
          if (index(devel_flag, 'final') .gt. 0) then
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '('//trim(out_string)//')') &
                  (sum(qe_tsm(1:nbands, n_eigen, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          else
            do N_k = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '('//trim(out_string)//')') &
                  (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N_k, atom)), n_eigen=1, nbands)
              end do
            end do
          end if
        elseif (index(photo_model, '1step') .gt. 0) then
          do N_k = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (matrix_unit, '('//trim(out_string)//')') &
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
      if (index(photo_model, '3step') .gt. 0) then
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
    use od_parameters, only: photo_model, photo_momentum

    implicit none
    real(kind=dp), dimension(:, :, :), allocatable :: fem_mat_temp
    ! real(kind=dp), dimension(:, :, :, :), allocatable :: tsm_reduced
    integer, intent(in)                         :: kpt_total
    character(len=99)                           :: filename
    character(len=100)                          :: out_string
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
      write (matrix_unit, '(a31,a12)') '## Transverse Momentum Model : ', trim(adjustl(photo_momentum))
      write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
      write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'
      write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
      write (matrix_unit, *) '## at a certain k-point, spin, and atom'
      write (matrix_unit, '(1x,a31,3(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, kpt_total, nspins, ')'
      allocate (fem_mat_temp(nbands, num_kpoints_on_node(0), nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate fem_mat_temp on root')
      token = -1
    end if
    write (out_string, '(I0,"(1x,",a,")")') nbands, 'ES16.8E3'

    ! On non root nodes
    if (.not. on_root) then
      ! - wait for the token
      call comms_recv(token, 1, 0)
      ! - send the respective qe_matrix for that specific atom
      call comms_send(foptical_matrix_weights(1, 1, 1), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
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
            write (matrix_unit, '('//trim(out_string)//')') (fem_mat_temp(n_eigen, N_k, N_spin), n_eigen=1, nbands)
          end do
        end do
        ! - receive the token from a node
        call comms_recv(token, 1, inode)
      end do
      ! - write root qe_matrix elements
      do N_k = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          write (matrix_unit, '('//trim(out_string)//')') (foptical_matrix_weights(n_eigen, N_spin, N_k), n_eigen=1, nbands)
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

    if (allocated(phi_accept_frac)) then
      deallocate (phi_accept_frac, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate phi_accept_frac')
    end if

    if (allocated(theta_arpes)) then
      deallocate (theta_arpes, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate theta_arpes')
    end if

    if (allocated(theta_internal)) then
      deallocate (theta_internal, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate theta_internal')
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

    if (allocated(photo_matrix_weights)) then
      deallocate (photo_matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate photo_matrix_weights')
    end if

    if (allocated(E_transverse)) then
      deallocate (E_transverse, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate E_transverse')
    end if

    if (allocated(absorp_photo)) then
      deallocate (absorp_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate absorp_photo')
    end if

    if (allocated(layer_gap)) then
      deallocate (layer_gap, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate layer_gap')
    end if

    if (allocated(atom_order)) then
      deallocate (atom_order, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate atom_order')
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

    if (allocated(E_kinetic)) then
      deallocate (E_kinetic, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate E_kinetic')
    end if

    if (allocated(field_emission)) then
      deallocate (field_emission, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate field_emission')
    end if

    if (allocated(qe_tsm)) then
      deallocate (qe_tsm, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate qe_tsm')
    end if

    if (allocated(ds_qe_den)) then
      deallocate (ds_qe_den, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate ds_qe_den')
      deallocate (ds_qe_num, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate ds_qe_num')
      deallocate (ds_mte_num, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate ds_mte_num')
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

    if (allocated(gkgrid_weight)) then
      deallocate (gkgrid_weight, stat=ierr)
      if (ierr /= 0) call io_error('Error: const_binding_energy_map - failed to deallocate gkgrid_weight')
    end if

  end subroutine photo_deallocate

end module od_photo
