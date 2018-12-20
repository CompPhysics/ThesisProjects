SUBROUTINE set_outputfiles(outputfile)
    USE constants
    CHARACTER (LEN=100), INTENT(IN) :: outputfile

    OPEN(unit=6,file=outputfile)

END SUBROUTINE set_outputfiles

SUBROUTINE set_inputfiles(renormfile, spdatafile)
    USE constants
    CHARACTER (LEN=100), INTENT(IN) :: renormfile, spdatafile

    OPEN(UNIT=7,FILE=renormfile)
    OPEN(UNIT=12,FILE=spdatafile)

END SUBROUTINE set_inputfiles

SUBROUTINE set_interaction(interaction, order)
    USE constants
    CHARACTER (LEN=100), INTENT(IN) :: interaction, order

    order_of_interaction = order
    type_of_interaction = interaction

END SUBROUTINE set_interaction

SUBROUTINE set_hf_iterations(hf)
    USE constants
    INTEGER, INTENT(IN) :: hf

    hf_iterations = hf

END SUBROUTINE set_hf_iterations

SUBROUTINE set_hffiles(hf_renorm, hf_spdata)
    CHARACTER (LEN=100), INTENT(IN) :: hf_renorm, hf_spdata

    OPEN(UNIT=8,FILE=hf_spdata)
    OPEN(UNIT=9,FILE=hf_renorm)

END SUBROUTINE set_hffiles
    
SUBROUTINE setup()
    USE constants
    USE ang_mom_functions
    
    keep_originalenergies = "no"
    CALL read_data
    CALL setup_ho_cutoff
    CALL commons_to_angmom

END SUBROUTINE setup

SUBROUTINE set_energies(num, energy)
    USE constants
    REAL(KIND(1.0D0)), INTENT(IN) :: energy
    INTEGER, INTENT(IN) :: num
    INTEGER :: i

    n_startenergy_veff = num
    starting_energy = energy

    IF ( MOD(n_startenergy_veff,2) == 0)n_startenergy_veff = &
        n_startenergy_veff+1

    ALLOCATE ( wcn (n_startenergy_veff) )
    IF ( n_startenergy_veff == 1) THEN
        wcn(1) = 0.0_dp
    ELSE
        DO i=1, n_startenergy_veff
            wcn(i)=-1.*i+(n_startenergy_veff/2+1)
        ENDDO
    ENDIF

END SUBROUTINE

SUBROUTINE set_hyperons(l, sp, s0, sm)
    USE constants
    INTEGER, INTENT(IN) :: l, sp, s0, sm

    n_lambda = l
    n_sigmap = sp
    n_sigma0 = s0
    n_sigmam = sm

    sp_mass(1) = proton_mass
    sp_mass(2) = neutron_mass
    sp_mass(3) = lambda_mass
    sp_mass(4) = sigmap_mass
    sp_mass(5) = sigma0_mass
    sp_mass(6) = sigmam_mass

    total_mass = l*lambda_mass + sp*sigmap_mass + s0*sigma0_mass &
        + sm*sigmam_mass
END SUBROUTINE set_hyperons

SUBROUTINE do_hartree_fock
  
  CALL setup_hfmatrix
END SUBROUTINE do_hartree_fock

SUBROUTINE do_onebody

    CALL onebody_contribution
END SUBROUTINE do_onebody

