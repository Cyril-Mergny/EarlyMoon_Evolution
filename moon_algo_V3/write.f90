SUBROUTINE write_structure
    !write structure profile
    !in profile.txt file
    use variables

    open(11,file='Outputs/structure.txt',position="append")
    do i =1 ,n
        write(11,"(I7,F15.3,F10.1,F10.3,F10.3,ES16.8,ES16.8)") i, r_in(i), rho_in(i), vp_in(i), vs_in(i), eta_in(i),D_in(i)
        enddo
	close(11)
END SUBROUTINE write_structure


SUBROUTINE write_Hprofile
    !write Htide profile
    use variables
    use param_y, only : htide
    integer :: i

    open(12,file='Outputs/H_profile.txt',position="append")
    do i=1,n
        write(12,*) real(htide(i))
        enddo
        close(12)
END SUBROUTINE write_Hprofile


SUBROUTINE write_Tcrustprofile
    !Write Temperature
    use variables
    integer :: i
    open(16,file="Outputs/Tcrust_profile.txt",position="append")
    do i = n_lim(4),n
        write(16,*) T(i)
        enddo
	write(16,*) "-"
	close(16)
        
END SUBROUTINE write_Tcrustprofile

SUBROUTINE write_radius
    !Write Temperature
    use variables, only : r_in,n,ncond,n_lim
    integer :: i
	
    open(17,file="Outputs/radius.txt",position="append")
    do i = n_lim(4),n
        write(17,*) r_in(i)
        enddo
	write(17,*) "-"
    close(17)
END SUBROUTINE write_radius


SUBROUTINE write_rlim
    !write
    use variables

    open(20,file='Outputs/rlim_evol.txt',position="append")
	write(20,*) Rlim(4),time,ecc,a,real(klove2),imag(klove2),mean_melt
	close(20)
        
END SUBROUTINE write_rlim



!SUBROUTINE write_power
!    !write
!    use variables
!
!    open(19,file='Outputs/power_e_a.txt',position="append")
!        write(19,*) ecc,2._dp*pi/w_in,Pcrust,Pmant,Pbase
!        close(19)
!END SUBROUTINE write_power
