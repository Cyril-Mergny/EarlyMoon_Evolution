
MODULE orbit_var
	use prec
	
	
	
	
	REAL(KIND=dp), DIMENSION(2) :: klove
	
	
END MODULE orbit_var

PROGRAM main_test
	use orbit
	use variables
	use tides, only : k2_out

	
	
	klove = k2_out
	
	
	Pdiss = 1e12_dp
	
	print *,e,a/Rearth,Pspin
	
	CALL calc_orbitalEvolution(ecc, w_in, beta, xincl, Pspin, timestep)
	
	print *,e,a/Rearth,Pspin
	
END PROGRAM main_test


