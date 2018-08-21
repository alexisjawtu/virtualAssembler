[K,F,num_fc,num_el] 	= assemble();
det(K)
U 		= K\F;
U 		= U(1:num_fc,1);

SEGUIR ACA: ver solsh2.m y ver errores2.m