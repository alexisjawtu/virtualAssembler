%% assmbl_tetr: function description
function [outputs] = assmbl_tetr(vertices,face_indices,face_types,face_pts,normalFacesE,measFacesE)
  %% for now, Vh(tetra) == Wh(tetra)
  %
  %% measFacesE in R(1x5)
  %
  %  TODO: hacer un benchmark de inicializar pts_of_faces cada vez vs. 
  %  pasarlo a la funcion.
  %
  %% face_quad_coef =    |meas(fi)/3| in case of triangle
  %                      |meas(fi)/3|
  %                      |meas(fi)/3|
  %

  dim_Vh 			= 4;
  n_vertices 		= 4;
  n_vol_pts 		= 4;
  vol_weights 		= .25*ones(n_vol_pts,1);
  face_quad_coef 	= repmat(measFacesE/3,3,1);
  M_Element    = [vertices(:,2)-vertices(:,1), vertices(:,3)-vertices(:,1), vertices(:,4)-vertices(:,1)];
  det_M_Element = det(M_Element);
  measE           = abs(det_M_Element)/3;
  quad_nrmlztn    = measE;   ***es este o meas/2 como en pyr?    
  we_basis          = zeros(3, dim_Vh, n_vol_pts);
  % coefficients tetrahedral cubature on 4 points GELLERT and HARBORD 91
  const_a     = .58541019662496852;
  const_b     = .1381966011250105;
  
  % tetrahedral cubature points from GELLERT AND HARBORD 1991
  % direct from the physical vertices      
  vol_pts(:,1) = const_a*vertices(:,1) + const_b*sum(vertices(:,[3,4,5]),2);
  vol_pts(:,2) = const_a*vertices(:,2) + const_b*sum(vertices(:,[2,4,5]),2);
  vol_pts(:,3) = const_a*vertices(:,3) + const_b*sum(vertices(:,[2,3,5]),2);
  vol_pts(:,4) = const_a*vertices(:,4) + const_b*sum(vertices(:,[2,3,4]),2);


  for l = 1:n_vol_pts
    we_basis(:,:,l)    = WE_basis(vol_pts(:,l),n_vertices);   
  end                                                                    
  
  int_E_w_w = zeros(dim_Vh);   %% int_E < w_k; w_r > dx
  for r = 1:dim_Vh
    for k = 1:dim_Vh
      %% ( (wr*wk)(p1), ..., (wr*wk)(p4) )
      vals            = reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts);
      int_E_w_w(r,k)  = vals*vol_weights;
    end
  end
  
  int_E_w_w *= quad_nrmlztn; 
  outputs   = int_E_w_w;
endfunction
