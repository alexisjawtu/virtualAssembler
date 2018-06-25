%% assmbl_tetr: function description
function [outputs] = assmbl_tetr(face_indices,face_types,Mdistances,det_M_Element,vol_pts,face_pts,P0,normalFacesE,measFacesE)
  %% for now, Vh(tetra) == Wh(tetra)
  %
  %% measFacesE in R(1x5)
  %
  %% face_pts is a cell. face_pts comes with five GLOBAL UNIQUE face indices,
  %% those of elements_by_faces.txt
  %
  % TODO: hacer un benchmark de inicializar pts_of_faces cada vez vs. pasarlo a la funcion
  %
  %%   face_quad_coef =    |meas(fi)/3| in case of triangle
  %                        |meas(fi)/3|
  %                        |meas(fi)/3|
  %
  

ACA: modificar este codigo. prisma ---> tetra

  face_quad_coef              = repmat([1;1;1;1;4;4;4;4;16],1,5);
  quadrilat                   = find(face_types == 4);
  triangles                   = find(face_types == 3);    
  %% TODO: instead of measFacesE(triangles) do: measFacesE(map(triangles)) 
  %% so that we don't care of "preserving orders"
  face_quad_coef(:,quadrilat) = face_quad_coef(:,quadrilat)*diag(measFacesE(quadrilat)/36);
  face_quad_coef(:,triangles) = [repmat([measFacesE(triangles)/3],3,1);zeros(6,2)];

  measE                       = abs(det_M_Element)/3;
  quad_nrmlztn                = measE/2;    

  we_potentials               = zeros(8, 4);
  we_basis                    = zeros(3, 4, 8);
  for l = 1:8
    we_potentials(l,:) = WE_potentials(vol_pts(:,l),5);
    we_basis(:,:,l)    = WE_basis(vol_pts(:,l),5);   
  end                                                                    
  
  int_E_w_w = zeros(5);   %% int_E < w_k; w_r > dx
  for r = 1:5
    for k = 1:5
      %% ( (wr*wk)(p1), ..., (wr*wk)(p8) )
      vals            = reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,8);
      int_E_w_w(r,k)  = vals*vol_wg;
    end
  end
  
  int_E_w_w *= quad_nrmlztn; 
  %% READ: -quad_nrmlztn*we_potentials.'*vol_wg/measE  
  %% -1/mu(E) * int_E q_j
  outputs   = int_E_w_w;
endfunction
