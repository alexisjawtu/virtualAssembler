%% assmbl_prism: function description
function [outputs] = assmbl_prism(face_indices,face_types,Mdistances,det_M_Element,vol_pts,face_pts,P0,normalFacesE,measFacesE)
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
  %%   face_quad_coef =    |1 |*(1/36) in case of rectangle
  %                        |1 |
  %                        |1 |
  %                        |1 |
  %                        |4 |
  %                        |4 |
  %                        |4 |
  %                        |4 |
  %                        |16|	

  clear('Mdistances');

  face_quad_coef              = repmat([1;1;1;1;4;4;4;4;16],1,5);
  quadrilat                   = find(face_types == 4);
  triangles                   = find(face_types == 3);    
  %% TODO: instead of measFacesE(triangles) do: measFacesE(map(triangles)) 
  %% so that we don't care of "preserving orders"
  face_quad_coef(:,quadrilat) = face_quad_coef(:,quadrilat)*diag(measFacesE(quadrilat)/36);
  face_quad_coef(:,triangles) = [repmat([measFacesE(triangles)/3],3,1);zeros(6,2)];

  measE                       = abs(det_M_Element)/3;

  quad_nrmlztn                = measE/2;    
    we_potentials       = zeros(8, 4);
    we_basis                = zeros(3, 4, 8);
    for l = 1:8
      we_potentials(l,:) = WE_potentials(vol_pts(:,l),5);
      we_basis(:,:,l)    = WE_basis(vol_pts(:,l),5);   
    end                                                                    
    
    int_E_w_w = zeros(4);   %% int_E < w_k; w_r > dx
    for r = 1:4
      for k = 1:4
        %% ( (wr*wk)(p1), ..., (wr*wk)(p8) )
        vals            = reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,8);
        int_E_w_w(r,k)  = vals*vol_wg;
      end
    end
    
    int_E_w_w *= quad_nrmlztn; 
    %% READ: -quad_nrmlztn*we_potentials.'*vol_wg/measE  
    %% -1/mu(E) * int_E q_j
    b1 = repmat(-we_potentials.'*vol_wg/2,1,5);
    % Boundary term. Formerly: b2_{i,j} = Int_{fi} qj dS. page 62 in the middle.
    b2 = zeros(4,5); 
    we_potentials_faces = zeros (4,9,5);
    for f = 1:5
      for pt = 1:size(face_pts,2) 
        we_potentials_faces(:,pt,f) = WE_potentials(face_pts{face_indices[f+1]}(:,pt),5).';
      end
    end
    %% now the evaluated arrays have 1 to 5 face numbers.
    for f = 1:5
      b2(:,f) = we_potentials_faces(:,:,f) * face_quad_coef(face_indices(f+1)); %% (4x1)
    end
    b2        = b2./repmat(measFacesE,4,1);
    b         = b1 + b2;
    PROJ      = int_E_w_w\b;
    K_comput  = PROJ.' * int_E_w_w * PROJ;
    outputs   = K_comput;
endfunction