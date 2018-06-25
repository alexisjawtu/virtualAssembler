%% assmbl_pyram: function description
%% normalFacesE: the last is the rectangle
function [outputs] = assmbl_pyram(face_indices,face_types,Mdistances,det_M_Element,vol_pts,face_pts,P0,normalFacesE,measFacesE)
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
    %%   face_quad_coef =    (1/36)*|1 | in case of rectangle
    %                               |1 |
    %                               |1 |
    %                               |1 |
    %                               |4 |
    %                               |4 |
    %                               |4 |
    %                               |4 |
    %                               |16|

    vol_wg 		   = .25*ones(8,1); % the weights for the volume quadrature points
    face_quad_coef = {};
    for index = 1:5
        face_quad_coef(face_indices(index)) = repmat(measFacesE(index)/3,3,1);
    end
    face_quad_coef(face_indices(find(face_types == 4))) = [1;1;1;1;4;4;4;4;16]*measFacesE(find(face_types == 4))/36;

    measE                       = abs(det_M_Element)/3;
    quad_nrmlztn                = measE/2;    
    rescale_factor              = 1/max(norm(Mdistances,2,'columns')); %% 1/diameter
    clear('Mdistances');

    we_potentials  			= zeros(8, 4);
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

    b2    = b2./repmat(measFacesE,4,1);
    b     = b1 + b2;
    PROJ  = int_E_w_w\b;

    K_comput  = PROJ.' * int_E_w_w * PROJ;

    %% for the stabilizating bilinear form

    % D_{j,i} = dof_j ( w_i ). See page 64 in the middle.
    D         = zeros(5,4);
        
    for face = 1:5
      for pts = 1:n_face_pts{n_VERT}(face)
        w_on_faces  = WE_basis (diameter, centroid, face_pts(:, pts, face), n_VERT);
        D(face,:)  += face_quad_coef(pts,face) * normalFacesE(:,face).' * w_on_faces;
      end
    end

    Proj_in_base_W  = D*PROJ;
    K_stab          = rescale_factor * (eye(5) - Proj_in_base_W).'*(eye(5) - Proj_in_base_W);
	outputs         = K_comput + K_stab;

endfunction