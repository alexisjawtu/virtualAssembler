%% assmble_pyram: function description
%% normalFacesE: the last is the rectangle
function [outputs] = assmble_pyram(vertices,faces_of_E,face_types,face_pts,normalFacesE,measFacesE)
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
    %
    % coefficients tetrahedral cubature on 4 points A. H. Stroud, Approximate 
    % Calculation of Multiple Integrals, Prentice-Hall, Englewood Cliffs, NJ, 1971
    % We are assuming that the pyramid is divided into two tetrahera of the same measure.
    % When this stops being the case, we have to write meas(T1)/4 * sum + meas(T2)/4 * sum.

    %%% old 
    %  face_quad_coef              = repmat([1;1;1;1;4;4;4;4;16],1,5);
    %  quadrilat                   = find(face_types == 4);
    %  triangles                   = find(face_types == 3);    
    %% TODO: instead of measFacesE(triangles) do: measFacesE(map(triangles)) 
    %% so that we don't care of "preserving orders"
    %  face_quad_coef(:,quadrilat) = face_quad_coef(:,quadrilat)*diag(measFacesE(quadrilat)/36);
    %  face_quad_coef(:,triangles) = [repmat(measFacesE(triangles)/3,3,1);zeros(6,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_vertices  = 5;
    nFaces      = 5;
    dim_Wh      = 4;
    n_vol_pts   = 8;
    vol_pts     = zeros(3,n_vol_pts); 
    const_a     = .58541019662496852;
    const_b     = .1381966011250105;
    vol_weights = .25*ones(8,1); % the weights for the volume quadrature points

    n_face_pts = [9,3,3,3,3]; %% Number of quad pts per face

    %% vertices(:,6) is the top of the pyramid
    if ~are_coplanar (vertices(:,1:4))
      fprintf('pyramid %d with vertices in different order', el);
      exit;
    end
    %% sum of 2 tetrahedral cubatures from A. H. Stroud, Approximate 
    %% Calculation of Multiple Integrals, Prentice-Hall, Englewood Cliffs, NJ, 1971
    %
    %% tetrahedron: p_0,p_1,p_2,p_4 == v1,v2,v3,v5
    % cubature points
    vol_pts(:,1) = const_a*vertices(:,1) + const_b*sum(vertices(:,[2,3,5]),2);
    vol_pts(:,2) = const_a*vertices(:,2) + const_b*sum(vertices(:,[1,3,5]),2);
    vol_pts(:,3) = const_a*vertices(:,3) + const_b*sum(vertices(:,[1,2,5]),2);
    vol_pts(:,4) = const_a*vertices(:,5) + const_b*sum(vertices(:,[1,2,3]),2);
    %% tetrahedron: vertices(:,1) vertices(:,3) vertices(:,4) vertices(:,5)
    % cubature points
    vol_pts(:,5) = const_a*vertices(:,1) + const_b*sum(vertices(:,[3,4,5]),2);
    vol_pts(:,6) = const_a*vertices(:,3) + const_b*sum(vertices(:,[1,4,5]),2);
    vol_pts(:,7) = const_a*vertices(:,4) + const_b*sum(vertices(:,[1,3,5]),2);
    vol_pts(:,8) = const_a*vertices(:,5) + const_b*sum(vertices(:,[1,3,4]),2);

    Mdistances = [vertices(:,2)-vertices(:,1), vertices(:,1)-vertices(:,3),...
                  vertices(:,4)-vertices(:,1), vertices(:,5)-vertices(:,1),...
                  vertices(:,2)-vertices(:,3), vertices(:,2)-vertices(:,4),...
                  vertices(:,5)-vertices(:,2), vertices(:,3)-vertices(:,4),...
                  vertices(:,3)-vertices(:,5), vertices(:,4)-vertices(:,5)];
   
    M_Element = Mdistances(:, [3, 1, 4]);

    face_quad_coef = {};
    %% after testing
    for index = 1:5
        face_quad_coef(faces_of_E(index)) = repmat(measFacesE(index)/3,3,1);
    end
    %% after testing
    face_quad_coef(faces_of_E(find(face_types == 4))) = [1;1;1;1;4;4;4;4;16]*measFacesE(find(face_types == 4))/36;

    quad_nrmlztn   = abs(det(M_Element))/24; %% <--- measE/8;    
    rescale_factor = 1/max(norm(Mdistances,2,'columns')); %% 1/h_E
    clear('Mdistances');

    we_potentials = zeros(n_vol_pts, dim_Wh);
    we_basis      = zeros(3, dim_Wh, n_vol_pts);

    for l = 1:n_vol_pts
      we_potentials(l,:) = WE_potentials(vol_pts(:,l),n_vertices);
      we_basis(:,:,l)    = WE_basis(vol_pts(:,l),n_vertices);   
    end                                                                    
    
    int_E_w_w = zeros(dim_Wh);   %% int_E < w_k; w_r > dx

    for r = 1:dim_Wh
      for k = 1:dim_Wh
        %% ( (wr*wk)(p1), ..., (wr*wk)(p8) )
        int_E_w_w(r,k) = sum(reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts));
      end
    end
    
    int_E_w_w *= quad_nrmlztn; 

	%%%%%%%%%%%%%%%%% for the COMPUTABILITY bilinear form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% after testing
    %% READ: -const * we_potentials.'*vol_weights/measE  
    %%                    -1/meas(E) * int_E q_j, because div(vk) == 1/meas(E), for every k = 1 .. 5
    b1 = repmat(-we_potentials.'*vol_weights/2,1,nFaces);
    % Boundary term. Formerly: b2_{i,j} = Int_{fi} qj dS. page 62 in the middle.
    %% 1/meas(fk) * ( iint_{fk} potential_i dS ) == iint_{partial E} ( (vk.n)*potential_i ) dS
    b2 = zeros(dim_Wh,nFaces); 
    we_potentials_faces = zeros (dim_Wh,max(n_face_pts),nFaces);
    for f = 1:nFaces
      for pt = 1:size(face_pts,2) 
        we_potentials_faces(:,pt,f) = WE_potentials(face_pts{faces_of_E(f)}(:,pt),n_vertices).';
      end
    end
    %% now the evaluated arrays have 1 to 5 face numbers.
    for f = 1:nFaces
      b2(:,f) = we_potentials_faces(:,:,f) * face_quad_coef(faces_of_E(f)); %% (4x1)
    end

    b2       = b2./repmat(measFacesE,dim_Wh,1);
    b        = b1 + b2;
    PROJ     = int_E_w_w\b;
    K_comput = PROJ.' * int_E_w_w * PROJ;

	%%%%%%%%%%%%%%%%% for the STABILIZING bilinear form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % after testing
    % D_{j,i} = dof_j ( w_i ).
    D         = zeros(nFaces,dim_Wh);
        
    for face = 1:nFaces
      for pts = 1:n_face_pts(face)
        w_on_faces  = WE_basis (diameter, centroid, face_pts(:, pts, face), n_VERT);
        D(face,:)  += face_quad_coef(pts,face) * normalFacesE(:,face).' * w_on_faces;
      end
    end
    Proj_in_base_W  = D*PROJ;

    K_stab          = rescale_factor * (eye(nFaces) - Proj_in_base_W).'*(eye(nFaces) - Proj_in_base_W);
	outputs         = K_comput + K_stab;

endfunction