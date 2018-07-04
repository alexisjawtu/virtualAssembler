%% assmbl_pyram: function description
%% normalFacesE: the last is the rectangle
function [outputs] = assmbl_pyram(vertices,faces_of_E,face_types,face_pts,normalFacesE,measFacesE)
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

    % coefficients tetrahedral cubature on 4 points A. H. Stroud, Approximate Calculation of Multiple Integrals, Prentice-Hall, Englewood Cliffs, NJ, 1971
    n_vol_pts = 8;
    vol_pts = zeros(3,n_vol_pts); 
    const_a     = .58541019662496852;
    const_b     = .1381966011250105;
    vol_weights = .25*ones(8,1); % the weights for the volume quadrature points

**** reducir esto solo para pyrs
n_face_pts        = {};    %% Number of quad pts per element type per face
n_face_pts(4) = [3,3,3,3];
n_face_pts(5) = [9,3,3,3,3];
n_face_pts(6) = [9,9,9,3,3];

    %% vertices(:,6) is the top of the pyramid
           
    if ~are_coplanar (vertices(:,1:4))
      fprintf('pyramid %d with vertices in different order', el);
      exit;
    end
    %% sum of 2 tetrahedral cubatures from A. H. Stroud, Approximate Calculation of Multiple Integrals, Prentice-Hall, Englewood Cliffs, NJ, 1971
    %% tetrahedron: p_0 p_1 p_2 p_4
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
    %Mcoords    = [vertices(:,1), vertices(:,2), vertices(:,3), vertices(:,4), vertices(:,5)];
    Mdistances = [vertices(:,2)-vertices(:,1), vertices(:,1)-vertices(:,3),...
                  vertices(:,4)-vertices(:,1), vertices(:,5)-vertices(:,1),...
                  vertices(:,2)-vertices(:,3), vertices(:,2)-vertices(:,4),...
                  vertices(:,5)-vertices(:,2), vertices(:,3)-vertices(:,4),...
                  vertices(:,3)-vertices(:,5), vertices(:,4)-vertices(:,5)];
    M_Element = Mdistances(:, [3, 1, 4]);

    face_quad_coef = {};
    for index = 1:5
        face_quad_coef(faces_of_E(index)) = repmat(measFacesE(index)/3,3,1);
    end
    face_quad_coef(faces_of_E(find(face_types == 4))) = [1;1;1;1;4;4;4;4;16]*measFacesE(find(face_types == 4))/36;

    quad_nrmlztn                = abs(det(M_Element))/24; %% <--- measE/8;    
    rescale_factor              = 1/max(norm(Mdistances,2,'columns')); %% 1/diameter
    clear('Mdistances');

    we_potentials  			= zeros(n_vol_pts, 4);
    we_basis                = zeros(3, 4, n_vol_pts);

    for l = 1:n_vol_pts
      we_potentials(l,:) = WE_potentials(vol_pts(:,l),5);
      we_basis(:,:,l)    = WE_basis(vol_pts(:,l),5);   
    end                                                                    
    
    int_E_w_w = zeros(4);   %% int_E < w_k; w_r > dx

    for r = 1:4
      for k = 1:4
        %% ( (wr*wk)(p1), ..., (wr*wk)(p8) )
        int_E_w_w(r,k) = sum(reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts));
      end
    end
    
    int_E_w_w *= quad_nrmlztn; 

    %% READ: -quad_nrmlztn*we_potentials.'*vol_weights/measE  
    %% -1/mu(E) * int_E q_j
    b1 = repmat(-we_potentials.'*vol_weights/2,1,5);
    % Boundary term. Formerly: b2_{i,j} = Int_{fi} qj dS. page 62 in the middle.
    b2 = zeros(4,5); 
    we_potentials_faces = zeros (4,9,5);
    for f = 1:5
      for pt = 1:size(face_pts,2) 
        we_potentials_faces(:,pt,f) = WE_potentials(face_pts{faces_of_E(f)}(:,pt),5).';
      end
    end
    %% now the evaluated arrays have 1 to 5 face numbers.
    for f = 1:5
      b2(:,f) = we_potentials_faces(:,:,f) * face_quad_coef(faces_of_E(f)); %% (4x1)
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