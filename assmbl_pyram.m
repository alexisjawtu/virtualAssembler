%% assmbl_pyram: function description
%% normalFacesE: the last is the rectangle
function [outputs] = assmbl_pyram(Mdistances,det_M_Element,vol_pts,face_pts,P0,normalFacesE,measFacesE)
    % TODO: hacer un benchmark de inicializar pts_of_faces cada vez vs. pasarlo a la funcion
    n_face_pts              = [3,3,3,3,9];
    

    %pts_of_faces          	= zeros(3,9,5);
    %pts_of_faces(:,:,5)   	= [0 .5 1  0 .5  1 0 .5 1;
    %                           0  0 0 .5 .5 .5 1  1 1; 
    %                           0  0 0  0  0  0 0  0 0];
%
%  	%pts_of_faces(:,1:3,1) 	= [0 0 0;
%    %                           .5 .5 0;
%    %                           0 .5 .5];
%
%	  %pts_of_faces(:,1:3,2) 	= [.5 .5 0;
%    %                           1 .5 .5;
%    %                           0 .5 .5];
%
%  	%pts_of_faces(:,1:3,3) 	= [1 .5 .5;
%    %                           .5 .5 0;
%    %                           0 .5 .5];
%
%  	%pts_of_faces(:,1:3,4) 	= [.5 .5 0;
%    %                           0 0   0;
    %                           0 .5 .5];

    vol_wg 					        = .25*ones(8,1);
    face_quad_coef          = [[repmat(measFacesE(2:5)/3,3,1);zeros(6,4)],measFacesE(1)/36*[1; 4; 1; 4; 16; 4; 1; 4; 1]];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                              /                            \
    %%                              | m(f2) m(f3) m(f4) m(f5)  1 | 
    %%                              | m(f2) m(f3) m(f4) m(f5)  4 |
    %%                              | m(f2) m(f3) m(f4) m(f5)  1 |
    %%                              |   0     0     0     0    4 |
    %%   face_quad_coef =           |   0     0     0     0   16 |
    %%                              |   0     0     0     0    4 |
    %%                              |   0     0     0     0    1 |
    %%                              |   0     0     0     0    4 |
    %%                              |   0     0     0     0    1 |
    %%                              \                            /
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    measE         			    = abs(det_M_Element)/3;
    quad_nrmlztn  			    = measE/2;    
    rescale_factor 				  = 1/max(norm(Mdistances,2,'columns')); %% 1/diameter
    clear('Mdistances');
%   rescale_factor  		    = 1/diameter;
    

    face_pts
    


    we_potentials  			    = zeros(8, 4);
    we_basis                = zeros(3, 4, 8);

    for l = 1:8
      we_potentials(l,:) = WE_potentials(vol_pts(:,l),5);
      we_basis(:,:,l)    = WE_basis(vol_pts(:,l),5);   
    end                                                                    
    
    int_E_w_w = zeros(4);   %% int_E < w_k; w_r > dx

    for r = 1:4
      for k = 1:4
        %% ( (wr*wk)(p1), ..., (wr*wk)(p8) )
        vals            = reshape(dot(we_basis(:,r,:), we_basis(:,k,:), 1), 1, 8);  
        int_E_w_w(r,k)  = vals*vol_wg;
      end
    end
    
    int_E_w_w *= quad_nrmlztn;  %% hasta aca esta  y chequeado a mano

    %% READ: -quad_nrmlztn*we_potentials.'*vol_wg/measE  
    %% -1/mu(E) * int_E q_j
    b1 = repmat(-we_potentials.'*vol_wg/2,1,5);
    b2 = zeros(4,5);
    we_potentials_faces = zeros (4,9,5);
    for face = 1:5
      for pt = 1:n_face_pts(face)  %% square face first
        we_potentials_faces(:,pt,face) = WE_potentials(face_pts(:,pt,face),5).';
      end
    end

    for f = 1:5
      b2(:,f) = we_potentials_faces(:,:,f) * face_quad_coef(:,f);   %% ( 4 x 1)
    end

    b2    = b2./repmat(measFacesE,4,1);
    b     = b1 + b2;
    PROJ  = int_E_w_w\b;

    K_comput  = PROJ.' * int_E_w_w * PROJ;

    %% for the stabilizating bilinear form
OJO CON LA NORMAL local:
en este archivo la primera cara es cuadrilat
    D(face,:)  += face_quad_coef(pts,face) * normalFacesE(:,face).' * w_on_faces;

    Proj_in_base_Vh = D*PROJ;
    K_stab = rescale_factor * (eye(5) - PiW).'*(eye(5) - PiW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    W         = ones(1,5);  %% size (1 x Ndof_E). eqref(45) page 62 and page 63.
    % Boundary term. B2_{i,j} = Int_{fi} qj dS. page 62 in the middle.
    B2        = zeros(4,5); 
    % D_{j,i} = dof_j ( w_i ). See page 64 in the middle.
    D         = zeros(5,4);
    %% loop F
    for face = 1:5
      for pts = 1:n_face_pts{n_VERT}(face)
        p           = WE_potentials(diameter, centrosid, face_pts(:, pts, face), n_VERT);
        w_on_faces  = WE_basis (diameter, centroid, face_pts(:, pts, face), n_VERT);

        %% Next line: we are assembling the quadrature by columns
        B2(:,face) += face_quad_coef(pts,face) * p.';
        %% Next line: item 13 page 68. the quadrature
        D(face,:)  += face_quad_coef(pts,face) * normalFacesE(:,face).' * w_on_faces;
      end
    end

    %format long;
    %'max difference between G and B*D. Page 19 on Drive'
    %max(max(abs(int_E_w_w - check_)))

    PiW       = D*Pi_zero_star; %% item 14, page 68.
    %% next lines: build local matrix K 11.3 page 69
    %% computability term.
%    K_comp_a  = Pi_zero_star.' * int_E_w_w * Pi_zero_star;
%    temp_mtrx = eye(size(PiW,1)) - PiW;
%    %% stability term.
%
%    K_stab_a(5)  = rescale_factor * measE * temp_mtrx.'*temp_mtrx; 
%
%    clear('temp_mtrx');
%
%    A = K_comp_a + K_stab_a{n_VERT};
%
	outputs = {int_E_w_w, measE}
endfunction1