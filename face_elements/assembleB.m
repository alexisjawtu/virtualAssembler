%% just for experiments and testing
function [res] = assembleB(num_faces)

  quad_pts_prism    = [.5 0 0; .5 .5 0; 0 .5 0; .5 0 .5; .5 .5 .5; 0 .5 .5; .5 0 1; .5 .5 1; 0 .5 1].';
  prism_face_coefs  = [1; 4; 1; 4; 16; 4; 1; 4; 1]; %% weights over faces

  % coefficients tetrahedral cubature on 4 points GELLERT and HARBORD 91
  const_a     = .58541019662496852;
  const_b     = .1381966011250105;

  vertices              = importdata ('vertices.txt');
  elements              = importdata ('elements_by_vertices.txt');
  elements_by_faces     = importdata ('elements_by_faces.txt');
  faces                 = importdata ('faces.txt');

 %vertices              = importdata ('vertices_by_elements.txt');
 %shared_faces          = importdata ('shared_faces.txt');
 %faces_local_to_global = importdata ('faces_local_to_global.txt');

  vertices              = vertices.';
  num_el                = size (elements,1);
  num_fc                = size (faces,1);

  %% read as columns in R3.

  %% points in the faces of the pyramid
  
  % En octave es por valor; hace una copia si modifico el argumento dentro
  % de la funcion. Una funcion llamante no modifica variables llamadas
  dim_WE            = {};  
  n_Faces           = {};
  meas_Refs         = {};
  n_vol_pts         = {};
  vol_wg            = {};
  n_face_pts        = {};    %% Number of quad pts per element type per face
  quad_aux_const    = {};    %% adjusting the volume quadrature
  rescale_factor    = {};
  K_stab_a          = {};
  
  K_stab_a(4) = zeros(4);
  K_stab_a(6) = zeros(5);

  rescale_factor(4) = 1;
  rescale_factor(6) = 1;

  n_vol_pts(4) = 4;
  n_vol_pts(5) = 8;
  n_vol_pts(6) = 9;

  vol_wg(4) = .25*ones(n_vol_pts{4},1);
  vol_wg(5) = .25*ones(n_vol_pts{5},1);
  vol_wg(6) = [1; 4; 1; 1;  4; 1; 1; 4; 1];

  dim_WE(4)   = 4;  
  dim_WE(5)   = 4;  
  dim_WE(6)   = 5;  

  n_Faces(4)  = 4;  
  n_Faces(5)  = 5;  
  n_Faces(6)  = 5;

  n_face_pts(4) = [3,3,3,3];
  n_face_pts(5) = [9,3,3,3,3];
  n_face_pts(6) = [9,9,9,3,3];

  meas_Refs(4) = 1/6;
  meas_Refs(5) = 1/3;
  meas_Refs(6) = 1/2;

  quad_aux_const(4) = 1;
  quad_aux_const(5) = .5;
  quad_aux_const(6) = 1/36; 

  %% quadrature points are requested only once:
  cell_index    = {};
  cell_index(4) = 1;
  cell_index(5) = 2;
  cell_index(6) = 3;
  Mdistances = 0;

  bit = 1;
  bit2 = 0;

  pyr_counter = 0;

  fprintf('Elements loop. %d\n',num_el);
  tic
  for el = 1:num_el
    clear ('p_0','p_1','p_2','p_3','p_4');
    n_VERT = elements(el,1);         %% to know which type of element

    if n_VERT == 5
      pyr_counter += 1; 
      %% TODO: ver si esto es igual!
      %% en ese caso, mejorar el Mcoords y el vol_pts
      P = vertices(:,elements(el,2:6));

      bit = bit && are_coplanar(P);
      bit2 = bit2 + counter_clockwise(P);

      p_0 = vertices(:,elements(el,2)); 
      p_1 = vertices(:,elements(el,3));
      p_2 = vertices(:,elements(el,4));
      p_3 = vertices(:,elements(el,5));
      %%% this is the top of the pyramid
      p_4 = vertices(:,elements(el,6));   
   
    elseif n_VERT == 4

      %% TODO: mismo que para n_VERT = 5

      %P = vertices(:,elements(el,2:5));

      p_0 = vertices(:,elements(el,2));
      p_1 = vertices(:,elements(el,3));
      p_2 = vertices(:,elements(el,4));
      p_3 = vertices(:,elements(el,5));
    
    elseif n_VERT == 6

      % TODO: misma cosa que para n_VERT = 4

      %P = vertices(:,elements(el,2:7));

      p_0 = vertices(:,elements(el,2));
      p_1 = vertices(:,elements(el,3));
      p_2 = vertices(:,elements(el,4));
      p_3 = vertices(:,elements(el,5));
      p_4 = vertices(:,elements(el,6));
      p_5 = vertices(:,elements(el,7));

      Mdistances  = [p_0-p_1, p_0-p_2, p_0-p_3, p_0-p_4, p_0-p_5, p_1-p_2, p_1-p_3, p_1-p_4, p_1-p_5, ...
                      p_2-p_3, p_2-p_4, p_2-p_5, ...
                      p_3-p_4, p_3-p_5, ...
                      p_4-p_5];

 
    else 
      error('invalid number of vertices. elements ' + num2str(el));
    
%%%% end of if 
    end
    


  end    %% end of main for
  toc
  res = 1;
  fprintf('coplanar: %d\ncounter: %d\npyr counter: %d\n', bit, bit2, pyr_counter);
endfunction
