  function [Kexp]= K_d(di_ii);
% --------------------------------------------------------------------------
% Input data of the cantilever beam
% --------------------------------------------------------------------------

 E=200*10^9;         % Young's modulus:N/m^2
 Wb=8.55/100;          % beam width:m
 Tb=0.75/100;           % beam thickness:m
 Lb=200/100;          % beam length:m
 pb=7850;            % beam density (per unit volume):kg/m^3
 I=Wb*Tb^3/12;       % moment of inertia:m^4
 A=Wb*Tb;            % cross section area of the beam:m^2

% --------------------------------------------------------------------------
% Defining the elements and nodes 
% --------------------------------------------------------------------------

N_elem=10;
node=zeros(N_elem+1,2);
for i=1:N_elem+1
   node(i,1)=i;
   node(i,2)=Lb/N_elem*(i-1);
end
 NUM_NODE=length(node);   % the size of nodes
 NUM_ELEM=length(node)-1; % the size of elements
 matrix_size=2*NUM_ELEM+2;
% --------------------------------------------------------------------------
% Assembling of K matrix
% --------------------------------------------------------------------------
 K=zeros(matrix_size,matrix_size);
 K_bc=zeros(matrix_size-2,matrix_size-2);

  

 ELNO=0; % ELNO:the ith element
 for ii=1:1:matrix_size-3  
    ELNO=ELNO+1;
    [KE]=K_elemu(node(ELNO,2),node(ELNO+1,2),Lb,E,I);
    K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)=K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)+KE;
 end % end of FOR loop -- assembly of K matrix
%
  K_bc=K(3:matrix_size,3:matrix_size);

  for j = 1 : 1 : (2*N_elem)
       if j==1
       Kexp(j,j) = (di_ii(1,j)+di_ii(1,j+2)) * K_bc(j,j);
        Kexp(j,j+1) = (di_ii(1,j+1)) * K_bc(j,j+1);
        elseif j==(2*N_elem)
   Kexp(j,j) = di_ii(1,0.5*j) * K_bc(j,j);
   Kexp(j,j-1) = di_ii(1,0.5*j) * K_bc(j,j-1);
   elseif j==(N_elem)
   Kexp(j,j) = di_ii(1,j) * K_bc(j,j);
   Kexp(j,j-1) = di_ii(1,j) * K_bc(j,j-1);
       elseif (j>1) && (j<(N_elem))
   Kexp(j,j) =(di_ii(1,j)+di_ii(1,j+1)) * K_bc(j,j);       
   Kexp(j,j-1) = di_ii(1,j) * K_bc(j,j-1);
   Kexp(j,j+1) = di_ii(1,j+1) * K_bc(j,j+1);
      elseif (j>(N_elem)) && (j<(2*N_elem))
   Kexp(j,j) =(di_ii(1,(floor(0.5*j)))+di_ii(1,(floor(0.5*j))+1)) * K_bc(j,j);       
   Kexp(j,j-1) = di_ii(1,(floor(0.5*j))) * K_bc(j,j-1);
   Kexp(j,j+1) = di_ii(1,(floor(0.5*j))+1) * K_bc(j,j+1);
   % Y = floor( X ) rounds each element of X to the nearest integer less than or equal to that element.
       end
  end
  end

  