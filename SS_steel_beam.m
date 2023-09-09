% Date: 3/9/2022
% By: Partha Sengupta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Provide all the related input in the sub-programm
% "SS_steel_beam.m" for each example problem


clear all;clf;
format long e;
% --------------------------------------------------------------------------
% Input data of the beam
% --------------------------------------------------------------------------

 E=200*10^9;         % Young's modulus:N/m^2
 Wb=10/100;          % beam width:m
 Tb=1/100;           % beam thickness:m
 Lb=35/100;          % beam length:m
 pb=7850;            % beam density (per unit volume):kg/m^3
 I=Wb*Tb^3/12;       % moment of inertia:m^4
 A=Wb*Tb;            % cross section area of the beam:m^2

% --------------------------------------------------------------------------
% Defining the elements and nodes 
% --------------------------------------------------------------------------

N_elem=10; % elements
node=zeros(N_elem+1,2); % nodes
for i=1:N_elem+1
   node(i,1)=i;
   node(i,2)=Lb/N_elem*(i-1);
end
 NUM_NODE=length(node);   % the size of nodes
 NUM_ELEM=length(node)-1; % the size of elements
 matrix_size=2*NUM_ELEM+2;
%
% --------------------------------------------------------------------------
% Initializes M & K matrix
% --------------------------------------------------------------------------
 M=zeros(matrix_size,matrix_size);
 K=zeros(matrix_size,matrix_size);
% --------------------------------------------------------------------------
% Assembling of M matrix
% --------------------------------------------------------------------------
%
 ELNO=0; % ELNO:the ith element
 for ii=1:2:matrix_size-3  
    ELNO=ELNO+1;
    [ME]=m_elemu(node(ELNO,2),node(ELNO+1,2),Lb,pb,A);
    M((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)=M((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)+ME;  
 end
 % end of FOR loop -- assembly of M matrix

   for i=1:2:matrix_size-3  
    ELNO=ELNO+1;
    [KE]=di_ii(ELNO)*(K_elemu(node(ELNO,2),node(ELNO+1,2),Lb,E,I));
    K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)=K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)+KE; % matrix assembly 
   end  % end of FOR loop -- assembly of K matrix
    K_bc=zeros(matrix_size-3,matrix_size-3);
    K_bc=K(3:matrix_size-1,3:matrix_size-1);
%
 M_bc=zeros(matrix_size-3,matrix_size-3);
 M_bc=M(3:matrix_size-1,3:matrix_size-1);
 
 for jj=1:N_elem
theta1(jj)=unifrnd(0.1,0.9);
 end

    theta1;

  
  for j = 1 : 1 : (2*N_elem)
       if j==1
       Kexp(j,j) = (theta1(1,j)+theta1(1,j+2)) * K_bc(j,j);
        Kexp(j,j+1) = (theta1(1,j+1)) * K_bc(j,j+1);
        elseif j==(2*N_elem)
   Kexp(j,j) = theta1(1,0.5*j) * K_bc(j,j);
   Kexp(j,j-1) = theta1(1,0.5*j) * K_bc(j,j-1);
   elseif j==(N_elem)
   Kexp(j,j) = theta1(1,j) * K_bc(j,j);
   Kexp(j,j-1) = theta1(1,j) * K_bc(j,j-1);
       elseif (j>1) && (j<(N_elem))
   Kexp(j,j) =(theta1(1,j)+theta1(1,j+1)) * K_bc(j,j);       
   Kexp(j,j-1) = theta1(1,j) * K_bc(j,j-1);
   Kexp(j,j+1) = theta1(1,j+1) * K_bc(j,j+1);
   elseif (j>(N_elem)) && (j<(2*N_elem))
   Kexp(j,j) =(theta1(1,(floor(0.5*j)))+theta1(1,(floor(0.5*j))+1)) * K_bc(j,j);       
   Kexp(j,j-1) = theta1(1,(floor(0.5*j))) * K_bc(j,j-1);
   Kexp(j,j+1) = theta1(1,(floor(0.5*j))+1) * K_bc(j,j+1);
   % Y = floor( X ) rounds each element of X to the nearest integer less than or equal to that element.
       end
end
Kexp;


 
% Calculate eigenvalues by the finite element formulation 
[p,lamda]=eig(Kexp,M_bc);
l=sort(diag(lamda));
w=(real(sort(diag(sqrt(lamda)))));

 