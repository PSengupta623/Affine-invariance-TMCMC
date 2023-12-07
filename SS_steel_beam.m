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
 Wb=8.55/100;          % beam width:m
 Tb=0.75/100;           % beam thickness:m
 Lb=200/100;          % beam length:m
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
 M_bc=zeros(matrix_size-2,matrix_size-2);
 M_bc=M(2:matrix_size-1,2:matrix_size-1);
 
 % for one element damage training dataset:
 dlevel=10;

 for jj=1:1:N_elem 
  for ii=1:1:dlevel
      if (1<=ii)&&(ii<=5) 
       di_ii= ones(1,N_elem);
       di_ii(jj)=normrnd(2,1);
      elseif (6<=ii)&&(ii<=10)
       di_ii= ones(1,N_elem);
       di_ii(jj)=normrnd(2,1);
      elseif (11<=ii)&&(ii<=15)
       di_ii= ones(1,N_elem);
       di_ii(jj)= normrnd(2,1);
      else
       di_ii= ones(1,N_elem);
       di_ii(jj)=normrnd(2,1);

      end
   %Kexp = K_d(di_ii);
   % (1). assembly of K matrix
 ELNO=0; % ELNO:the ith element
   for i=1:2:matrix_size-3  
    ELNO=ELNO+1;
    [KE]=di_ii(ELNO)*(K_elemu(node(ELNO,2),node(ELNO+1,2),Lb,E,I));
    K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)=K((ELNO*2-1):(ELNO+1)*2,(ELNO*2-1):(ELNO+1)*2)+KE; % matrix assembly 
   end  % end of FOR loop -- assembly of K matrix
    K_bc=zeros(matrix_size-2,matrix_size-2);
    K_bc=K(2:matrix_size-1,2:matrix_size-1);
%
% Calculate eigenvalues by the finite element formulation 
[p,lamda]=eig(K_bc,M_bc);
l=sort(diag(lamda));
w=(real(sort(diag(sqrt(lamda)))));

% Considering only the translational dof:
p1=zeros(N_elem,N_elem);
p1((1:1:N_elem),(1:1:N_elem))= p((1:2:(2*N_elem)),(1:2:(2*N_elem)));

% Normalization:
for k=1:1:length(p1)
    for kk=1:1:length(p1)
    pn(kk,k)= p1(kk,k)/max(abs(p1(:,k))); %normalized p1
    end
end
  pnd((1+((ii-1)*N_elem):ii*N_elem),1:N_elem) = pn;
  di_jj(1+(ii-1)*N_elem:ii*N_elem) = di_ii;
 end
    phi(1+(jj-1)*N_elem*dlevel:jj*N_elem*dlevel,1:N_elem) = pnd;
    di(1+(jj-1)*N_elem*dlevel:jj*N_elem*dlevel) = di_jj;
 end
    di=di.'; % damage parameter: 1= undamaged, 0= fully damaged
    phi=(-1)*phi; % Normalized Mode shape matrix
    
% xnode=1:1:N_elem;
% figure
% plot(xnode,phi(1:10,1));
