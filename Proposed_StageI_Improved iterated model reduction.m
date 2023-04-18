clc;
close all;
clear all;

%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage-I: Improved iterated model reduction technique based on Sengupta
% and Chakraborty (2023) (https://doi.org/10.1016/j.ymssp.2022.109586)
% Date: 17/03/2022
% By: Partha Sengupta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global M_mm M_ms M_sm M_ss K_mm K_ms K_sm K_ss K_ssinv phi_mmr l_mrmr;

ndof = 10;
n=ndof;
N=n;
M=zeros(n);
    for i = 1 : 1 : ndof
        M(i,i)=10^5;
    end

ki = 2 *(10^8); 
Kexp=zeros(n);
count = 0;
   
%block containng equal stiffness
    for j = 1 : 1 : ndof
       Kexp(j,j) =  (ki+ki);
       if j<ndof
       Kexp(j,j+1) = (-ki);  
       Kexp(j+1,j) = (-ki);
       end
    end
    Kexp(ndof,ndof) = (ki);
   K=Kexp
   
% eigen values and eigen vector
[p,lamda]=eig(K,M);
disp("The natural frequencies of the system are");
l=sort(diag(lamda));
wf=sort(diag(sqrt(lamda)))

p


%% Matrix partitioning
m=2;mr=2;
N=n;

M_mm = M(1:m,1:m);
M_ms = M(1:m,m+1:N);
M_sm = M(m+1:N,1:m);
M_ss = M(m+1:N,m+1:N);

K_mm = K(1:m,1:m);
K_ms = K(1:m,m+1:N);
K_sm = K(m+1:N,1:m);
K_ss = K(m+1:N,m+1:N);



K_ssinv=inv(K_ss);

phi_mmr=p(1:m,1:mr);

l_mrmr=lamda(1:mr,1:mr);

t=-K_ssinv*K_sm ;

  for i=1:N
a1 =t*phi_mmr*inv(l_mrmr)*phi_mmr'*t';

b1=t*phi_mmr*inv(l_mrmr)*phi_mmr';

c1=inv(phi_mmr*inv(l_mrmr)*phi_mmr');

d1=phi_mmr*inv(l_mrmr)*phi_mmr'*t';

 e1=inv(phi_mmr*inv(l_mrmr)*phi_mmr');
% 
% g1=phi_mmr*inv(l_mrmr)*phi_mmr'*t';

h1=t*phi_mmr*inv(l_mrmr)*phi_mmr'*t';

k1=t*phi_mmr*inv(l_mrmr)*phi_mmr';

l1=inv(phi_mmr*inv(l_mrmr)*phi_mmr');

m1=phi_mmr*inv(l_mrmr)*phi_mmr'*t';

m11 = inv(phi_mmr*inv(l_mrmr)*phi_mmr');

mm1=(pinv(h1-(k1*l1*m1)));

Kss_mod = mm1;

mmm1 = h1-(k1*l1*m1);

Kss_invmod = mmm1;

Kms_mod = -m11*m1*mm1;
 
%a11 = Kss_invmod*Kms_mod'
 % a = -((a1)-((b1)*(c1)*(d1))*(e1)*(g1)*(inv(((h1)-((k1)*(l1))*(m1)))))' 



   
    
     a = Kss_invmod*Kms_mod';
     b = -Kss_invmod*(M_sm + M_ss*t)*phi_mmr*l_mrmr*inv(phi_mmr);
     
%     b = -((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')))*(M_sm + M_ss*t)*phi_mmr*l_mrmr*l_mrmr*inv(phi_mmr)
%     c = -((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')))*(C_sm + C_ss*t)*phi_mmr*l_mrmr*inv(phi_mmr)
%     %c = -K_ssinv*(ab(1)*M_sm + ab(2)*((-inv(phi_mmr*inv(l_mrmr)*phi_mmr')*(phi_mmr*inv(l_mrmr)*phi_mmr*t'))*(inv((t*phi_mmr*inv(l_mrmr)*phi_mmr*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*inv(phi_mmr*inv(l_mrmr)*phi_mmr')*(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))))'+ab(1)*M_ss*t + ab(2)*inv(((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))))*t))*phi_mmr*l_mrmr*inv(phi_mmr) 
% 
     t= a+b;
  end
 
 %((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')))
  
    %(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')*(inv(((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*(inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')))))'

 
    t = (((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')* ...
      (inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))) * ...
      (inv(phi_mmr*inv(l_mrmr)*phi_mmr')*(phi_mmr*inv(l_mrmr)*phi_mmr*t')*...
      ((t*phi_mmr*inv(l_mrmr)*phi_mmr*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*inv(phi_mmr*inv(l_mrmr)*phi_mmr')*inv(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))))'...
      -(((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')* ...
      (inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))) * ...
      (M_sm + M_ss*t)*phi_mmr*l_mrmr*l_mrmr*inv(phi_mmr))...
      -((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')* ...
      (inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))) * ...
      (ab(1)*M_sm + ab(2)*(-inv(phi_mmr*inv(l_mrmr)*phi_mmr')*(phi_mmr*inv(l_mrmr)*phi_mmr*t')*...
      ((t*phi_mmr*inv(l_mrmr)*phi_mmr*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')*inv(phi_mmr*inv(l_mrmr)*phi_mmr')*inv(phi_mmr*inv(l_mrmr)*phi_mmr'*t'))))'...
      + ab(1)*M_ss*t + ab(2)*inv(((t*phi_mmr*inv(l_mrmr)*phi_mmr'*t')-((t*phi_mmr*inv(l_mrmr)*phi_mmr')* ...
      (inv(phi_mmr*inv(l_mrmr)*phi_mmr'))*(phi_mmr*inv(l_mrmr)*phi_mmr'*t')))) *t)*phi_mmr*l_mrmr*inv(phi_mmr)) 

         [L,U,P] = lu(t);

[U1,S1,V1] = svd(L);
s1 = diag(S1);   % vector of singular values
tolerance = max(size(L))*eps(max(s1));
p1 = sum(s1>tolerance);
% Define spaces
Up1 = p1+U1(:,1:p1);
%U0 = U(:,p+1:Nx);
Vp1 = p1+(V1(:,1:p1));
%V0 = V(:,p+1:Nx);
%Sp = spdiags( s(1:p), 0, p, p );
%Sp1 = spdiags(1.0.*s1(1:p1), 0, p1, p1 );
Sp1 = spdiags(max(s1).*s1(1:p1), 0, p1, p1 );
% Calc AInv such that x = AInv * b
Lnv = Up1 * Sp1 * Vp1';

[U2,S2,V2] = svd(U);




s2 = diag(S2);   % vector of singular values
tolerance = max(size(U))*eps(max(s1));
p2 = sum(s2>tolerance);
% Define spaces
Up2 = zeros(m,m);
for i =1:m
for j =1:m
Up2(i,j)  = phi_mmr(i,j) + V2(i,j) ;
end
end

Up22 = (p1+Up2(:,1:p1));

%P = 10^-(N/2) ;
%Up22 = Up2+[8.71 ; 5 ; 7 ; 9 ; 5];
%U0 = U(:,p+1:Nx);
Vp2 =  (p2+V1(:,1:p2));
%V0 = V(:,p+1:Nx);
%Sp = spdiags( s(1:p), 0, p, p );
%Sp2 = spdiags(1.0./s1(1:p2), 0, p2, p2 );
Sp2 = spdiags(max(s1).*s1(1:p1), 0, p1, p1 );
%Sp21 = % Calc AInv such that x = AInv * b
Unv = Vp1 * S2 * Up2;
t1 =   Lnv*Unv;
 
t1;

epstol = [-3.44436578288908e+16,1.75201530508875e+16;-10071473191869.8,5122967852217.86];

t2 = t1 * epstol ;
  
 



T=[eye(m);t2];


M_R=T'*M*T;
K_R=T'*K*T;





%solving eigen vaulue problem for reduced mass and stiffness matrix
[v,lam]=eig(K_R,M_R);

wr=sort(diag(sqrt(lam)))




%wr=wr./wr(1);

%wr=wr*wf(1) %from reduced model

wf(1:m)     %from full model

phi=(T*p(1:m,:))

phi_n = zeros(n);

for i=1:n
x=max(phi(:,i));
y=min(phi(:,i));
if abs(x)>abs(y)
    phi_n(:,i)=phi(:,i)./x;
elseif abs(x)<abs(y)
    phi_n(:,i)=phi(:,i)./y;
end
end

phi_n

p(1:m,:)

pr=[0 0 0 0 0 0 0 0;phi];
%% plot of mode shapes
Y=[0;0;0;0;0;0;0;0;0];
floor=[0;1;2;3;4;5;6;7;8];

figure,
for i=1:N
    subplot(1,8,i);
    plot([pr(:,i)],floor,'b-o','LineWidth',1.5)
    yticks([0 1 2 3 4 5 6 7 8]);
    ylabel('Mode Numbers');
    title(['Mode Shape ?',num2str(i)],'Fontsize',12);
end

 