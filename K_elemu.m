% k_elemu
function [KE]=K_elemu(prenode,postnode,Lb,E,I)
% PURPOSE : This is a subprogram as for Stiffness matrice
%
%  K=[ 12   6*lb    -12    6*lb
%           4*lb^2  -6*lb  2*lb^2  
%                    12   -6**lb
%      symmetric           4*lb^2  ];      
%
  l=postnode-prenode;  
  lb=l/Lb;
     % the length of present element per unit length, lb:length_bar
  KE(1,1)=12     ;                
  KE(2,1)=-6*lb  ; KE(2,2)=4*lb^2 ; 
  KE(3,1)=-12    ; KE(3,2)=6*lb   ; KE(3,3)=12     ; 
  KE(4,1)=-6*lb  ; KE(4,2)=2*lb^2 ; KE(4,3)=6*lb   ; KE(4,4)=4*lb^2;
  KE=KE/lb^3*E*I/Lb^3;
for i=2:4
    for j=1:i-1
        KE(j,i)=KE(i,j);
    end
end