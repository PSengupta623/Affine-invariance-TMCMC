% m_elemu
function [ME]=m_elemu(prenode,postnode,Lb,pb,A)
% PURPOSE : This is a subprogram as for Mass matrice
%
%  M=1/420[ 156  22*lb   54     -13*lb
%                4*lb^2  13*lb  -3*lb^2  
%                        156    -22*lb
%           symmetric            4*lb^2];      
%
    l=postnode-prenode; 
    lb=l/Lb;
       % the length of the present element per unit length, lb:length_bar
    ME(1,1)=156    ;
    ME(2,1)=22*lb  ; ME(2,2)=4*lb^2 ; 
    ME(3,1)=54     ; ME(3,2)=13*lb  ; ME(3,3)=156    ; 
    ME(4,1)=-13*lb ; ME(4,2)=-3*lb^2; ME(4,3)=-22*lb ; ME(4,4)=4*lb^2 ;
    ME=ME*lb/420*pb*A*Lb;
    for i=2:4
      for j=1:i-1
          ME(j,i)=ME(i,j);
      end
  end