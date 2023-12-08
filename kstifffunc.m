function [Kmat] = kstifffunc(Nd,theta)

% Formation of stiffness matrix
       elementmatrix    = [1 -1;-1 1];     % Element Matrix for Stiffness Matrix
       Kmat        = zeros(Nd,Nd);  
       Kstiffsub        = zeros(Nd,Nd); 
       Kstiffsub(:,:,1) = theta(1) * [1 zeros(1,Nd-1);zeros(Nd-1,1) zeros(Nd-1,Nd-1)];
       for ii = 2:Nd  %% Definition of substructure stiffness matrix
           Kstiffsub((ii-1):ii,(ii-1):ii,ii) = theta(ii) * elementmatrix;
       end
       for jj = 1:Nd
           Kmat = Kmat + Kstiffsub(:,:,jj); % Distribution of Subtructure stiffness to Global Matrix
       end