function [phi0, f0] = priorformodal(theta,modem,Mmat,N)
    Nd = N(1); Nm = N(4);
    [Kmat] = kstifffunc(Nd,theta);
    [modesh0,w2r]     = eig(Kmat,Mmat);
    modesh            = modesh0(:,modem);
    w2r               = diag(w2r);
    f0(:,1) = w2r(modem);
    phimodel = zeros(Nd,Nm);
    for m = 1:Nm
        modesh(:,m)   = modesh(:,m)/norm(modesh(:,m)); % to normalize all modeshapes to have unit norm
        phimodel(:,m) = modesh(:,m)/sign(modesh(end,m));
    end

    phi0 = reshape(phimodel,Nd*Nm,1);