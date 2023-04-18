function [theta] =  drawtheta(f, phi, sigthetasqr, Massmat, invB0, theta0, N)

    P = N(1); No = N(2); J = N(3); M = N(4); K = N(5);
    invsigthetasqr = 1/sigthetasqr;
    H = zeros(P*M,K);
    b = zeros(M*P,1);
    [Ksub] = substiffnessmatrix(P);
    for i = 1 : M
        phip = phi((i-1)*P+1:i*P);
            for k = 1 : K
                H((i-1)*P+1:i*P,k) = Ksub(:,:,k) * phip;
            end
        b((i-1)*P+1:i*P,1)=(f(i)*Massmat)*phip;
    end  

    invB1 = invB0 + invsigthetasqr* (H)'* H;
    B1 = invB1\eye(size(invB1));
    thetatilde = B1 * (invB0 * theta0 + invsigthetasqr * (H)'* b);
    [NorRand] = mvn_rand(B1);
    theta = thetatilde + NorRand; % Drawing samples for Stiffness Parameters 
    