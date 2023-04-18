% A simulated data is generated by adding noise to the modal data
% (frequency and modeshapes)
clear; 
damagecase   = 0;   % 0: undamaged, 1: damage case

% Information for the Model
N          = 10;                   %  number of DOFs
J          = 300;                  % No. of data sets or experimental data
if damagecase == 0
    modem    = [1:10];                 % modes to be observed for the undamaged structure
else
    modem    = [1 2 3];               % modes to be observed for the damaged structure
end
M = length(modem);
sen         = [1:10];                 % measured DOFs
S = length(sen);
elemass       = (100/1000) * ones(1,N);
theta = 176.729*ones(N,1);
indexdamage = zeros(N,1); % Put 1 in the desired damage location
damageFactor = zeros(N,1); % the degree of damage in decimal value
[w2r, phi] = FEmodelSolver(N, damagecase, elemass, theta, sen, indexdamage, damageFactor);
w2m = w2r(modem);
phim = phi(:,modem);

% ***************************
%  Generation of modal data
%  ***************************
wm      = w2m.^0.5; % circular frequency in rad/s
fr      = 1/(2*pi)*wm; % natural frequency in Hz
phinoise = zeros(S*M,J);
w2noise = zeros(M,J);
noisyphi = zeros(S,M);
nu = 5; % degrees of freedom, put nu>200 for normal distribution
np = 5; % noise percentage in frequency
psh = (5/100); % noise percentage in modeshapes

 posterior(i,s) = normrnd(p(s,i),((np/100)*(p(s,i))));
for j = 1 : J
    for m = 1 : M
        noise           = np*phim(:,m).*trnd(nu,[S,1]); 
        noisyphi(:,m)  = phim(:,m) + noise;
    end
    phinoise(:,j)   = reshape(noisyphi,S*M,1);
    clear noise
    noise = pw2 * w2m.* trnd(nu,[M,1]); % nw2*w2r.* randn
    w2noise(:,j)            = w2m + noise;
end

truesig2phi = zeros(S,M);
for m = 1 : M
    truesig2phi(:, m) = (psh*phim(:,m)).^2;
end
truesig2z = (pw2 * w2m).^2;

save('truevariance.mat', 'truesig2phi', 'truesig2z')
save('modaldata.mat','phinoise','phim','w2noise','w2m')

% y = [linspace(15,70,100); linspace(200,500,100); linspace(500,1500,100)];
% figure()
% for i = 1:3
%     subplot(3,1,i)
%     h = histogram(w2noise(i,:),'Normalization', 'pdf');
% %     histogram(w2noise(i,:))
%     hold on 
%     plot(y(i,:), normpdf(y(i,:),w2r(i),truesigz(i)));
% %     hold on
% %     plot(y(i,:), tpdf(y(i,:),5));
% end


