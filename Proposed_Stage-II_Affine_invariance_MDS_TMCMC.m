
load modaldata.mat % loading the modal data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage-II: Affine-invariance based Transitional Markov Chain Monte Carlo sampler
% Date: 3/9/2022
% By: Partha Sengupta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Provide all the related input in the sub-programm
% "input_parameter.m" for each example problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ General information of the structure ---
N          = 10;                   %  # DOFs
sen         = [1:2];         % measured DOFs 1 3 5 7 8 10
N0          = length(sen);          % total number of measured DOFs
Nset          = 10;                   % No. of data set
modem       = [1:1];              % modes to be observed
M          = length(modem);        % # modes considered
P          = 10;                   % No. of parameters
Nset1           = [N,S,J,M,P];
thetaj   = normrnd(0,1) * ones(P,1); %  stiffness parameters
Mmat        = (100/1000) * eye(N); % Mass mastrix of uniform mass
% d1 = 0.0; d4 = 0.0; % Damage factor

% phicap = zeros(S*M,J);
% phiallmode = zeros(N,M);
% for i = 1:M
%     for j = 1:J
%         phiset = phinoise(:,j); % per set of experimental data
%         for n = 1:M
%             phiallmode(:,n) = phiset((n-1)*N+1:n*N);
%         end
%         phimode = phiallmode(:,modem(i)); % per mode
%         phicap((i-1)*S+1:i*S,j) = phimode(sen); % only where sensor is placed
%     end
% end

% observed response
phicap      = phinoise(:,1:Nset);
zcap   = w2noise(modem,1:Nset);

% Observation Matrix that picks the components of phi corresponding to No-measured DOFs.
Lo0     = zeros(N0,N); % Matrix of size No x Nd
for i = 1 : N0
    Lo0(i,sen(i)) = 1;
end
% for m = 1:M
%     Lo((m-1)*No+1:m*No,(m-1)*P+1:m*P) = Lo0;
% end
Gamma   = Lo0;

% ----- Prior  Specification ----
%   For Stiffness
theta0      = 200 * ones(P,1);
B0   = (100)^(2)*eye(P);
invB0 = B0\eye(size(B0));

% for phi
[phi0, z0] = priorformodal(theta0,modem,Mmat,Nset);

%   uncertainties for sigma2z
alpha0 = 5;
delta0 = 1;
sig2z0 = 1/gamrnd(alpha0/2,2/delta0);

%   uncertainties for sigma2phi
p0 = 5;
q0 = 1e-3;
sig2phi0 = 1/gamrnd(p0/2,2/q0);

%   uncertainties for sigmathetasqr
a0 = 5;
b0 = 10;
sig2theta0 = 1/gamrnd(a0/2,2/b0);

% hyper-parameters
nu = [0.01 1 50]; % for lambda
eta = 2; % for h

% quantities calculated outside the loop
alpha1 = alpha0 + M*Nset;
p1 = p0 + S*M*Nset;
a1 = a0 + N*M;

%------------ Affine-invariance to estimate MDS factor
% Define default values: 
dflts =  {[],[],[],[],[],0,2,3}; % define default values
[nsamples,loglikelihood,w2noise,priorformodal,burnin,lastBurnin,stepsize,thinchain] = ...
       internal.stats.parseArgs(pnames, dflts, modeshape_w2freq{:});
   
persistent isoctave;  
if isempty(isoctave)
	isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end

if nargin<3
    error('GWMCMC:toofewinputs','GWMCMC requires atleast 3 inputs.')
end
M=size(minit,2);
if size(minit,1)==1
    minit=bsxfun(@plus,minit,randn(M*5,M));
end


p=inputParser;
if isoctave
    p=p.addParamValue('StepSize',2,@isnumeric); %addParamValue is chosen for compatibility with octave. Still Untested.
    p=p.addParamValue('ThinChain',10,@isnumeric);
    p=p.addParamValue('ProgressBar',false,@islogical);
    p=p.addParamValue('Parallel',false,@islogical);
    p=p.addParamValue('BurnIn',0,@isnumeric);
    p=p.parse(varargin{:});
else
    p.addParameter('StepSize',2,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
    p.addParameter('ThinChain',10,@isnumeric);
    p.addParameter('ProgressBar',false,@islogical);
    p.addParameter('Parallel',false,@islogical);
    p.addParameter('BurnIn',0,@isnumeric);
    p.parse(varargin{:});
end
p=p.Results;

Nwalkers=size(minit,1);

if size(minit,2)*2>size(minit,1)
    warning('Affinedimensions','Check minit dimensions.\nIt is recommended that there be atleast twice as many walkers in the ensemble as there are model dimension.')
end

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end



Nkeep = Nsamples + p.BurnIn; % number of samples drawn per walker

models=nan(Nwalkers,M,Nkeep); % pre-allocate output matrix

models(:,:,1)=minit; % models: A WxMxT matrix, minit: A Mx(W*T) matrix

if ~iscell(logPfuns)
    logPfuns={logPfuns};
end

NPfun=numel(logPfuns);

%calculate logP state initial pos of walkers
logP=nan(Nwalkers,NPfun,Nkeep); %logP = WxPxT
for wix=1:Nwalkers
    for fix=1:NPfun
        v=logPfuns{fix}(minit(wix,:));
        if islogical(v) %reformulate function so that false=-inf for logical constraints.
            v=-1/v;logPfuns{fix}=@(m)-1/logPfuns{fix}(m); %experimental implementation of experimental feature
        end
        logP(wix,fix,1)=v;
    end
end

if ~all(all(isfinite(logP(:,:,1))))
    error('Starting points for all walkers must have finite logP')
end


reject=zeros(Nwalkers,1);

% models: A WxMxT matrix; logP: WxPxT matrix
curm = models(:,:,1);  %curm: W x M matrix
curlogP = logP(:,:,1); %curlogP: W x P matrix
progress(0,0,0)
totcount=Nwalkers;
for row = 1:Nkeep % number of samples drawn per walker
    for jj=1:p.ThinChain
        %generate proposals for all walkers
        rix = mod((1:Nwalkers)+floor(rand*(Nwalkers-1)),Nwalkers)+1; % pick a random partner (Nwalker x 1 vector)
        
        proposedm = zeros(Nwalkers, size(minit,2));  % Nwalkers x dim matrix
        zz = zeros(Nwalkers, 1);                     % Nwalkers x 1 vector
        for i = 1:Nwalkers
        while true
        zz(i) = ((p.StepSize - 1)*rand(1,1) + 1).^2/p.StepSize;  % scalar
        proposedm(i,:) = curm(rix(i),:) - bsxfun(@times,(curm(rix(i),:)-curm(i,:)),zz(i)); % Nwalkers x dim matrix
        if box(proposedm(i,:)) % The box function is the Prior PDF in the feasible region.
        % Note: If a point is out of bounds, this function will return 0 = false.
        break;
        end
        end
        end
        
        logrand=log(rand(Nwalkers,NPfun+1)); %moved outside because rand is slow inside parfor
        if p.Parallel
            %parallel/non-parallel code is currently mirrored in
            %order to enable experimentation with separate optimization
            %techniques for each branch. Parallel is not really great yet.
            %TODO: use SPMD instead of parfor.

            parfor wix=1:Nwalkers
                cp=curlogP(wix,:);
                lr=logrand(wix,:);
                acceptfullstep=true;
                proposedlogP=nan(1,NPfun);
                if lr(1)<(numel(proposedm(wix,:))-1)*log(zz(wix))
                    for fix=1:NPfun
                        proposedlogP(fix)=logPfuns{fix}(proposedm(wix,:)); 
                        if lr(fix+1)>proposedlogP(fix)-cp(fix) || ~isreal(proposedlogP(fix)) || isnan( proposedlogP(fix) )
                            acceptfullstep=false;
                            break
                        end
                    end
                else
                    acceptfullstep=false;
                end
                if acceptfullstep
                    curm(wix,:)=proposedm(wix,:); curlogP(wix,:)=proposedlogP;
                else
                    reject(wix)=reject(wix)+1;
                end
            end
        else %NON-PARALLEL
            for wix=1:Nwalkers
                acceptfullstep=true;
                proposedlogP=nan(1,NPfun);
                if logrand(wix,1)<(numel(proposedm(wix,:))-1)*log(zz(wix))
                    for fix=1:NPfun
                        proposedlogP(fix)=logPfuns{fix}(proposedm(wix,:));
                        if logrand(wix,fix+1)>proposedlogP(fix)-curlogP(wix,fix) || ~isreal(proposedlogP(fix)) || isnan(proposedlogP(fix))
                            acceptfullstep=false;
                            break
                        end
                    end
                else
                    acceptfullstep=false;
                end
                if acceptfullstep
                    curm(wix,:)=proposedm(wix,:); curlogP(wix,:)=proposedlogP;
                else
                    reject(wix)=reject(wix)+1;
                end
            end

        end
        totcount=totcount+Nwalkers;
        progress((row-1+jj/p.ThinChain)/Nkeep,curm,sum(reject)/totcount)
    end
    models(:,:,row)=curm;
    logP(:,:,row)=curlogP;
end
progress(1,0,0);

acceptance = 1 - (sum(reject)/totcount);

if p.BurnIn>0
    crop=p.BurnIn;
    models(:,:,1:crop)=[]; 
    logP(:,:,1:crop)=[];
end

 me(j) = mean(modeshape_w2freq(j));
sd(j) = std(modeshape_w2freq(j));
MDS(j) = max(((sd(j)/me(j))-1));

j      = 0;                   % Initialise loop for the transitional likelihood
thetaj = prior_rnd(nsamples); % theta0 = N x D
pj     = 0;                   % p0 = 0 (initial tempering parameter)
Dimensions = size(thetaj, 2); % size of the vector theta

count = 1; % Counter
samps(:,:,count) = thetaj;
beta_j(count) = pj;

%% Initialization of matrices and vectors
thetaj1   = zeros(nsamples, Dimensions);
step(count) = stepsize;
%% Main loop
while pj < 1    
    j = j+1;
    
    %% Calculate the tempering parameter p(j+1):
    for l = 1:nsamples
        log_fD_T_thetaj(l) = loglikelihood(thetaj(l,:));
    end
    if any(isinf(log_fD_T_thetaj))
        error('The prior distribution is too far from the true region');
    end
    pj1 = calculate_pj1(log_fD_T_thetaj, pj);
    fprintf('TEMCMC: Iteration j = %2d, pj1 = %f\n', j, pj1);
    
    %% Compute the plausibility weight for each sample wrt f_{j+1}
    fprintf('Computing the weights ...\n');
    a       = (pj1-pj)*log_fD_T_thetaj;
    wj      = exp(a);
    wj_norm = wj./sum(wj);                % normalization of the weights
    
    %% Compute S(j) = E[w{j}] (eq 15)
    S(j) = mean(wj);
    
    log_posterior = @(t) log(priorpdf(t)) + pj1*loglikelihood(t);
    
    %% During the last iteration we require to do a better burnin in order
    % to guarantee the quality of the samples:
    if pj1 == 1
        burnin = lastBurnin;
    end
    
    %% Start Nc different Markov chains
    fprintf('Markov chains ...\n\n');
    idx = randsample(nsamples, nsamples, true, wj_norm);
    
    start = zeros(nsamples, Dimensions);
    for i = 1:nsamples
            start(i,:) = thetaj(idx(i), :);
    end

% Preallocating

[samples,logp,acceptance_rate] = EMCMCsampler(start, log_posterior, 1, priorpdf, ...
                            'StepSize', stepsize,...
                            'BurnIn', burnin,...
                            'ThinChain', thinchain); 
        
        samples_nominal = permute(samples, [2 1 3]);
        
        % To compress thetaj1 into a nsamples x Dimension vector
        thetaj1 = samples_nominal(:,:)';
       
        % According to Cheung and Beck (2009) - Bayesian model updating ...,
        % the initial samples from reweighting and the resample of samples of
        % fj, in general, do not exactly follow fj1, so that the Markov
        % chains must "burn-in" before samples follow fj1, requiring a large
        % amount of samples to be generated for each level.
        
        %% Adjust the acceptance rate (optimal = 23%)
        % See: http://www.dms.umontreal.ca/~bedard/Beyond_234.pdf
        %{
      if acceptance_rate < 0.3
         % Many rejections means an inefficient chain (wasted computation
         %time), decrease the variance
         beta = 0.99*beta;
      elseif acceptance_rate > 0.5
         % High acceptance rate: Proposed jumps are very close to current
         % location, increase the variance
         beta = 1.01*beta;
      end
        %}
    
    fprintf('\n');
    acceptance(count) = acceptance_rate;
    
    %% Prepare for the next iteration
    c_a = (acceptance_rate - ((0.21./Dimensions) + 0.23));
    stepsize_nominal = stepsize.*exp(c_a);
    
    if stepsize_nominal <= 1
    stepsize = 1.01;
    else
    stepsize = stepsize_nominal;
    end
    
    count = count+1;
    samps(:,:,count) = thetaj1;
    step(count) = stepsize;
    thetaj = thetaj1;
    pj     = pj1;
    beta_j(count) = pj;
end

% estimation of f(D) -- this is the normalization constant in Bayes
log_fD = sum(log(S(1:j)));

%% Description of outputs:

output.samples = thetaj;         % To only show samples from the final posterior
output.allsamples = samps;       % To show all samples from the initial prior to the final posterior
output.log_evidence = log_fD;    % To generate the logarithmic of the evidence
output.acceptance = acceptance;  % To show the mean acceptance rates for all iterations
output.beta = qsj;            % To show the values of temepring parameters, qsj 
output.step = step;              % To show the values of step-size

return; % End


%% Calculate the plusibility q(j+1)
function qj1 = calculate_qj1(log_fD_T_thetaj, qj)

thres = 1; 

wj = @(e) exp(abs(e)*log_fD_T_thetaj); % N x 1
fmin = @(e) std(wj(e)) - threshold*mean(wj(e)) + realmin;
e = abs(fzero(fmin, 0)); % e is >= 0, and fmin is an even function
if isnan(e)
    error('There is an error finding e');
end

qj1 = ((MDS(j))^(N0-1))*( min(1, qj + e);

return; % End

storez  = zeros(M,Nsim); % phi1
storephi  = zeros(N*M,Nsim); % phi1
storetheta = zeros(P,Nsim);
storesig2z = zeros(1,Nsim);
storesig2phi = zeros(1,Nsim);
storesig2theta = zeros(1,Nsim);
storelambda= zeros(M,Nsim);
storeh = zeros(M,Nsim);

% Store Value
storetheta(:,1) = theta0(:,1); % theta
storesig2z(1) = sig2z0; % sigma1
storesig2phi(1)  = sig2phi0; % sigma2
storesig2theta(1) = sig2theta0;
storephi(:,1)  = phi0; % phi1
for i=1:M
    storelambda(i,1) = gamrnd(nu(i)/2,nu(i)/2);
end
storeh(:,1) = gamrnd(eta/2,eta/2,[M,1]);
storez(:,1) = z0;

% Initial values
sig2z = sig2z0;
sig2phi = sig2phi0;
sig2theta = sig2theta0;
theta = truetheta + randi([10 100],N,1);
lambda = storelambda(:,1);
h = storeh(:,1);
[~, z] = priorformodal(theta,modem,Mmat,Nset);

tic
hwait = waitbar(0,'Simulation in Progress');
for n = 1:Nsim-1
    % draw phi conditional on phicap, f, sigphisqr, h, theta, sigthetasqr
    [phi] = drawphi(phicap, z, sig2phi, h, theta, sig2theta, Nset, Gamma, Mmat);
    
    % draw f conditional on fcap, phi, sigfsqr, lambda, theta, sigthetasqr
    [z] = drawz(zcap, phi, sig2z, lambda, theta, sig2theta, Nset, Mmat);
    
    % draw theta conditional on f, phi, sigthetasqr
    [theta] = drawtheta(z, phi, sig2theta, Mmat, invB0, theta0, Nset);
    
    % draw sigfsqr conditional on fcap, f, lambda
    [sig2z] = drawsig2z(zcap, z, lambda, alpha1, delta0, Nset);
    
    % draw sigphisqr conditional on phicap, phi, h
    [sig2phi] = drawsig2phi(phicap, phi, Gamma, h, p1, q0, Nset);
    
    % draw sigthetasqr conditional on f, , phi, theta
    [sig2theta] = drawsig2theta(phi,z,theta,Mmat,Nset,b0,a1);
    
    % draw lambda conditional on fcap, f, sigfsqr
    [lambda] = drawlambda(zcap, z, sig2z, nu, Nset);
    
    % draw h conditional on phicap, phi, sigphisqr
    [h] = drawh(phicap, phi, sig2phi, Gamma, eta, Nset);
    
    storetheta(:,n+1) = theta;
    storephi(:,n+1) = phi;
    storez(:,n+1) = z;
    storesig2phi(n+1) = sig2phi;
    storesig2z(n+1) = sig2z;
    storelambda(:,n+1) = lambda;
    storesig2theta(n+1) = sig2theta;
    storeh(:,n+1) = h;
    waitbar(n/Nsim);
end
close(hwait)
toc

%%% Post Burn Value
postphi = storephi(:,burnin+1:Nsim);
postz = storez(:,burnin+1:Nsim);
posttheta = storetheta(:,burnin+1:Nsim);
postsig2z = storesig2z(burnin+1:Nsim);
postsig2phi = storesig2phi(burnin+1:Nsim);
postsig2theta = storesig2theta(burnin+1:Nsim);
postlambda = storelambda(:,burnin+1:Nsim);
posth = storeh(:,burnin+1:Nsim);

clear results
results(:,1) = truetheta;
results(:,2) = mean(posttheta,2);
results(:,3) = std(posttheta,0,2);
results(:,4) = results(:,3)./results(:,2);
results(:,5) = abs(results(:,1) - results(:,2))./results(:,1);
disp(['  --------------------------------------------------'])
disp(['      True      Mean     s.d.      c.o.v.      Error'])
disp(['  --------------------------------------------------'])
disp([results])
atheta = results;

clear results
zmodel = w2r;
results(:,1) = zmodel(modem);
results(:,2) = mean(postz,2);
results(:,3) = std(postz,0,2);
results(:,4) = results(:,3)./results(:,2);
results(:,5) = abs(results(:,1) - results(:,2))./results(:,1);
disp(['  --------------------------------------------------'])
disp([' True_z   Mean_z    s.d.      c.o.v.    error'])
disp(['  --------------------------------------------------'])
disp([results])
az = results;

clear results
results(1) = mean(postsig2z);
results(2) = std(postsig2z);
results(:,3) = results(:,2)./results(:,1);
disp(['  --------------------------------------------------'])
disp(['    sig2z    s.d.     c.o.v.'])
disp(['  --------------------------------------------------'])
disp([results])


clear results
results(1) = mean(postsig2phi);
results(2) = std(postsig2phi);
results(:,3) = results(:,2)./results(:,1);
disp(['  --------------------------------------------------'])
disp(['    sig2phi    s.d.     c.o.v.'])
disp(['  --------------------------------------------------'])
disp([results])

clear results
results(1) = mean(postsig2theta);
results(2) = std(postsig2theta);
results(:,3) = results(:,2)./results(:,1);
disp(['  --------------------------------------------------'])
disp(['    sig2theta    s.d.     c.o.v.'])
disp(['  --------------------------------------------------'])
disp([results])
asig2theta= results;

clear results
results(:,1) = mean(postlambda,2);
results(:,2) = std(postlambda,0,2);
results(:,3) = results(:,2)./results(:,1);
disp(['  --------------------------------------------------'])
disp(['    lambda    s.d.    c.o.v.'])
disp(['  --------------------------------------------------'])
disp([results])

clear results
results(:,1) = mean(posth,2);
results(:,2) = std(posth,0,2);
results(:,3) = results(:,2)./results(:,1);
disp(['  --------------------------------------------------'])
disp(['    h    s.d.    c.o.v.'])
disp(['  --------------------------------------------------'])
disp([results])

clear results
results= (1./mean(posth,2))*(mean(postsig2phi,2));
disp(['  --------------------------------------------------'])
disp(['   sig2phi  '])
disp(['  --------------------------------------------------'])
disp([results])
asig2phi = results;

clear results
results = (1./mean(postlambda,2))*(mean(postsig2z,2));
disp(['  --------------------------------------------------'])
disp(['   sig2z  '])
disp(['  --------------------------------------------------'])
disp([results])
asig2z = results;

% clear results
% results = sqrt(sig2z)./mean(postz,2);
% disp(['  --------------------------------------------------'])
% disp(['   c.o.v of z'])
% disp(['  --------------------------------------------------'])
% disp([results])

% Calculation of MAC value
U = zeros(S,M);
meanphi = mean(postphi,2);
for i = 1:M
    phip = meanphi((i-1)*N+1:i*N);
    U(:,i) = Gamma * phip;
end

meancap = mean(phicap,2);
Mac = zeros(M,1);
for m = 1:M
    phiexp = meancap((m-1)*S+1:m*S);
    Mac(m) = ((phiexp'*(U(:,m))))^2/(((U(:,m))'*U(:,m))*((phiexp)'*phiexp));
end
disp(Mac')

---- plots for modeshapes ------%

figure()
for ii=1:M
    phiexp = meancap((ii-1)*S+1:ii*S);
    subplot(1,M,ii)
    plot([0;U(:,ii)],([0 sen]),'-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])%,'LineWidth',2.,'Marker','.',...
        %'MarkerSize',20,'Color',[0 1 0],'markeredgecolor','k')
    hold on
    plot([0;phiexp],([0 sen]),'k--o')%,'LineWidth',2.,'Marker','.',...
        %'MarkerSize',20,'Color',[0 1 0],'markeredgecolor','b')
    xlim([-max(abs(U(:,ii))) max(abs(U(:,ii)))])
    grid on
    xlabel(['\phi_',num2str(modem(ii))]);
    title(['Mode ',num2str(modem(ii))],'FontSize',18,'FontWeight','normal')
    legend('Identified \phi','Actual \phi')
    set(gca,'FontSize',25)
end

--------- Trace Plots --------------- %
s ='a':'z';
figure()
for jj = 1:9
    subplot(3,3,jj)
%     plot(jj,'Parent',subplot(3,3,jj),'LineWidth',2,'Color',[0 0 0]);
    plot(storetheta(jj,:),'LineWidth',0.5,'Color',[0 0 0]); %line([0 Nsim],[truepara(jj) truepara(jj)])
    %     subplot(Npara,2,(jj-1)*2+2); hist(thetastorage(jj,:),100)
    %     title(['Traceplot for \theta_',num2str(jj)]);
%     xlabel('simulations')
    set(subplot(3,3,jj),'FontSize',7,'FontWeight','normal','LineWidth',1);
    ylabel([])
    xlim([0 Nsim])
    xlabel(['(',s(jj),')'])
    title(['\theta_',num2str(jj)],'FontWeight','normal')
    set(gca,'FontWeight','normal')
end
print -deps partial_traceplot_shear_model1
% print -deps full_traceplot_shear_model1
% 
%---- Histogram of theta ------%
data_set = [ postlambda; posth];
figure()
xlabs = {'z_1', 'z_2', 'z_3', '\lambda_1', '\lambda_2', '\lambda_3'...
    , 'h_1', 'h_2', 'h_3', '\sigma^2_z', '\sigma^2_\phi', '\sigma^2_\theta'};
for k = 1:(3+3)
    subplot(2,M,k)
    hist(data_set(k,:))
    xlabel(['(',s(k),')']);
    set(subplot(2,M,k),'FontSize',7,'FontWeight','normal','LineWidth',1)
    h = findobj(gca,'Type','patch');
%     set(gca,'FontWeight','bold')
    h.FaceColor = [1 1 1];
    h.EdgeColor = 'k';
    h.LineWidth = 0.5;
end
% print -deps partial_model1
print -deps full_hist_shear_model1

save('AllParameters.mat','postphi','postz','posttheta','postsig2z',...
    'postsig2phi','postlambda','posth','storetheta')

