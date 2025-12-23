% System Parameters
sigma2_n = 1e-3;         % Noise variance  
sigma2_x = 1e-1;         % Signal variance
sampleInterval = 0.8e-6; % Sampling frequency is 1.25MHz
N = 500;                % Number of samples
dopplerFrequency = 1000; % Doppler frequency
m = 5;                   % Channel length
ENS = 500;               % Ensemble

%Kalman parameters
alpha = 0.9999;
q1 = .1;
q2 = 1;
lambda = .9;

% Initialization
MSEexrls = zeros(N,1);
MSEhexrls = zeros(N,1);
MSEgexrls = zeros(N,1);
MSErls = zeros(N,1);
mseEXRLS = zeros(N,1);
mseHEXRLS = zeros(N,1);
mseGEXRLS = zeros(N,1);
mseRLS = zeros(N,1);
epsilon = 1/sigma2_x;
P_R = zeros(m,m,N);
P_H = zeros(m,m,N);
P_G = zeros(m,m,N);
channel  = randn(m,N);

%  Nonlinearity
typeNonlinear = 1;
paramNonlinear = 2;

for ens = 1:ENS
    disp(ens);

    %% Data Formatting
    % Rayleigh Fading Channel
    for i=1:m
        channel(i,:) = rayleigh(sampleInterval,N,dopplerFrequency);
    end

    % Type of signal 
    complex = 1; % 0 for a real signal

    if complex == 0
        channel = real(channel);
        inputSignal = randn(1,N+m);% inputSignal = inputSignal - mean(inputSignal); inputSignal = sqrt(sigma2_x)*inputSignal/std(inputSignal); 
        noise = randn(N,1); noise = noise - mean(noise); noise = sqrt(sigma2_n)*noise/std(noise);
    else
        inputSignal = randn(1,N+m) + 1j*randn(1,N+m); 
        noise = (randn(N,1) + 1j*randn(N,1)); noise = noise - mean(noise); noise = sqrt(sigma2_n)*noise/std(noise);
    end

    U = zeros(m,N);      
    for kk = 1:N
        U(:,kk) = inputSignal(kk:kk+m-1);
    end

    %Desired training signal
    d = zeros(N,1);
    for ii=1:N
        d(ii) = U(:,ii)'*channel(:,ii);
    end
    d = d + noise;

    %Pass through the nonlinearity
    % d = nlG(d,paramNonlinear,typeNonlinear);

    % type = 0 => e*(n) = y*(n) - u^H(n)w(n)
    % type = 1 => e(n)  = y(n)  - w^H(n)u(n)
    type = 0;


    %% RLS
    [~,mseRLS,~] = RLS(U, d, lambda, epsilon, type);
    MSErls = MSErls + mseRLS;
 
    %% EX-RLS
    [~,mseEXRLS,~,P_R] = EXRLS(U, d, lambda, q1, q2, alpha, epsilon, type);
    MSEexrls = MSEexrls + mseEXRLS;

    %% G-EX-RLS
    % aprox = 0 => Cholesky
    % aprox = 1 => Square-root aproximation
    % aprox = 2 => Taylor series aproximation
    aprox = 0;
    [~,mseGEXRLS,~,P_G] = G_EX_RLS(U, d, lambda, q1, q2, alpha, epsilon, aprox, type);
    MSEgexrls = MSEgexrls + mseGEXRLS;

    %% H-EX-RLS
    [~,mseHEXRLS,~,P_H] = H_EX_RLS(U, d, lambda, q1, q2, alpha, epsilon, aprox, type);
    MSEhexrls = MSEhexrls + mseHEXRLS;
end

MSEdBexrls = 10*log10(MSEexrls/ENS);
MSEdBhexrls = 10*log10(MSEhexrls/ENS);
MSEdBgexrls = 10*log10(MSEgexrls/ENS);
MSEdBrls = 10*log10(MSErls/ENS);

figure
plot(10 * log10(sigma2_n) * ones(N, 1), '--r', 'linewidth', 1)
hold on
plot(1:N,MSEdBrls,'m','LineWidth',2)
plot(1:N,MSEdBexrls,'k','LineWidth',2)
plot(1:N,MSEdBgexrls,'--g','LineWidth',2)
plot(1:N,MSEdBhexrls,'-.b','LineWidth',1)
legend('Noise floor','RLS','EX-RLS','G-EX-RLS','H-EX-RLS')
grid
axis tight
set(gca, 'FontSize', 11);
set(gca, 'FontName', 'Times-Roman');
xlabel('n'),ylabel('MSE (dB)')

energyR = zeros(N,1);
energyH = zeros(N,1);
energyG = zeros(N,1);
for n = 1:N
    % Energy = Squared Frobenius norm
    energyR(n) = norm(P_R(:,:,n),'fro')^2;
    energyH(n) = norm(P_H(:,:,n),'fro')^2;
    energyG(n) = norm(P_G(:,:,n),'fro')^2;
end

figure;
plot(10*log10(energyR),'k','LineWidth',2); 
hold on;
plot(10*log10(energyG),'b','LineWidth',2);
plot(10*log10(energyH),'r--','LineWidth',2);
xlabel('n'); ylabel('Energy (dB)');
legend('EX-RLS','G-EX-RLS','H-EX-RLS');
grid on;
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14);

