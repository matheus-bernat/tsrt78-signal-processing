%% Settings
M=1000;               % Number of simulations
N=100;              % Number of samples per realization
T=1;                % Sample period [s]


% True model parameters
A=[1 T; 0 1];       % System matrix
G=[T^2/2; T];       % Noise gain matrix
C=[1 0];            % Observation matrix

Q=0.1^2;            % Process noise [(m/s^2)^2]
R=1^2;              % Measurement noise [m^2]

PI0=diag([1^2 0.1^2]);   % Initial state covariance


% Filter parameters
A_kf=[1 T; 0 1];       % System matrix
G_kf=[T^2/2; T];       % Noise gain matrix
C_kf=[1 0];            % Observation matrix

Q_kf=0.1^2;            % Process noise [(m/s^2)^2]
R_kf=1^2;              % Measurement noise [m^2]

PI0_kf=diag([1^2 0.1^2]);   % Initial state covariance



% Allocate memory for Monte Carlo loop
MSE=zeros(2,N);
for m=1:M
    
    % Display iteration number
    disp(['Simulation ' num2str(m)])
    
    %% Simulate the system
    x=zeros(2,N);   % Allocate memory for the states
    y=zeros(1,N);   % Allocate memory for the observations
    
    x(:,1)=chol(PI0)*randn(2,1); % Draw a random initial state estimate
    
    % Loop
    for k=1:N-1
        
        % Generate a process noise sample
        w=sqrt(Q)*randn;
        
        % Propagate the state
        x(:,k+1)=A*x(:,k)+G*w;
        
        % Generate a measurement noise sample
        v=sqrt(R)*randn;
        
        % Observations
        y(k+1)=C*x(:,k+1)+v;
    end
    
    %% Kalman filter
    x_hat=zeros(2,N);       % Allocate memory for the state estimates
    x_hat(:,1)=zeros(2,1);  % Initial state estimate
    
    P=zeros(2,2,N);         % Allocate the memory for the state error covariance
    P(:,:,1)=PI0_kf;           % Initial covariance
    
    
    for k=1:N-1
        
        % Time update
        x_hat(:,k+1)=A_kf*x_hat(:,k);
        P(:,:,k+1)=A_kf*P(:,:,k)*A_kf'+G_kf*Q_kf*G_kf';
        
        % Measurement update
        K=P(:,:,k+1)*C_kf'/(C_kf*P(:,:,k+1)*C_kf'+R_kf);
        x_hat(:,k+1)=x_hat(:,k+1)+K*(y(k+1)-C*x_hat(:,k+1));
        P(:,:,k+1)=P(:,:,k+1)-(P(:,:,k+1)*C_kf'/(C_kf*P(:,:,k+1)*C_kf'+R_kf))*C_kf*P(:,:,k+1);
    end
    
    
    %% Statistics
    MSE=MSE+(x-x_hat).^2;
end

MSE=MSE./M;

figure(3)
clf
subplot(2,1,1)
plot(squeeze(P(1,1,:)),'k')
hold on
plot(MSE(1,:),'r')
legend({'Kalman filter calculate $[P_k]_{1,1}$','Simulated (true) $E[(x^{(1)}_k-\hat{x}^{(1)}_k)^2]$'},'Interpreter','latex','FontSize',12)
box on;
grid on;
ylabel('Position error [m^2]')
title('Position error vs time')
xlim([1 N])

subplot(2,1,2)
plot(squeeze(P(2,2,:)),'k')
hold on
plot(MSE(2,:),'r')
legend({'Kalman filter calculate $[P_k]_{2,2}$','Simulated (true) $E[(x^{(2)}_k-\hat{x}^{(2)}_k)^2]$'},'Interpreter','latex','FontSize',12)
box on;
grid on;
ylabel('Speed error [(m/s)^2]')
title('Speed error vs time')
xlim([1 N])
