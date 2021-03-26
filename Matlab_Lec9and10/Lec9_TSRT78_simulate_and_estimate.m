%% Settings
N=100;              % Number of samples per realization
T=1;                % Sample period [s]


% True model parameters
A=[1 T; 0 1];       % System matrix
G=[T^2/2; T];       % Noise gain matrix
C=[1 0];            % Observation matrix

Q=0.01^2;            % Process noise [(m/s^2)^2]
R=1^2;              % Measurement noise [m^2]

PI0=diag([1^2 0.1^2]);   % Initial state covariance


% Filter parameters
A_kf=[1 T; 0 1];       % System matrix
G_kf=[T^2/2; T];       % Noise gain matrix
C_kf=[1 0];            % Observation matrix

Q_kf=0.01^2;            % Process noise [(m/s^2)^2]
R_kf=1^2;              % Measurement noise [m^2]

PI0_kf=diag([1^2 0.1^2]);   % Initial state covariance



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

P=zeros(2,2,N);             % Allocate the memory for the state error covariance
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


figure(2)
clf


subplot(2,1,1)
sigma=3*sqrt(squeeze(P(1,1,:)))';
t=(1:length(sigma));
fill([t fliplr(t)],[x_hat(1,:)+sigma fliplr(x_hat(1,:)-sigma)],'k','LineStyle','none')
alpha(0.25)
hold on
plot(x(1,:),'k')
plot(x_hat(1,:),'r--')
plot(y(1,:),'b*')
legend({'$\pm 3\sqrt{[P_k]_{1,1}}$','$x^{(1)}_k$','$\hat{x}^{(1)}_k$','$y_k$'},'Interpreter','latex','FontSize',12)
box on;
grid on;
ylabel('Position [m]')
title('Position vs time')
xlim([1 N])

subplot(2,1,2)
sigma=3*sqrt(squeeze(P(2,2,:)))';
t=(1:length(sigma));
fill([t fliplr(t)],[x_hat(2,:)+sigma fliplr(x_hat(2,:)-sigma)],'k','LineStyle','none')
alpha(0.25)
hold on
plot(x(2,:),'k')
plot(x_hat(2,:),'r--')
legend({'$\pm 3\sqrt{[P_k]_{2,2}}$','$x^{(2)}_k$','$\hat{x}^{(2)}_k$'},'Interpreter','latex','FontSize',12)
box on;
grid on;
xlabel('Time [s]')
ylabel('Speed [m/s]')
title('Speed vs time')
xlim([1 N])
