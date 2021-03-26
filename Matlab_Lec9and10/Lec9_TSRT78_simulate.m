%% Settings
N=100;              % Number of samples per realization
T=1;                % Sample period [s]

A=[1 T; 0 1];       % System matrix
G=[T^2/2; T];       % Noise gain matrix
C=[1 0];            % Observation matrix

Q=0.01^2;            % Process noise [(m/s^2)^2]
R=1^2;              % Measurement noise [m^2]
    
PI0=diag([1^2 0.1^2]);   % Initial state covariance         



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

figure(1)

subplot(3,1,1)
plot(y,'k*')
hold on;
box on;
grid on;
ylabel('Position [m]')
xlabel('Sample $k$','Interpreter','latex')
title('Measurement: $y_k$','Interpreter','latex','FontSize',12)
xlim([1 N])

subplot(3,1,2)
plot(x(1,:),'k')
hold on;
box on;
grid on;
ylabel('Position [m]')
xlabel('Sample $k$','Interpreter','latex')
title('State: $x^{(1)}_k$','Interpreter','latex','FontSize',12)
xlim([1 N])

subplot(3,1,3)
plot(x(2,:),'k')
hold on;
box on;
grid on;
ylabel('Speed [m]')
xlabel('Sample $k$','Interpreter','latex')
title('State: $x^{(2)}_k$','Interpreter','latex','FontSize',12)
xlim([1 N])


