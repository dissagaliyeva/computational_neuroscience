%% ordinary differential equations manual
x0 = 2; % first value of x 
t0 = 0; % first time point 
dt = 0.001; % time step for simulations 
tmax = 20;  % max time length
tvec = t0:dt:tmax; % vector of time points 
disp(tvec)

x = zeros(size(tvec)); 
x(1) = x0; % value of x at t=0

for i = 2:length(tvec)
    tval = tvec(i);
    dxdt = 2 * sin(0.1 * tval * tval); % rate of change of x, dx/dt 
    x(i) = x(i - 1) + dt * dxdt; % forward Euler
end

plot(tvec, x)

%% ordinary differential equations in-built functionality
[t, x] = ode45(@test_ode, tvec, x0); 
figure(2)
plot(t, x)


%% coupled (=connected) ODEs with multiple variables (p 150)

tmax = 100;
x0 = 1.0;
v0 = 0.0;                   % initial velocity at t=0
omega = 3;                  % angular frequency
omega_sqr = omega^2;        % square of angular frequency
dt = 0.00001;               % small ft for forward euler
tvector = 0:dt:tmax;        % vector of time points
x1 = zeros(size(tvector));  % vector of position values
x2 = zeros(size(tvector));  % vector of velocity values

x1(1) = x0;
x2(1) = v0;

% for loop
for i = 2:length(tvector)
    x1(i) = x1(i - 1) + x2(i - 1) * dt; 
    x2(i) = x2(i - 1) + omega_sqr * x1(i - 1) * dt;
end

figure(1)
plot(tvector, x1)
figure(2)
plot(tvector, x2)


%%
















