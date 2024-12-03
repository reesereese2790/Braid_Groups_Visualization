clear 
clc 
%addpath('./braidlab-master/+braidlab');
addpath('C:\Users\arm20dh\OneDrive\Documents\GitHub\Braid_Groups_Visualization\braidlab-master');
warning('off','BRAIDLAB:braid:entropy:noconv'); 

% GRAPHS TO REPLICATE
% create FTBE vs. sigma 
% create topological entropy vs. sigma
% create number of generators vs. sigma
% create number of generators vs. number of strands
% create number of generators vs. time 

% FTBE vs. sigma
% initialize variables for braid diagram
N = 10;
x0= [0:1/(N-1):1;zeros(1,N)]; 

%model for brownian motion, mu=drift, sigma=volatility
mu=zeros(2,1); 
sigma = 10;
sigma_matrix=sigma*eye(2);    
T0=0; dt=0.01; Tf=1; t=T0:dt:Tf;
num_steps = length(t);
trajectories = zeros(N,2,length(t));

% Simulate Brownian motion with constraints
for i = 1:N
    % Brownian motion initialization (start at x0b(i))
    x = zeros(2, num_steps);
    x(:, 1) = x0(:, i);  % Set the initial position as x0b(i)
    for j = 2:num_steps
        % Brownian motion increment
        dB = sqrt(dt) * randn(2, 1); 
        
        % Brownian motion update with drift and volatility
        x(:, j) = x(:, j-1) + mu * dt + sigma_matrix * dB;
    end
    
    % Enforce the constraint of the endpoint being equal to x0b(i)
    final_point = x0(:, i); % The endpoint should be x0b(i)
    x(:, end) = final_point;  % Set the final point directly to x0b(i)
    
    % Store the trajectory for the i-th path
    trajectories(i, :, :) = x;
end

% Plot the generated random Brownian motion trajectories
figure;
hold on;
for i = 1:N
    plot(t, squeeze(trajectories(i, 1, :)));
end
xlabel('Time');
ylabel('Position');
title('Brownian Motion Trajectories with Constraints');
hold off;

% encode braid, make sure trajectories is the right input
% b = braidlab.braid(trajectories);



% initialize variables for ftbe vs. sigma graph
% sigma_slider = 1:1:50;
% ftbe_array = zeros(1,length(sigma_slider));

% find the FTBE

% plot
% figure 
% xlabel('Volatility')
% ylabel('Finite Time Braiding Exponent')
% title('FTBE vs. Volatility')
% average over several iterations





% OLD CODE
% 
% % iterations of generating different random entropies 
% M=20; 
%
% % iterating number of strands
% N_initial=2; 
% N_final=20; 
% N_increment=1; 
% N_slider=N_initial:N_increment:N_final; 
% 
% % initialize entropy values and averages
% e_array=zeros(M,length(N_slider)); 
% average_line=zeros(1,length(N_slider)); 
