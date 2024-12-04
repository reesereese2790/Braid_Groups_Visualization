clear 
clc 
addpath('C:\Users\arm20dh\OneDrive\Documents\GitHub\Braid_Groups_Visualization\braidlab-master');
warning('off','BRAIDLAB:braid:entropy:noconv'); 

% initialize variables for ftbe vs. sigma graph
M = 100;
sigma_slider = 1:1:10;
ftbe_array = zeros(M,length(sigma_slider));

for v = 1:M
    for p=1:length(sigma_slider)
        % initialize variables for braid diagram
        N = 10;
        x0= [0:1/(N-1):1;zeros(1,N)]; 
        mu=zeros(2,1); 
        sigma = sigma_slider(p);
        sigma_matrix=sigma*eye(2);    
        T0=0; dt=0.01; Tf=1; t=T0:dt:Tf;
        num_steps = length(t);
        trajectories = zeros(N,2,length(t));
        
        % simulate Brownian motion
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
            
            % store the trajectory for the i-th path
            trajectories(i, :, :) = x;
        end
        
        % permute the dimensions to be compatible with a braid
        trajectories = permute(trajectories,[3,2,1]);
        
        % create databraid
        warning ('off ','BRAIDLAB : braid : colorbraiding : notclosed ')
        b = braidlab.databraid(trajectories);
        b = compact(b);
        ftbe_array(v,p) = ftbe(b);
    end
end

figure
plot(sigma_slider,mean(ftbe_array,1))
xlabel('Volatility')
ylabel('Finite Time Braiding Exponent')
title('FTBE vs. Volatility')