clear 
clc 
addpath('C:\Users\arm20dh\OneDrive\Documents\GitHub\Braid_Groups_Visualization\braidlab-master');
warning('off','BRAIDLAB:braid:entropy:noconv'); 

% M is the number of iterations we are averaging over, N is strands
% sigma slider is the x-axis for the volatility we are varying
function [ftbe_matrix,top_entropy_matrix,count_matrix] = vary_strands(M,N,sigma_slider)
    ftbe_matrix = zeros(M,length(sigma_slider));
    top_entropy_matrix = zeros(M,length(sigma_slider));
    count_matrix = zeros(M,length(sigma_slider));
    
    for v = 1:M
        for p=1:length(sigma_slider)
            % initialize variables for braid diagram
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
            b_data = braidlab.databraid(trajectories);
            b_braid = b_data.braid;
            count_matrix(v,p) = b_braid.length;
            ftbe_matrix(v,p) = ftbe(compact(b_data));
            top_entropy_matrix(v,p) = entropy(b_braid);
            
            
            
        end
    end
end

sigma_slider = 1:1:20;
[a,b,c] = vary_strands(10,10,1:1:20);
[d,e,f] = vary_strands(10,20,1:1:20);
[g,h,i] = vary_strands(10,30,1:1:20);

figure
plot(sigma_slider,mean(a,1),'r','LineWidth', 2)
hold on;
plot(sigma_slider,mean(d,1),'b','LineWidth', 2)
hold on;
plot(sigma_slider,mean(g,1),'g','LineWidth', 2)
hold on;
legend('N=10','N=20','N=30');
xlabel('Volatility of Trajectory')
ylabel('Finite Time Braiding Exponent')
title('FTBE vs. Volatility')

figure
plot(sigma_slider,mean(b,1),'r','LineWidth', 2)
hold on;
plot(sigma_slider,mean(e,1),'b','LineWidth', 2)
hold on;
plot(sigma_slider,mean(h,1),'g','LineWidth', 2)
hold on;
legend('N=10','N=20','N=30');
xlabel('Volatility of Trajectory')
ylabel('Topological Entropy')
title('Topological Entropy vs. Volatility')

figure
plot(sigma_slider,mean(c,1),'r','LineWidth', 2)
hold on;
plot(sigma_slider,mean(f,1),'b','LineWidth', 2)
hold on;
plot(sigma_slider,mean(i,1),'g','LineWidth', 2)
hold on;
legend('N=10','N=20','N=30');
xlabel('Volatility of Trajectory')
ylabel('Braid Length')
title('Braid Length vs. Volatility')

% to show all points and not only average
% plot(sigma_slider,vary_strands(10,10),'o', 'MarkerSize', 1)