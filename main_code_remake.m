clear 
clc 
addpath('./braidlab-master/+braidlab');
warning('off','BRAIDLAB:braid:entropy:noconv'); 

% iterations of generating different random entropies 
M=20; 

% volatility
sigma=3; 

% iterating number of strands
N_initial=2; 
N_final=20; 
N_increment=1; 
N_slider=N_initial:N_increment:N_final; 

% initialize entropy values and averages
e_array=zeros(M,length(N_slider)); 
average_line=zeros(1,length(N_slider)); 

%initialize plot 
figure 
xlabel('Number of Strands') 
ylabel('Entropy of Braid (1/s)') 
title('entropy vs. number of strands') 
hold on 


N_index=1; 
for N=N_slider 
    
    x0b= [0:1/(N 1):1;zeros(1,N)]; 
    %initialize starting points 
    x0 = permute(x0b, [3,1,2]); %rearrange dimensions, instead of spreadsheet boxes use a 3 rd dimension to better represent the slices of which is starting points 
    constraints=[x0;x0]; %start and end at same point 
    
    %model for brownian motion, mu=drift, sigma=volatility
    Mu=zeros(2,1); 
    Sigma=sigma*eye(2);  
    brown=bm(Mu,Sigma); 
    
    %generate 1000 points between 0 and 1
    T0=0; dt=0.01; Tf=1; t=T0:dt:Tf; 
end