close all;
clear all;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reffered paper: [A] "Computation of Information Rates by Particle Methods", J
% Dauwels, Hans Andrea Loeliger.
%There is another paper: [B] "Simulation based computation of Information Rates
%for Channels with Memory", by D M Arnold, Pascal Vontoebel, Hans Andrea Loeliger. 
%The reffered paper is an extension of this latter work.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE 1: This program computes the entropy rate of the  output  sequence of
%an AWGN channel. The channel is powered by an energy harvester with an
%infinite energy buffer. The channel input is sampled from Gauusian
%distribution with mean 0 and variance rho. The sampled vector is projected back to the
%constrained vector according to some rule viz. the causal neighbour
%projection rule.%%%%

%NOTE 2: Causal Projection Rule: Project the components of the vector to the 
%extreme value along the axis by retaining the polarity.




%%%Initialization%%%

rho_vector=0:1:5   % harvested energy Y_k = rho a.s.
eps1=0.001; %Epsilon factor which helps in the positive drift. Randomly chosen.
eps2=.3; %Epsilon factor deciding the "effective" particles. Chosen as given in the reffered paper.
nv=1; %Noise variance
N=10^3; %Number of particles. Value chosen as in the reffered paper.
num_iter=20 %This is the n in \frac{1}{n}h(W^n). W^n is the channel output vector. 
average_n=500; % \frac{1}{n}h(W^n) is averaged for average_n times.
R_causal_gaussian=zeros(1,length(rho_vector));

%%%%%%%%%%%%%%%%%%%

for loop_r=1:length(rho_vector)

    rho=rho_vector(loop_r);
    
    entrpy_Wn=0;

    
for loop_a=1:average_n
    
%Initialization%%%%%%%%%%%%%

mu_k=zeros(2,N);
mu_k(2,:)=1/N;  %Intial weights issued
extnd_mu_k=zeros(3,N); % Extended particle. See the reference.
emp_entrpy_W=0; %Empirical output entropy viz. -log(P(W^n). Note that it is NOT sclaed by \frac{1}{n}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

x_kv=normrnd(0,rho,1,num_iter);%Sampling an "num_iter" length input vector from Gaussian distribution with mean zero and variance rho.
x_kv=causal_projection(x_kv,rho );%Projecting the vector back to the constraint set using nearest neighbour projection rule.


for loop=1:num_iter
    
    %%%Channel output and energy buffer updation%%%%%
    
    x_k=x_kv(loop);%Channel input at time k
    n_k=normrnd(0,nv);%Additive noise at time k
    w_k=x_k+n_k;%Channel output at time k.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Message passing algorithm part%%%%%%%%%%%%%%%%
    
    %The algorith is given in [A]
    
    extnd_mu_k(1,:)=mu_k(1,:);
    x_k_prime=normrnd(0,rho,1,N);
    y_k=rho*ones(1,N);
    eta_k=y_k-(x_k_prime).^2;
    temp_mu=mu_k(1,:)+eta_k;
    extnd_mu_k(2,:)=max(0,temp_mu);%Update of energy buffer
    temp_mu=mu_k(1,:)+y_k;
    extnd_mu_k(3,:)=sign(x_k_prime).*min(sqrt(temp_mu),abs(x_k_prime));
    temp_lambda=(1/sqrt(2*pi*nv))*exp((-.5/nv)*(w_k-extnd_mu_k(3,:)).^2);% To see full expression of 1/lambda_k see the reffered paper. 
    lambda_inv=dot(temp_lambda,mu_k(2,:));
    lambda_k=1/lambda_inv;
    mu_k(1,:)=extnd_mu_k(2,:);
    temp_weight=mu_k(2,:);
    mu_k(2,:)=lambda_k*(temp_lambda.*temp_weight);
    temp_mu=mu_k(2,:);
    
    if 1/dot(temp_mu,temp_mu) < eps2*N
    temp_mu=rand_gen(mu_k(1,:),mu_k(2,:),N); %Generates a random vector of length N  with support mu_k(1,:) and probabilities mu_k(2,:).
    mu_k(1,:)=temp_mu;
    mu_k(2,:)=1/N;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    emp_entrpy_W=emp_entrpy_W+log(lambda_k);
end

h_W_n=emp_entrpy_W/num_iter;
entrpy_Wn=entrpy_Wn+h_W_n;


end


h_W_n=entrpy_Wn/average_n;
R_causal_gaussian(loop_r)=h_W_n-.5*log(2*pi*exp(1)*nv);%In nats. 


end

R_causal_gaussian
C_infty=.5*log(1+rho_vector/nv)                   %Infinite_buffer capacity. 

toc   
