clc; clear; close all

A = load('adjacency_277_reorder.mat');
A = A.mydata;
A = double(A);

%%

deg = sum(A,2);
[N_k,k] = groupcounts(deg);
% figure
% scatter(k,N_k,'red','filled')
% 
% set(gca,'Xscale','log')
 
[alpha, xmin, ~]=plfit(deg);
plplot(deg, xmin, alpha)

%%

c1 = 1:43;
c2 = 44:112;
c3 = 113:132;
c4 = 133:240;
c5 = 241:258;
c6 = 259:277;

N = size(A,2);

E = zeros(N); %intra-community
T = zeros(N); %inter-community

% l = [length(c1) length(c2) length(c3) length(c4) length(c5)];

l = [ones(1,length(c1)) 2*ones(1,length(c2)) 3*ones(1,length(c3)) 4*ones(1,length(c4)) 5*ones(1,length(c5)) 6*ones(1,length(c6))];


for i = 1:N
    for j = 1:N
        if (A(i,j) == 1)
            if (l(i) == l(j))    
                E(i,j) = 1;
            else
                T(i,j) = 1;
            end
        end
    end
end

% isequal(T+E,A)

%%

L = E - diag(sum(E,2)');


figure(1)

imagesc(E)


save("E.mat",'E')
save("T.mat",'T')

%% RK4 method for solving HR equations 

n_tot = 277; % number of neurons
t_tot = 700; % number of time steps 
h = 0.05;

t = 0:h:t_tot;

a = 1;
b = 3;
c = 1;
d = 5;
s = 4;
p_o = -1.6;
I_ext = 3.25;
r = 0.005;
g_ch = 0.015;
g_el = 0.5;
V_syn = 2;
lambda = 10;
theta_syn = -0.25;
p = zeros(length(t),N);
q = zeros(length(t),N);
n = zeros(length(t),N);

p(1,:) = -1.30784489 + 0.5*rand(1,N);
q(1,:) = -7.3218313 + 0.5*rand(1,N);
n(1,:) = 3.35299859 + 0.5*rand(1,N);


dp_dt = @(t,p,q,n,L,T) q - a * p^3 + b * p^2 - n + I_ext + g_el * L - g_ch * (p - V_syn) * T;
dq_dt = @(t,p,q,n) c - d*p^2 - q;
dn_dt = @(t,p,q,n) r*(s*(p - p_o) - n);

%theta_syn = -2.5;

p_temp1 = p(1,:);
q_temp1 = q(1,:);
n_temp1 = n(1,:);

p_temp2 = p_temp1;
q_temp2 = q_temp1;
n_temp2 = n_temp1;
%%

for i = 2:length(t)

parfor k = 1:n_tot
    
    sig = 1./(1 + exp(-1.* lambda .* (p_temp2 - theta_syn)));
    p_k1 = h*dp_dt(t,p_temp2(k),q_temp2(k),n_temp2(k),L(k,:)*p_temp2',T(k,:)*sig');
    q_k1 = h*dq_dt(t,p_temp2(k),q_temp2(k),n_temp2(k));
    n_k1 = h*dn_dt(t,p_temp2(k),q_temp2(k),n_temp2(k));

    p_k2 = h*(dp_dt(t + h/2, p_k1/2 + p_temp2(k), q_temp2(k) + q_k1/2, n_temp2(k) + n_k1/2, L(k,:)*p_temp2',T(k,:)*sig'));
    q_k2 = h*(dq_dt(t + h/2, p_k1/2 + p_temp2(k), q_temp2(k) + q_k1/2, n_temp2(k) + n_k1/2));
    n_k2 = h*dn_dt(t + h/2,p_k1/2 + p_temp2(k), q_temp2(k) + q_k1/2, n_temp2(k) + n_k1/2);

    p_k3 = h*(dp_dt(t + h/2, p_k2/2 + p_temp2(k), q_temp2(k) + q_k2/2, n_temp2(k) + n_k2/2, L(k,:)*p_temp2',T(k,:)*sig'));
    q_k3 = h*(dq_dt(t + h/2, p_k2/2 + p_temp2(k), q_temp2(k) + q_k2/2, n_temp2(k) + n_k2/2));
    n_k3 = h*dn_dt(t + h/2,p_k2/2 + p_temp2(k), q_temp2(k) + q_k2/2, n_temp2(k) + n_k2/2);


    p_k4 = h*(dp_dt(t, p_k3 + p_temp2(k), q_temp2(k) + q_k3, n_temp2(k) + n_k3, L(k,:)*p_temp2',T(k,:)*sig'));
    q_k4 = h*(dq_dt(t, p_k3 + p_temp2(k), q_temp2(k) + q_k3, n_temp2(k) + n_k3));
    n_k4 = h*dn_dt(t,p_k3 + p_temp2(k), q_temp2(k) + q_k3, n_temp2(k) + n_k3);


    p(i,k) = p_temp2(k) + 1/6 * (p_k1 + 2*p_k2 + 2*p_k3 + p_k4);
    q(i,k) = q_temp2(k) + 1/6 * (q_k1 + 2*q_k2 + 2*q_k3 + q_k4);
    n(i,k) = n_temp2(k) + 1/6 * (n_k1 + 2*n_k2 + 2*n_k3 + n_k4);

end

p_temp2 = p(i,:);
q_temp2 = q(i,:);
n_temp2 = n(i,:);


end

%%


phi = zeros(length(t),N);

for i=1:length(t)

    for j = 1:N
    sig = 1./(1 + exp(-1.* lambda .* (p(i,:) - theta_syn)));
    phi(i,j) = (dq_dt(t,p(i,j),q(i,j),n(i,j)) * p(i,j) - q(i,j) * dp_dt(t,p(i,j),q(i,j),n(i,j),L(j,:)*p(i,:)',T(j,:)*sig'))./(q(i,j)^2 + p(i,j)^2);

    end

end


%%










