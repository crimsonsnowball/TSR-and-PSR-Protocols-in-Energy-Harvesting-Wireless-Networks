% system model with P,X as primary network
% Secondary network with direct link b/w S and D

clear
clc


N = 1e4;
mu1 = 1;
mu2 = 1;
mu3 = 1;
Nd = 4;
Ns = 2;
Pt_p = 5;
ita = 0.8;
Ith = 2;
No = 1;
R = 1;
T = 1;
i=1;
sigma_n_sq = 1;
% alpha = 0.3;
alpha = 0:0.1:1;
P_out= zeros(1,length(alpha));

for alpha = 0:0.1:1
    
    hpsj_jx1 = exprnd(mu1, Ns, 1, N);  
    Pt_Sj_max = ita *Pt_p* (alpha/(1-alpha)) .* sum (hpsj_jx1,1);  
    Pt_Sj_max_1xN = reshape(Pt_Sj_max,[],1);

    hsjx = exprnd(mu3, 1, Ns, N);
    hsjx_Nxj = reshape(hsjx,[],Ns);
    X2_Nxj = (Ith./hsjx_Nxj);
    
    X_Nxj = min(Pt_Sj_max_1xN,X2_Nxj);
    
    hsjdk = exprnd(mu2, Nd, Ns, N);
    Y_1xjxN = sum(hsjdk,1);
    Y_Nxj = reshape(Y_1xjxN,[],Ns); % Re think this 
     
    psi_Nxj = Y_Nxj.*X_Nxj;
    psi_Nx1 = max(psi_Nxj,[],2);
     
    Ro = (No + Pt_p .* sigma_n_sq)*(2.^(R./(T.*(1-alpha)))-1);

    P_out(i) = sum(psi_Nx1' < Ro)/N;
    
    i=i+1;

 end

[P_out_min, index] = min(P_out);
alpha = 0:0.1:1;
alpha_min = alpha(index);
plot(alpha,P_out,alpha_min,P_out_min,'o')
hold on