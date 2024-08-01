clear
clc
L=1; % Number of Relays

N = 1e5;    % Samples   
Ns = 5;


NRi = zeros(L,1)'+ 4 ;  % Assuming all Realy have same number of antennas
Nd = 4;

mu_1i = [5 1]; 
mu_2i = [1 1];
mu_3i = [1 4];
mu_4i = [2 1];

I_sx = 2;
I_RiX = 2;
lambda2 = 1;
lambda3 = 1;
mu2=1/lambda2;
mu3=1/lambda2;
Pp = 1;
ita = 1;

R = 1;
T = 1;

% No = 1;
% sigma1sq = 1;
% C1 = sigma1sq*Pp + No;
% sigma2sq = 1;
% C2 = sigma2sq*Pp + No;

C1=1;C2=1;

% alpha = 0.2;
o=1;
for alpha = 0:0.045:1
                                     % Channel Initialization
    
                        % Ui
                        % i=1 for all term with suffix 'i'
    
    % Energy Harvesting Phase of Source
    
    hpsj_Nxj = exprnd(mu2, N, Ns);    % [ P --> Sj ] Nxj 
    Pt_Sj_max_Nx1 = (2*ita .*Pp.* (alpha./(1-alpha))).* sum (hpsj_Nxj,2);
    
    % End of Energy Harvesting Phase of Source
    
    hsjx_Nxj = exprnd(mu3, N, Ns); % [Sj --> X] Nxj
    X2_Nxj_U = (I_sx./hsjx_Nxj); % X2 = Ith/hsjx
    X_Nxj_U = min(Pt_Sj_max_Nx1,X2_Nxj_U); % X = min{Pt sj max , Ith/hsjx}
    
    hsjrik_kxjxN = exprnd(mu_1i(1),NRi(1),Ns,N); % [Sj --> R1k] NR1xNs
    Yi_1xjxN_U = sum(hsjrik_kxjxN,1); % Sum k = 1 --> k = NR1
    
    if Ns>1                                 % Logic make sure concatenate is correct for Ns=1 (it has to do with how squeeze works) 
        Yi_Nxj_U = squeeze(Yi_1xjxN_U)';   % Logic to reduce one dimention
    else
         Yi_Nxj_U = squeeze(Yi_1xjxN_U) ;   % Logic to reduce one dimention
    end
    
    Zi_U = (Yi_Nxj_U.*X_Nxj_U)./C1;
    [Ui_Nx1 , index_j_Nx1 ] = max(Zi_U,[],2) ;
    
    %                 % Vi
    
    % Energy Harvesting Phase of Relay
    
    hprik_Nxk = exprnd(mu_3i(1),N,NRi(1));
    Pt_Rik_max_Nx1 = (2*ita*Pp*(alpha./(1-alpha)) ).* (sum (hprik_Nxk,2));
    
    % End of Energy Harvesting Phase of Relay
    
    hrikx_Nxk = exprnd(mu_4i(1),N,NRi(1));  
    X2_Nxk_V = I_RiX./hrikx_Nxk;
    X_Nxk_V = min(Pt_Rik_max_Nx1,X2_Nxk_V);
    
    hrikdn_nxkxN = exprnd(mu_2i(1),Nd,NRi(1),N); % [R1k --> Dn] Nd x NR1
    Yi_1xkxN_V = sum(hrikdn_nxkxN,1);% Sum n = 1 --> n = Nd
    
    if NRi(1)>1                            
        Yi_Nxj_V = squeeze(Yi_1xkxN_V)';
    else
        Yi_Nxj_V = squeeze(Yi_1xkxN_V);
    end
    
    Zi_V = (Yi_Nxj_V.*X_Nxk_V)./C2;
    [Vi_Nx1 , index_k_Nx1 ] = max(Zi_V,[],2) ;
    
    
    
    gamma_e2e_Nx1 = min(Ui_Nx1,Vi_Nx1);
    Ro = (2.^((2*R)./(T.*(1-alpha)))-1);
    
    P_out(o) = sum(gamma_e2e_Nx1'<Ro)/N;
    o=o+1;
end

alpha = 0:0.045:1;
plot(alpha,P_out)
xlabel('alpha')
ylabel('P_out')
hold on