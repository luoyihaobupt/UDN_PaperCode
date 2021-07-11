clear all;
close all;
clc;

	
W = 1;         	
alpha = 4;
lambda_M = 1*10^(-6);     	

P_M_dB = 30;            	
P_Su_dB = 15;           	
P_Sd_dB = 20;           	
	
P_M = 10^(P_M_dB/10);
P_Su = 10^(P_Su_dB/10);
P_Sd = 10^(P_Sd_dB/10);

xi_M_dB = -5;       	
xi_Sd_dB = 0;       	
xi_Su_dB = 0;       	
	
xi_M = 10^(xi_M_dB/10);
xi_Sd = 10^(xi_Sd_dB/10);
xi_Su = 10^(xi_Su_dB/10);

R_M0 = 100;
R_Su0 = 50;
R_Sd0 = 60;

eta_M = pi*(xi_M^(2/alpha))*(R_M0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Su = pi*(xi_Su^(2/alpha))*(R_Su0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Sd = pi*(xi_Sd^(2/alpha))*(R_Sd0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);

	
lambda_S_ini = 0;
lambda_S_end = 1.5*10^(-4);
lambda_S_num = 39;
delta_lambda_S = (lambda_S_end-lambda_S_ini)/lambda_S_num;
lambda_S_array = lambda_S_ini:1*delta_lambda_S:1*lambda_S_end;		

	
Pr_d_ini = 0;
Pr_d_end = 1;
Pr_d_num = lambda_S_num;
delta_Pr_d = (Pr_d_end-Pr_d_ini)/Pr_d_num;
Pr_d_array = Pr_d_ini:delta_Pr_d:Pr_d_end;
Pr_u_array = 1- Pr_d_array;

ASR_SC = zeros(lambda_S_num, Pr_d_num);
ASR_uplink_SC = zeros(lambda_S_num, Pr_d_num);
ASR_downlink_SC = zeros(lambda_S_num, Pr_d_num);
	
for lambda_S_ite = 1:1:(lambda_S_num+1);
	lambda_S = lambda_S_array(1,lambda_S_ite);
    for Pr_d_ite = 1:1:(Pr_d_num+1);
        Pr_d = Pr_d_array(1,Pr_d_ite);
        Pr_u = Pr_u_array(1,Pr_d_ite);    
        	
        ASR_uplink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_u*lambda_S*log2(1+xi_Su)*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
            lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
        	
        ASR_downlink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_d*lambda_S*log2(1+xi_Sd)*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
            lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
        	
        ASR_SC(lambda_S_ite,Pr_d_ite) = ASR_uplink_SC(lambda_S_ite,Pr_d_ite) + ASR_downlink_SC(lambda_S_ite,Pr_d_ite);
    end
end

ASR_SC_min = min(min(ASR_SC));
ASR_SC_max = max(max(ASR_SC));

figure;
	
surf(Pr_d_array(1,:),lambda_S_array(1,:),ASR_SC);
axis([Pr_d_ini,Pr_d_end,lambda_S_ini,lambda_S_end,ASR_SC_min,ASR_SC_max]);
xlabel({'Pr_d','(1-Pr_u)'});
ylabel({'Small cell BS density','(Small cell BS/m^{2})'});
zlabel('ASR of small cells (bit/S/Hz/m^{2})')
hold on;


	
	
W = 1;         	
alpha = 4;
lambda_M = 1*10^(-6);     	

P_M_dB = 40;            	
P_Su_dB = 5;           	
P_Sd_dB = 10;           	
	
P_M = 10^(P_M_dB/10);
P_Su = 10^(P_Su_dB/10);
P_Sd = 10^(P_Sd_dB/10);

xi_M_dB = -5;       	
xi_Sd_dB = 0;       	
xi_Su_dB = 0;       	
	
xi_M = 10^(xi_M_dB/10);
xi_Sd = 10^(xi_Sd_dB/10);
xi_Su = 10^(xi_Su_dB/10);

R_M0 = 100;
R_Su0 = 50;
R_Sd0 = 60;

eta_M = pi*(xi_M^(2/alpha))*(R_M0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Su = pi*(xi_Su^(2/alpha))*(R_Su0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Sd = pi*(xi_Sd^(2/alpha))*(R_Sd0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);

	
lambda_S_ini = 0;
lambda_S_end = 1.5*10^(-4);
lambda_S_num = 39;
delta_lambda_S = (lambda_S_end-lambda_S_ini)/lambda_S_num;
lambda_S_array = lambda_S_ini:1*delta_lambda_S:1*lambda_S_end;		

	
Pr_d_ini = 0;
Pr_d_end = 1;
Pr_d_num = lambda_S_num;
delta_Pr_d = (Pr_d_end-Pr_d_ini)/Pr_d_num;
Pr_d_array = Pr_d_ini:delta_Pr_d:Pr_d_end;
Pr_u_array = 1- Pr_d_array;

ASR_SC = zeros(lambda_S_num, Pr_d_num);
ASR_uplink_SC = zeros(lambda_S_num, Pr_d_num);
ASR_downlink_SC = zeros(lambda_S_num, Pr_d_num);
	
for lambda_S_ite = 1:1:(lambda_S_num+1);
	lambda_S = lambda_S_array(1,lambda_S_ite);
    for Pr_d_ite = 1:1:(Pr_d_num+1);
        Pr_d = Pr_d_array(1,Pr_d_ite);
        Pr_u = Pr_u_array(1,Pr_d_ite);    
        	
        ASR_uplink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_u*lambda_S*log2(1+xi_Su)*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
            lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
        	
        ASR_downlink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_d*lambda_S*log2(1+xi_Sd)*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
            lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
        	
        ASR_SC(lambda_S_ite,Pr_d_ite) = ASR_uplink_SC(lambda_S_ite,Pr_d_ite) + ASR_downlink_SC(lambda_S_ite,Pr_d_ite);
    end
end

ASR_SC_min = min(min(ASR_SC));
ASR_SC_max = max(max(ASR_SC));

surf(Pr_d_array(1,:),lambda_S_array(1,:),ASR_SC);

	
	
W = 1;         	
alpha = 4;
lambda_M = 1*10^(-6);     	

P_M_dB = 30;            	
P_Su_dB = 15;           	
P_Sd_dB = 20;           	
	
P_M = 10^(P_M_dB/10);
P_Su = 10^(P_Su_dB/10);
P_Sd = 10^(P_Sd_dB/10);

xi_M_dB = -5;       	
xi_Sd_dB = 0;       	
xi_Su_dB = 0;       	
	
xi_M = 10^(xi_M_dB/10);
xi_Sd = 10^(xi_Sd_dB/10);
xi_Su = 10^(xi_Su_dB/10);

R_M0 = 100;
R_Su0 = 70;
R_Sd0 = 100;

eta_M = pi*(xi_M^(2/alpha))*(R_M0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Su = pi*(xi_Su^(2/alpha))*(R_Su0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Sd = pi*(xi_Sd^(2/alpha))*(R_Sd0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);

	
lambda_S_ini = 0;
lambda_S_end = 1.5*10^(-4);
lambda_S_num = 39;
delta_lambda_S = (lambda_S_end-lambda_S_ini)/lambda_S_num;
lambda_S_array = lambda_S_ini:1*delta_lambda_S:1*lambda_S_end;		

	
Pr_d_ini = 0;
Pr_d_end = 1;
Pr_d_num = lambda_S_num;
delta_Pr_d = (Pr_d_end-Pr_d_ini)/Pr_d_num;
Pr_d_array = Pr_d_ini:delta_Pr_d:Pr_d_end;
Pr_u_array = 1- Pr_d_array;

ASR_SC = zeros(lambda_S_num, Pr_d_num);
ASR_uplink_SC = zeros(lambda_S_num, Pr_d_num);
ASR_downlink_SC = zeros(lambda_S_num, Pr_d_num);
	
for lambda_S_ite = 1:1:(lambda_S_num+1);
	lambda_S = lambda_S_array(1,lambda_S_ite);
    for Pr_d_ite = 1:1:(Pr_d_num+1);
        Pr_d = Pr_d_array(1,Pr_d_ite);
        Pr_u = Pr_u_array(1,Pr_d_ite);    
        	
        ASR_uplink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_u*lambda_S*log2(1+xi_Su)*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
            lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
        	
        ASR_downlink_SC(lambda_S_ite,Pr_d_ite) = W*Pr_d*lambda_S*log2(1+xi_Sd)*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
            lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
        	
        ASR_SC(lambda_S_ite,Pr_d_ite) = ASR_uplink_SC(lambda_S_ite,Pr_d_ite) + ASR_downlink_SC(lambda_S_ite,Pr_d_ite);
    end
end

ASR_SC_min = min(min(ASR_SC));
ASR_SC_max = max(max(ASR_SC));

surf(Pr_d_array(1,:),lambda_S_array(1,:),ASR_SC);
