clear all;
close all;
clc;

	
W = 1;         	
alpha = 4;
lambda_M = 1*10^(-6);     	

theta_M = 0.1;          	
theta_S = 0.1;          	

P_M_dB = 36;            	
P_Su_dB = 14;           	
P_Sd_dB = 11;           	
	
P_M = 10^(P_M_dB/10);
P_Su = 10^(P_Su_dB/10);
P_Sd = 10^(P_Sd_dB/10);

Pr_u = 0.6;                 	
Pr_d = 1-Pr_u;              	

xi_M_dB = -5;        	
xi_Sd_dB = 0;       	
xi_Su_dB = 0;       	
	
xi_M = 10^(xi_M_dB/10);
xi_Sd = 10^(xi_Sd_dB/10);
xi_Su = 10^(xi_Su_dB/10);

R_M0 = 80;
R_Su0 = 60;
R_Sd0 = 70;

eta_M = pi*(xi_M^(2/alpha))*(R_M0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Su = pi*(xi_Su^(2/alpha))*(R_Su0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Sd = pi*(xi_Sd^(2/alpha))*(R_Sd0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);

	
lambda_S_ini = 0;
lambda_S_end = 1.5*10^(-4);
lambda_S_num = 200;
delta_lambda_S = (lambda_S_end-lambda_S_ini)/lambda_S_num;
lambda_S = lambda_S_ini:delta_lambda_S:lambda_S_end;		

	
ASR_uplink_SC = W*Pr_u*lambda_S.*log2(1+xi_Su).*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
    lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
	
ASR_downlink_SC = W*Pr_d*lambda_S.*log2(1+xi_Sd).*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
    lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
	
ASR_SC = ASR_uplink_SC + ASR_downlink_SC;

	
[ASR_SC_max, lambda_S_max_index] = max(ASR_SC);

	
A3 = eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha)));
B3 = eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d);
ASR_uplink_SC_max = W*Pr_u*(1/A3)*log2(1+xi_Su)*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
    (1/A3)*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
ASR_downlink_SC_max = W*Pr_d*(1/B3)*log2(1+xi_Sd).*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
    (1/B3)*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));

	
lambda_sup1 = (-lambda_M*(P_M^(2/alpha))-((P_M^(2/alpha)/eta_M)*log(1-theta_M)))/...
    (Pr_u*(P_Su^(2/alpha))+Pr_d*(P_Sd^(2/alpha)));
ASR_uplink_SC_sup1 = W*Pr_u*(lambda_sup1)*log2(1+xi_Su)*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
    (lambda_sup1)*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
ASR_downlink_SC_sup1 = W*Pr_d*(lambda_sup1)*log2(1+xi_Sd)*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
    (lambda_sup1)*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
ASR_SC_sup1 = ASR_uplink_SC_sup1+ASR_downlink_SC_sup1;

	
delta_lambda_S_feasible_num = 40;
delta_lambda_S_feasible = lambda_sup1/delta_lambda_S_feasible_num;
lambda_S1 = 0:delta_lambda_S_feasible:lambda_sup1;
lambda_S2 = lambda_sup1:-delta_lambda_S_feasible:0;
lambda_S_feasible = [lambda_S1, lambda_S2];
y1 = zeros(1,size(lambda_S1,2));
ASR_uplink_SC_feasible = W*Pr_u*(lambda_S2).*log2(1+xi_Su).*exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
    (lambda_S2)*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
ASR_downlink_SC_feasible = W*Pr_d*(lambda_S2).*log2(1+xi_Sd).*exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
    (lambda_S2)*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
ASR_SC_feasible = ASR_uplink_SC_feasible + ASR_downlink_SC_feasible;
y = [y1, ASR_SC_feasible];
pic1 = fill(lambda_S_feasible,y,[0.878 0.878 0.878]);
hold on;
	
pic2 = plot(lambda_S(1,:),ASR_uplink_SC(1,:),'-m');
hold on;
pic3 = plot(lambda_S(1,:), ASR_downlink_SC(1,:),'-r');
pic4 = plot(lambda_S(1,:),ASR_SC(1,:),'-k');
	
pic5 = plot((1/A3),ASR_uplink_SC_max,'mv');
pic6 = plot((1/B3),ASR_downlink_SC_max,'r^');
pic7 = plot(lambda_S(1,lambda_S_max_index),ASR_SC_max,'ko');
pic8 = plot(lambda_sup1,ASR_SC_sup1,'k*');

xlabel('Small cell BS density (Small cell BS/m^{2})');
ylabel('ASR of Small cells (bit/s/Hz/{m^2})');
legend([pic2, pic3, pic4, pic5, pic6, pic7],'uplink ASR',...
    'downlink ASR','total ASR','maximum uplink ASR',...
    'maximum downlink ASR','maximum total ASR');
axis([lambda_S_ini,lambda_S_end,0,max(ASR_SC)*1.1]);


