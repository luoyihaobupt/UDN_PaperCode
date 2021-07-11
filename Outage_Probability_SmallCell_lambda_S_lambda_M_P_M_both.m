clear all;
close all;
clc;

	
alpha = 4;
	
	

	
lambda_M_array = [1 3]*(10^(-6));     	
lambda_M_num = size(lambda_M_array,2);

	
P_Su_dB = 22;           	
P_Sd_dB = 30;           	

P_M_array = [35 40];                	
P_M_num = size(P_M_array,2);

	
P_M_array = 10.^(P_M_array/10);
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

R_M0 = 100;
R_Su0 = 60;
R_Sd0 = 70;

eta_M = pi*(xi_M^(2/alpha))*(R_M0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Su = pi*(xi_Su^(2/alpha))*(R_Su0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);
eta_Sd = pi*(xi_Sd^(2/alpha))*(R_Sd0^2)*gamma(1+2/alpha)*gamma(1-2/alpha);

line_style = char('-k*','-rv','-ms','-mo','-r^','-k');

	
lambda_S_ini = 0;
lambda_S_end = 2*10^(-4);
lambda_S_num = 20;
delta_lambda_S = (lambda_S_end-lambda_S_ini)/lambda_S_num;
lambda_S = lambda_S_ini:1*delta_lambda_S:1*lambda_S_end;		

	
P_M = P_M_array(1,1);
for lambda_M_ite = 1:lambda_M_num;
    lambda_M = lambda_M_array(1,lambda_M_ite);
		
	Outage_SmallCell_uplink = 1 - exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
	lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
	plot(lambda_S(1,:), Outage_SmallCell_uplink(1,:), line_style(lambda_M_ite,:))
	hold on;
end
P_M = P_M_array(1,2);
	
Outage_SmallCell_uplink = 1 - exp(-lambda_M*eta_Su*((P_M/P_Su)^(2/alpha))-...
    lambda_S*eta_Su*(Pr_u+Pr_d*((P_Sd/P_Su)^(2/alpha))));
plot(lambda_S(1,:), Outage_SmallCell_uplink(1,:), line_style(lambda_M_num+1,:))
hold on;

P_M = P_M_array(1,1);
	
for lambda_M_ite = 1:lambda_M_num;
    lambda_M = lambda_M_array(1,lambda_M_ite);
	Outage_SmallCell_downlink = 1 - exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
        lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
    plot(lambda_S(1,:), Outage_SmallCell_downlink(1,:), line_style(3+lambda_M_ite,:))
end

P_M = P_M_array(1,2);
	
Outage_SmallCell_downlink = 1 - exp(-lambda_M*eta_Sd*((P_M/P_Sd)^(2/alpha))-...
    lambda_S*eta_Sd*(Pr_u*((P_Su/P_Sd)^(2/alpha))+Pr_d));
plot(lambda_S(1,:), Outage_SmallCell_downlink(1,:), line_style(6,:))
hold on;

	
	
xlabel('Small cell BS density (Small cell BS/m^2)');
ylabel('Outage probability');
legend('\lambda_{M}=1x10^{-6} Macro Cell UE/m^2, P_{M}=35dBm',...
    '\lambda_{M}=3x10^{-6} Macro Cell UE/m^2, P_{M}=35dBm',...
    '\lambda_{M}=3x10^{-6} Macro Cell UE/m^2, P_{M}=40dBm',...
    '\lambda_{M}=1x10^{-6} Macro Cell UE/m^2, P_{M}=35dBm',...
    '\lambda_{M}=3x10^{-6} Macro Cell UE/m^2, P_{M}=35dBm',...
    '\lambda_{M}=3x10^{-6} Macro Cell UE/m^2, P_{M}=40dBm',...
        'Location','SouthEast');
axis([lambda_S_ini,lambda_S_end,0,1])
	
