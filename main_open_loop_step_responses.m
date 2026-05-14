clear all
close all
clc

load initial_conditions.mat

fs = 12;

Tsim = 200;
ST = 0.01;
step_time = 0;

step_Tpo = 0;
step_Tci = 0;
step_Fp = 0;
step_Tpi = 0;

U = 35; % [Btu/hr °F ft^2]
A = 20; % [ft^2]
Cpp = 0.38; % [Btu/lbm]
Mp = 15; % [lbm]
Fp = 1500/7936.64; % [kg/s]
Fc = 2500/7936.64; % [kg/s]
Fc_min = 500/7936.64; % [kg/s]
Tpi = 366.48; % [K]
Cpc = 0.456; % [Btu/lbm]
Mc = 40; % [lbm]
Tci = 294.26; % [K]
N_stages = 20;

Tco_init1 = Tco_init1.Data(end);
Tco_init2 = Tco_init2.Data(end);
Tco_init3 = Tco_init3.Data(end);
Tco_init4 = Tco_init4.Data(end);
Tco_init5 = Tco_init5.Data(end);
Tco_init6 = Tco_init6.Data(end);
Tco_init7 = Tco_init7.Data(end);
Tco_init8 = Tco_init8.Data(end);
Tco_init9 = Tco_init9.Data(end);
Tco_init10 = Tco_init10.Data(end);
Tco_init11 = Tco_init11.Data(end);
Tco_init12 = Tco_init12.Data(end);
Tco_init13 = Tco_init13.Data(end);
Tco_init14 = Tco_init14.Data(end);
Tco_init15 = Tco_init15.Data(end);
Tco_init16 = Tco_init16.Data(end);
Tco_init17 = Tco_init17.Data(end);
Tco_init18 = Tco_init18.Data(end);
Tco_init19 = Tco_init19.Data(end);
Tco_init20 = Tco_init20.Data(end);

Tpo_init1 = Tpo_init1.Data(end);
Tpo_init2 = Tpo_init2.Data(end);
Tpo_init3 = Tpo_init3.Data(end);
Tpo_init4 = Tpo_init4.Data(end);
Tpo_init5 = Tpo_init5.Data(end);
Tpo_init6 = Tpo_init6.Data(end);
Tpo_init7 = Tpo_init7.Data(end);
Tpo_init8 = Tpo_init8.Data(end);
Tpo_init9 = Tpo_init9.Data(end);
Tpo_init10 = Tpo_init10.Data(end);
Tpo_init11 = Tpo_init11.Data(end);
Tpo_init12 = Tpo_init12.Data(end);
Tpo_init13 = Tpo_init13.Data(end);
Tpo_init14 = Tpo_init14.Data(end);
Tpo_init15 = Tpo_init15.Data(end);
Tpo_init16 = Tpo_init16.Data(end);
Tpo_init17 = Tpo_init17.Data(end);
Tpo_init18 = Tpo_init18.Data(end);
Tpo_init19 = Tpo_init19.Data(end);
Tpo_init20 = Tpo_init20.Data(end);

Tpo_init = (Tpo_init1 - 32)*(5/9) + 273.15;

a1 = N_stages/(3600*Mc);
a2 = (U*A)/(3600*Mc*Cpc);
a3 = N_stages/(3600*Mp);
a4 = (U*A)/(3600*Mp*Cpp);


%%
step_Fc = 0.025;
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = -0.025; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = 0.05;
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,':k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = -0.05; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,':k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = 0.1; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'-.k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = -0.1; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'-.k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = 0.15; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'--k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 

%%
step_Fc = -0.15; 
[t_PIDA,~,y_PIDA,Tco_PIDA,u_PIDA] = sim('hx_model_alsop.slx');
plot(t_PIDA,y_PIDA,'--k','linewidth',1.5)
hold on
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
set(gca, 'FontSize', fs); 