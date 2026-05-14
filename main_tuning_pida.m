clear all
close all
clc

load initial_conditions.mat

s = tf('s');

Tsim = 800;
ST = 0.01;
step_time = 0;

U = 35; % [Btu/hr °F ft^2]
A = 20; % [ft^2]
Cpp = 0.38; % [Btu/lbm]
Mp = 15; % [lbm]
Fp = 1500/7936.64; % [kg/s]
Fc = 2500/7936.64; % [kg/s]
Tpi = 366.48; % [K]
Cpc = 0.456; % [Btu/lbm]
Mc = 40; % [lbm]
Tci = 294.26; % [K]
N = 20;

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

a1 = N/(3600*Mc);
a2 = (U*A)/(3600*Mc*Cpc);
a3 = N/(3600*Mp);
a4 = (U*A)/(3600*Mp*Cpp);


step_Fc = -0.05;

[t,~,y,Tco,u] = sim('hx_model_alsop_20stages.slx');

% Gain estimation
Dy = y(end)-y(1);
Du = step_Fc;
K = Dy/Du;

% Time instants to reach 5%, 35.3%, 85.3% of the final value
i5 = find(y > y(1)+(0.05*Dy),1);
t5 = t(i5);
i35 = find(y > y(1)+(0.353*Dy),1);
t35 = t(i35);
i85 = find(y > y(1)+(0.853*Dy),1);
t85 = t(i85);

figure
plot(t,y,'k','linewidth',1.5)
hold on 
plot(t5,y(i5),'*k')
plot([0 t5],[y(i5) y(i5)],'--k','linewidth',1.5)
plot([t5 t5],[y(1) y(i5)],'--k','linewidth',1.5)
plot(t35,y(i35),'*k')
plot([0 t35],[y(i35) y(i35)],'--k','linewidth',1.5)
plot([t35 t35],[y(1) y(i35)],'--k','linewidth',1.5)
plot(t85,y(i85),'*k')
plot([0 t85],[y(i85) y(i85)],'--k','linewidth',1.5)
plot([t85 t85],[y(1) y(i85)],'--k','linewidth',1.5)
grid on
ylabel('Tpo [K]')
xlabel('Time [s]')
ylim([y(1) y(end)+0.2])

% Transfer function parameters of the approximate model of the original 
% system
L = 1.3*t35-0.29*t85;
tau = 0.67*(t85-t35);
alfa = 0.598 + 0.4799*(t5/L) - (0.41/((t5/tau)^0.6));

% Approximate model of the original system, resulting in a third-order
% system
Gs3 = K/((1+tau*s)*(1+(1-alfa)/2*L*s)^2)*exp(-alfa*L*s);

y = y-y(1);
opt = stepDataOptions('StepAmplitude',step_Fc);
[y_TOPDT,t_TOPDT] = step(Gs3,Tsim,opt);
figure
plot(t,y)
hold on
plot(t_TOPDT,y_TOPDT)

% Tuning the PIDA controller without using filters, via pole-zero
% cancellation
denominatore = Gs3.den{1};
a1 = denominatore(3);
a2 = denominatore(2);
a3 = denominatore(1);
TiPIDA = a1;
TdPIDA = a2/a1;
TaPIDA = a3/a1;

Ms=1.4;

for KpPIDA=-0.001:-0.0001:-0.5
    PIDA = KpPIDA*(1 + 1/(TiPIDA*s) + TdPIDA*s/(TdPIDA/10*s+1) +...
        TaPIDA*s^2/(TaPIDA/10*s+1)^2);
    Ls_PIDA = minreal(PIDA*Gs3, 0.001);
    [MAG,~] =bode(feedback(1,Ls_PIDA)); 
    MS_PIDA = max(MAG);
    if MS_PIDA>Ms
        break
    end
end
MS_PIDA


[GmPIDA,PmPIDA,WcgPIDA,WcpPIDA]=margin(Ls_PIDA); % 51.8 with Haalman

FsPIDA=1/(0.05/WcpPIDA*s+1)^2;
PIDA=PIDA*FsPIDA;    

% save('Tuning_PIDA_14','KpPIDA','TiPIDA','TdPIDA','TaPIDA','WcpPIDA')