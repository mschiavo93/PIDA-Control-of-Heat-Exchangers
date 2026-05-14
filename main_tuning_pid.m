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
mu = Dy/Du;

erroremin=10000;
T1_est_opt=0;
T2_est_opt=0;
delay_opt=0;
delaymax=100;

y = y-y(1);
opt = stepDataOptions('StepAmplitude',step_Fc);

for delay=0:0.01:delaymax
    delay
    A1= (y(end)*Tsim-sum(y)*ST)/y(end)-delay;
    T1piuT2piudelay = A1+delay;
    g2 = y(1:round(T1piuT2piudelay/ST)); 
    A2suA1 = (1/(A1))*sum(g2)*ST; 
    p1 =     0.01043 ; 
    q1 =     -0.2594 ;
    r_est = (p1) / (A2suA1/y(end) + q1);       
    T1_est = A1/(1+r_est);
    T2_est = r_est*T1_est;
    G_est=tf([0 mu],[T1_est*T2_est T1_est+T2_est 1]);
    G_est.outputdelay=delay;
    g_est = step(G_est,0:ST:Tsim,opt);
    errore=sum((y-g_est).^2);
    if (errore<erroremin)
        erroremin=errore;
        delay_opt=delay;
        T1_est_opt=T1_est;
        T2_est_opt=T2_est;
    end
end
Gs2=tf([0 mu],[T1_est_opt*T2_est_opt T1_est_opt+T2_est_opt 1]);
Gs2.outputdelay=delay_opt;
[y_SOPDT,t_SOPDT] = step(Gs2,Tsim,opt);
figure
plot(t,y)
hold on
plot(t_SOPDT,y_SOPDT)
denominatore = Gs2.den{1};
a1 = denominatore(2);
a2 = denominatore(1);
TiPID = a1;
TdPID = a2/a1;

Ms=1.4;

for KpPID=-0.001:-0.0001:-0.5
    PID = KpPID*(1 + 1/(TiPID*s) + TdPID*s/(TdPID/10*s+1));
    Ls_PID = minreal(PID*Gs2, 0.001);
    [MAG,~] =bode(feedback(1,Ls_PID)); 
    MS_PID = max(MAG);
    if MS_PID>Ms
        break
    end
end
MS_PID

[GmPID,PmPID,WcgPID,WcpPID]=margin(Ls_PID);

FsPID=1/(0.05/WcpPID*s+1)^2;   

% save('Tuning_PID_14_2','KpPID','TiPID','TdPID','WcpPID')

