%% Integrate all the pulses Z1
load('Z1.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 1;
[meanmm1,suma1,I1] = promig(ECG, Z, type, fs);
%% Integrate all the pulses Z2
load('Z2.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 3;
[meanmm2,suma2,I2] = promig(ECG, Z, type, fs);
%% Integrate all the pulses Z3
load('Z3.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 2;
[meanmm3,suma3,I3] = promig(ECG, Z, type, fs);
%% Integrate all the pulses Z4
load('Z4.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 2;
[meanmm4,suma4,I4] = promig(ECG, Z, type, fs);
%% Integrate all the pulses Z5
load('Z5.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 2;
[meanmm5,suma5,I5] = promig(ECG, Z, type, fs);
%% Integrate all the pulses PPG
load('Z1.mat'); PPG=data(:,4); ECG = data(:,3); fs = 2000; type = 2;
[meanmm6,suma6,IPPG] = promig(ECG, PPG, type, fs);
%% Integrate all the pulses Z7
load('armtoarm.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 2;
[meanmm7,suma7,I7] = promig(ECG, Z, type, fs);
%% Integrate all the pulses Z8
load('legtoleg.mat'); Z=-data(:,1); ECG = data(:,3); fs = 2000; type = 2;
[meanmm8,suma8,I8] = promig(ECG, Z, type, fs);

suma1 = suma1(1:meanmm1); suma2 = suma2(1:meanmm2); suma3 = suma3(1:meanmm3);
suma4 = suma4(1:meanmm4); suma5 = suma5(1:meanmm5); suma6 = suma6(1:meanmm6);
suma7 = suma7(1:meanmm7); suma8 = suma8(1:meanmm8);

suma1 = suma1/max(suma1); suma2 = suma2/max(suma2); suma3 = suma3/max(suma3);
suma4 = suma4/max(suma4); suma5 = suma5/max(suma5); suma6 = suma6/max(suma6);
suma7 = suma7/max(suma7); suma8 = suma8/max(suma8);

I1
I2
I3
I4
I5
IPPG
I7
I8
figure(1)
plot(suma1)
hold on
plot(suma2)
plot(suma6)
grid on
legend('ICG','Neck','PPG')

plot(suma3)
hold on
plot(suma4)
plot(suma5)
grid on
legend('Uparm','Midarm','Downarm')

figure(3)
plot(suma7)
hold on
plot(suma8)
legend('legtoleg','armtoarm')
