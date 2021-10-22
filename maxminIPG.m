load('Z1.mat'); fs = 2000; type = 1; ECG = data(:,3); Z = -data(:,1); fs=2000; ks=0;
%% Filtratge de les senyals
%% Low Filter
kmax = 0.3*fs;
Wn = 5*2/fs;
N = 5;                                                                 
[a,b] = butter(N,Wn,'low');                                            
Z_l = filtfilt(a,b,Z); 
%% Calculem les derivades
h=0.001; x=[1:length(ECG)]./fs;
ECG_D = diff(ECG)/h;
Z_D = diff(Z_l)/h;
%% Cridem a Pan-Tompkins
[IBP,MeanRR] = PTfun(ECG,fs,ks);
%% Caracteristiques per als mínims:
if type == 1
    kmin =(MeanRR/4)*fs
end
if type == 2
    IBP = IBP + floor(MeanRR*fs/8);
    kmin = (MeanRR/4)*fs
end
%% Trobem màxims i mínims de la IPG

%% ***********************SENYAL FILTRADA**********************************
IPGmax=[]; IPGmax_index=[]; IPGmin=[]; IPGmin_index=[];

for ii = 1:length(IBP)
    ECGpeak=IBP(ii);
    
    bottommin = max(1,ECGpeak+kmin/2); 
    bottommax = max(1,ECGpeak); 
    topmin = min(ECGpeak+kmin, length(ECG)-1);
    topmax = min(ECGpeak+kmax, length(ECG)-1);

    peakIPG = max(Z_l(bottommax:topmax)); peakIPG_index = find(Z_l(bottommax:topmax)==peakIPG); peakIPG_index = peakIPG_index + bottommax;
    lowIPG = min(Z_l(bottommin:topmin)); lowIPG_index = find(Z_l(bottommin:topmin)==lowIPG); lowIPG_index = lowIPG_index + bottommin;
    IPGmax=[IPGmax,peakIPG]; IPGmin=[IPGmin,lowIPG];
    IPGmax_index=[IPGmax_index,peakIPG_index]; IPGmin_index=[IPGmin_index,lowIPG_index];
end
figure(2)
plot([1:length(Z)]./fs,Z,'b');
hold on
plot([1:length(Z_l)]./fs,Z_l,'--k');
hold on
plot(IPGmax_index./fs,IPGmax,'or'); 
hold on
plot(IPGmin_index./fs,IPGmin,'og'); 
axis([10 40 -0.3 0.3])
grid on
hold off

% % ********************SENYAL SENSE FILTRAR********************************
IPGmax=[]; IPGmax_index=[]; IPGmin=[]; IPGmin_index=[];

for ii = 1:length(IBP)
    ECGpeak=IBP(ii);
    
    bottommin = max(1,ECGpeak+kmin/2); 
    bottommax = max(1,ECGpeak); 
    topmin = min(ECGpeak+kmin, length(ECG)-1);
    topmax = min(ECGpeak+kmax, length(ECG)-1);

    peakIPG = max(Z(bottommax:topmax)); peakIPG_index = find(Z(bottommax:topmax)==peakIPG); peakIPG_index = peakIPG_index + bottommax;
    lowIPG = min(Z(bottommin:topmin)); lowIPG_index = find(Z(bottommin:topmin)==lowIPG); lowIPG_index = lowIPG_index + bottommin;
    IPGmax=[IPGmax,peakIPG]; IPGmin=[IPGmin,lowIPG]; %peakIPG_index(ceil(length(peakIPG_index)/2))
    IPGmax_index=[IPGmax_index,peakIPG_index(ceil(length(peakIPG_index)/2))];
    IPGmin_index=[IPGmin_index,lowIPG_index(ceil(length(lowIPG_index)/2))];
end
figure(3)
plot([1:length(Z)]./fs,Z,'b');
hold on
plot([1:length(Z_l)]./fs,Z_l,'--k');
hold on
plot(IPGmax_index./fs,IPGmax,'or'); 
hold on
plot(IPGmin_index./fs,IPGmin,'og'); 
axis([10 40 -0.3 0.3])
grid on
hold off

% % ******************SENYAL SENSE FILTRAR AMB PICS ECG******************************
IPGmax=[]; IPGmax_index=[]; IPGmin=[]; IPGmin_index=[];
for ii = 1:length(IBP)
    ECGpeak=IBP(ii);
    
    bottom = max(1,ECGpeak); 
    topmin = min(ECGpeak+kmin, length(ECG)-1);
    topmax = min(ECGpeak+kmax, length(ECG)-1);

    peakIPG = max(Z(bottom:topmax)); peakIPG_index = find(Z(bottom:topmax)==peakIPG); peakIPG_index = peakIPG_index + bottom;
    lowIPG = min(Z(bottom:topmin)); lowIPG_index = find(Z(bottom:topmin)==lowIPG); lowIPG_index = lowIPG_index + bottom;
    IPGmax=[IPGmax,peakIPG]; IPGmin=[IPGmin,lowIPG]; %peakIPG_index(ceil(length(peakIPG_index)/2))
    IPGmax_index=[IPGmax_index,peakIPG_index(ceil(length(peakIPG_index)/2))];
    IPGmin_index=[IPGmin_index,lowIPG_index(ceil(length(lowIPG_index)/2))];
end

[IBP,MeanRR] = PTfun(ECG,fs,ks);
IBP8 = IBP + floor(MeanRR*fs/8);
IBP_new = IBP + floor(MeanRR*fs/4);


figure(4)
plot([1:length(Z)]./fs,Z,'b');
hold on
plot(IPGmax_index./fs,IPGmax,'or'); 
hold on
plot(IPGmin_index./fs,IPGmin,'og'); 
hold on
line(repmat(IBP/fs,[2 1]),...
repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(IBP8)),...
'LineWidth',1.0,'LineStyle','-.','Color','k');
hold on
line(repmat((IBP+kmin)/fs,[2 1]),...
repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(IBP_new)),...
'LineWidth',1.0,'LineStyle','-.','Color','y');
axis([10 40 -0.3 0.3])
grid on
hold off