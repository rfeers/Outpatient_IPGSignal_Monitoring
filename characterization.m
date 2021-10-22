function[IBP,IPGmax, IPGmin, IPGmax_index,IPGmin_index]=characterization3_0(ECG,Z,type,fs)
fs=2000; ks=0;
%% Filtratge de les senyals
%% Low Filter
kmax = 0.4*fs;
Wn = 5*2/fs;
N = 5;                                                                 
[a,b] = butter(N,Wn,'low');                                            
Z_l = filtfilt(a,b,Z); 
%% Calculem les derivades
h=0.001; x=[1:length(ECG)]./fs;
ECG_D = diff(ECG)/h;
Z_D = diff(Z_l)/h;
%% Cridem a Pan-Tompkins
[IBP,MeanRR] = PTfun1_0(ECG,fs,ks);
%% Caracteristiques per als mínims:
if type == 1
    kmin =(MeanRR/4)*fs; kmin = floor(kmin);
end
if type == 2
    IBP = IBP + floor(MeanRR*fs/8);
    kmin = (MeanRR/4)*fs; kmin = floor(kmin);
end
%% Trobem màxims i mínims de la IPG

% % ********************SENYAL SENSE FILTRAR********************************
IPGmax=[]; IPGmax_index=[]; IPGmin=[]; IPGmin_index=[];

for ii = 1:length(IBP)
    ECGpeak=IBP(ii);
    bottommin = max(1,ECGpeak);  
    if type == 2
    bottommin = max(1,ECGpeak+kmin/2); 
    end
    bottommax = max(1,ECGpeak); 
    topmin = min(ECGpeak+kmin, length(ECG)-1);
    topmax = min(ECGpeak+kmax, length(ECG)-1);

    peakIPG = max(Z(bottommax:topmax)); peakIPG_index = find(Z(bottommax:topmax)==peakIPG); peakIPG_index = peakIPG_index + bottommax;
    lowIPG = min(Z(bottommin:topmin)); lowIPG_index = find(Z(bottommin:topmin)==lowIPG); lowIPG_index = lowIPG_index + bottommin;
    IPGmax=[IPGmax,peakIPG]; IPGmin=[IPGmin,lowIPG]; %peakIPG_index(ceil(length(peakIPG_index)/2))
    IPGmax_index=[IPGmax_index,peakIPG_index(ceil(length(peakIPG_index)/2))];
    IPGmin_index=[IPGmin_index,lowIPG_index(ceil(length(lowIPG_index)/2))];
end
end



