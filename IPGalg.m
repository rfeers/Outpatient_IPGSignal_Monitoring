kplot=1; 
%% Filtratge de les senyals
%% Low Filter
% load('S3_2.mat');fs=1000; ks=0;
load('ECG-IPGH-Miren.txt');ECG = ECG_IPGH_Miren(1:12000,1); impHH = ECG_IPGH_Miren(1:12000,2); fs=210; ks=0;
load('ECG-IPGH-JosepM.txt');ECG = ECG_IPGH_Miren(1:12000,1); impHH = ECG_IPGH_Miren(1:12000,2); fs=210; ks=0;
if fs==1000
   k=round(400);
   heigth=3;
else
   k=round(100);
   heigth=1;
end
    
Wn = 5*2/fs;
N = 5;                                                                 
[a,b] = butter(N,Wn,'low');                                            
impHH_l = filtfilt(a,b,impHH); 

%% Calculem les derivades
h=0.001; x=[1:length(ECG)]./fs;
ECG_D = diff(ECG)/h;
impHH_D = diff(impHH_l)/h;

%% Cridem a Pan-Tompkins
[IBP,MeanRR] = PTfun(ECG,fs,ks);

%% Trobem màxims i mínims de la IPG
IPGmax=[]; IPGmax_index=[]; IPGmin=[]; IPGmin_index=[];

for ii = 1:length(IBP)
    ECGpeak=IBP(ii);
    
    bottom = max(1,ECGpeak); top = min(ECGpeak+k, length(ECG)-1);

    peakIPG = max(impHH_l(bottom:top)); peakIPG_index = find(impHH_l==peakIPG);
    lowIPG = min(impHH_l(bottom:top)); lowIPG_index = find(impHH_l==lowIPG);
    IPGmax=[IPGmax,peakIPG]; IPGmin=[IPGmin,lowIPG];
    IPGmax_index=[IPGmax_index,peakIPG_index]; IPGmin_index=[IPGmin_index,lowIPG_index];
    
%     figure(5)
%     plot([1:length(impHH_l)]./fs,impHH_l,'b');
%     hold on
%     line([bottom/fs bottom/fs], [-1 heigth]);
%     line([top/fs top/fs], [-1 heigth]);
%     hold on
%     plot(IPGmax_index/fs,IPGmax,'or')
%     hold on
%     plot(IPGmin_index/fs,IPGmin,'og')
%     drawnow
%     hold off
%     pause
end
%% Plots
if kplot
    ax1 = subplot(3,1,1);
    plot(x,ECG)
    title('ECG')
    line(repmat(IBP/fs,[2 1]),...
    repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(IBP)),...
    'LineWidth',1.0,'LineStyle','-.','Color','r');
    title('ECG with peaks corresponding to QRS')
    axis tight
    grid on
    
    ax2 = subplot(3,1,2);
    plot(x,impHH_l)
    title('Impedance')
    hold on
    plot(IPGmax_index/fs,IPGmax,'or')
    hold on
    plot(IPGmin_index/fs,IPGmin,'og')
%     line(repmat(IBP/fs,[2 1]),...
%     repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(IBP)),...
%     'LineWidth',1.0,'LineStyle','-.','Color','r');
    axis tight
    grid on
    
    ax3 = subplot(3,1,3);
    plot([1:length(ECG_D)]./fs,impHH_D)
    linkaxes([ax1,ax2,ax3],'x')
    title('Derivative of Impedance')
    xlabel('Time [s]')
    axis tight
    grid on
    hold off
end