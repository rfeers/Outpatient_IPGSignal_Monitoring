%% ======================= Pan-Thompskins Meu ==========================%%
%Sampling Freq; 
ks=1;k=0;
% load('ECG-IPGH-JosepM.txt');ECG = ECG_IPGH_JosepM(1:12000,1); fs=210;
load('ECG-IPGH-Miren.txt');ECG = ECG_IPGH_Miren(1:12000,1); fs=210; 
% load('S16_2.mat'); fs=1000;
%% ======================= BandPass [5,12] Hz ==========================%%
%% ============ Noise cancelation(Filtering)( 5-15 Hz) =============== %%
if fs == 210
% ------------------ remove the mean of Signal -----------------------%
%% ==== Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2 ==== %%
%%It has come to my attention the original filter doesnt achieve 12 Hz
%    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
%    a = [1 -2 1];
%    ecg_l = filter(b,a,ecg); 
%    delay = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Wn = 12*2/fs;
   N = 3;                                                                  % order of 3 less processing
   [a,b] = butter(N,Wn,'low');                                             % bandpass filtering
   ECG_L = filtfilt(a,b,ECG); 
 %% ======================= start figure ============================= %%
%% ==== High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1)) ==== %%
%%It has come to my attention the original filter doesn achieve 5 Hz
%    b = zeros(1,33);
%    b(1) = -1; b(17) = 32; b(33) = 1;
%    a = [1 1];
%    ecg_h = filter(b,a,ecg_l);    % Without Delay
%    delay = delay + 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Wn = 5*2/fs;
   N = 3;                                                                  % order of 3 less processing
   [a,b] = butter(N,Wn,'high');                                            % bandpass filtering
   ECG_BP = filtfilt(a,b,ECG_L); 
else
f1=5;                                                                      % cuttoff low frequency to get rid of baseline wander
f2=15;                                                                     % cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
N = 3;                                                                     % order of 3 less processing
[a,b] = butter(N,Wn);                                                      % bandpass filtering
ECG_BP= filtfilt(a,b,ECG);
end
%% ======================= Derivative Filter ==========================%%
% --------- H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) ------------- %
if fs ~= 210
 int_c = (5-1)/(fs*1/40);
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
 b = [1 2 0 -2 -1].*(1/8)*fs;   
end
ECG_D = filtfilt(b,1,ECG_BP);

%% ========================== Squaring  ================================%%
ECG_S = ECG_D.^2;
%% ======================= Moving Windows  =============================%%
% N = 300; %Number of samples to convolate
% ECG_W = conv(ECG_S, ones(1,N));
ECG_W = conv(ECG_S ,ones(1 ,round(0.150*fs))/round(0.150*fs));
%Generamos un vector de muestras y lo pasamos a tiempo:
x = [1:length(ECG)]; tvec = x/fs;
if k
    figure(1)
    subplot(211)
    plot(tvec, ECG_BP)
    grid on
    xlabel('Time [s]')
    ylabel('Amplitude')
    title('Band Pass Filtered Signal')
    axis tight
    subplot(212)
    plot([1:length(ECG_W)]./10^4,ECG_W);
    grid on
    xlabel('Nº of Samples x 10^4')
    ylabel('Amplitude')
    axis tight
    title('Window Integrated Filtered Signal')
end
%% ========================== ALGORITHM ================================%%
%% ======================== FIDUCIAL MARK ==============================%%
%Los picos tienen que estar separados por 200ms como mínimo. Los picos de 
%la señal integrada coinciden con los QRS. Haremos una primera selección de
%picos y continuaremos para seleccionar solo los picos de verdad.
[peaksW,indexW] = findpeaks(ECG_W,'MINPEAKDISTANCE',round(0.2*fs)); 
LpeaksW = length(peaksW);
if k
    figure(2)
    plot([1:length(ECG_W)]./10^4,ECG_W,'b');
    hold on
    plot(indexW./10^4,peaksW,'or');
end
%% ======================== INICIALIZATION =============================%%
%% Windows Integrated Signal
SigPeakW = mean(ECG_W(1:fs))/2; SigPeakW_vec = zeros(1,LpeaksW);
NoiPeakW = mean(ECG_W(1:fs))/3; NoiPeakW_vec = zeros(1,LpeaksW);
THRW = NoiPeakW + 0.25*(SigPeakW-NoiPeakW);THRW_vec = zeros(1,LpeaksW);
THRW_Noise = 0.5*THRW;
%% BandPass Filtered Signal
SigPeakBP = mean(ECG_BP(1:fs))/2; SigPeakBP_vec = zeros(1,LpeaksW);
NoiPeakBP = mean(ECG_BP(1:fs))/3; NoiPeakBP_vec = zeros(1,LpeaksW);
THRBP = NoiPeakBP + 0.25*(SigPeakBP-NoiPeakBP);THRBP_vec = zeros(1,LpeaksW);
THRBP_Noise = 0.5*THRBP;
if k
    figure(3)
    subplot(211)
    title('Windows Integrated Signal with Thresholds')
    plot([1:length(ECG_W)]./10^4,ECG_W,'b');
    hold on
    plot(indexW./10^4,SigPeakW_vec,'g');
    hold on
    plot(indexW./10^4,NoiPeakW_vec,'k');
    hold on
    plot(indexW./10^4,THRW_vec,'r');
    axis tight
    hold off
    subplot(212)
    title('Band Pass Filtered Signal with Thresholds')
    plot([1:length(ECG_BP)]./fs,ECG_BP,'b');
    hold on
    plot([1:LpeaksW]./fs,SigPeakBP_vec,'g');
    hold on
    plot([1:LpeaksW]./fs,NoiPeakBP_vec,'k');
    hold on
    plot([1:LpeaksW]./fs,THRBP_vec,'r');
    axis tight
    hold off
end
%% ======================== THRESHOLDS ==============================%%
AW=[]; IW=[]; ABP=[]; IBP=[]; Aprova=[];Iprova=[]; %Void vectors are generated to save the possible peaks
for ii = 1:LpeaksW
%     if length(IW)>2
%         if (IW(end)-IW(end-1))<350
%             if length(IW)>length(IBP)
%                 IBP(end) = []; ABP(end) = [];  
%             end
%             if (IW(end-1)-IW(end-2))<350
%                IW(end-1)=[]; AW(end-1)=[];
%             else
%             IW(end) = []; AW(end) = [];
%             end
%         end
%     end
    if peaksW(ii) > THRW %The peak can be QRS
            pos_peak = find(ECG_W==peaksW(ii)); %The index of the peak is found
            AW = [AW, peaksW(ii)]; IW = [IW, pos_peak]; %Vectors are updated
            SigPeakW =0.125*peaksW(ii) + 0.875*SigPeakW;
            
            %% Now we search for peaks in the Band Pass Filteres Signal
            pos_peak_BP = round(pos_peak*length(ECG)/length(ECG_W));
            bottom = max(1,(pos_peak_BP-100)); top = min(pos_peak_BP+100, length(ECG_BP)-1);
            peakBP = max(ECG_BP(bottom:top));
%             [pksBP,indexBP] = findpeaks(ECG_BP(bottom:top),'MINPEAKDISTANCE',round(0.2*fs));
%             for jj=1:length(pksBP)
%                 peakBP=pksBP(jj);
                if peakBP > THRBP
                    pos_peak = find(ECG_BP==peakBP);
                    ABP = [ABP, peakBP];IBP = [IBP, pos_peak];
                    Aprova=[Aprova, peakBP];Iprova=[Iprova,pos_peak];
                    SigPeakBP =0.125*peakBP + 0.875*SigPeakBP;
%                     figure(5)
%                     plot([1:length(ECG_BP)]./fs,ECG_BP,'b');
%                     hold on
%                     plot(IBP./fs,ABP,'ob'); 
%                     line([bottom/fs bottom/fs], [0 1]);
%                     line([top/fs top/fs], [0 1]);
%                     drawnow
%                     hold off
%                     pause
                    
                else
                    NoiPeakBP =0.125*peakBP + 0.875*NoiPeakBP;
                end
%             end
            elseif THRW_Noise<peaksW(ii) && THRW_Noise<peaksW(ii)%The peak is noise
                NoiPeakW =0.125*peaksW(ii) + 0.875*NoiPeakW;
                %NoiPeakBP =0.125*peakBP + 0.875*NoiPeakBP;
    end
%% ======================== SEARCH BACK ==============================%%
    if length(AW)>9
        MeanRR = mean(diff(IW));
        if (IW(end)-IW(end-1))>1.66*MeanRR
            THRW = THRW/2;SigPeakW = SigPeakW/2; NoiPeakW = NoiPeakW/2;
           peakW = max(ECG_W(IW(end-1)-500:IW(end)+500));
           if peakW > THRW/2
               pos_peak = find(ECG_W==peakW);
               pk_temp = AW(end); pos_temp = IW(end);
               AW(end) = []; IW(end) = [];
               AW = [AW, peakW, pk_temp]; IW = [IW, pos_peak, pos_temp];
           end
               
        end
        if (IBP(end)-IBP(end-1))>1.66*MeanRR
               THRBP = THRBP/2;SigPeakBP = SigPeakBP/2; NoiPeakBP = NoiPeakBP/2;
               pos_peak = (IBP(end-2)+IBP(end))/2;
               pos_peak_BP = round(pos_peak*length(ECG)/length(ECG_W));
               peakBP = max(ECG_BP(pos_peak_BP-400:pos_peak_BP+400));
               if peakBP > THRBP/2
                    pos_peak = find(ECG_BP==peakBP);
                    pk_temp = ABP(end); pos_temp = IBP(end);
                    ABP(end) = []; IBP(end) = [];
                    ABP = [ABP, peakBP, pk_temp]; IBP = [IBP, pos_peak, pos_temp];
               end
        end
    end
%% ==================== Detection for T Waves ==========================%%
if length(IBP)>3
    if (IBP(end)-IBP(end-1)) <= round(0.3600*fs)
        Slope1 = mean(diff(ECG_BP(IBP(end)-round(0.075*fs):IBP(end))));   
        Slope2 = mean(diff(ECG_BP(IBP(end-1)-round(0.075*fs):IBP(end-1))));
        if abs(Slope1) <= abs(0.5*(Slope2))  
            IW(end) = []; AW(end) = [];
            IBP(end) = []; ABP(end) = [];
        end
    end
end

THRW = NoiPeakW + 0.25*(SigPeakW-NoiPeakW);
THRW_Noise = 0.5*THRW;  
THRBP = NoiPeakBP + 0.25*(SigPeakBP-NoiPeakBP);
THRBP_Noise = 0.5*THRBP;  
%% ======================== UPDATE VECTORS ==============================%%
SigPeakW_vec(ii) = SigPeakW;
NoiPeakW_vec(ii) = NoiPeakW;
THRW_vec(ii) = THRW;
SigPeakBP_vec(ii) = SigPeakBP;
NoiPeakBP_vec(ii) = NoiPeakBP;
THRBP_vec(ii) = THRBP;

% figure(5)
% plot([1:length(ECG_W)]./fs,ECG_W,'b');
% hold on
% plot(IW./fs,AW,'ob');    
% drawnow
end
if ks
    figure(4)
    ax1 = subplot(3,1,1);
    plot([1:length(ECG_BP)]./fs,ECG_BP,'b');
    hold on
    plot(IBP./fs,ABP,'ob');
    hold on
    plot(Iprova./fs,Aprova,'oy');
    hold on
    plot(indexW./fs,SigPeakBP_vec,'g');
    hold on
    plot(indexW./fs,NoiPeakBP_vec,'k');
    hold on
    plot(indexW./fs,THRBP_vec,'r');
    hold off
    axis tight
    title('Band Pass Filtered Signals with peaks corresponding to QRS')
    
    ax2 = subplot(3,1,2);
    plot([1:length(ECG_W)]./fs,ECG_W,'b');
    hold on
    plot(IW./fs,AW,'ob');    
    hold on
    plot(indexW./fs,SigPeakW_vec,'g');
    hold on
    plot(indexW./fs,NoiPeakW_vec,'k');
    hold on
    plot(indexW./fs,THRW_vec,'r');
    axis tight
    hold off
    title('Window Integrated Signals with peaks corresponding to QRS')
    
    ax3 = subplot(3,1,3);
    plot([1:length(ECG)]./fs,ECG_BP,'b');
    hold on
    line(repmat(IBP/fs,[2 1]),...
    repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(IBP)),...
    'LineWidth',1.0,'LineStyle','-.','Color','r');
    title('ECG with peaks corresponding to QRS')
    axis tight
    linkaxes([ax1,ax2,ax3],'x')
    hold off
end