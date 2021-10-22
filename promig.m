function[meanmm,suma,I] = promig(ECG, Z, type, fs)
%% Ens dóna la integral de la senyal i un promig de la forma de la ona de pols
%The variables we are going to work with are called
[IBP,IPGmax, IPGmin, IPGmax_index,IPGmin_index]=characterization4_0(ECG,Z,type,fs);
%We generate both vectors and matrixes we are goint to use to store data
Ivec = [];
mat = zeros(10000,100);
%The average minimum to minimum distance is computed 
meanmm = mean(diff(IPGmin_index));
%For each minimum, take the pulse wave
for i=1:length(IPGmin)-1

    xmin = IPGmin_index(i);
    xmax = IPGmin_index(i+1);
    if abs(xmin-xmax) < 1.5*meanmm &&  abs(xmin-xmax)>0.5*meanmm
        Z_temp = Z(xmin:xmax); 
        mat(i,1:length(Z_temp))=Z_temp;
        if min(Z_temp)<0
            Z_temp = Z_temp + abs(min(Z_temp));
        end
        I = trapz(Z_temp);
        Ivec = [Ivec, I];
    end
end
Ivec1 = Ivec;
suma = 0;
for i=1:100
    suma = suma + mat(i,:);
end
hold off
suma = suma/length(IPGmin);
I = sum(Ivec1)/length(Ivec1);
end