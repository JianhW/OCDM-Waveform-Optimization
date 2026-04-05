function Bitsre=Receiver(M,Block_Num,N,C,Equal,Symbols1,H0,SNR,methodCfg)
%1=ZF, 2=MMSE
methodCfg = get_method_cfg(methodCfg);
P=N+C;
R=[zeros(N,P-N),eye(N)];
Symbols2=zeros(N,1,Block_Num);
%% Remove Cyclic Prefix
for a=1:Block_Num
    Symbols2(:,:,a)=R*Symbols1(:,:,a);
end
%% Construct and Apply DFT matrix
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);
Symbols3=zeros(size(Symbols2));
for count=1:Block_Num
    Symbols3(:,:,count)=FFT*Symbols2(:,:,count);
end
%% Construct and Compensate Zadoff-Chu Sequence
Y=zeros(N);
if mod(N,2)==0
    for k=1:N
        Y(N-k+1,N-k+1)=exp(-1i*pi*k^2/N);
    end
end
if mod(N,2)==1
    for k=1:N
        Y(N-k+1,N-k+1)=exp(-1i*pi*k*(k+1)/N);
    end
end
%% Construct Equalization matrix
% 均衡器类型由主脚本 OCDM.m 中的 equalizerList 控制：Equal=1 为 ZF，Equal=2 为 MMSE
S=eye(N);
T=[S(2*N-P+1:N,:);S];
H=R*H0*T;
D=FFT*H*IFFT; %%Diagonal Matrix
if Equal==1
    G=inv(D);
else
    G=D'/(D*D'+eye(size(D))/SNR);
end
%% Equalization
Symbols4=zeros(size(Symbols3));
Symbols5=zeros(size(Symbols3));

for count=1:Block_Num
    if methodCfg.rxUseY
        Symbols4(:,:,count)=G*Y*Symbols3(:,:,count);
    else
        Symbols4(:,:,count)=G*Symbols3(:,:,count);
    end

    if methodCfg.rxPostIFFT
        Symbols5(:,:,count)=IFFT*Symbols4(:,:,count);
    else
        Symbols5(:,:,count)=Symbols4(:,:,count);
    end
end
%% Demodulation
if M==4
    Symbols6=qamdemod(Symbols5/sqrt(1/2),M);
end
if M==16
    Symbols6=qamdemod(Symbols5/sqrt(1/10),M);
end
if M==64
    Symbols6=qamdemod(Symbols5/sqrt(1/42),M);
end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbols6(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
end















