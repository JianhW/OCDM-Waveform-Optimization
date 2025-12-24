%% ===========1.Parameter==============

lambda=0.2;
Length=260;
c_length=4;
o_length=256;
c=1/8;
Ec=0.3;

max_interation=100;
threshould=1e-7;


%% ==========2.Function==============
%AISL
function aisl = AISL(x)

    P = length(x);
    
    aisl = 0;
    
    for p = -(P-1):(P-1)
        if p == 0
            continue;
        end
        if p > 0
            r = x(1:P-p).' * conj(x(1+p:P));
        else
            r = x(1-p:P).' * conj(x(1:P+p));
        end
        aisl = aisl + abs(r)^2;
    end
end
% PISL
function pil = PISL(x)

    P = length(x);
    pil = 0;
    
    for k = 1:(P-1)  % 正旁瓣
        r = x.' * conj(circshift(x, k));  % 周期性循环
        pil = pil + abs(r)^2;
    end
    
    for k = -(P-1):-1  % 负旁瓣
        r = x.' * conj(circshift(x, k));  % circshift自动支持负移
        pil = pil + abs(r)^2;
    end
end


%PAPR  
function papr=PAPR(x)
    power=abs(x).^2;
    papr=max(power)/mean(power);
end
%DFnT matrix
function DFnT_matrix=DFnT(P)
    m=(0:P-1)'*ones(1,P);
    n=ones(P,1)*(0:P-1);

    fresnel_kernel=exp(1j*pi/P*(m-n).^2);
    normalization=(1/sqrt(P))*exp((pi/4)*1j);
    DFnT_matrix=normalization*fresnel_kernel;
end

%DFT matrix
function DFT_matrix=DFT(P)
    m=(0:P-1)'*ones(1,P);
    n=ones(P,1)*(0:P-1);
    phase_term=exp(-1j*2*pi/P*m.*n);
    DFT_matrix=phase_term/sqrt(P);
end

%% ===========3.Generate code===========
bits=randi([0,1],1,c_length);
communication=(2*bits-1);
communication=communication.'/norm(communication)*sqrt(Ec);

optimization=randn(o_length,1)+1j*randn(o_length,1);
optimization=optimization/norm(optimization)*sqrt(1-Ec);

s=[communication;optimization];

disp(norm(s));

s_timedomain=DFnT(Length)'*s;
%disp("======Signal=======");

%disp("======PISL=======");
%disp(AISL(s_timedomain));
%disp("======PAPR=======");
%disp(PAPR(s_timedomain));


%% ===========4.MM optimize===========
%扩充矩阵T
matrix2=DFnT(2*Length);
matrix2=matrix2(Length+1:2*Length,1:Length);
T=[eye(Length);matrix2*DFnT(Length)';eye(Length);zeros(Length)];


%两个矩阵B、A
B=[];
A=[];
for k=0:2*Length-1
    a_k = (1/sqrt(Length))*exp(1i * 2*pi*k/(2*Length) * (0:(2*Length-1)))';
    b_k=DFnT(2*Length)*a_k;
    fit_0=zeros(2*Length,1);
    a_k=[fit_0;a_k];
    b_k=[b_k;fit_0];
    A=[A,a_k];
    B=[B,b_k];
end


%最大特征值
lambda_max=max(1,lambda);

data_AISL=[];
data_PAPR=[];

for k=1:max_interation
    %fprintf("==========当前迭代第 %d 次==========\n",k);

    %hat(s)'
    ss=T*s;

    %mu_Phi and mu_Gamma
    mu_Phi=B'*ss;
    mu_Phi=abs(mu_Phi).^2;
    mu_Gamma=A'*ss;
    mu_Gamma=abs(mu_Gamma).^2;

    lambda1=max(mu_Phi);
    lambda2=max(mu_Gamma);
    lambda_max1=max(lambda1,lambda2);
    %W阵
    W=B*(diag(mu_Phi)-lambda1/2*eye(2*Length))*B'+lambda*A*(diag(mu_Gamma)-lambda2/2*eye(2*Length))*A'-lambda_max/2*eye(4*Length);

    %-Upsilion
    upsilion=-T'*W*ss;
    o = upsilion(c_length+1:Length);
    optimize = sqrt(1 - Ec) * o / norm(o);
    s=[communication;optimize];

    s_timedomain=DFnT(Length)'*s;

    AISL_k=AISL(s_timedomain);
    PAPR_k=PAPR(s_timedomain);

    data_AISL=[data_AISL,AISL_k];
    data_PAPR=[data_PAPR,PAPR_k];
    


end

data_obj=data_AISL+lambda*data_PAPR;

subplot(1,3,1);
plot(data_AISL);
subplot(1,3,2);
plot(data_PAPR);
subplot(1,3,3);
plot(data_obj);
