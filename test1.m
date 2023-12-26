% 清理环境
clear all; close all; clc;

% 输入参数-1（与第一段代码相匹配）
Tu = 200e-6; % 有效OFDM符号周期
T = Tu/2048; % 基带基本周期
G = 1/4; % 循环前缀因子，选择1/4、1/8、1/16和1/32
delta = G*Tu; % 保护带时长
Ts = delta+Tu; % 总OFDM符号周期
Kmax = 1705; % 子载波数量
FS = 4096; % IFFT/FFT长度
q = 10; % 载波周期与基本周期的比率
fc = q*1/T; % 载波频率
Rs = 4*fc; % 仿真周期

% 数据生成器-1（与第一段代码相匹配）
M = Kmax + 1;
a = -1 + 2*round(rand(M,1)) + 1i*(-1 + 2*round(rand(M,1)));

% IFFT操作准备-1（与第一段代码相匹配）
info = zeros(FS,1); % 初始化IFFT向量
info(1:(M/2)) = a(1:(M/2)); % 前一半子载波
info(end-(M/2-1):end) = a((M/2+1):end); % 后一半子载波

% 子载波生成-1（与第一段代码相匹配）
carriers = FS.*ifft(info, FS); % IFFT
carriers = carriers.*sqrt(FS); % 归一化功率

% 生成OFDM信号-1（与第一段代码相匹配）
ofdm_signal = [carriers(end-FS*G+1:end); carriers]; % 添加循环前缀

% 绘制时域信号-1（与第一段代码相匹配）
figure;
subplot(2,1,1);
plot(real(ofdm_signal(1:100)));
title('Time Domain Signal (Real Part)');
xlabel('Sample');
ylabel('Amplitude');

% 信号的频谱-1（与第一段代码相匹配）
N = length(ofdm_signal); % 信号长度
f = Rs*(-N/2:N/2-1)/N; % 频率向量
Y = fftshift(fft(ofdm_signal, N))/N; % 频谱归一化

% 绘制频谱-1（与第一段代码相匹配）
subplot(2,1,2);
plot(f,abs(Y));
title('Magnitude Spectrum of the OFDM Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 绘制功率谱密度-1（与第一段代码相匹配）
figure;
pwelch(ofdm_signal,[],[],[],Rs);
title('Power Spectral Density of the OFDM Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% 输入参数-2（B）
Tu1 = 200e-6; % 有效OFDM符号周期
T1 = Tu1/2048; % 基带基本周期
G1 = 0; % 选择1/4、1/8、1/16和1/32
delta1 = G1*Tu1; % 保护带时长
Ts1 = delta1+Tu1; % 总OFDM符号周期
Kmax1 = 1705; % 子载波数量
Kmin1 = 0;
FS1 = 4096; % IFFT/FFT 长度
q1 = 10; % 载波周期与基本周期的比率
fc1 = q1*1/T1; % 载波频率
Rs1 = 4*fc1; % 仿真周期
t1 = 0:1/Rs1:Tu1;

% 数据生成器-2（B）
M1 = Kmax1 + 1;
rand('state',0);
a1 = -1 + 2*round(rand(M1,1)).' + i*(-1 + 2*round(rand(M1,1))).';
A1 = length(a1);
info1 = zeros(FS1,1);
info1(1:(A1/2)) = [a1(1:(A1/2)).']; % 零填充
info1((FS1-((A1/2)-1)):FS1) = [a1(((A1/2)+1):A1).'];

% 子载波生成-2（B）
carriers1 = FS1.*ifft(info1, FS1);
tt1 = 0:T1/2:Tu1;
figure(11);
subplot(211);
stem(tt1(1:20),real(carriers1(1:20)));

subplot(212);
stem(tt1(1:20),imag(carriers1(1:20)));
figure(12);
f1 = (2/T1)*(1:(FS1))/(FS1);
subplot(211);
plot(f1,abs(fft(carriers1, FS1))/FS1);
subplot(212);
pwelch(carriers1,[],[],[],2/T1);
