% 清理环境
clear all; close all; clc;

% 参数设置
Tu = 200e-6; % 有效OFDM符号周期
T = Tu/2048; % 基带基本周期
G = 1/4; % 循环前缀因子，选择1/4、1/8、1/16和1/32
delta = G*Tu; % 保护带时长
Ts = delta+Tu; % 总OFDM符号周期
Kmax = 1705; % 子载波数量
FS = 4096; % IFFT/FFT长度
q = 10; % 载波周期与基本周期的比率
fc = q*1/T; % 载波频率
Rs = 4*fc; % 仿真采样率

% 数据生成器
M = Kmax+1;
a = -1 + 2*round(rand(M,1)) + 1i*(-1 + 2*round(rand(M,1))); % QPSK调制

% IFFT操作准备
info = zeros(FS,1); % 初始化IFFT向量
info(1:(M/2)) = a(1:(M/2)); % 前一半子载波
info(end-(M/2-1):end) = a((M/2+1):end); % 后一半子载波

% 子载波生成
carriers = ifft(info, FS); % IFFT
carriers = carriers.*sqrt(FS); % 归一化功率

% 生成OFDM信号
ofdm_signal = [carriers(end-FS*G+1:end); carriers]; % 添加循环前缀

% 信号聚合（这里只有一个信号，故简单复制）
combined_signal = ofdm_signal; % 在实际中这里可以是多个信号的叠加

% 绘制时域信号
figure;
subplot(2,1,1);
plot(real(combined_signal(1:100)));
title('Time Domain Signal (Real Part)');
xlabel('Sample');
ylabel('Amplitude');

% 信号的频谱
N = length(combined_signal); % 信号长度
f = Rs*(-N/2:N/2-1)/N; % 频率向量
Y = fftshift(fft(combined_signal, N))/N; % 频谱归一化

% 绘制频谱
subplot(2,1,2);
plot(f,abs(Y));
title('Magnitude Spectrum of the Combined Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 绘制功率谱密度
figure;
pwelch(combined_signal,[],[],[],Rs);
title('Power Spectral Density of Combined Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
