
%% データの読み込み

folder = '/Users/guregoa/Documents/MATLAB/snirf/snirf/2/RS4_SL_4216.snirf';
raw = nirs.io.loadSNIRF(folder);
raw_data_HbO = raw.data(:, 2:2:end);

% 図を作る
figure;
% 秒単位の時間ベクトルを作成
temps_total = 560; % 合計秒数
intervalles = temps_total / (length(raw_data_HbO) - 1); % 点間の間隔
temps = 0:intervalles:temps_total; %　時間ベクトル


plot(temps,raw_data_HbO);
title('波長850nmのデータ');
xlabel('時間（秒）');
ylabel('データ');


%% 前処理

% ステップ １：データを５０００サンプルに切り出す
raw.data = raw.data(end - 5000 + 1:end, :);
raw_data_HbO = raw.data(:, 2:2:end);

%　図を作る
figure;
temps_total = 560; 
intervalles = temps_total / (length(raw_data_HbO) - 1);
temps = 0:intervalles:temps_total; 
plot(temps,raw_data_HbO);
title('波長850nmのデータ');
xlabel('時間（秒）');
ylabel('データ');

% ステップ２：光度への変換
j = nirs.modules.OpticalDensity();
hb = j.run(raw);
data_HbO = hb.data(:, 2:2:end);

% 図を作る
figure;
temps_total = 560; 
intervalles = temps_total / (length(data_HbO) - 1); 
temps = 0:intervalles:temps_total; 
plot(temps,data_HbO);
title('波長850nmの光度データ');
xlabel('時間（秒）');
ylabel('データ');


% ステップ３：ウェーブレットを用いてノイズを減少する
j = nirs.modules.WaveletFilter();
j.sthresh = 3; %標準偏差値
hb = j.run(hb);
data_HbOwavelet = hb.data(:, 2:2:end);

% 図を作る
figure;
temps_total = 560; 
intervalles = temps_total / (length(data_HbOwavelet) - 1); 
temps = 0:intervalles:temps_total; 
plot(temps,data_HbOwavelet);
title('ウェーブレット後のデータ');
xlabel('時間（秒）');
ylabel('データ');

% 図を作る

figure;
temps_total = 560; 
intervalles = temps_total / (length(data_HbO) - 1); 
temps = 0:intervalles:temps_total;

% サブプロット一番目
subplot(2, 1, 1);
plot(temps,data_HbO);
title('');
xlabel('時間（秒）');
ylabel('データ');

% サブプロット二番目
subplot(2, 1, 2);
plot(temps,data_HbOwavelet);
title('Plot de data\_HbOwavelet');
xlabel('時間（秒）');
ylabel('データ');


% ステップ ４：ＭＬＢ則によるヘモグロビンへの変換
j = nirs.modules.BeerLambertLaw();
j.PPF = 4.2;
hb = j.run(hb);
data_HbO = hb.data(:, 2:2:end);
data_HbO_2 = hb.data(:, 2:2:end);

% 図を作る
figure;
plot(temps,data_HbO);
title('MLB後のデータ');
xlabel('時間（秒）');
ylabel('データ');

% ステップ ５: グローバル信号回帰
global_signal = mean(hb.data, 2);
lpf = 0.09;
Fs = hb.Fs;

% レジェンド多項式の作成 (ハイパスフィルタ)
n_data = size(hb.data, 1);
s_data = n_data / Fs;
k = 1 + floor(s_data / 150);
n = linspace(-1, 1, n_data)';
L = zeros(n_data, k + 1);

for i = 1:k + 1
  tmp = legendre(i - 1, n);
  tmp = tmp(1, :);
  L(:, i) = tmp / max(abs(tmp));  
end


% ローパス・フィルター用のサインとコサインの行列の作成
dft_matrix = dftmtx(n_data);
idx = floor((lpf / Fs) * n_data);
dft_matrix_lpf = dft_matrix(idx:n_data - idx + 1, :);
sin_lpf_mtx = imag(dft_matrix_lpf);
cos_lpf_mtx = real(dft_matrix_lpf);
lpf_mtx = [cos_lpf_mtx' sin_lpf_mtx'];

% 回帰行列の作成（ハイパスフィルター、ローパスフィルター、グローバル信号）
reg_mat = [L lpf_mtx global_signal];
beta_data = pinv(reg_mat) * hb.data;
hb.data = hb.data - reg_mat * beta_data;
data_HbO = hb.data(:, 2:2:end);

% 図を作る
figure;
plot(temps,data_HbO);
title('フィルターの後');
xlabel('時間（秒）');
ylabel('データ');


fs = hb.Fs;
t = 0:1/fs:1-1/fs;
x = data_HbO_2;
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

plot(freq,pow2db(psdx));
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")

fs = hb.Fs;
t = 0:1/fs:1-1/fs;
x = data_HbO;
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;
gain_factor = 100;

plot(freq,10*log10(gain_factor * psdx));
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")

