clc;
close all;
clear;

fs = 5;
Amp = 1;
t = 0:1/fs:2*pi; % time vector
sine_wave = Amp * sin(t);

% Add noise
a = 0.1; % upper limit
b = 0; % lower limit
noise = (b - a) .* rand(length(sine_wave), 1) + a; 
noise = noise';
sine_noise = sine_wave + noise;

% Convert from real to integers
total_wordlength = 16;
scaling = 7;
sine_noise_integers = round(sine_noise .* (2^scaling));

% Generate a Gaussian impulse signal
n = length(t);
mu = pi; % Mean
sigma = 0.5; % Standard deviation
impulse = exp(-0.5 * ((t - mu) / sigma).^2);

% Convert from real to integers
impulse_integers = round(impulse .* (2^scaling));

% Generate a direct impulse signal
direct_impulse = zeros(1, n);
mid_index = round(n / 2);
direct_impulse(mid_index) = 1;

% Convert from real to integers
direct_impulse_integers = round(direct_impulse .* (2^scaling));

% Function to write data to a file with a maximum word limit
function write_data_to_file(data, file_name, max_words)
    num_words = length(data);
    if num_words > max_words
        disp(['Data exceeds ', num2str(max_words), ' words. Truncating to ', num2str(max_words), ' words.']);
        data = data(1:max_words); % Truncate data to max_words
    end
    yy = cellstr(data);
    fid = fopen(file_name, 'wt');
    fprintf(fid, '%8s \n', yy{:});
    fclose(fid);
end

% Convert from integers to binary and write to files with max word limit
max_words = 32;
sine_noise_in_binary_signed = dec2bin(mod(sine_noise_integers, 2^total_wordlength), total_wordlength);
write_data_to_file(sine_noise_in_binary_signed, '/opt/Xilinx/Vivado/2023.1/Projects/FIR_Filter/Moving_Average_FIR_Filter/Moving_Average_FIR_Filter.srcs/sim_1/sine.data', max_words);

impulse_in_binary_signed = dec2bin(mod(impulse_integers, 2^total_wordlength), total_wordlength);
write_data_to_file(impulse_in_binary_signed, '/opt/Xilinx/Vivado/2023.1/Projects/FIR_Filter/Moving_Average_FIR_Filter/Moving_Average_FIR_Filter.srcs/sim_1/gaussian_impulse.data', max_words);

direct_impulse_in_binary_signed = dec2bin(mod(direct_impulse_integers, 2^total_wordlength), total_wordlength);
write_data_to_file(direct_impulse_in_binary_signed, '/opt/Xilinx/Vivado/2023.1/Projects/FIR_Filter/Moving_Average_FIR_Filter/Moving_Average_FIR_Filter.srcs/sim_1/direct_impulse.data', max_words);

disp('Text files representing input signals were generated with a maximum of 32 words each');

% Moving Average FIR Filter Model
% Parameters
N = 16;           % Data width
WINDOW_SIZE = 4;  % Window size

% Moving Average FIR Filter Function
function y = moving_average_fir(x, window_size)
    coeff = ones(1, window_size) / window_size;  % Equal coefficients
    y = filter(coeff, 1, x);
end

% Apply the filter to the Gaussian impulse signal
impulse_filtered = moving_average_fir(impulse_integers, WINDOW_SIZE);

% Apply the filter to the centered direct impulse signal
direct_impulse_filtered = moving_average_fir(direct_impulse_integers, WINDOW_SIZE);

% Apply the filter to the sine wave + noise signal
sine_noise_filtered = moving_average_fir(sine_noise_integers, WINDOW_SIZE);

% Interpolation factor
interp_factor = 10;  % Interpolation factor for smoother plots

% Interpolate the signals for smoother plots
t_interp = linspace(min(t), max(t), length(t) * interp_factor);
sine_noise_interp = interp1(1:length(sine_noise_integers), sine_noise_integers, linspace(1, length(sine_noise_integers), length(sine_noise_integers) * interp_factor), 'spline');
sine_noise_filtered_interp = interp1(1:length(sine_noise_filtered), sine_noise_filtered, linspace(1, length(sine_noise_filtered), length(sine_noise_filtered) * interp_factor), 'spline');
impulse_interp = impulse_integers; % No smoothing for Gaussian impulse
impulse_filtered_interp = impulse_filtered; % No smoothing for filtered Gaussian impulse
direct_impulse_interp = interp1(1:length(direct_impulse_integers), direct_impulse_integers, linspace(1, length(direct_impulse_integers), length(direct_impulse_integers) * interp_factor), 'spline');
direct_impulse_filtered_interp = interp1(1:length(direct_impulse_filtered), direct_impulse_filtered, linspace(1, length(direct_impulse_filtered), length(direct_impulse_filtered) * interp_factor), 'spline');

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plot for Sine Wave + Noise (Figure 1)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure(1);
subplot(3, 1, 1);
plot(t, sine_wave);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Original Sine Wave');
grid on;
xlim([min(t) 6.2]); % Limit x-axis to 6.2
ylim([-1.5 1.5]); % y-axis limits from -1.5 to 1.5

subplot(3, 1, 2);
plot(t_interp, sine_noise_interp);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Sine Wave + Noise');
grid on;
xlim([min(t) 6.2]); % Limit x-axis to 6.2
ylim([-150 150]); % y-axis limits from -150 to 150

subplot(3, 1, 3);
plot(t_interp, sine_noise_filtered_interp);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Filtered Sine Wave + Noise');
grid on;
xlim([min(t) 6.2]); % Limit x-axis to 6.2
ylim([-150 150]); % y-axis limits from -150 to 150

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plot for Gaussian Impulse (Figure 2)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure(2);
subplot(2, 1, 1);
plot(1:length(impulse_integers), impulse_integers);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Gaussian Impulse Signal');
grid on;
xlim([1 32]); % x-axis limits from 1 to 32

subplot(2, 1, 2);
plot(1:length(impulse_filtered), impulse_filtered);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Filtered Gaussian Impulse Signal');
grid on;
xlim([1 32]); % x-axis limits from 1 to 32

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plot for Direct Impulse without interpolation (Figure 3)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
direct_impulse_amplified = direct_impulse * 128;

figure(3);
subplot(2, 1, 1);
plot(t, direct_impulse_amplified);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Direct Impulse Signal');
grid on;
xlim([min(t) 6.2]); % Limit x-axis to 6.2
ylim([-0.1*128 1.1*128]);

% Adding marker for peak value
[direct_impulse_amplified, peak_index_amplified] = max(direct_impulse_amplified);
hold on;
plot(t(peak_index_amplified), direct_impulse_amplified, 'ro', 'MarkerSize', 10, 'LineWidth', 1);
legend('Signal', 'Peak');
text(t(peak_index_amplified), direct_impulse_amplified, ['Peak: ', num2str(direct_impulse_amplified)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
hold off;

subplot(2, 1, 2);
plot(t, direct_impulse_filtered);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Filtered Direct Impulse Signal');
grid on;
xlim([min(t) 6.2]); % Limit x-axis to 6.2
ylim([-5 35]); % y-axis limits from -5 to 35

% Ensure that the filter coefficients are plotted correctly with the filtered signal
hold on;
coeff_positions = mid_index:mid_index + 3; % Adjust positions for the filter coefficients
coeff_values = direct_impulse_filtered(coeff_positions);
plot(t(coeff_positions), coeff_values, 'bo', 'MarkerSize', 10, 'LineWidth', 1);
legend('Filtered Signal', 'Filter Coefficients');
text(t(coeff_positions), coeff_values, arrayfun(@num2str, coeff_values, 'UniformOutput', false), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plot for Direct Impulse with interpolation (Figure 4)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure(4);
subplot(2, 1, 1);
plot(t_interp, direct_impulse_interp);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Interpolated Direct Impulse Signal');
grid on;
xlim([min(t_interp) 6.2]); % Limit x-axis to 6.2

[peak_value_interp, peak_index_interp] = max(direct_impulse_interp); 
hold on; 
plot(t_interp(peak_index_interp), peak_value_interp, 'ro', 'MarkerSize', 10, 'LineWidth', 1); 
legend('Signal', 'Peak');
text(t_interp(peak_index_interp), peak_value_interp, arrayfun(@num2str, peak_value_interp, 'UniformOutput', false), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'); 
hold off;

subplot(2, 1, 2);
plot(t_interp, direct_impulse_filtered_interp);
xlabel('\bf Time');
ylabel('\bf Amplitude');
title('\bf Filtered Interpolated Direct Impulse Signal');
grid on;
xlim([min(t_interp) 6.2]);
ylim([-5 40]); % y-axis limits from -5 to 40

hold on;
coeff_positions = round(linspace(mid_index, mid_index + 3, 4) * interp_factor); % Adjust positions for interpolation factor
coeff_values = direct_impulse_filtered_interp(coeff_positions);
plot(t_interp(coeff_positions), coeff_values, 'bo', 'MarkerSize', 10, 'LineWidth', 1);
legend('Filtered Signal', 'Filter Coefficients');
text(t_interp(coeff_positions), coeff_values, arrayfun(@num2str, coeff_values, 'UniformOutput', false), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off;

%% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
% FFT Plots
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_fft(sine_noise_interp, fs, '\bf FFT of Sine Wave + Noise', fs/10);
plot_fft(sine_noise_filtered_interp, fs, '\bf FFT of Filtered Sine Wave + Noise', fs/10);
plot_fft(impulse_interp, fs, '\bf FFT of Gaussian Impulse Signal', fs/10);
plot_fft(impulse_filtered_interp, fs, '\bf FFT of Filtered Gaussian Impulse Signal', fs/10);
plot_fft(direct_impulse_interp, fs, '\bf FFT of Direct Impulse Signal', fs/10);
plot_fft(direct_impulse_filtered_interp, fs, '\bf FFT of Filtered Direct Impulse Signal', fs/10);

function plot_fft(signal, fs, title_str, noise_freq)
    L = length(signal);
    N = 2^nextpow2(L); % Zero-padding to next power of 2
    Y = fft(signal, N); % FFT without shift to keep zero-frequency at beginning
    f = fs * (0:(N/2)) / N; % Adjust frequency axis for positive frequencies only
    P2 = abs(Y / L);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    % Interpolate for smoother plot
    f_interp = linspace(min(f), max(f), 1000);
    P1_interp = interp1(f, P1, f_interp, 'spline');

    [~, fundamental_idx] = max(P1_interp); % Find the index of the fundamental frequency
    fundamental_freq = f_interp(fundamental_idx); % Find the fundamental frequency

    % Define padding for the text annotations
    padding_factor_fundamental = -0.035; % percentage of the maximum amplitude
    padding_factor_noise = 0.010; % percentage of the maximum amplitude
    text_offset_fundamental = padding_factor_fundamental * max(P1_interp); % Calculate text offset
    text_offset_noise = padding_factor_noise * max(P1_interp); % Calculate text offset 

    figure;
    plot(f_interp, P1_interp, 'LineWidth', 1); % Smoother lines with increased thickness
    title(title_str, 'FontWeight', 'bold', 'FontSize', 14);
    xlabel('\bf Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('\bf Amplitude', 'FontWeight', 'bold');
    xlim([0 1.0]); 
    ylim([min(-0.25, min(P1_interp)) max(P1_interp) * 1.2]); % Ensure y-axis starts at 0 or below if negative amplitudes
    grid on;

    % Add marker for fundamental frequency
    hold on;
    plot(fundamental_freq, P1_interp(fundamental_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 1);
    % Position text within the plotting area
    text(fundamental_freq, P1_interp(fundamental_idx) + text_offset_fundamental, num2str(fundamental_freq), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);

    % Add marker for noise frequency
    noise_idx = find(f_interp >= noise_freq, 1); % Find the index of the noise frequency
    if ~isempty(noise_idx)
        plot(noise_freq, P1_interp(noise_idx), 'bo', 'MarkerSize', 10, 'LineWidth', 1);
        % Position text within the plotting area
        text(noise_freq, P1_interp(noise_idx) + text_offset_noise, num2str(noise_freq), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
    end

    % Add legend
    legend({'Signal', 'Fundamental freq.', 'Noise freq.'}, 'Location', 'northeast', 'FontSize', 10);
    hold off;
end




