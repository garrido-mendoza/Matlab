% Define filter specifications
FILTER_TAPS = 60;              % Number of filter taps (59th-order + 1 for symmetry)
INPUT_WIDTH = 24;              % Bit-width of the input data
OUTPUT_WIDTH = 24;             % Bit-width of the output data
MAC_WIDTH = INPUT_WIDTH;       % Multiply-Accumulate width (ignoring COEFF_WIDTH)

% Design a 59th-order Low Pass Filter with a cutoff frequency of 1 kHz
Fs = 44.1e3;  % Sampling frequency
Fc = 1e3;     % Cutoff frequency
coeff_s_float = fir1(FILTER_TAPS - 1, Fc / (Fs / 2));  % Normalized cutoff frequency

% Convert floating-point coefficients to fixed-point (2's complement, 16-bit)
coeff_s_fixed = fi(coeff_s_float, true, 16, 15);

% Display the converted fixed-point coefficients in hexadecimal
coeff_s_hex = arrayfun(@(x) dec2hex(bin2dec(x.bin)), coeff_s_fixed, 'UniformOutput', false);
% disp('Converted Fixed-Point Coefficients (Hexadecimal):');
% disp(coeff_s_hex);

% Assuming the original floating-point coefficients for comparison
original_coefficients = coeff_s_float;  % Use the designed floating-point coefficients

% Compare the converted coefficients with the original coefficients
% disp('Original Floating-Point Coefficients:');
% disp(original_coefficients);

% Calculate the difference
coeff_difference = double(coeff_s_fixed) - original_coefficients;
% disp('Difference between Converted and Original Coefficients:');
% disp(coeff_difference);

% Define fixed-point types
T_input = sfix(INPUT_WIDTH);
T_output = sfix(OUTPUT_WIDTH);
T_mac = sfix(MAC_WIDTH);

% Define simulation parameters
num_samples = 1000;  % Increase number of samples for higher resolution
clk_period = 1 / Fs;

% Step input
step_input = fi([zeros(1, 500), ones(1, 500)], T_input, INPUT_WIDTH, INPUT_WIDTH-1);

% Sine wave input (1 kHz sine wave)
t = (0:num_samples-1) / Fs;
sine_wave_input = fi(sin(2 * pi * 1e3 * t), T_input, INPUT_WIDTH, INPUT_WIDTH-1);

% Random noise input
noise_amplitude = 0.1;
rand_noise_input = fi(noise_amplitude * (rand(1, num_samples) - 0.5) * 2, T_input, INPUT_WIDTH, INPUT_WIDTH-1);

% Impulse input
impulse_input = fi([1, zeros(1, num_samples-1)], T_input, INPUT_WIDTH, INPUT_WIDTH-1);

% Create a function to simulate the pipelined transposed FIR filter
function data_o = simulate_fir_filter(input_data, FILTER_TAPS, T_input, T_mac, coeff_fp)
    % Initialize registers
    areg_s = fi(zeros(FILTER_TAPS, 1), T_input);
    mreg_s = fi(zeros(FILTER_TAPS, 1), T_mac);
    preg_s = fi(zeros(FILTER_TAPS, 1), T_mac);
    data_o = fi(zeros(length(input_data), 1), T_mac);
    
    % Process each sample individually
    for sample_idx = 1:length(input_data)
        % Shift the input registers
        areg_s(2:end) = areg_s(1:end-1);
        areg_s(1) = input_data(sample_idx);
        
        % Perform multiplication and accumulate with pipelining
        for i = 1:FILTER_TAPS
            % Perform multiplication
            mreg_s(i) = areg_s(i) * coeff_fp(i);
            
            % Perform accumulation
            if i == 1
                preg_s(i) = mreg_s(i);
            else
                preg_s(i) = mreg_s(i) + preg_s(i-1);
            end
        end
        
        % Output data
        data_o(sample_idx) = preg_s(FILTER_TAPS);
    end
end

% Run simulations for each input type
step_output = simulate_fir_filter(step_input, FILTER_TAPS, T_input, T_mac, coeff_s_fixed);
sine_wave_output = simulate_fir_filter(sine_wave_input, FILTER_TAPS, T_input, T_mac, coeff_s_fixed);
rand_noise_output = simulate_fir_filter(rand_noise_input, FILTER_TAPS, T_input, T_mac, coeff_s_fixed);
impulse_output = simulate_fir_filter(impulse_input, FILTER_TAPS, T_input, T_mac, coeff_s_fixed);

% Plot the outputs with higher resolution
figure;
plot(double(step_output));
title('Step Output');
ylabel('Amplitude');
xlabel('Sample');
grid on;
axis auto;

figure;
plot(double(sine_wave_output));
title('Sine Wave Output');
ylabel('Amplitude');
xlabel('Sample');
grid on;
axis auto;

figure;
plot(double(rand_noise_output));
title('Random Noise Output');
ylabel('Amplitude');
xlabel('Sample');
grid on;
axis auto;

figure;
plot(double(impulse_output));
title('Impulse Output');
ylabel('Amplitude');
xlabel('Sample');
grid on;
axis auto;

% % Display the outputs
% disp('Step input output:');
% disp(double(step_output));
% disp('Sine wave input output:');
% disp(double(sine_wave_output));
% disp('Random noise input output:');
% disp(double(rand_noise_output));
% disp('Impulse input output:');
% disp(double(impulse_output));
