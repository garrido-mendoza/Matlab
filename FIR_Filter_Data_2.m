% MATLAB Script to Generate Direct Impulse Data for Verilog Testbench

% Parameters
total_wordlength = 16; % Data width
fs = 5; % Manageable sampling frequency for MATLAB (100 Hz)
t = 0:1/fs:2*pi; % Time vector

% Generate a direct impulse signal (single high value in the center)
n = length(t);
direct_impulse = zeros(1, n);
mid_index = round(n / 2);
direct_impulse(mid_index) = 1;

% Convert the impulse signal to 16-bit integer values
direct_impulse_integers = round(direct_impulse * (2^7)); % Scale appropriately

% Upscale the signal for Verilog model
scaling_factor = 500000; % Scaling factor to approximate 50 MHz
direct_impulse_upscaled = repmat(direct_impulse_integers, 1, scaling_factor);

% Convert to binary representation
direct_impulse_in_binary_signed = dec2bin(mod((direct_impulse_upscaled), 2^total_wordlength), total_wordlength);
yy_direct_impulse = cellstr(direct_impulse_in_binary_signed);

% Save the binary data to a text file
fid_direct_impulse = fopen('direct_impulse2.txt', 'wt');
fprintf(fid_direct_impulse, '%16s\n', yy_direct_impulse{:});
fclose(fid_direct_impulse);

disp('Direct impulse data was saved to direct_impulse2.txt');
