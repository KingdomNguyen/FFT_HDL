function test_fft_radix2_hdl()
% Testbench for radix2_fft_hdl

% Test input (length 8)
x_real = [1 1 1 1 0 0 0 0];
x_imag = zeros(1,8);
x = complex(x_real, x_imag);

% HDL FFT
X_hdl = fft_radix2_hdl(x);

% MATLAB reference FFT
X_ref = fft(x);

% Hiển thị
disp('HDL FFT Output:');
disp(X_hdl);
disp('MATLAB FFT:');
disp(X_ref);
disp('Max Error:');
disp(max(abs(X_hdl - X_ref)));

% Plot
figure;
subplot(2,1,1); stem(abs(X_hdl)); title('HDL FFT');
subplot(2,1,2); stem(abs(X_ref)); title('Reference FFT');
end
