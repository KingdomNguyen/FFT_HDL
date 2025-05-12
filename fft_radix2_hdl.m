function [fft_out] = fft_radix2_hdl(input_data)
% RADIX2_FFT_HDL Radix-2 FFT implementation for HDL Coder (N = 8 fixed)
%   [fft_out] = radix2_fft_hdl(input_data) computes the FFT of input_data
%   using radix-2 algorithm optimized for HDL generation

% Ensure fixed size and format
N = 8;
input_data = input_data(:);  % ensure column

% Validate input length
assert(length(input_data) == N, 'Input length must be 8');

% Bit-reverse the input
bit_reversed_data = bit_reverse(input_data);

% Compute FFT (statically unrolled)
fft_out = fft_stage(bit_reversed_data);
end

function [reversed_data] = bit_reverse(input_data)
% BIT_REVERSE Perform bit reversal permutation for N = 8
reversed_order = [1, 5, 3, 7, 2, 6, 4, 8];  % MATLAB indexing (bit-reverse of 0:7)
reversed_data = complex(zeros(8,1));
for i = 1:8
    reversed_data(i) = input_data(reversed_order(i));
end
end

function [stage_out] = fft_stage(stage_in)
% FFT_STAGE Explicitly unrolled FFT stages for N = 8
stage_out = stage_in;

% === Stage 1 (distance = 1) ===
for k = 0:2:6
    idx1 = k + 1;
    idx2 = idx1 + 1;
    W = get_twiddle(0, 2);
    t = stage_out(idx2) * W;
    u = stage_out(idx1);
    stage_out(idx1) = u + t;
    stage_out(idx2) = u - t;
end

% === Stage 2 (distance = 2) ===
for k = 0:4:4  % chỉ chạy 2 block
    for j = 0:1
        idx1 = k + j + 1;
        idx2 = idx1 + 2;
        W = get_twiddle(j, 4);
        t = stage_out(idx2) * W;
        u = stage_out(idx1);
        stage_out(idx1) = u + t;
        stage_out(idx2) = u - t;
    end
end

% === Stage 3 (distance = 4) ===
for j = 0:3
    idx1 = j + 1;
    idx2 = idx1 + 4;
    W = get_twiddle(j, 8);
    t = stage_out(idx2) * W;
    u = stage_out(idx1);
    stage_out(idx1) = u + t;
    stage_out(idx2) = u - t;
end
end


function W = get_twiddle(k, N)
% Return twiddle factor W_N^k from static table (N = 8 fixed)

% Precomputed twiddle factors for N=8
persistent W_ROM_8;

if isempty(W_ROM_8)
    angles = -2*pi*(0:7)/8;
    W_ROM_8 = complex(cos(angles), sin(angles));
end

W = W_ROM_8(k + 1); % MATLAB indexing
end
