function cc = extract_cc(response)
    % response size should be (num_points, num_DOFs)
    cc = ifft(2*real(log(fft(response))));
end