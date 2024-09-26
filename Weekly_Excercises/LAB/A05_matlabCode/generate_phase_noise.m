function theta_n = generate_phase_noise(length_of_noise, sigmaDeltaTheta)
    % Create phase noise
    theta_n = zeros(length_of_noise,1);
    %% TODO
    %generate white gaussina noise
    theta_n(1) = rand(1);
    for i = 2:length_of_noise
        theta_n(i) = theta_n(i-1) + sigmaDeltaTheta * randn(1);
    end
    
end