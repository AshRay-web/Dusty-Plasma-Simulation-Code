load experimentaloutput.mat
% Given parameters
% Qd = ...   % your array of dust charge
% epsilon = ... % permittivity
% rd = ...      % dust radius
% lambda_D = ...% Debye length
% beta = ...    % some parameter

% Compute phi_d
phi_d = Qd ./ (4 * pi * epsilon * rd .* (1 + (rd ./ lambda_D) .* sqrt(1 + beta)));


