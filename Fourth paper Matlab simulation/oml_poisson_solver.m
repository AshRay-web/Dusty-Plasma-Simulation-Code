% oml_poisson_solver.m
% Solves OML locally for each gridpoint and returns Qd_array and Qd0
% Inputs: ne(x), ni(x) (arrays across grid)
% Output: Qd_array (Coulombs), Qd0 (reference at presheath entrance)

function [Qd_array, Qd0,phi_d] = oml_poisson_solver(ne, ni)

    % Load constants and plasma parameters
    Input   % <-- your file that defines: Tpse, Tpsi, rd, epsilon, e, Me, Mi, nps, n_dust, etc.

    N = length(ne);
    Qd_array = zeros(1, N);
    phi_d    = zeros(1, N);

    % ---- Reference presheath entrance values ----
    ne0 = nps;  % from Input
    ni0 = nps;  % quasi-neutral entrance

    % Compute reference dust charge at presheath entrance
    Qd0 = compute_Qd_OML(nps, nps, Tpse, Tpsi, rd, epsilon, e, Me, Mi);

    % ---- Local Qd(x) at each gridpoint ----
    for i = 1:N
        [Qd_array(i),phi_d(i)] = compute_Qd_OML(ne(i), ni(i), Tpse, Tpsi, rd, epsilon, e, Me, Mi);
    end

end

