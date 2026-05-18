% compute_Qd_OML.m
% Computes dust grain charge Qd using OML current balance
% Inputs: ne0, ni0, Tpse, Tpsi, rd, epsilon, e, Me, Mi
% Output: Qd (Coulombs)

function [Qd,phi_d] = compute_Qd_OML(ne0, ni0, Tpse, Tpsi, rd, epsilon, e, Me, Mi)

    % Boltzmann constant (if not already global)
    kB = 1.380649e-23;

    % Convert temperatures to Joules
    Te_J = Tpse * kB;
    Ti_J = Tpsi * kB;

    % Debye length (electron)
    lambda_D = sqrt(epsilon * Te_J / (ne0 * e^2));

    % Collection area of dust
    A_surf = 4 * pi * rd^2;

    % Electron and ion thermal flux velocities
    vth_e_flux = sqrt(Te_J / (2 * pi * Me));
    vth_i_flux = sqrt(Ti_J / (2 * pi * Mi));

    % Electron current to grain (OML, for negative phi)
    Ie_fun = @(phi) - e * A_surf * ne0 * vth_e_flux * exp(e * phi / Te_J);

    % Ion current to grain (OML linearized form)
    Ii_fun = @(phi)   e * A_surf * ni0 * vth_i_flux * (1 - (e * phi) / Ti_J);

    % Net current
    I_net = @(phi) Ii_fun(phi) + Ie_fun(phi);

    % Solve for floating potential phi_d (negative root, start ~ -3Te/e)
    phi_d = fzero(I_net, -3*Te_J/e);

    % Temperature ratio beta = Te/Ti (in eV)
    beta = Tpse / Tpsi;

    % Yukawa-corrected capacitance for dust charge
    Qd = 4 * pi * epsilon * rd * (1 + (rd / lambda_D) * sqrt(1 + beta)) * phi_d;

end

