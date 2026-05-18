%This is a new OML SOLVER from the thesis


function rho_d_array = oml_poisson_solver(ne, ni)
% OML_SOLVER  Solves dust grain floating potential using OML theory.
%
%   phi_d = oml_solver(ne, Te, ni, Ti, mi, Z, rd)
%
% Inputs:
%   ne  - electron density [m^-3]
%   Te  - electron temperature [eV]
%   ni  - ion density [m^-3]  (can be vector for multiple ion species)
%   Ti  - ion temperature [eV] (same size as ni)
%   mi  - ion masses [kg]      (same size as ni)
%   Z   - ion charge state (e.g., 1 for singly charged) (same size as ni)
%   rd  - dust grain radius [m]
%
% Output:
%   phi_d - floating potential of dust grain [V]


Input  % Loads constants: Me, Mi, epsilon, kB, e, Tpse, Tpsi, rd, ndust, etc.

N = length(ne);              % number of points
rho_d_array = zeros(1, N);   % preallocate output

for i = 1:N
      % Current densities for this iteration
      ne0 = ne(i);
      ni0 = ni(i);

      me = Me;
      mi = Mi;

      % convert temps to Joules
      Te_J = Tpse * kB;
      Ti_J = Tpsi * kB;

      % Debye length (electron)
      lambda_D = sqrt(epsilon * Te_J / (ne0 * e^2));

      %Compute phi_d by OML current balance
      A_surf = 4 * pi * rd^2;
      vth_e_flux = sqrt(Te_J / (2 * pi * me));
      vth_i_flux = sqrt(Ti_J / (2 * pi * mi));

  % Define electron current function (depends on phi_d sign)
  Ie = @(phi) -pi*rd^2 * e * ne * ve .* exp(e*phi/(kB*Te_J));

  % Ion current function
  Ii = @(phi) sum( pi*rd^2 .* ni .* e .* Z .* vi .* (1 - 2*Z*e*phi./(rd*mi.*vi.^2)) .* ...
                    (phi < 0) ...  % valid for negative phi_d
                   + pi*rd^2 .* ni .* e .* Z .* vi .* exp(-2*Z*e*phi./(rd*mi.*vi.^2)) .* (phi > 0) );

  % Net current balance
  I_net = @(phi) Ie(phi) + Ii(phi);

  % Solve for root: I_net(phi_d) = 0
  phi_guess = -Te; % good initial guess (negative few Te)
  phi_d = fzero(I_net, phi_guess);

end
end

