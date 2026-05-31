function [aR, bT] = nearfield_array_response(params, pl)
% =========================================================================
% nearfield_array_response.m
% =========================================================================
% * Mathematical Background:
%   Computes the exact near-field UPA array response vectors aR and bT
%   for the l-th scattering path, using the spherical wavefront model.
%   Eq. (12) from paper:
%     a_R(nR) = exp(-j*2*pi/lambda * (d^nR_R,l - d^c_R,l))
%     b_T(nT) = exp(-j*2*pi/lambda * (d^nT_T,l - d^c_T,l))
%
%   The UPA antenna index is:
%     nR = (NRy*NRz+1)/2 + nz_R*NRy + ny_R
%   with ny_R in [-(NRy-1)/2 ... (NRy-1)/2]
%        nz_R in [-(NRz-1)/2 ... (NRz-1)/2]
%
% * Inputs:
%   params - struct with fields: NRy, NRz, NTy, NTz, d, lambda, pT, pR
%   pl     - 3x1 position of l-th scattering point [m]
%
% * Outputs:
%   aR - NR x 1 array response at UT (receiver)
%   bT - NT x 1 array response at BS (transmitter)
%
% * MATLAB Implementation:
%   Uses meshgrid to create antenna positions and computes distances
%   using vectorized Euclidean norm.
%
% * Complexity Analysis:
%   O(NR + NT) (fully vectorized)
% =========================================================================

    lambda = params.lambda;
    d      = params.d;
    pT     = params.pT;   % 3x1 BS position
    pR     = params.pR;   % 3x1 UT position (may be true or estimated)
    NRy    = params.NRy;  NRz = params.NRz;
    NTy    = params.NTy;  NTz = params.NTz;

    % --- Receiver (UT) array response ---
    % UPA on yoz plane; antenna indices ny_R, nz_R centered at 0
    ny_R_vec = (-(NRy-1)/2 : (NRy-1)/2);   % 1 x NRy
    nz_R_vec = (-(NRz-1)/2 : (NRz-1)/2);   % 1 x NRz

    % Grid all (ny,nz) pairs: nR index = (NRy*NRz+1)/2 + nz*NRy + ny
    % Build NR x 2 matrix of [ny, nz] pairs (row-major: nz outer, ny inner)
    [NY_R, NZ_R] = meshgrid(ny_R_vec, nz_R_vec);   % NRz x NRy each
    ny_R = NY_R(:);   % NR x 1
    nz_R = NZ_R(:);   % NR x 1

    % Antenna positions relative to UT center (on yoz plane)
    % Antenna at (0, ny*d, nz*d) + pR  in world coords
    % pR = [xR; yR; zR], UT array center
    ant_pos_R = pR + [zeros(1,length(ny_R)); ny_R'*d; nz_R'*d];  % 3 x NR

    % Distances from SP to each UT antenna
    diff_R = ant_pos_R - pl;   % 3 x NR
    d_nR   = sqrt(sum(diff_R.^2, 1))';  % NR x 1

    % Distance from SP to UT array center
    d_c_R  = norm(pl - pR);

    % Array response (Eq. 12)
    aR = exp(-1j * 2*pi/lambda * (d_nR - d_c_R));   % NR x 1

    % --- Transmitter (BS) array response ---
    [NY_T, NZ_T] = meshgrid((-(NTy-1)/2:(NTy-1)/2), (-(NTz-1)/2:(NTz-1)/2));
    ny_T = NY_T(:);
    nz_T = NZ_T(:);

    ant_pos_T = pT + [zeros(1,length(ny_T)); ny_T'*d; nz_T'*d];  % 3 x NT
    diff_T = ant_pos_T - pl;   % 3 x NT
    d_nT   = sqrt(sum(diff_T.^2, 1))';  % NT x 1
    d_c_T  = norm(pl - pT);

    bT = exp(-1j * 2*pi/lambda * (d_nT - d_c_T));   % NT x 1
end
