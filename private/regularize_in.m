function [in_removes, out_removes] = regularize_in(int_order, ext_order, S, ismag, extended_remove)

if nargin<4 || isempty(ismag)
  ismag = true(size(S,1),1);
end
if nargin<5 || isempty(extended_remove)
  extended_remove = [];
end

n_in  = (int_order+2).*int_order;
n_out = (ext_order+2).*ext_order;

[a_lm_sq, rho_i]  = compute_sphere_activation_in(0:int_order);
[degrees, orders] = get_degrees_orders(int_order);
a_lm_sq           = a_lm_sq(degrees);

I_tots     = zeros(1, n_in);
in_keepers = 1:n_in; 
out_removes = regularize_out(int_order, ext_order, S, ismag, extended_remove);
out_keepers = (n_in+1):size(S, 2);

normS = sqrt(sum(S.^2));
S     = S*diag(1./normS);
eigv  = zeros(n_in,2);

remove_order = zeros(1,0);

noise_lev = 5e-13.^2; % T/m, is this for grads?
for ii = 1:n_in
  this_S  = S(:, [in_keepers out_keepers]);
  [u,s,v] = svd(this_S, 'econ');
  diags   = diag(s);
  eigv(ii,:) = [diags(1) diags(end)];
  v = v(:,1:numel(in_keepers))';
  v = diag(1./normS(in_keepers))*v;
  eta_lm_sq = v*diag(1./diags)*u';
  eta_lm_sq = sum(eta_lm_sq.^2, 2) .* noise_lev;

  % this is some mysterious scaling, that is done in MNE-python to match
  % the elekta implementation
  eta_lm_sq(orders(in_keepers)==0) = eta_lm_sq(orders(in_keepers)==0).*2;
  eta_lm_sq = eta_lm_sq.*0.0025;

  snr = a_lm_sq(in_keepers)./eta_lm_sq';
  
  I_tots(ii) = 0.5 * sum(log(snr+1)./log(2));

  [minsnr, m]  = min(snr);
  remove_order = cat(2, remove_order, in_keepers(m));
  in_keepers   = setdiff(in_keepers, remove_order);
  
  maxI = max(I_tots);
  if ii > 10 && all(I_tots(ii-1:ii) < 0.95 * maxI)
    break
  end
end

if n_in>0
  lim_idx = find(I_tots >= 0.98.*maxI, 1, 'first');
  in_removes = remove_order(1:lim_idx);
else
  in_removes = [];
end

% def _regularize_in(int_order, ext_order, S_deco
% mp, mag_or_fine,
%                    extended_remove):
%     """Regularize basis set using idealized SNR measure."""
%     n_in, n_out = _get_n_moments([int_order, ext_order])
% 
%     # The "signal" terms depend only on the inner expansion order
%     # (i.e., not sensor geometry or head position / expansion origin)
%     a_lm_sq, rho_i = _compute_sphere_activation_in(
%         np.arange(int_order + 1))
%     degrees, orders = _get_degrees_orders(int_order)
%     a_lm_sq = a_lm_sq[degrees]
% 
%     I_tots = np.zeros(n_in)  # we might not traverse all, so use np.zeros
%     in_keepers = list(range(n_in))
%     out_removes = _regularize_out(int_order, ext_order, mag_or_fine,
%                                   extended_remove)
%     out_keepers = list(np.setdiff1d(np.arange(n_in, S_decomp.shape[1]),
%                                     out_removes))
%     remove_order = []
%     S_decomp = S_decomp.copy()
%     use_norm = np.sqrt(np.sum(S_decomp * S_decomp, axis=0))
%     S_decomp /= use_norm
%     eigs = np.zeros((n_in, 2))
% 
%     noise_lev = 5e-13  # noise level in T/m
%     noise_lev *= noise_lev  # effectively what would happen by earlier multiply
%     for ii in range(n_in):
%         this_S = S_decomp.take(in_keepers + out_keepers, axis=1)
%         u, s, v = _safe_svd(this_S, full_matrices=False, **check_disable)
%         del this_S
%         eigs[ii] = s[[0, -1]]
%         v = v.T[:len(in_keepers)]
%         v /= use_norm[in_keepers][:, np.newaxis]
%         eta_lm_sq = np.dot(v * 1. / s, u.T)
%         del u, s, v
%         eta_lm_sq *= eta_lm_sq
%         eta_lm_sq = eta_lm_sq.sum(axis=1)
%         eta_lm_sq *= noise_lev
% 
%         # Mysterious scale factors to match MF, likely due to differences
%         # in the basis normalizations...
%         eta_lm_sq[orders[in_keepers] == 0] *= 2
%         eta_lm_sq *= 0.0025
%         snr = a_lm_sq[in_keepers] / eta_lm_sq
%         I_tots[ii] = 0.5 * np.log2(snr + 1.).sum()
%         remove_order.append(in_keepers[np.argmin(snr)])
%         in_keepers.pop(in_keepers.index(remove_order[-1]))
%         # heuristic to quit if we're past the peak to save cycles
%         if ii > 10 and (I_tots[ii - 1:ii + 1] < 0.95 * I_tots.max()).all():
%             break
%         # if plot and ii == 0:
%         #     axs[0].semilogy(snr[plot_ord[in_keepers]], color='k')
%     if n_in:
%         max_info = np.max(I_tots)
%         lim_idx = np.where(I_tots >= 0.98 * max_info)[0][0]
%         in_removes = remove_order[:lim_idx]
%         for ii, ri in enumerate(in_removes):
%             eig = eigs[ii]
%             logger.debug(
%                 f'            Condition {eig[0]:0.3f} / {eig[1]:0.3f} = '
%                 f'{eig[0] / eig[1]:03.1f}, Removing in component '
%                 f'{ri}: l={degrees[ri]}, m={orders[ri]:+0.0f}'
%             )
%         logger.debug(
%             f'        Resulting information: {I_tots[lim_idx]:0.1f} '
%             f'bits/sample ({100 * I_tots[lim_idx] / max_info:0.1f}% of peak '
%             f'{max_info:0.1f})'
%         )
%     else:
%         in_removes = remove_order[:0]
%     return in_removes, out_removes
% 
% 


function [a_power, rho_i] = compute_sphere_activation_in(degrees)

% def _compute_sphere_activation_in(degrees):
%     """Compute the "in" power from random currents in a sphere.
%     Parameters
%     ----------
%     degrees : ndarray
%         The degrees to evaluate.
%     Returns
%     -------
%     a_power : ndarray
%         The a_lm associated for the associated degrees (see
%         :footcite:`KnuutilaEtAl1993`).
%     rho_i : float
%         The current density.
%     References
%     ----------
%     .. footbibliography::
%     """
%     r_in = 0.080  # radius of the randomly-activated sphere
% 
%     # set the observation point r=r_s, az=el=0, so we can just look at m=0 term
%     # compute the resulting current density rho_i
% 
%     # This is the "surface" version of the equation:
%     # b_r_in = 100e-15  # fixed radial field amplitude at distance r_s = 100 fT
%     # r_s = 0.13  # 5 cm from the surface
%     # rho_degrees = np.arange(1, 100)
%     # in_sum = (rho_degrees * (rho_degrees + 1.) /
%     #           ((2. * rho_degrees + 1.)) *
%     #           (r_in / r_s) ** (2 * rho_degrees + 2)).sum() * 4. * np.pi
%     # rho_i = b_r_in * 1e7 / np.sqrt(in_sum)
%     # rho_i = 5.21334885574e-07  # value for r_s = 0.125
%     rho_i = 5.91107375632e-07  # deterministic from above, so just store it
%     a_power = _sq(rho_i) * (degrees * r_in ** (2 * degrees + 4) /
%                             (_sq(2. * degrees + 1.) *
%                             (degrees + 1.)))
%     return a_power, rho_i

b_r_in = 100e-15;
r_in   = 0.08; % in m
r_s    = 0.05 + r_in;
rho_degrees = 1:100;
in_sum = sum( ( rho_degrees.*(rho_degrees+1)./(2.*rho_degrees+1) ) .* (r_in/r_s) .^ (2.*rho_degrees+2)) .* 4 .* pi;
rho_i  = b_r_in * 1e7 / sqrt(in_sum);
%rho_i = 5.91107375632e-7; % value for r_s at 0.13?
a_power = (rho_i^2) .* (degrees .* (r_in.^(2.*degrees+4)) ./ ((degrees+1).*(2.*degrees+1).^2));

