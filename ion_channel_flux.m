function [J, J_rhs, J_factor] = ion_channel_flux(c, phi, G, channel, ion_idx, t)
%J = ion_channel_flux(c, phi, G, channel, ion_idx, t)

% Compute Nernst equilibrium potential
E = (G.kB*G.T/(G.z(channel.ion(ion_idx))*G.e))*log(c(channel.m_e,channel.ion(ion_idx))./c(channel.m_i,channel.ion(ion_idx)));

% Compute transmembrane potential
if isfield(G, 'use_diff') && G.use_diff
    V = G.current_V;
else
    V = phi(channel.m_i)-phi(channel.m_e);
end

% Compute channel open probability
if channel.num_states > 0
    O = 1 - sum(channel.states,2);
else
    O = channel.open(t)*ones(size(channel.m_e));
end

% Compute flux (factor 10^9 is used to achive the required unit of mMnm/ms)
J = (1e9/(G.z(channel.ion(ion_idx))*G.F*channel.A)).*channel.g(ion_idx).*channel.nc.*(V-E).*O;

% Compute flux conductance
J_factor = (1e9/(G.z(channel.ion(ion_idx))*G.F*channel.A)).*channel.g(ion_idx).*channel.nc.*O;

% Compute rhs part of the flux
J_rhs = -E.*(1e9/(G.z(channel.ion(ion_idx))*G.F*channel.A)).*channel.g(ion_idx).*channel.nc.*O;



end