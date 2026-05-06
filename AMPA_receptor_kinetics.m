function ds = AMPA_receptor_kinetics(states, c)
%ds = AMPA_receptor_kinetics(states, c)
% AMPA receptor model from Jonas, P., Major, G., & Sakmann, B. (1993). 
% Quantal components of unitary EPSCs at the mossy fibre synapse on CA3 
% pyramidal cells of rat hippocampus. The Journal of physiology, 472(1), 
% 615-663.

% Extract glutamate concentration
glu = c(:, 5);

% Define parameters
k1p = 4.59;   % 1/(mM*ms)
k1m = 4.26;   % 1/ms
k2p = 28.4;   % 1/(mM*ms)
k2m = 3.26;   % 1/ms
k3p = 1.27;   % 1/(mM*ms)
k3m = 0.0457; % 1/ms
a = 4.24;     % 1/ms
b = 0.9;      % 1/ms
a1 = 2.89;    % 1/ms
b1 = 0.0392;  % 1/ms
a2 = 0.172;   % 1/ms
b2 = 7.27e-4; % 1/ms
a3 = 0.0177;  % 1/ms
b3 = 4.0e-3;  % 1/ms
a4 = 0.0168;  % 1/ms
b4 = 0.1904;  % 1/ms

% Load states
C0 = states(:,1);
C1 = states(:,2);
C2 = states(:,3);
C3 = states(:,4);
C4 = states(:,5);
C5 = states(:,6);
O = 1 - (C0 + C1 + C2 +C3 + C4 + C5);

% Set up temporal derivatives
ds = zeros(size(states));
ds(:,1) = -k1p*glu.*C0 + k1m*C1;                                   % [C0]
ds(:,2) = k1p*glu.*C0 - (k1m + k2p*glu + a1).*C1 + k2m*C2 + b1*C3; % [C1]
ds(:,3) = k2p*glu.*C1 - (k2m + a + a2)*C2 + b*O + b2*C4;           % [C2]
ds(:,4) = a1*C1 - (b1 + k3p*glu).*C3 + k3m*C4;                     % [C3]
ds(:,5) = k3p*glu.*C3 + a2*C2 - (k3m + b2 + a4)*C4 + b4*C5;        % [C4]
ds(:,6) = a3*O + a4*C4 - (b3 + b4)*C5;                             % [C5]

end