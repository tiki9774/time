%% KF 

x0 = [0 0 0]';
t = 1:100;

xhat_plus = zeros(3,length(t));
pplus = zeros(3,3,length(t));


I = eye(3);

dt = 0.1;

F = [1 dt; 0 1];
G = [0.5 * dt^2; dt];



%% Propagate

for ii  = 1:length(t)
    % Prediction Step
    xhat_minus = F*xhat_plus(:,ii);
    pminus = F * Pplus * F' + Q;
    K = Pminus * H'*inv(H*pminus*H' + R);

    % Correction Step
    xhat_plus(:,ii) = xhat_minus + K * (Y - H * xhat_minus);
    pplus(:,:,ii) = (I - K*H)*pminus;


end