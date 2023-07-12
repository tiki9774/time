%Lec22_1Drobotstatefilter.m
clc,clear,close all

%%fix random number seed to get repeatable results (useful for debugging):
rng(100) %%comment out to get different random results each time


m0 = [1,1e-12,0];
P0 = eye(3);

tau = 1;
%%CT Model spec
A = [0 1 0;
     0 0 1;
     0 0 0];
B = [0;1];
Gamma = [0;0;1]; %look into vanloan method, and psd
C = [1 0];
Qtilde = eye(3); %AWG process noise intensity: fix later,  Q= psd
Rtilde = 50e-9; %AWG measurement noise intensity: seconds, double check
% CTsys = ss(A,B,C,0);  %%

%%DT model conversion
deltaT = 1; %second, maybe 1000?
tau = deltaT;
%%Use c2d to cheat/convert parts of model

% DTsys = c2d(CTsys,deltaT,'zoh'); %%
F = [1 tau tau^2/2; 0 1 tau; 0 0 1]; %DT space matrix
G = eye(3)
H = [-1 0 0 1 0 0];
H = [1 0 0];
%%Use Van Loan method to find covariance of process noise
% M = deltaT*[-A, Gamma*Qtilde*Gamma';
%             zeros(size(A)), A']
% matrexp_M = expm(M)        
% invFQblock = matrexp_M(1:2,3:4);
%use process noise from paper
q1 = 1 %look at this more
q2 = 1
q3 = 1
Q = q1*[tau 0 0; 0 0 0; 0 0 0] + q2*[tau^3/3 tau^2/2 0; tau^2/2 tau 0; 0 0 0] + q3 * [tau^5/20 tau^4/8 tau^2/6; tau^4/8 tau^3/3 tau^2/2; tau^2/6 tau^2/2 tau];
%%Find the measurement noise covariance
R = Rtilde / deltaT

tvec = 0:deltaT:12;
u = 2*cos(0.75*tvec); %acceleration input

%% 1. Pure forward prediction
%Show recursion for prediction of mean and covariance
mk = m0;
Pk = P0;
mk_pred_hist = zeros(3,length(tvec));
Pk_pred_hist = zeros(3,3,length(tvec));
for k=1:length(tvec)
  
    mkp1 = F*mk';% + G*u(k);
    Pkp1 = F*Pk*F' + Q;
    
    %%%pxyk = reshape( mvnpdf(XY,mkp1',Pkp1) , size(X)); %%update mvnpdf
    
    mk = mkp1';
    mk_pred_hist(:,k)=mkp1;
    Pk = Pkp1;
    Pk_pred_hist(:,:,k)=Pkp1;
end

figure(3), 
subtitle(['1D Robot State Prediction Results for \Delta t = ',num2str(deltaT)])
subplot(211)
plot(tvec,mk_pred_hist(1,:),'b','LineWidth',3), hold on
plot(tvec,mk_pred_hist(1,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'b--','LineWidth',2)
plot(tvec,mk_pred_hist(1,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'b--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('\xi_k (m)','FontSize',18),
legend('mean \xi(k)','2\sigma bounds')

subplot(212)
plot(tvec,mk_pred_hist(2,:),'r','LineWidth',3), hold on
plot(tvec,mk_pred_hist(2,:)+2*sqrt(squeeze(Pk_pred_hist(2,2,:))'),'r--','LineWidth',2)
plot(tvec,mk_pred_hist(2,:)-2*sqrt(squeeze(Pk_pred_hist(2,2,:))'),'r--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('dot \xi_k (m/s)','FontSize',18),
legend('mean dot \xi(k)','2\sigma bounds')

%% 2. Simulate ground truth trajectory and measurements from DT LTI model
xk_truehist = zeros(3,length(tvec));
ykhist = zeros(1,length(tvec));
xk = [0 0 0]'%mvnrnd(m0,P0)'; %sample initial robot state
for k=1:length(tvec)
  
    %%simulate process noise and add to actual state
%    wk = mvnrnd(zeros(1,2),Q)';
    a = Q(1,1)*randn;
    b = Q(2,2)*randn;
    c = Q(3,3)*randn;
    wk = [a;b;c];
    xkp1 = F*xk + wk;     
    
    %%simulate measurement noise and add to sensor data
    %vkp1 = mvnrnd(zeros(1,1),R)';
    vkp1 = R(1,1)*randn;
    % b = R(2,2)*randn;
    % c = R(3,3)*randn;
    % vkp1 = [a 0 0; 0 b 0;0 0 c];
    H_tmp = [1 0 0];
    ykp1 = H_tmp*xkp1 + vkp1;
    
    %%store and iterate
    xk_truehist(:,k) = xkp1;
    ykhist(:,k) = ykp1; 
    xk = xkp1;
end


%% 3. Kalman Filter

%%Initialize
mk = m0;
Pk = P0;
mk_filt_hist = zeros(3,length(tvec));
Pk_filt_hist = zeros(3,3,length(tvec));

%%Specify KF's Q and R matrices: 
%%uncomment the appropriate lines below to try out different cases:

%%Case 1: KF uses the correct/true Q and R values for DT system:
Qkf = Q;
Rkf = R;

%%Case 2a: KF assumes smaller Q for DT system, but correct R:
% Qkf = Q*0.001;
% Rkf = R;

%%Case 2b: KF assumes smaller Q for DT system, but correct R:
% Qkf = Q*100;
% Rkf = R;

%%Case 3: KF *thinks* sensors are not as good, but KF uses correct Q
% Qkf = Q;
% Rkf = 10*R; 

%%Apply Kalman filter updates
for k=1:length(tvec)
  
    %%Perform prediction step
    mkp1_minus = F*mk' ;%+ G*u(k);
    Pkp1_minus = F*Pk*F' + Qkf;
         
    %%Compute Kalman gain
    Kkp1 = Pkp1_minus*H'/(H*Pkp1_minus*H' + Rkf);
    
    %%Perform measurement update step
    ykp1_report = ykhist(:,k); %pull report of actual data from sensor
    ykp1_pred = H*mkp1_minus; %predicted measurement
    innov_kp1 = ykp1_report - ykp1_pred; %compute innovation
    mkp1_plus = mkp1_minus + Kkp1*innov_kp1; %compute update to state mean
    Pkp1_plus = (eye(3) - Kkp1*H)*Pkp1_minus; %compute update to covar
    
    %%store results and cycle for next iteration
    mk = mkp1_plus'; 
    mk_filt_hist(:,k) = mkp1_plus;
    Pk = Pkp1_plus;
    Pk_filt_hist(:,:,k)=Pkp1_plus;
end

figure(4), 
subtitle(['1D Robot State Filtering Results for \Delta t = ',num2str(deltaT)])
subplot(311)
plot(tvec,mk_filt_hist(1,:),'b','LineWidth',3), hold on
plot(tvec,mk_filt_hist(1,:)+2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'b--','LineWidth',2)
plot(tvec,mk_filt_hist(1,:)-2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'b--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('\xi_k (m)','FontSize',18),
plot(tvec,xk_truehist(1,:),'k-o','LineWidth',1.5)
%%%compare to prediction results
plot(tvec,mk_pred_hist(1,:),'c','LineWidth',1.5), hold on
plot(tvec,mk_pred_hist(1,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'c--','LineWidth',1.5)
plot(tvec,mk_pred_hist(1,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'c--','LineWidth',1.5)
ylim([-10 10])
%%plot raw sensor data
plot(tvec,ykhist(1,:),'kx','LineWidth',1.5,'MarkerSize',7)
legend('filter mean \xi(k)','filter 2\sigma bounds','','truth',...
        'mean pred \xi(k)','pred 2\sigma bounds','','y_k sensor data')

subplot(312)
plot(tvec,mk_filt_hist(2,:),'r','LineWidth',3), hold on
plot(tvec,mk_filt_hist(2,:)+2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'r--','LineWidth',2)
plot(tvec,mk_filt_hist(2,:)-2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'r--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('dot \xi_k (m/s)','FontSize',18),
legend('mean dot \xi(k)','2\sigma bounds')
plot(tvec,xk_truehist(2,:),'k-o','LineWidth',1.5)
%%%compare to prediction results
plot(tvec,mk_pred_hist(2,:),'m','LineWidth',1.5), hold on
plot(tvec,mk_pred_hist(2,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
plot(tvec,mk_pred_hist(2,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
legend('filter mean dot \xi(k)','filter 2\sigma bounds','','truth',...
        'mean pred dot \xi(k)','pred 2\sigma bounds','')
ylim([-10 10])

subplot(313)
plot(tvec,mk_filt_hist(3,:),'r','LineWidth',3), hold on
plot(tvec,mk_filt_hist(3,:)+2*sqrt(squeeze(Pk_filt_hist(3,2,:))'),'r--','LineWidth',2)
plot(tvec,mk_filt_hist(3,:)-2*sqrt(squeeze(Pk_filt_hist(3,2,:))'),'r--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('dot \xi_k (m/s)','FontSize',18),
legend('mean dot \xi(k)','2\sigma bounds')
plot(tvec,xk_truehist(2,:),'k-o','LineWidth',1.5)
%%%compare to prediction results
plot(tvec,mk_pred_hist(3,:),'m','LineWidth',1.5), hold on
plot(tvec,mk_pred_hist(3,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
plot(tvec,mk_pred_hist(3,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
legend('filter mean dot \xi(k)','filter 2\sigma bounds','','truth',...
        'mean pred dot \xi(k)','pred 2\sigma bounds','')
ylim([-10 10])

%%Plot state estimation errors versus time
figure(5), 
subtitle(['1D Robot State Estimation Errors for \Delta t = ',num2str(deltaT)])
subplot(211)
plot(tvec,xk_truehist(1,:)-mk_filt_hist(1,:),'b','LineWidth',3), hold on
plot(tvec,xk_truehist(1,:)-mk_filt_hist(1,:)+2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'b--','LineWidth',2)
plot(tvec,xk_truehist(1,:)-mk_filt_hist(1,:)-2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'b--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('\xi_k error (m)','FontSize',18),
%%%compare to prediction results
plot(tvec,xk_truehist(1,:)-mk_pred_hist(1,:),'c','LineWidth',1.5), hold on
plot(tvec,xk_truehist(1,:)-mk_pred_hist(1,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'c--','LineWidth',1.5)
plot(tvec,xk_truehist(1,:)-mk_pred_hist(1,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'c--','LineWidth',1.5)
plot(tvec,zeros(length(xk_truehist(1,:))),'k-o')
legend('filter error \xi(k)','filter 2\sigma bounds','',...
        'pred error \xi(k)','pred 2\sigma bounds','','truth')
ylim([-12 12])

subplot(212)
plot(tvec,xk_truehist(2,:)-mk_filt_hist(2,:),'r','LineWidth',3), hold on
plot(tvec,xk_truehist(2,:)-mk_filt_hist(2,:)+2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'r--','LineWidth',2)
plot(tvec,xk_truehist(2,:)-mk_filt_hist(2,:)-2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'r--','LineWidth',2)
xlabel('time step, k','FontSize',18), ylabel('dot \xi_k error (m/s)','FontSize',18),
%%%compare to prediction results
plot(tvec,xk_truehist(2,:)-mk_pred_hist(2,:),'m','LineWidth',1.5), hold on
plot(tvec,xk_truehist(2,:)-mk_pred_hist(2,:)+2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
plot(tvec,xk_truehist(2,:)-mk_pred_hist(2,:)-2*sqrt(squeeze(Pk_pred_hist(1,1,:))'),'m--','LineWidth',1.5)
plot(tvec,zeros(length(xk_truehist(2,:))),'k-o')
legend('filter error dot \xi(k)','filter 2\sigma bounds','',...
        'pred error dot \xi(k)','pred 2\sigma bounds','','truth')
ylim([-12 12])

figure(6)
subplot(211)
plot(tvec,xk_truehist(1,:),'k-o','LineWidth',1.5), hold on
plot(tvec,ykhist(1,:),'rx','LineWidth',1.5,'MarkerSize',7)
legend('true pos','y meas')
xlabel('time step, k','FontSize',18), ylabel('\xi_k (m/s)','FontSize',18),
subplot(212)
plot(tvec,xk_truehist(2,:),'k-o','LineWidth',1.5), hold on
legend('true vel')
xlabel('time step, k','FontSize',18), ylabel('dot \xi_k (m/s)','FontSize',18),
