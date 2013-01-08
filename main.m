clear
clc
close all

Fd = 44.2e6 * 1e-3; %MHz
Td = 1/Fd;

pic_Init = 0; 
pic_Ext = 0;
pic_Est = 1;

stdn = 8e-3;
dstdn2 = 2*stdn^2;
invdstdn2 = 1/2*stdn^2;

T = 0.02;
% F = [1 T 0;
%        0 1 T;
%        0 0 1];
F = [1 T;
       0 1];
G = [0; 1];
nx = 2;

Sksi = 1e-3;
Dksi = Sksi/T;
Dxx = Dksi; %G'*Dksi*G;
detDxx = det(Dxx);
Dxxm1 = Dxx^-1;
Da = Dksi * (1/T) * (1/(2*pi*1.6e9)*3e8)^2; % D[a], (m/s^2)^2

K = 400;
L = round(T/Td);
M = 2;

qcno_dB = 45;
qcno = 10.^(qcno_dB/10);
A = 2*stdn*sqrt(qcno*Td);
invstdn2 = 1/stdn^2;

lbase = 2; lambda0 = 0.19;

% Initial state

% Initial PD of PhD
% H_psi1 = 2 * lbase/lambda0 * 2*pi;
H_psi1 = pi;
D_extr_psi1 = (H_psi1^2)/12;
psi1s = (rand(1,1) - 0.5)*H_psi1; % psi_0

% Initial PD of PhD rate
H_psi2 = 2*pi*lbase/lambda0 / 20;
D_extr_psi2 = (H_psi2^2)/12; 
psi2s = (rand(1,1) - 0.5)*H_psi2; % (diff psi)_0

% Initial PD of rate of PhD rate
H_psi3 = 6*sqrt(12 * Dksi);
D_extr_psi3 = (H_psi3^2)/12; 
psi3s = (rand(1,1) - 0.5)*H_psi3; % (diff diff psi)_0

Xs = [psi1s; psi2s];

Npsi = [50; 40; 10]; % Number of points by axes
maxpsi = 6*[sqrt(D_extr_psi1); 80*sqrt(D_extr_psi2); sqrt(D_extr_psi3)]; % Argument's area
minpsi = -maxpsi;
dpsi = (maxpsi-minpsi) ./ Npsi; % differential step
psi1 = minpsi(1):dpsi(1):maxpsi(1); % Argument's vectors
psi2 = minpsi(2):dpsi(2):maxpsi(2);
psi3 = minpsi(3):dpsi(3):maxpsi(3);

pest_psi1 = 1/H_psi1 .* ( (psi1 >= (-H_psi1/2))&(psi1 <= (H_psi1/2)) ); % Initial PDs
pest_psi2 = 1/H_psi2 .* ( (psi2 >= (-H_psi2/2))&(psi2 <= (H_psi2/2)) );
pest_psi3 = 1/H_psi3 .* ( (psi3 >= (-H_psi3/2))&(psi3 <= (H_psi3/2)) );

pest = pest_psi1'*pest_psi2; % Common PD

if pic_Init
    hF = figure(1);
    surf(psi2/2/pi, psi1/2/pi, pest)
    title('Initial');
    ylabel('\psi, cycles');
    xlabel('\psi'', Hz');
    % zlabel('p(x_k|Y_k)')
    zlabel('p(x_k)')
    saveas(hF, sprintf('%04.0f.png', 0), 'png')
    drawnow
    pause(0.1)
end

% Cycles
max_arg_exp = 100; max_arg_exp_stdn = max_arg_exp*stdn^2;
phi0_int = 0:(2*pi/10):2*pi; % \phi_1,0 axe for integration
dphi = phi0_int(2) - phi0_int(1); % differential step of \phi_1_0
for k = 1:K

    pextr = zeros(length(psi1), length(psi2));
    for j1 = 1:length(psi1)
        for i1 = 1:length(psi1)
            for i2 = 1:length(psi2)
                delta_to_psi1 = -(psi1(j1) - (psi1(i1) + psi2(i2)*T)); % = (psi1(i1) - (psi1(j1) - psi2(i2)*T))
                if (abs(delta_to_psi1)< dpsi(1)/2) || (delta_to_psi1 == dpsi(1)/2)  % psi1(j1) - psi2(i2)*T  \approx psi1(i1)
%                     t = psi1(j1) - psi2(i2)*T;
%                     p(t) =  p(x) + (p(x+1) - p(x))/dx*(t-x)                    
                    if delta_to_psi1 < 0 % t < x
                        if i1 > 1 % x > xmin
                            pest_t = (pest(i1,i2)-pest(i1-1,i2))/dpsi(1) * delta_to_psi1 + pest(i1,i2);
                        elseif i1 == 1 % x = xmin
                            pest_t = (pest(i1+1,i2)-pest(i1,i2))/dpsi(1) * delta_to_psi1 + pest(i1,i2); % df/dx approx const on ends
                        end
                    else % t > x
                        if i1 == length(psi1)
                            pest_t = (pest(i1,i2)-pest(i1-1,i2))/dpsi(1) * delta_to_psi1 + pest(i1,i2); % df/dx approx const on ends
                        else
                            pest_t = (pest(i1+1,i2)-pest(i1,i2))/dpsi(1) * delta_to_psi1 + pest(i1,i2); 
                        end
                    end
                    for j2 = 1:length(psi2)
                        pextr(j1, j2) = pextr(j1, j2) + pest_t*exp( -0.5*(psi2(j2) -  psi2(i2))^2*Dxxm1 ); %  * 1 / sqrt( (2*pi)^nx *detDxx) * dpsi(2)
                    end
                end
            end
        end
    end

    pextr = pextr * 1 / sqrt( (2*pi)^nx *detDxx) * dpsi(2); % multiplicators from cycle
    
    if pic_Ext
        hF = figure(2);
        surf(psi2/2/pi, psi1/2/pi, pextr);
        title(['Extrapolation, t = ' sprintf('%.3f s', k*T)])
        %     title(['Progress of PD, t = ' sprintf('%.3f s', k*T)])
        ylabel('\psi, cycles');
        xlabel('\psi'', Hz');
        zlabel('p(x_{k+1}|Y_k)')
        %     zlabel('p(x_k)')
        drawnow
        %     saveas(hF, sprintf('%04.0f.png', k), 'png')
        pause(0.1);
    end

    Xs = 0*F*Xs + 0*G*randn(1,1)*sqrt(Dksi);
    PW = 2*pi*Fd/3.3712*(0:(L-1))*Td; % Intermediate freq phase
    
    phi0 =  0*rand(1,1)*2*pi;
    phi0_int = phi0;
    S2 = A*cos(Xs(1) + Xs(2)*(0:(L-1))*Td +phi0 + PW);
    S1 = A*cos(phi0 + PW);
    y2 = S2 + 0*stdn*randn(1,L);
    y1 = S1 + 0*stdn*randn(1,L);

    j = 0;
    fprintf('Measurements are captured, my general, at t = %.3f s\n', k*T);
    Likehood = ones(length(psi1), length(psi2));
    lnLikehood = zeros(length(psi1), length(psi2));
    for j2 = 1:length(psi2)
        Spsi2base = exp(1i*(psi2(j2)*(0:(L-1))*Td + PW));
        for ip  = 1:length(phi0_int)
            Sphibase = Spsi2base*exp(1i*phi0_int(ip));
            Sop1 = A*cos(phi0_int(ip) + PW);
            for j1 = 1:length(psi1)
                Sop2 =A*real(Sphibase*exp(1i*psi1(j1)));
                Lstep = L;
                l_ok = 0;
                while l_ok < L
                    
                    step_nook = 1;
                    while step_nook
                        Lstep = min([Lstep, L-l_ok]); 
                        argarr = (l_ok+1):(l_ok+Lstep);
                        argexp =  (Sop1(argarr))*(y1(argarr))'+(Sop2(argarr))*(y2(argarr))' - 0.5*((Sop1(argarr))*(Sop1(argarr))' + (Sop2(argarr))*(Sop2(argarr))');
                        if abs(argexp)<=max_arg_exp_stdn
                            step_nook = 0;
                        else
                            if Lstep > 1
                                Lstep = Lstep - 1;
                            else
                                fprintf('Bad dynamic range \n');
                                step_nook = 0;
                            end
                        end
                    end
%                     Likehood(j1, j2) = Likehood(j1, j2) + exp( (  -0.5*(y1*y1'+y2*y2') + Sop1*y1'+Sop2*y2' - 0.5*(Sop1*Sop1' + Sop2*Sop2')  ) /  stdn^2  ); % / (2*pi*stdn^2) * dphi / (2*pi);
                    expmult = exp( argexp * invstdn2 );
                    j = j+1;
                    eargexp_j(j) = argexp * invstdn2;
                    expmult_j(j) = expmult;
                    Likehood(j1, j2) = Likehood(j1, j2) * ( expmult ); 
                    lnLikehood(j1, j2) = lnLikehood(j1, j2) + argexp * invstdn2;
                    if expmult == 0
                        fprintf('All is lost!\n');
                    end
                    
                    l_ok = l_ok + Lstep;
                    if 3*abs(argexp) < max_arg_exp_stdn
                        Lstep = (fix(max_arg_exp_stdn/abs(argexp)) - 2) + Lstep;
                    end

                end
            end
        end
    
    if pic_Est
        figure(3)
        surf(psi2/2/pi, psi1/2/pi, lnLikehood)
        title('Ln of Likehood function');
        ylabel('\psi, cycles');
        xlabel('\psi'', Hz');
        zlabel('ln p(y|x)')
        drawnow
        pause(0.1)
    end
    return
    pest = pextr;
    pest = pest/(sum(sum(pest))*dpsi(1)*dpsi(2));
 
end



%             Phi = psi1(j1) + psi2(j2)*(0:(L-1))*Td;
%             I2 = y2 * (cos(Phi))';
%             Q2 = y2 * (sin(Phi))';
%             
%             I1 = y1 * ones(L,1);
%             Q1 = y1 * zeros(L,1);
%             X = sqrt( (I1+I2).^2 + (Q1+Q2).^2 );
%             Likehood(j1, j2) = ...
%             besseli(0, A/stdn^2*X)*... % Infinity
%                 exp(0);


%     pextr = nan(length(psi1), length(psi2));
%     for j1 = 1:length(psi1)
%         for j2 = 1:length(psi2)
%             pextr(j1, j2) = 0;
%             xk = [psi1(j1); psi2(j2)];
%             for i1 = 1:length(psi1)
%                 xkm1 = [psi1(i1); 0];
%                 for i2 = 1:length(psi2)
%                     xkm1(2) = psi2(i2);
%                     pxx_norm = exp(-0.5*(xk-F*xkm1)'*Dxxm1*(xk-F*xkm1)); % * 1 / sqrt( (2*pi)^nx *detDxx)
%                     pextr(j1, j2) = pextr(j1, j2) + pest(i1, i2)*pxx_norm*dpsi(1)*dpsi(2);
%                 end
%             end
%         end
%     end
