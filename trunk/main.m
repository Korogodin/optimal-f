clear
clc
close all

Fd = 44.2e6 * 1e-3; %MHz
Td = 1/Fd;

pic_Init = 1; 
pic_Ext = 1;
pic_Bess = 1;
pic_Est = 1;

stdn = 8e-2;
dstdn2 = 2*stdn^2;
invdstdn2 = 1/2*stdn^2;

dBtoNp = 8.685889638; % = 1 Np / 1 dB
T = 0.02;
% F = [1 T 0;
%        0 1 T;
%        0 0 1];
F = [1 T;
       0 1];
G = [0; 1];
nx = 2;

% Sksi = 1e-3;
Sksi = 1e-0;
Dksi = Sksi/T;
Dxx = Dksi; %G'*Dksi*G;
detDxx = det(Dxx);
Dxxm1 = Dxx^-1;
Da = Dksi * (1/T) * (1/(2*pi*1.6e9)*3e8)^2; % D[a], (m/s^2)^2

K = 200;
L = round(T/Td);
M = 2;

qcno_dB = 20;
qcno = 10.^(qcno_dB/10);
A = 2*stdn*sqrt(qcno*Td);
invstdn2 = 1/stdn^2;
if qcno_dB > 30
    lnmode = 1;
else
    lnmode = 0;
end

lbase = 2; lambda0 = 0.19;

% Initial state

% Initial PD of PhD
% H_psi1 = 2 * lbase/lambda0 * 2*pi;
H_psi1 = pi;
D_extr_psi1 = (H_psi1^2)/12;
psi1s = (rand(1,1) - 0.5)*H_psi1; % psi_0

% Initial PD of PhD rate
H_psi2 = 2*pi*lbase/lambda0 / 3 ;
D_extr_psi2 = (H_psi2^2)/12; 
psi2s = (rand(1,1) - 0.5)*H_psi2; % (diff psi)_0

% Initial PD of rate of PhD rate
H_psi3 = 6*sqrt(12 * Dksi);
D_extr_psi3 = (H_psi3^2)/12; 
psi3s = (rand(1,1) - 0.5)*H_psi3; % (diff diff psi)_0

Xs = [psi1s; psi2s];

Npsi = [100; 120; 10]; % Number of points by axes
maxpsi = 6*[2*sqrt(D_extr_psi1); 6*sqrt(D_extr_psi2); sqrt(D_extr_psi3)]; % Argument's area
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
    subplot(2,2,1);
    surf(psi2/2/pi, psi1/2/pi, pest)
    title('Initial');
    ylabel('\psi, cycles');
    xlabel('\psi'', Hz');
    zlabel('p(x_k)')
%     saveas(hF, sprintf('%04.0f.png', 0), 'png')
    drawnow
    pause(0.1)
end

% Cycles
tint = (0:(L-1))*Td;
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

%     pextr = pextr * 1 / sqrt( (2*pi)^nx *detDxx) * dpsi(2); % multiplicators from cycle
    pextr = pextr/(sum(sum(pextr))*dpsi(1)*dpsi(2));
    
    if pic_Ext
        hF = figure(1);
        subplot(2,2,2);
        surf(psi2/2/pi, psi1/2/pi, pextr);
        title(['Extrapolation, ' sprintf('k = %.0f, t = %.3f s', k, (k-1)*T)])
        ylabel('\psi, cycles');
        xlabel('\psi'', Hz');
        zlabel('p(x_{k}|Y_{k-1})')
        drawnow
        %     saveas(hF, sprintf('%04.0f.png', k), 'png')
        pause(0.1);
    end

%     Xs = F*Xs + G*randn(1,1)*sqrt(Dksi);
    cfreq = 0.3;
    phi_Xs = 2*pi*cfreq*k*T + pi/2;
    Xs(1) = psi1(end)/2*cos(phi_Xs);
    Xs(2) = -psi1(end)/2*sin(phi_Xs)*2*pi*cfreq;

    PW = 2*pi*Fd/3.3712*tint; % Intermediate freq phase
    argSop1 = PW';
    Sop1cos = cos(argSop1);
    Sop1sin = sin(argSop1);
    
    phi0 =  rand(1,1)*2*pi;
    S2 = A*cos(Xs(1) + Xs(2)*tint +phi0 + PW);
    S1 = A*cos(phi0 + PW);
    y2 = S2 + 1*stdn*randn(1,L);
    y1 = S1 + 1*stdn*randn(1,L);

    fprintf('Measurements are captured, my general, at t = %.3f s\n', k*T);
    
    calc_lnLike = 0;
    if calc_lnLike
        phi0_int = 0:(pi/4):pi;
        for ip  = 1:length(phi0_int)
            lnLikehood = zeros(length(psi1), length(psi2));
            for j2 = 1:length(psi2)
                Spsi2base = exp(1i*(PW + psi2(j2)*tint+ phi0_int(ip)));
                Sop1 = A*cos(PW + phi0_int(ip));
                for j1 = 1:length(psi1)
                    Sop2 =A*real(Spsi2base*exp(1i*psi1(j1)));
                    lnLikehood(j1, j2) = lnLikehood(j1, j2) +   Sop1*y1'+Sop2*y2' - 0.5*(Sop1*Sop1' + Sop2*Sop2')    ; % / (2*pi*stdn^2) * dphi / (2*pi);
                end
            end
            lnLikehood = lnLikehood * invstdn2  ;
            
%             lnmode = 0;
            figure(2)
            if lnmode 
                surf(psi2/2/pi, psi1/2/pi, dBtoNp*lnLikehood)
                if phi0_int(ip) == 0
                    title('ln L(x|\phi_1) for \phi_1 = 0');
                else
                    title(['ln L(x|\phi_1) for \phi_1 = \pi/' sprintf('%.0f', pi/phi0_int(ip)) ]);
                end
                zlabel('L(x|\mu), dB')
            else
                surf(psi2/2/pi, psi1/2/pi, exp(lnLikehood))
                if phi0_int(ip) == 0
                    title('L(x|\phi_1) for \phi_1 = 0');
                else
                    title(['L(x|\phi_1) for \phi_1 = \pi/' sprintf('%.0f', pi/phi0_int(ip)) ]);
                end
                zlabel('L(x|\mu)')
            end
            ylabel('\psi, cycles');
            xlabel('\psi'', Hz');
            drawnow
            pause(0.1)
%             saveas(hF, sprintf('%04.0f_%.0f.png', k, ip), 'png')
        end
    end
    
    argBessel = zeros(length(psi1), length(psi2));
    for j2 = 1:length(psi2)
        for j1 = 1:length(psi1)
            argSop2 = (PW + psi2(j2)*tint + psi1(j1))';
            Sop2cos = cos(argSop2);
            Sop2sin = sin(argSop2);

            I2 = y2 * Sop2cos;
            Q2 = y2 * Sop2sin;
            I1 = y1 * Sop1cos;
            Q1 = y1 * Sop1sin;
            argBessel(j1, j2) = sqrt( (I1+I2).^2 + (Q1+Q2).^2 );
        end
    end
    argBessel = A*argBessel*invstdn2;
    L_in_Np = argBessel - 0.5*log(argBessel);    
   
    if pic_Bess
%         lnmode = 0;
        figure(1)
        subplot(2,2,3);
        if lnmode
%             L_in_Np = argBessel - 0.5*(log(argBessel) + log(2*pi))
            surf(psi2/2/pi, psi1/2/pi, dBtoNp*L_in_Np) 
            title('Mean Likehood L(x) in dB scale');            
            zlabel('L(x), dB')
        else
            surf(psi2/2/pi, psi1/2/pi, besseli(0, argBessel))
            title('Mean Likehood L(x)');
            zlabel('L(x)')
        end

        ylabel('\psi, cycles');
        xlabel('\psi'', Hz');
        drawnow
        pause(0.1)
    end

%     return

    dB_extr = 120;
    minpextr = max(max(pextr))*exp(-dB_extr);
    fixed_pextr = (pextr > minpextr).*pextr  + (pextr <= minpextr)*minpextr;
    fixed_pextr = fixed_pextr / minpextr;
    
    pest = log(fixed_pextr) + L_in_Np;

    dB_est = 120;
    minpest = max(max(pest))-dB_est;
    fixed_pest = (pest > minpest).*pest  + (pest <= minpest)*minpest;
    fixed_pest = fixed_pest - minpest;
    
    pest = exp(fixed_pest).*(pest > minpest); % *0, if  arg = argmax - dB_est
    
    pesr = pextr;
    pest = pest/(sum(sum(pest))*dpsi(1)*dpsi(2));
 
    if pic_Est
        figure(1)
        subplot(2,2,4);
        surf(psi2/2/pi, psi1/2/pi, pest)
        title(['Aposteriori probability density and True value, ' sprintf('k = %.0f, t = %.3f s', k, k*T)]);
        zlabel('p(x_k|Y_k)')
        ylabel('\psi, cycles');
        xlabel('\psi'', Hz');
        hold on
        plot3(  [Xs(2)/2/pi, Xs(2)/2/pi], [Xs(1)/2/pi, Xs(1)/2/pi], [0, 1.2*max(max(pest))], 'r', ...
                  Xs(2)/2/pi, Xs(1)/2/pi, 1.2*max(max(pest)), 'r*');
        hold off        
        drawnow
        pause(0.1)        
    end
    
end






