
pest1 = zeros(1, length(psi1));
pest2 = zeros(1, length(psi2));
for j1 = 1:length(psi1)
    for j2 = 1:length(psi2)
        pest1(j1) = pest1(j1) + pest(j1, j2)*dpsi(2);
        pest2(j2) = pest2(j2) + pest(j1, j2)*dpsi(1);
    end
end

figure(1)
subplot(2,2,1);
surf(psi2, psi1, pest)
xlabel('\psi'', Hz');
ylabel('\psi, cycles');
zlabel('p(\psi, \psi'')')
subplot(2,2,2);
surf(psi2, psi1, pest)
xlabel('\psi'', Hz');
ylabel('\psi, cycles');
zlabel('p(\psi, \psi'')')
subplot(2,2,3);
plot(psi1, pest1);
xlabel('\psi, cycles');
ylabel('p(\psi)')
subplot(2,2,4);
plot(psi2, pest2);
xlabel('\psi'', Hz');
ylabel('p(\psi'')')

m1 = 0;
for j1 = 1:length(psi1)
    m1 = m1 + psi1(j1)*pest1(j1);
end
m1 = m1 * dpsi(1);

m2 = 0;
for j2 = 1:length(psi2)
    m2 = m2 + psi2(j2)*pest2(j2);
end
m2 = m2 * dpsi(2);

D1 = 0;
for j1 = 1:length(psi1)
    D1 = D1 + (psi1(j1) - m1)^2*pest1(j1);
end
D1 = D1 * dpsi(1);

D2 = 0;
for j2 = 1:length(psi2)
    D2 = D2 + (psi2(j2) - m2)^2*pest2(j2);
end
D2 = D2 * dpsi(2);

Cova = 0;
MM = 0;
for j1 = 1:length(psi1)
    for j2 = 1:length(psi2)
        Cova = Cova +  (psi1(j1) - m1)*(psi2(j2) - m2)*pest(j1, j2)*dpsi(2)*dpsi(1);
        MM = MM +  psi1(j1)*psi2(j2)*pest(j1, j2)*dpsi(2)*dpsi(1);
    end
end
ro = Cova / D1 / D2;


fprintf('Norm pest1 = %f, m1 = %f, D1 = %f\n', sum(pest1)*dpsi(1), m1, D1);
fprintf('Norm pest2 = %f, m2 = %f, D2 = %f\n', sum(pest2)*dpsi(2), m2, D2);
fprintf('Cov = %f, corr = Cov/(D1*D2) = %f, M[PhikOmk] = %f\n', Cova, ro, MM);