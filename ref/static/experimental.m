J = linspace(0,1.5,50);

for i=1:length(J)
    CT(i) = 0.0966*J(i)^3 - 0.4296*J(i)^2 + 0.0946*J(i) + 0.42;

    CP(i) = -0.0320 *J(i)^4 +0.4146*J(i)^3 - 1.4525*J(i)^2 + 1.3564*J(i) - 0.0453;
    eff(i) = J(i)*CT(i)/CP(i);
end
figure
plot(J, CT);
hold on;
plot(J, CP);
% hold on;
% plot(J(40:end), eff(40:end));