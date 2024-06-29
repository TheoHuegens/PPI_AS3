close all

x = tests(1,2:end)*1e-9;
y = tests(2,2:end);

figure(1)
contourf(x,y,log10(ResultTestsK));
c = colorbar;
c.Label.String = 'log10(k2) [-]';
c.Ticks = -3:1;
c.TickLabels = compose('10^{%d}',c.Ticks);
xlabel("mu, Shear Modulus [GPa]");
ylabel("eta, Viscosity [Pa s]");
movegui(figure(1), [600 25]);

figure(2)
contourf(x,y,log10(ResultTestsH));
c = colorbar;
c.Label.String = 'log10(h2) [-]';
c.Ticks = -3:1;
c.TickLabels = compose('10^{%d}',c.Ticks);
xlabel("mu, Shear Modulus [GPa]");
ylabel("eta, Viscosity [Pa s]");
movegui(figure(2), [0 25]);

figure(3)
pcolor(x,y,ResultTestsK);
p.EdgeColor = 'none';
set(gca,'ColorScale','log');

c = colorbar;
c.Label.String = 'log10(k2) [-]';
xlabel("mu, Shear Modulus [GPa]");
ylabel("eta, Viscosity [Pa s]");
movegui(figure(3), [600 450]);

figure(4)
pcolor(x,y,ResultTestsH);
p.EdgeColor = 'none';
set(gca,'ColorScale','log');

c = colorbar;
c.Label.String = 'log10(h2) [-]';
xlabel("mu, Shear Modulus [GPa]");
ylabel("eta, Viscosity [Pa s]");
movegui(figure(4), [0 450]);

