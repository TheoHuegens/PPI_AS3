close all

param_legend = {'density [kg/m3]','bulk modulus [GPa]','shear modulus [GPa]','viscosity [Pa s]'};
aa = 20;
bb = 10;

figure(1)
for v = 1:length(param_legend)
    semilogx(MagVar,ResultTestsH(:,v),'LineWidth',2);
    hold on;
end
axis([MagVar(1) MagVar(end) -2 2])
legend(param_legend, 'Location', 'southwest','Fontsize',bb);
ylabel("h2 [-]",'Fontsize',aa);
xlabel("parameter change factor [-]",'Fontsize',aa)
movegui(figure(1), [600 25]);
hold off;

figure(2)
for v = 1:length(param_legend)
    semilogx(MagVar,ResultTestsK(:,v),'LineWidth',2);
    hold on;
end
axis([MagVar(1) MagVar(end) -2 2])
legend(param_legend, 'Location', 'southwest','Fontsize',bb);
ylabel("k2 [-]",'Fontsize',aa);
xlabel("parameter chaneg factor [-]",'Fontsize',aa)
movegui(figure(2), [0 25]);
hold off;