%% clear
clearvars
clc
set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name=matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))

%% Layer Model

% multi layer
MercuryModel = [
    [427.8042,7225,127e9,100e9,1e20],
    [1113.349,7019,85.9e9,0,0],
    [2326.946,3307.6,129.9e9,65e9,1e11],
    [2439.4,3100,120e9,55e9,1e23]
    ];

% single layer
MercuryModelBase = [
    [2439.4,3758.61,57.23e9,8.55e9,1.03e9]
    ];

% 1D parameters variations
Nvars = 10;
MagVar = logspace(-1,1,Nvars);
ModelTests = cell(1);
for t = 1:length(MercuryModel(1,:))-1 % for each variable, except boundaries
    for v = 1:Nvars % for each variable change
        TestModel = MercuryModel; % copy of base model
        TestModel(:,t+1) = MercuryModel(:,t+1)*MagVar(v); % replace with test value
        ModelTests{t,v} = TestModel; % save model as test
    end
end
%ModelTests = reshape(ModelTests,[],1);
[szx,szy] = size(ModelTests);
ResultTests = zeros(szx,szy,2);

%% test variations
for t = 1:szx
    for v = 1:szy
        MercuryLayers = ModelTests(t,v);
        [h2,k2] = Multi_Layer_Icy(MercuryLayers{1});
        ResultTests(t,v,:) = [h2,k2];
    end
end
ResultTestsH = transpose(real(ResultTests(:,:,1)));
ResultTestsK = transpose(real(ResultTests(:,:,2)));

%% plot
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