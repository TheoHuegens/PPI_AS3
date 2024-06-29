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
        [h2,k2] = Multi_Layer_Icy_Eval(MercuryLayers{1});
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

function [h2,k2] = Multi_Layer_Icy_Eval(MercuryLayers)

    %% planetary model
    % Fe,sol Inner-Core, dummy liquid layer (1)
    Interior_Model_Mercury(1).R0= 5;
    Interior_Model_Mercury(1).rho0= MercuryLayers(1,2);
    Interior_Model_Mercury(1).Ks0= MercuryLayers(1,3);
    Interior_Model_Mercury(1).mu0= MercuryLayers(1,4);
    Interior_Model_Mercury(1).eta0= MercuryLayers(1,5);
    % defaults from bulk
    % Interior_Model_Mercury(1).rho0= 3758.609;
    % Interior_Model_Mercury(1).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(1).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(1).eta0= 1.03e19;  %viscosity
    
    % Fe,sol Inner-Core (2)
    Interior_Model_Mercury(2).R0= MercuryLayers(1,1);
    Interior_Model_Mercury(2).rho0= MercuryLayers(1,2);
    Interior_Model_Mercury(2).Ks0= MercuryLayers(1,3);
    Interior_Model_Mercury(2).mu0= MercuryLayers(1,4);
    Interior_Model_Mercury(2).eta0= MercuryLayers(1,5);
    % defaults from bulk
    % Interior_Model_Mercury(2).rho0= 3758.609;
    % Interior_Model_Mercury(2).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(2).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(2).eta0= 1.03e19;  %viscosity
    
    % FeS,liq Outer-Core (3)
    Interior_Model_Mercury(3).R0= MercuryLayers(2,1);
    Interior_Model_Mercury(3).rho0= MercuryLayers(2,2);
    Interior_Model_Mercury(3).Ks0= MercuryLayers(2,3);
    Interior_Model_Mercury(3).mu0= MercuryLayers(2,4);
    Interior_Model_Mercury(3).eta0= MercuryLayers(2,5);
    Interior_Model_Mercury(3).ocean= 1;
    % defaults from bulk
    % Interior_Model_Mercury(3).rho0=(4, 3758.609;
    % Interior_Model_Mercury(3).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(3).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(3).eta0= 1.03e19;  %viscosity
    
    % MA mantle (4) 
    Interior_Model_Mercury(4).R0= MercuryLayers(3,1);
    Interior_Model_Mercury(4).rho0= MercuryLayers(3,2);
    Interior_Model_Mercury(4).Ks0= MercuryLayers(3,3);
    Interior_Model_Mercury(4).mu0= MercuryLayers(3,4);
    Interior_Model_Mercury(4).eta0= MercuryLayers(3,5);
    % defaults from bulk
    % Interior_Model_Mercury(4).rho0= 3758.609;
    % Interior_Model_Mercury(4).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(4).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(4).eta0= 1.03e19;  %viscosity
    
    %  crust (5) 
    Interior_Model_Mercury(5).R0= MercuryLayers(4,1);
    Interior_Model_Mercury(5).rho0= MercuryLayers(4,2);
    Interior_Model_Mercury(5).Ks0= MercuryLayers(4,3);
    Interior_Model_Mercury(5).mu0= MercuryLayers(4,4);
    Interior_Model_Mercury(5).eta0= MercuryLayers(4,5);
    % defaults from bulk
    % Interior_Model_Mercury(5).rho0= 3758.609;
    % Interior_Model_Mercury(5).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(5).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(5).eta0= 1.03e19;  %viscosity
    
    %% forcing
    Forcing_Mercury(1).Td=87.969*24*3600; 
    Forcing_Mercury(1).n=2; 
    Forcing_Mercury(1).m=0; 
    Forcing_Mercury(1).F=1;
    
    %% Numerics
    %radial discretization
    Numerics_Mercury.Nlayers = length(Interior_Model_Mercury); % number of concentric layers. Including the core!
    Numerics_Mercury.method = 'variable'; % method of setting the radial points per layer
    Numerics_Mercury.Nrbase = 2000;
    %code parallelization
    Numerics_Mercury.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
    Numerics_Mercury.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
    % lateral variations
    Numerics_Mercury.perturbation_order = 2; %maximum order to which couplings are considered
    Numerics_Mercury.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
    Numerics_Mercury.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
    Numerics_Mercury.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
    Numerics_Mercury.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
    [Numerics_Mercury, Interior_Model_Mercury] = set_boundary_indices(Numerics_Mercury, Interior_Model_Mercury,'verbose');
    
    %% tidal response
    Interior_Model_Mercury_U = get_rheology(Interior_Model_Mercury,Numerics_Mercury,Forcing_Mercury);
    [Love_Spectra_Mercury,y_rad_Mercury]=get_Love(Interior_Model_Mercury_U,Forcing_Mercury,Numerics_Mercury,'verbose');
    
    iforcing=find(Love_Spectra_Mercury.n==Forcing_Mercury.n & Love_Spectra_Mercury.m==Forcing_Mercury.m);
    
    k2=Love_Spectra_Mercury.k(iforcing); 
    h2=Love_Spectra_Mercury.h(iforcing);

end