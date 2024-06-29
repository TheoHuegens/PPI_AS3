%% clear
clearvars
close all
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

%% choose evals
eval_1V = true; % verif
eval_2D = false; % eval

%% Layer Model

% single layer
MercuryModel = [
    [2439.4,3758.61,57.23e9,8.55e9,1.03e9]
    ];

% parameter baselines
mu0 = 70e9; %rock
eta0 = 1e20;

%% 2D parameters variations
if eval_2D == true
    Nvars = 10;
    tests = [ % 1=boundary, 2=density, 3=bulk modulus, 4=shear modulus, 5=viscosity
        [4,linspace(0.1*mu0,10*mu0,Nvars)],
        [5,linspace(eta0*1e-5,eta0*1e5,Nvars)]
        ];
    
    ModelTests = cell(1);
    for t = 1:length(tests(1,:))-1
        for v = 1:length(tests(2,:))-1
            TestModel = MercuryModel; % copy of base model
            TestModel(tests(1,1)) = tests(1,t+1); % replace with test value
            TestModel(tests(2,1)) = tests(2,v+1); % replace with test value
            ModelTests{t,v} = TestModel; % save model as test
        end
    end
    %ModelTests = reshape(ModelTests,[],1);
    [szx,szy] = size(ModelTests);
    ResultTests = zeros(szx,szy,2);
    
    % test variations
    for t = 1:szx
        for v = 1:szy
            MercuryLayers = ModelTests(t,v);
            [h2,k2] = Single_Layer_Eval(MercuryLayers{1});
            ResultTests(t,v,:) = [h2,k2];
        end
    end
    ResultTestsH = real(ResultTests(:,:,1));
    ResultTestsK = real(ResultTests(:,:,2));
    
    % plot    
    aa = 20;
    bb = 10;
    
    x = tests(1,2:end)*1e-9;
    y = tests(2,2:end);
    
    figure(1)
    contourf(x,y,log10(ResultTestsK));
    c = colorbar;
    c.Label.String = 'log10(k2) [-]';
    c.Ticks = -3:1;
    c.TickLabels = compose('10^{%d}',c.Ticks);
    xlabel("mu, Shear Modulus [GPa]",'Fontsize',aa);
    ylabel("eta, Viscosity [Pa s]",'Fontsize',aa);
    movegui(figure(1), [600 25]);
    
    figure(2)
    contourf(x,y,log10(ResultTestsH));
    c = colorbar;
    c.Label.String = 'log10(h2) [-]';
    c.Ticks = -3:1;
    c.TickLabels = compose('10^{%d}',c.Ticks);
    xlabel("mu, Shear Modulus [GPa]",'Fontsize',aa);
    ylabel("eta, Viscosity [Pa s]",'Fontsize',aa);
    movegui(figure(2), [0 25]);
    
    figure(3)
    pcolor(x,y,ResultTestsK);
    p.EdgeColor = 'none';
    set(gca,'ColorScale','log');
    
    c = colorbar;
    c.Label.String = 'log10(k2) [-]';
    xlabel("mu, Shear Modulus [GPa]",'Fontsize',aa);
    ylabel("eta, Viscosity [Pa s]",'Fontsize',aa);
    movegui(figure(3), [600 400]);
    
    figure(4)
    pcolor(x,y,ResultTestsH);
    p.EdgeColor = 'none';
    set(gca,'ColorScale','log');
    
    c = colorbar;
    c.Label.String = 'log10(h2) [-]';
    xlabel("mu, Shear Modulus [GPa]",'Fontsize',aa);
    ylabel("eta, Viscosity [Pa s]",'Fontsize',aa);
    movegui(figure(4), [0 400]);
end

%% elasticity variations
if eval_1V == true

    ModelTests = cell(1);
    
    % an elastic model (remove or comment Interior_Model.eta);
    TestModel = MercuryModel;
    TestModel(5) = 0; % replace with test value
    ModelTests{1} = TestModel; % save model as test
    
    % a viscoelastic model;
    TestModel = MercuryModel;
    ModelTests{2} = TestModel; % save model as test
    
    % an ideal fluid (this can be approximated by considering ω → 0).
    TestModel = MercuryModel;
    TestModel(4) = 5; % replace with test value
    ModelTests{3} = TestModel; % save model as test

    % test variations
    ResultTests = zeros(3,2);
    param_legend = {'elastic','viscoelastic','ideal fluid'};
    for t = 1:3
        MercuryLayers = ModelTests(t);
        [h2,k2] = Single_Layer_Eval(MercuryLayers{1});
        ResultTests(t,:) = [h2,k2];
    end
    ResultTestsH = real(ResultTests(:,1));
    ResultTestsK = real(ResultTests(:,2));

end

%% default model
[h2,k2] = Single_Layer_Eval(MercuryModel);
ResultTests = [h2,k2];

function [h2,k2] = Single_Layer_Eval(MercuryLayers)

    %% planetary model
    
    % Fe,sol Inner-Core, dummy liquid layer (1)
    Interior_Model_Mercury(1).R0=5;
    Interior_Model_Mercury(1).rho0= MercuryLayers(2);
    
    % Mantle layer (2)
    Interior_Model_Mercury(2).R0= MercuryLayers(1);  %surface radius
    Interior_Model_Mercury(2).rho0= MercuryLayers(2);
    Interior_Model_Mercury(2).Ks0= MercuryLayers(3); % Bulk modulus
    Interior_Model_Mercury(2).mu0= MercuryLayers(4);  %shear modulus
    if MercuryLayers(5) ~= 0
        Interior_Model_Mercury(2).eta0= MercuryLayers(5);  %viscosity
    end
    
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