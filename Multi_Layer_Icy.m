function [h2,k2] = Multi_Layer_Icy(MercuryLayers)

    %% planetary model
    % Fe,sol Inner-Core, dummy liquid layer (1)
    Interior_Model_Mercury(1).R0= 5;
    Interior_Model_Mercury(1).rho0= 7225;
    Interior_Model_Mercury(1).Ks0= 127e9;
    Interior_Model_Mercury(1).mu0= 100e9;
    Interior_Model_Mercury(1).eta0= 1e20;
    % defaults from bulk
    % Interior_Model_Mercury(1).rho0= 3758.609;
    % Interior_Model_Mercury(1).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(1).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(1).eta0= 1.03e19;  %viscosity
    
    % Fe,sol Inner-Core (2)
    Interior_Model_Mercury(2).R0= 427.8042;
    Interior_Model_Mercury(2).rho0= 7225;
    Interior_Model_Mercury(2).Ks0= 127e9;
    Interior_Model_Mercury(2).mu0= 100e9;
    Interior_Model_Mercury(2).eta0= 1e20;
    % defaults from bulk
    % Interior_Model_Mercury(2).rho0= 3758.609;
    % Interior_Model_Mercury(2).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(2).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(2).eta0= 1.03e19;  %viscosity
    
    % FeS,liq Outer-Core (3)
    Interior_Model_Mercury(3).R0= 1113.349;
    Interior_Model_Mercury(3).rho0= 7009;
    Interior_Model_Mercury(3).Ks0= 85.9e9;
    Interior_Model_Mercury(3).mu0= 70e9;
    Interior_Model_Mercury(3).eta0= 1e20;
    Interior_Model_Mercury(3).ocean= 1;
    % defaults from bulk
    % Interior_Model_Mercury(3).rho0= 3758.609;
    % Interior_Model_Mercury(3).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(3).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(3).eta0= 1.03e19;  %viscosity
    
    % MA mantle (4) 
    Interior_Model_Mercury(4).R0= 2326.946; 
    Interior_Model_Mercury(4).rho0= 3307.6; 
    Interior_Model_Mercury(4).Ks0= 129.9e9;
    Interior_Model_Mercury(4).mu0= 8.548402e9;  %65e6; 
    Interior_Model_Mercury(4).eta0= 1.03e19;  %1e11;
    % defaults from bulk
    % Interior_Model_Mercury(4).rho0= 3758.609;
    % Interior_Model_Mercury(4).Ks0= 57.22543e9; % Bulk modulus
    % Interior_Model_Mercury(4).mu0= 8.548402e9;  %shear modulus
    % Interior_Model_Mercury(4).eta0= 1.03e19;  %viscosity
    
    %  crust (5) 
    Interior_Model_Mercury(5).R0= 2439.4; 
    Interior_Model_Mercury(5).rho0= 3100; 
    Interior_Model_Mercury(5).Ks0= 120e9;
    Interior_Model_Mercury(5).mu0= 55e9; 
    Interior_Model_Mercury(5).eta0= 1e23;
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
    
    k2_Mercury=Love_Spectra_Mercury.k(iforcing); 
    h2_Mercury=Love_Spectra_Mercury.h(iforcing);
    
    disp('MERCURY');
    disp(['k_2 ' num2str(k2_Mercury)]);
    disp(['h_2 ' num2str(h2_Mercury)]);

end