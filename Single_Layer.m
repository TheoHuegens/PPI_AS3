function [h2,k2] = Single_Layer(MercuryLayers)

    %% planetary model
    
    % Fe,sol Inner-Core, dummy liquid layer (1)
    Interior_Model_Mercury(1).R0=5;
    Interior_Model_Mercury(1).rho0= MercuryLayers(2);
    
    % Mantle layer (2)
    Interior_Model_Mercury(2).R0= MercuryLayers(1);  %surface radius
    Interior_Model_Mercury(2).rho0= MercuryLayers(2);
    Interior_Model_Mercury(2).Ks0= MercuryLayers(3); % Bulk modulus
    Interior_Model_Mercury(2).mu0= MercuryLayers(4);  %shear modulus
    Interior_Model_Mercury(2).eta0= MercuryLayers(5);  %viscosity
    
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