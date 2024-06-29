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

% parameter baselines
mu0 = 70e9; %rock
eta0 = 1e20;

% multi layer
% MercuryModel = [
%     [427.8042,7225,127e9,100e9,1e20],
%     [1113.349,7019,85.9e9,0,0],
%     [2326.946,3307.6,129.9e9,65e9,1e11],
%     [2439.4,3100,120e9,55e9,1e23]
%     ];

% single layer
MercuryModel = [
    [2439.4,3758.61,57.23e9,8.55e9,1.03e9]
    ];

% parameters variations
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

%% test variations
for t = 1:szx
    for v = 1:szy
        MercuryLayers = ModelTests(t,v);
        [h2,k2] = Single_Layer(MercuryLayers{1});
        ResultTests(t,v,:) = [h2,k2];
    end
end
ResultTestsH = real(ResultTests(:,:,1));
ResultTestsK = real(ResultTests(:,:,2));

%% plot
% see plot HK2 to change labels accordingly