%% requirements:
% 1. gcc for compiling c codes. Load before starting Matlab:
%    module load gcc/6.3.0
% 2. Use Matlab 2020a
%% Setup the directory where the membrane object is located and add the directory to Matlab's function pool
dir_mod = './';
addpath(genpath(dir_mod));

% Create 'unit' u using the unit module, and 'membrane' m using the membrane module
u = ComUnit('erg', ComUnit.nm_to_cm(1000), 300, ComUnit.kBT_to_erg(10, 300));
m = ModMembrane(2, 'unit', u);
m.pm.Vdh.V0 = 0.1; % Setting internal force's strength using recommended value
m.pm.k_c = 5;
m.pm.k_V = 8;
m.pm.k_A = 16;
m.pm.k_P = 0;
m.pm.dt = 0.005;
m.pm.Vdh.V0 = 0.2;
m.pm.f_const_std_std = 0.001;
m.pm.kBT = 0.1;
% m.pm.V0 = 150 * 0.6; % 60% of spherical volume
% m.pm.A0 = 137;
m.pm.nAVmean = 4;
m_original = m;

% 1. Getting internal force:
Fi = Finternal(m, 'plot_or_not', false);

% Other parameter recommendations:
mu = 200;

% Parameters for volume and surface area forces
targetVolume = 0.6 * sum(Volume(m_original)); % Set your target volume here
targetSurfaceArea = 1.0 * sum(Area(m_original)); % Set your target surface area here
fprintf('Target Volume: %f, Surface Area: %f\n', targetVolume, targetSurfaceArea);

kv = 16; % Set appropriate value for kv
ks = 32; % Set appropriate value for ks

n_iter = 3000;
stds = zeros(n_iter, 1);

% Perturbation amount for finite difference calculation
epsilon = 1e-4;

for iter = 1:n_iter
    % 2. Spatial range for remesh later
    rLim = [min(Fi.rg), max(Fi.rg)];
    
    % Pre-compute volumes and areas for perturbation
    volumes = sum(Volume(m));
    areas = sum(Area(m));
    
    % Compute bending, volume, and surface forces using finite difference method
    Fb = zeros(size(m.var.coord));
    Fv = zeros(size(m.var.coord));
    Fs = zeros(size(m.var.coord));
    
    parfor i = 1:length(m.var.coord)
        for dim = 1:3
            m_update = m;
            r_orig = m_update.var.coord(i, dim);
            
            % Perturb positively
            m_update.var.coord(i, dim) = r_orig + epsilon;
            Ev_before = kv * (sum(Volume(m_update)) - targetVolume)^2 / targetVolume;
            Es_before = ks * (sum(Area(m_update)) - targetSurfaceArea)^2 / targetSurfaceArea;
            Hp = Helfrich(m_update); % Compute Helfrich free energy
            
            % Perturb negatively
            m_update.var.coord(i, dim) = r_orig - epsilon;
            Ev_after = kv * (sum(Volume(m_update)) - targetVolume)^2 / targetVolume;
            Es_after = ks * (sum(Area(m_update)) - targetSurfaceArea)^2 / targetSurfaceArea;
            Hm = Helfrich(m_update); % Compute Helfrich free energy
            
            % Reset the coordinate
            m_update.var.coord(i, dim) = r_orig;
            
            % Finite difference approximation of the gradient
            Fb(i, dim) = -(sum(Hp) - sum(Hm)) / (2 * epsilon);
            Fv(i, dim) = -(Ev_before - Ev_after) / (2 * epsilon);
            Fs(i, dim) = -(Es_before - Es_after) / (2 * epsilon);
        end
    end
    
    % Smoothing forces
    [Fv_smooth, Fs_smooth] = smoothForces(Fv, Fs, m.var.edge_all, m.var.coord);
    
    if mod(iter, 100) == 0
        fig = figure;
        col = m.Helfrich();
        plot(m, 'f', fig, 'col', col, 'col_min', 0, 'col_max', 1, 'colBar', true);
        fprintf('Volume: %f, Surface Area: %f\n', sum(Volume(m)), sum(Area(m)));
    end
    
    % 3. Getting the adaptive time step
    [dt, Ftotal, l] = varDt(m, Fi, Fb + Fv_smooth + Fs_smooth, mu);
    m.var.coord = m.var.coord + m.pm.mu * Ftotal * dt;
    
    % 4. Remeshing
    [m, ~] = RemeshCtrl(m, Fi, rLim, 'l0', l, 'print_or_not', false);
end

% Plot the results
plot(1:n_iter, stds);

% Plot the membrane 'm'
fig = figure;
subplot(1, 2, 1);
col = m_original.Helfrich();
plot(m_original, 'f', fig, 'col', col, 'col_min', 0, 'col_max', 1, 'colBar', true);
subplot(1, 2, 2);
col = m.Helfrich();
plot(m, 'f', fig, 'col', col, 'col_min', 0, 'col_max', 1, 'colBar', true);

% Smoothing function
function [Fv_smooth, Fs_smooth] = smoothForces(Fv, Fs, edge_all, coord)
    Fv_smooth = zeros(size(coord));
    Fs_smooth = zeros(size(coord));
    for i = 1:size(coord, 1)
        edges = edge_all(edge_all(:, 1) == i | edge_all(:, 2) == i, :);
        neighbors = unique(edges(edges ~= i));
        Fv_smooth(i,:) = Fv(i,:);
        Fs_smooth(i,:) = Fs(i,:);
        for j = 1:size(neighbors,1)
            Fv_smooth(i,:) = Fv_smooth(i,:) + Fv(neighbors(j),:);
            Fs_smooth(i,:) = Fs_smooth(i,:) + Fs(neighbors(j),:);
        end
        Fv_smooth(i,:) = Fv_smooth(i,:) / (size(neighbors, 1) + 1);
        Fs_smooth(i,:) = Fs_smooth(i,:) / (size(neighbors, 1) + 1);
    end
end
