%% requirements:
% 1. gcc for compiling c codes. Load before starting Matlab:
%    module load gcc/6.3.0
% 2. Use Matlab 2020a
%% Setup the directory where the membrane object is located and add the directory to Matlab's function pool
%dir_mod = '/home2/s171152/codes/matlab/mine/git/memCompCourse/memcompcourse';
dir_mod = './';
addpath(dir_mod);

% Parameter ranges
k_c_values = [1, 5, 10]; % example values
kv_values = [10, 20, 30]; % example values
ks_values = [8, 16, 32]; % example values

% Define the directory to save the results
results_dir = './results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Loop through all combinations of k_c, kv, and ks
for k_c = k_c_values
    for kv = kv_values
        for ks = ks_values
            % create 'unit' u using the unit module, and 'membrane' m using the membrane module
            u = ComUnit('erg', ComUnit.nm_to_cm(1000), 300, ComUnit.kBT_to_erg(10, 300));
            m = ModMembrane(2, 'unit', u);
            m.pm.Vdh.V0 = 0.1; % Setting internal force's strength using recommended value
            m.pm.k_c = k_c;
            m.pm.k_V = 4;
            m.pm.k_A = 8;
            m.pm.k_P = 0;
            m_original = m;
            
            % 1. getting internal force:
            Fi = Finternal(m, 'plot_or_not', false);
            
            % other parameter recommendations:
            mu = 100;
            
            % Parameters for volume and surface area forces
            targetVolume = 0.6 * sum(Volume(m_original)); % Set your target volume here
            targetSurfaceArea = 1.0 * sum(Area(m_original)); % Set your target surface area here
            % print targetVolume and targetSurfaceArea
            fprintf('Target Volume: %f, Target Surface Area: %f\n', targetVolume, targetSurfaceArea);
            
            n_iter = 5000;
            stds = zeros(n_iter, 1);
            
            % Perturbation amount for finite difference calculation
            epsilon = 1e-4;
            
            for iter = 1:n_iter
                % 2. spatial range for remesh later
                rLim = [min(Fi.rg), max(Fi.rg)];
                % Compute bending force (Fb) using finite difference method
                Fb = zeros(size(m.var.coord));
                Fv = zeros(size(m.var.coord));
                Fs = zeros(size(m.var.coord));
                parfor i = 1:length(m.var.coord)
                    for dim = 1:3
                        m_update = m;
                        r_orig = m_update.var.coord(i, dim);
                        % Perturb positively
                        m_update.var.coord(i, dim) = r_orig + epsilon;
                        Ev_before = kv * ((sum(Volume(m_update)) - targetVolume)^2) / targetVolume;
                        Es_before = ks * ((sum(Area(m_update)) - targetSurfaceArea)^2) / targetSurfaceArea;
                        Hp = Helfrich(m_update); % Compute Helfrich free energy
                        % Perturb negatively
                        m_update.var.coord(i, dim) = r_orig - epsilon;
                        Ev_after = kv * ((sum(Volume(m_update)) - targetVolume)^2) / targetVolume;
                        Es_after = ks * ((sum(Area(m_update)) - targetSurfaceArea)^2) / targetSurfaceArea;
                        Hm = Helfrich(m_update); % Compute Helfrich free energy
                        % Reset the coordinate
                        m_update.var.coord(i, dim) = r_orig;
                        % Finite difference approximation of the gradient
                        Fb(i, dim) = -(sum(Hp) - sum(Hm)) / (2 * epsilon);
                        Fv(i, dim) = -(Ev_before - Ev_after) / (2 * epsilon);
                        Fs(i, dim) = -(Es_before - Es_after) / (2 * epsilon);
                    end
                end
                Fb_smooth = zeros(size(m.var.coord));
                Fv_smooth = zeros(size(m.var.coord));
                Fs_smooth = zeros(size(m.var.coord));
                for i = 1:size(m.var.coord, 1)
                    edges = m.var.edge_all(m.var.edge_all(:, 1) == i | m.var.edge_all(:, 2) == i, :);
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

                % 3. getting the adaptive time step
                [dt, Ftotal, l] = varDt(m, Fi, Fb+Fv_smooth+Fs_smooth, mu);
                m.var.coord = m.var.coord + m.pm.mu * Ftotal * dt;
                % 4. remeshing
                [m, ~] = RemeshCtrl(m, Fi, rLim, 'l0', l, 'print_or_not', false);

            end

            % Plot and save the membrane 'm'
            fig = figure;
            subplot(1, 2, 1);
            col = rand(m_original.var.n_coord, 3);
            plot(m_original, 'f', fig, 'col', col, 'col_min', 0, 'col_max', 1, 'colBar', true);
            subplot(1, 2, 2);
            col = rand(m.var.n_coord, 3);
            plot(m, 'f', fig, 'col', col, 'col_min', 0, 'col_max', 1, 'colBar', true);

            % Save the figure
            saveas(fig, fullfile(results_dir, sprintf('membrane_kc_%d_kv_%d_ks_%d.png', k_c, kv, ks)));
            close(fig);

        end
    end
end
