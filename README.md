## README for RBC Membrane Shape Simulation
### Overview
This repository contains Matlab and C code to simulate the shape of a red blood cell using a membrane object. The code implements the Helfrich energy model, finite difference methods for calculating forces, and adaptive time steps for numerical stability.
### Requirements
1. **GCC**: Ensure GCC is installed for compiling the C codes. Load GCC before starting Matlab with the following command:
   ```sh
   module load gcc/6.3.0
   ```
2. **Matlab**: Use Matlab 2023a.
### Setup
1. Clone the repository to your local machine:
   ```sh
   git clone https://github.com/yaober/SE_RBC.git
   cd SE_RBC
   ```
2. Setup the directory where the membrane object is located and add the directory to Matlab's function pool. Update the directory path if needed:
   ```matlab
   dir_mod = './';
   addpath(dir_mod);
   ```
### Parameters
- **k_c_values**: [1, 5, 10] (Bending rigidity)
- **kv_values**: [10, 20, 30] (Volume constraint)
- **ks_values**: [8, 16, 32] (Surface area constraint)
Directory Structure
- **results_dir**: Directory to save the results. Create the directory if it doesn't exist:
  ```matlab
  results_dir = './results';
  if ~exist(results_dir, 'dir')
      mkdir(results_dir);
  end
  ```
### Running the Simulation
The simulation runs through all combinations of `k_c`, `kv`, and `ks` values. For each combination, it performs the following steps:
1. **Create unit and membrane objects**:
   ```matlab
   u = ComUnit('erg', ComUnit.nm_to_cm(1000), 300, ComUnit.kBT_to_erg(10, 300));
   m = ModMembrane(2, 'unit', u);
   m.pm.Vdh.V0 = 0.1;
   m.pm.k_c = 5;
   m.pm.k_V = 16;
   m.pm.k_A = 32;
   m.pm.k_P = 0;
   ```
2. **Calculate internal force**:
   ```matlab
   Fi = Finternal(m, 'plot_or_not', false);
   ```
3. **Set target volume and surface area**:
   ```matlab
   targetVolume = 0.6 * sum(Volume(m));
   targetSurfaceArea = 1.0 * sum(Area(m));
   ```
4. **Iterate through time steps**:
   ```matlab
   n_iter = 5000;
   epsilon = 1e-4;
   for iter = 1:n_iter
       % Compute forces, update membrane coordinates, and remesh
   end
   ```
5. **Plot and save results**:
   ```matlab
   fig = figure;
   subplot(1, 2, 1);
   plot(m_original, 'f', fig);
   subplot(1, 2, 2);
   plot(m, 'f', fig);
   saveas(fig, fullfile(results_dir, sprintf('membrane_kc_%d_kv_%d_ks_%d.png', k_c, kv, ks)));
   close(fig);
   ```
