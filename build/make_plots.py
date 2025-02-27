import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import argparse

def check_file(file_path, threshold = 1e5):
  """
  Check if the file exists and if the highest "DoFs" entry is above the threshold.
  Returns True if the file exists and passes the check, False otherwise.
  """
  if not os.path.exists(file_path):
    print(f"File not found: {file_path}")
    return False
  try:
    df = pd.read_csv(file_path)
  except Exception as e:
    print(f"Error reading {file_path}: {e}")
    return False

  if "DoFs" not in df.columns:
    print(f"'DoFs' column missing in {file_path}")
    return False

  max_dofs = df["DoFs"].max()
  if max_dofs > threshold:
    print(f"{file_path}: max DoFs = {max_dofs} (OK, threshold: {threshold})")
    return True
  else:
    print(f"{file_path}: max DoFs = {max_dofs} (insufficient, threshold: {threshold})")
    return False


def make_step_14_plots(file_name, output_file, reference_file_name=None, exact_value=None):
    data = pd.read_csv(file_name)
    pd.options.display.float_format = '{:.10e}'.format  # 10 digits of precision

    if reference_file_name is not None:
        reference_data = pd.read_csv(reference_file_name)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

    axes[0].loglog(data['DoFs'], data['ex POINT err']/exact_value, "r-+", linewidth=1.0, label='Exact Error')
    axes[0].loglog(data['DoFs'], data['est err']/exact_value, "g-+", linewidth=1.0, label='Estimated Error')
    if reference_file_name is not None:
        axes[0].loglog(reference_data['DoFs'], reference_data['ex POINT err']/exact_value, "b:+", linewidth=0.5, label='GlobRef exact error (Reference)')
    axes[0].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
    axes[0].set_ylabel('Relative error', fontsize=12)
    axes[0].set_title('Exact vs Estimated POINT error', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, which="both", linestyle='--', linewidth=0.5)

    axes[1].loglog(data['DoFs'], data['ex POINT err']/exact_value, "r-+", linewidth=0.5, label='Exact Error')
    axes[1].loglog(data['DoFs'], abs(data['ex POINT err'] - data["est err"])/exact_value, "y-+", linewidth=0.5, label=r'$|e_{\mathrm{ex}} - e_{\mathrm{est}}|$')
    axes[1].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
    axes[1].set_title('Precision of the estimate', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_classic_pair_flux(file_name, output_file, reference_file_name=None, exact_value= None):
  data = pd.read_csv(file_name)

  if reference_file_name is not None:
    reference_data = pd.read_csv(reference_file_name)

  fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

  axes[0].loglog(data['DoFs'], data['std FLUX err']/exact_value,"r-+", linewidth=1.0, label='Exact Error')
  axes[0].loglog(data['DoFs'], data['est err']/exact_value,"g-+", linewidth=1.0, label='Estimated Error')
  if reference_file_name is not None:
    axes[0].loglog(reference_data['DoFs'], reference_data['std FLUX err']/exact_value,"b:+", linewidth=0.5, label='GlobRef exact error (Reference)')
  axes[0].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[0].set_ylabel('Relative error', fontsize=12)
  axes[0].set_title('Exact vs Estimated error for the FLUX at the emitter', fontsize=14)
  axes[0].legend()
  axes[0].grid(True, which="both", linestyle='--', linewidth=0.5)

  axes[1].loglog(data['DoFs'], data['std FLUX err']/exact_value, "r-+", linewidth=0.5, label='Exact Error')
  axes[1].loglog(data['DoFs'], abs(data['std FLUX err']-data["est err"])/exact_value,"y-+", linewidth=0.5, label=r'$|e_{\mathrm{ex}} - e_{\mathrm{est}}|$')
  axes[1].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[1].set_title('Precision of the estimate', fontsize=14)
  axes[1].legend()
  axes[1].grid(True, which="both", linestyle='--', linewidth=0.5)

  # Adjust layout and show the plots
  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()

def plot_a_sei(glob_file_name1,fGO_file_name1,concentric_file_name1,glob_file_name2,fGO_file_name2,concentric_file_name2, output_file, name1, name2, title, exact_value):
  glob_data1 = pd.read_csv(glob_file_name1)
  fGO_data1 = pd.read_csv(fGO_file_name1)
  concentric_data1 = pd.read_csv(concentric_file_name1)

  glob_data2 = pd.read_csv(glob_file_name2)
  fGO_data2 = pd.read_csv(fGO_file_name2)
  concentric_data2 = pd.read_csv(concentric_file_name2)

  plt.figure(figsize=(7, 6))  # Adjust the figure size for two plots

  plt.subplot(1, 1, 1)  # Second subplot (bottom)
  plt.loglog(glob_data1['DoFs'], glob_data1['std FLUX err']/exact_value,"b:+", linewidth=1.0, label='['+name1+'] Global refinement')
  plt.loglog(fGO_data1['DoFs'], fGO_data1['std FLUX err']/exact_value,"r:+", linewidth=1.0, label='['+name1+'] Flux goal-oriented')
  plt.loglog(concentric_data1['DoFs'], concentric_data1['std FLUX err']/exact_value,":+", color='tab:brown', linewidth=1.0, label='['+name1+'] Concentric refinement')

  plt.loglog(glob_data2['DoFs'], glob_data2['std FLUX err']/exact_value,"b-+", linewidth=1.0, label='['+name2+'] Global refinement')
  plt.loglog(fGO_data2['DoFs'], fGO_data2['std FLUX err']/exact_value,"r-+", linewidth=1.0, label='['+name2+'] Flux goal-oriented')
  plt.loglog(concentric_data2['DoFs'], concentric_data2['std FLUX err']/exact_value,"-+", color='tab:brown', linewidth=1.0, label='['+name2+'] Concentric refinement')

  plt.xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  plt.ylabel('Relative flux error', fontsize=12)
  plt.title(title, fontsize=14)
  plt.legend()
  plt.grid(True, which="both", linestyle='--', linewidth=0.5)

  # Adjust layout and show the plots
  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()


def linea_gialla_accanto(file_name1,file_name2,output_file,name1,name2,exact_value1,exact_value2):

  data1 = pd.read_csv(file_name1)
  data2 = pd.read_csv(file_name2)

  data1['std FLUX err'] = data1['std FLUX err'] / exact_value1
  data1['est err'] = data1['est err'] / exact_value1
  data2['std FLUX err'] = data2['std FLUX err'] / exact_value2
  data2['est err'] = data2['est err'] / exact_value2


  # Creazione della figura con due subplot
  fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

  axes[0].loglog(data1['DoFs'], data1['std FLUX err'], "r-+", linewidth=0.5, label='Exact Error')
  axes[0].loglog(data1['DoFs'], abs(data1['std FLUX err'] - data1["est err"]), "y-+", linewidth=0.5, label=r'$|e_{\mathrm{ex}} - e_{\mathrm{est}}|$')
  axes[0].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[0].set_ylabel('Relative error', fontsize=12)
  axes[0].set_title(name1, fontsize=14)
  axes[0].legend()
  axes[0].grid(True, which="both", linestyle='--', linewidth=0.5)

  axes[1].loglog(data2['DoFs'], data2['std FLUX err'], "r-+", linewidth=0.5, label='Exact Error')
  axes[1].loglog(data2['DoFs'], abs(data2['std FLUX err'] - data2["est err"]), "y-+", linewidth=0.5, label=r'$|e_{\mathrm{ex}} - e_{\mathrm{est}}|$')
  axes[1].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[1].set_title(name2, fontsize=14)
  axes[1].legend()
  axes[1].grid(True, which="both", linestyle='--', linewidth=0.5)

  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()

def plot_a_cinque_L2_H1(glob_file_name, fGO_file_name, concentric_file_name, kelly_file_name, weighted_kelly_file_name, output_file, exact_value):
  glob_data = pd.read_csv(glob_file_name)
  fGO_data = pd.read_csv(fGO_file_name)
  concentric_data = pd.read_csv(concentric_file_name)
  kelly_data = pd.read_csv(kelly_file_name)
  weighted_kelly_data = pd.read_csv(weighted_kelly_file_name)

  fig, ax = plt.subplots(1, 1, figsize=(7, 6))

  ax.loglog(glob_data['DoFs'], glob_data['H1'],"b-.+", linewidth=1.0, label='H1 Global refinement')
  ax.loglog(fGO_data['DoFs'], fGO_data['H1'],"r-.+", linewidth=1.5, label='H1 Flux dual-weighted goal-oriented')
  ax.loglog(concentric_data['DoFs'], concentric_data['H1'],"-.+", color='tab:brown', linewidth=1.0, label='H1 Concentric refinement')
  ax.loglog(kelly_data['DoFs'], kelly_data['H1'],"m-.+", linewidth=1.0, label='H1 KellyErrorEstimator')
  ax.loglog(weighted_kelly_data['DoFs'], weighted_kelly_data['H1'],"c-.+", linewidth=1.0, label='H1 Weighted Kelly')


  ax.loglog(glob_data['DoFs'], glob_data['L2'],"b:+", linewidth=1.0, label='L2 Global refinement')
  ax.loglog(fGO_data['DoFs'], fGO_data['L2'],"r:+", linewidth=1.5, label='L2 Flux dual-weighted goal-oriented')
  ax.loglog(concentric_data['DoFs'], concentric_data['L2'],":+", color='tab:brown', linewidth=1.0, label='L2 Concentric refinement')
  ax.loglog(kelly_data['DoFs'], kelly_data['L2'],"m:+", linewidth=1.0, label='L2 KellyErrorEstimator')
  ax.loglog(weighted_kelly_data['DoFs'], weighted_kelly_data['L2'],"c:+", linewidth=1.0, label='L2 Weighted Kelly')


  ax.set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  ax.set_ylabel('error', fontsize=12)
  ax.set_title(r'Solution error in norms $L^2$ and $H^1$', fontsize=12)
  ax.legend(fontsize=8, loc='lower left')
  ax.grid(True, which="both", linestyle='--', linewidth=0.5)

  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()

def plot_a_cinque(glob_file_name, fGO_file_name, concentric_file_name, kelly_file_name, weighted_kelly_file_name, output_file, exact_value):
  glob_data = pd.read_csv(glob_file_name)
  fGO_data = pd.read_csv(fGO_file_name)
  concentric_data = pd.read_csv(concentric_file_name)
  kelly_data = pd.read_csv(kelly_file_name)
  weighted_kelly_data = pd.read_csv(weighted_kelly_file_name)

  fig, ax = plt.subplots(1, 1, figsize=(7, 6))
  
  ax.loglog(glob_data['DoFs'], glob_data['std FLUX err']/exact_value,"b-+", linewidth=0.5, label='Global refinement')
  ax.loglog(fGO_data['DoFs'], fGO_data['std FLUX err']/exact_value,"r-+", linewidth=1.0, label='Flux dual-weighted goal-oriented')
  ax.loglog(concentric_data['DoFs'], concentric_data['std FLUX err']/exact_value,"-+", color='tab:brown', linewidth=0.5, label='Concentric refinement')
  ax.loglog(kelly_data['DoFs'], kelly_data['std FLUX err']/exact_value,"m-+", linewidth=0.5, label='KellyErrorEstimator')
  ax.loglog(weighted_kelly_data['DoFs'], weighted_kelly_data['std FLUX err']/exact_value,"c-+", linewidth=0.5, label='Weighted Kelly')

  ax.set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  ax.set_ylabel('Relative flux error', fontsize=12)
  ax.set_title('Flux at the emitter error', fontsize=12)
  ax.legend()
  ax.grid(True, which="both", linestyle='--', linewidth=0.5)

  # Adjust layout and show the plots
  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()

def plot_classic_pair_real_application(file_name, output_file, reference_file_name=None, concentric_file_name = None):
  data = pd.read_csv(file_name)

  if reference_file_name is not None:
    reference_data = pd.read_csv(reference_file_name)

  if concentric_file_name is not None:
    concentric_data = pd.read_csv(concentric_file_name) 
  
  dofs = int(list(data['DoFs'])[-1])
  reference_value = list(data['std FLUX value'])[-1]
  max_dofs = dofs/50
  data = data[data['DoFs']<max_dofs]
  data['std FLUX err'] = abs(data['std FLUX value']-reference_value)
  if reference_file_name is not None:
    reference_data = reference_data[reference_data['DoFs']<max_dofs]
    reference_data['std FLUX err'] = abs(reference_data['std FLUX value']-reference_value)
  if concentric_file_name is not None:
    concentric_data = concentric_data[concentric_data['DoFs']<max_dofs]
    concentric_data['std FLUX err'] = abs(concentric_data['std FLUX value']-reference_value)

  fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

  axes[0].loglog(data['DoFs'], data['std FLUX err']/reference_value,"r-+", linewidth=1.0, label='Exact Error')
  axes[0].loglog(data['DoFs'], data['est err']/reference_value,"g-+", linewidth=1.0, label='Estimated Error')
  if reference_file_name is not None:
    axes[0].loglog(reference_data['DoFs'], reference_data['std FLUX err']/reference_value,"b-+", linewidth=0.5, label='GlobRef exact error (Reference)')
  if concentric_file_name is not None:
    axes[0].loglog(concentric_data['DoFs'], concentric_data['std FLUX err']/reference_value,"c-+", linewidth=0.5, label='Exact error refining only around emitters')
  axes[0].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[0].set_ylabel('Relative error', fontsize=12)
  axes[0].set_title('Exact vs Estimated error for the FLUX at the emitter', fontsize=14)
  axes[0].legend()
  axes[0].grid(True, which="both", linestyle='--', linewidth=0.5)

  axes[1].loglog(data['DoFs'], data['std FLUX err']/reference_value, "r-+", linewidth=0.5, label='Exact Error')
  axes[1].loglog(data['DoFs'], abs(data['std FLUX err']-data["est err"])/reference_value,"y-+", linewidth=0.5, label=r'$|e_{\mathrm{ex}} - e_{\mathrm{est}}|$')
  axes[1].set_xlabel('Degrees of Freedom (DoFs)', fontsize=12)
  axes[1].set_title('Precision of the estimate', fontsize=14)
  axes[1].legend()
  axes[1].grid(True, which="both", linestyle='--', linewidth=0.5)

  # Adjust layout and show the plots
  plt.tight_layout()
  plt.savefig(output_file, dpi=300)
  plt.close()


def main():
  parser = argparse.ArgumentParser(
        description="Check required results files and run creation commands if needed."
    )
  parser.add_argument(
      "-n", "--num_procs", type=int, required=True,
      help="Number of processors to use for mpirun"
  )
  args = parser.parse_args()

  required_files = {
      # From test-1
      "../tests/test-1_original-step14/results/config-1/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-1_original-step14/config-1-GlobRef.yaml",
          "threshold": 3e4
      },
      "../tests/test-1_original-step14/results/config-2/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-1_original-step14/config-2-GO.yaml",
          "threshold": 3e4
      },
      # NOTE: The reference case is ONLY available in the current branch, i.e. the serial code. (You are in the right place)

      # From test-16
      "../tests/test-16-LogCircular-1-2/results/config-5/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-16-LogCircular-1-2/config-5-basis-fGO.yaml",
          "threshold": 1.9e5
      },
      "../tests/test-16-LogCircular-1-2/results/config-4/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-16-LogCircular-1-2/config-4-basis-global.yaml",
          "threshold": 1.9e5
      },
      "../tests/test-16-LogCircular-1-2/results/config-6/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-16-LogCircular-1-2/config-6-FlatManif-global.yaml",
          "threshold": 1.9e5
      },
      "../tests/test-16-LogCircular-1-2/results/config-7/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-16-LogCircular-1-2/config-7-FlatManif-fGO.yaml",
          "threshold": 1.9e5
      },

      # From test-8
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-11/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-11-FlatManif-fGO.yaml",
          "threshold": 4e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-10/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-10-FlatManif-glob.yaml",
          "threshold": 4e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-12/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-12-FlatManif-only_concentric.yaml",
          "threshold": 3e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-19/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-19-MappingQ2-FlatManif-global.yaml",
          "threshold": 3e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-20-MappingQ2-FlatManif-fGO.yaml",
          "threshold": 3e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-21/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-21-MappingQ2-FlatManif-only_concentric.yaml",
          "threshold": 3e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20b/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-20b-MappingQ2-FlatManif-Kelly.yaml",
          "threshold": 6e5
      },
      "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20c/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-8_LogCircular-1-100-fluxGO/config-20c-MappingQ2-FlatManif-KellyWeight.yaml",
          "threshold": 6e5
      },

      # From test-18
      "../tests/test-18-WireWire/results/config-12/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-18-WireWire/config-12-myWWdelaunay3-MappingQ2-FlatM-global.yaml",
          "threshold": 2e4
      },
      "../tests/test-18-WireWire/results/config-13/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-18-WireWire/config-13-myWWdelaunay3-MappingQ2-FlatM-fGO.yaml",
          "threshold": 1e6
      },
      "../tests/test-18-WireWire/results/config-14/convergence_results.csv": {
          "command": "./ion-propulsion ../tests/test-18-WireWire/config-14-myWWdelaunay3-MappingQ2-FlatM-concentric.yaml",
           "threshold": 2e4
      },
      # NOTE: Limited to test-18: to obtain the SAME image of the manuscript, use the parallel DISTRIBUTED code instead.
      #       Since a precise enough estimate needs to be obtained to become the reference value, the computational cost of this simulation is high.
      #       The max number of DoFs has been REDUCED here.
      
  }

  for file_path, config in required_files.items():
      print(f"Processing: {file_path}")
      threshold = config.get("threshold", 1e5)  # Use default threshold if not specified
      if check_file(file_path, threshold):
          print("  -> File exists and condition is satisfied.")
      else:
          # File missing or condition not met: trigger the command.
          command = config["command"]
          print(f"  -> Running: {command}")
          try:
              subprocess.run(command, shell=True, check=True)
          except subprocess.CalledProcessError as err:
              print(f"Error while executing command for {file_path}: {err}")
  
  # Figure 4: Reference case (step-14)
  make_step_14_plots("../tests/test-1_original-step14/results/config-2/convergence_results.csv", 
                    "Figure_4.pdf",
                    "../tests/test-1_original-step14/results/config-1/convergence_results.csv", 
                    exact_value=0.0334473)
  # NOTE: The reference case is ONLY available in the current branch, i.e. the serial code. (You are in the right place)

  # Figure 6: Flux, Annulus 1:2
  plot_classic_pair_flux("../tests/test-16-LogCircular-1-2/results/config-5/convergence_results.csv",
                        "Figure_6.pdf",
                        "../tests/test-16-LogCircular-1-2/results/config-4/convergence_results.csv",
                        exact_value=1.812944056731e+05)

  # Figure 7: Flux, Annulus 1:2, after introduction manifold
  plot_classic_pair_flux("../tests/test-16-LogCircular-1-2/results/config-7/convergence_results.csv",
                        "Figure_7.pdf",
                        "../tests/test-16-LogCircular-1-2/results/config-6/convergence_results.csv",
                        exact_value=1.812944056731e+05)

  # Figure 8: Flux, Annulus 1:100
  plot_classic_pair_flux("../tests/test-8_LogCircular-1-100-fluxGO/results/config-11/convergence_results.csv",
                        "Figure_8.pdf",
                        "../tests/test-8_LogCircular-1-100-fluxGO/results/config-10/convergence_results.csv",
                        exact_value=2.728752707684e+04)

  # Figure 9: Comparison MappingQ1 vs MappingQ2 [Errros' plot]
  plot_a_sei(glob_file_name1 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-10/convergence_results.csv",
            fGO_file_name1 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-11/convergence_results.csv",
            concentric_file_name1 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-12/convergence_results.csv", 
            glob_file_name2 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-19/convergence_results.csv",
            fGO_file_name2 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20/convergence_results.csv", 
            concentric_file_name2 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-21/convergence_results.csv",
            output_file = "Figure_9.pdf", 
            name1 = 'MappingQ1', 
            name2 = 'MappingQ2', 
            title = 'Fixed: Flat Manifold', 
            exact_value = 2.728752707684e+04)

  # Figure 10: Comparison MappingQ1 vs MappingQ2 [Precision's plot]
  linea_gialla_accanto(file_name1 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-11/convergence_results.csv",
                      file_name2 = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20/convergence_results.csv",
                      output_file = "Figure_10.pdf", 
                      name1 = "MappingQ1",
                      name2 = "MappingQ2",
                      exact_value1 = 2.728752707684e+04,
                      exact_value2 = 2.728752707684e+04)

  # Figure 11: Comparison different refinement algorithms [L2 and H1 errors]
  plot_a_cinque_L2_H1(glob_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-19/convergence_results.csv",
                      fGO_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20/convergence_results.csv",
                      concentric_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-21/convergence_results.csv",
                      kelly_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20b/convergence_results.csv", 
                      weighted_kelly_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20c/convergence_results.csv",
                      output_file = "Figure_11.pdf", 
                      exact_value = 2.728752707684e+04)

  # Figure 12: Comparison different refinement algorithms [FLUX error]
  plot_a_cinque(glob_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-19/convergence_results.csv",
                fGO_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20/convergence_results.csv",
                concentric_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-21/convergence_results.csv",
                kelly_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20b/convergence_results.csv", 
                weighted_kelly_file_name = "../tests/test-8_LogCircular-1-100-fluxGO/results/config-20c/convergence_results.csv",
                output_file = "Figure_12.pdf", 
                exact_value = 2.728752707684e+04)

  # Figure 15: Real ion thruster mesh
  plot_classic_pair_real_application(file_name = "../tests/test-18-WireWire/results/config-13/convergence_results.csv",
                                    reference_file_name = "../tests/test-18-WireWire/results/config-12/convergence_results.csv",
                                    concentric_file_name = "../tests/test-18-WireWire/results/config-14/convergence_results.csv",
                                    output_file = "Figure_15.pdf")
  # NOTE: Limited to Figure 15: to obtain the SAME image of the manuscript, use the parallel DISTRIBUTED code instead.
  #       Since a precise enough estimate needs to be obtained to become the reference value, the computational cost of this simulation is high.
  #       The max number of DoFs has been REDUCED here.

if __name__ == "__main__":
    main()