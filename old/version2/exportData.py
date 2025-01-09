import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# List of J values
J_values = [0, 0.0001, 0.001, 0.01, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5]

# Directory containing the CSV files
results_folder = "Results_Running"

# Initialize dictionaries to store data for exporting
export_data_vcm = {}
export_data_Fvcm = {}
export_data_Zmax = {}
export_data_Rmax = {}
export_data_dt_Zmax = {}
export_data_dt_Rmax = {}
# Find the maximum length of all time series
max_length = 0
for J in J_values:
    filename = f"Bo0.0-We10-J{J}-Oh0.01-MAXlevel10-epsilon0.01.csv"
    filepath = os.path.join(results_folder, filename)
    df = pd.read_csv(filepath)
    max_length = max(max_length, len(df["t"]))

for J in J_values:
    # Construct the CSV filename
    filename = f"Bo0.0-We10-J{J}-Oh0.01-MAXlevel10-epsilon0.01.csv"
    filepath = os.path.join(results_folder, filename)
    
    # Read the data
    df = pd.read_csv(filepath)
    t = df["t"]
    vcm = df["vcm"]
    Rmax = df["Rmax"]
    Zmax = np.where(df["Zmax"] < 0, 0, df["Zmax"])  # Set values less than 0 to 0

    # Calculate derivatives
    dvcm_dt = np.gradient(vcm, t)
    dZmax_dt = np.gradient(Zmax, t)
    dRmax_dt = np.gradient(Rmax, t)

    t = np.pad(t, (0, max_length - len(t)), constant_values=np.nan)
    vcm = np.pad(vcm, (0, max_length - len(vcm)), constant_values=np.nan)
    Zmax = np.pad(Zmax, (0, max_length - len(Zmax)), constant_values=np.nan)
    Rmax = np.pad(Rmax, (0, max_length - len(Rmax)), constant_values=np.nan)
    dvcm_dt = np.pad(dvcm_dt, (0, max_length - len(dvcm_dt)), constant_values=np.nan)
    dZmax_dt = np.pad(dZmax_dt, (0, max_length - len(dZmax_dt)), constant_values=np.nan)
    dRmax_dt = np.pad(dRmax_dt, (0, max_length - len(dRmax_dt)), constant_values=np.nan)

    # Store the data in the dictionaries
    export_data_vcm[f"t_J={J}"] = t
    export_data_vcm[f"vcm_J={J}"] = vcm
    export_data_Fvcm[f"t_J={J}"] = t
    export_data_Fvcm[f"F_vcm_J={J}"] = dvcm_dt
    export_data_Zmax[f"t_J={J}"] = t
    export_data_Zmax[f"Zmax_J={J}"] = Zmax
    export_data_Rmax[f"t_J={J}"] = t
    export_data_Rmax[f"Rmax_J={J}"] = Rmax
    export_data_dt_Zmax[f"t_J={J}"] = t
    export_data_dt_Zmax[f"dZmax_dt_J={J}"] = dZmax_dt
    export_data_dt_Rmax[f"t_J={J}"] = t
    export_data_dt_Rmax[f"dRmax_dt_J={J}"] = dRmax_dt

# Convert dictionaries to DataFrames
df_vcm = pd.DataFrame(export_data_vcm)
df_Fvcm = pd.DataFrame(export_data_Fvcm)
df_Zmax = pd.DataFrame(export_data_Zmax)
df_Rmax = pd.DataFrame(export_data_Rmax)
df_dt_Zmax = pd.DataFrame(export_data_dt_Zmax)
df_dt_Rmax = pd.DataFrame(export_data_dt_Rmax)

# Export DataFrames to CSV
df_vcm.to_csv("export_vcm.csv", index=False)
df_Fvcm.to_csv("export_Fvcm.csv", index=False)
df_Zmax.to_csv("export_Zmax.csv", index=False)
df_Rmax.to_csv("export_Rmax.csv", index=False)
df_dt_Zmax.to_csv("export_dt_Zmax.csv", index=False)
df_dt_Rmax.to_csv("export_dt_Rmax.csv", index=False)
