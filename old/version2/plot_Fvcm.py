import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# List of J values
J_values = [0,0.0001,0.001,0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.4,0.5]

# Directory containing the CSV files
results_folder = "Results_Running"

# Create a new figure
plt.figure(figsize=(8, 6))

for J in J_values:
    # Construct the CSV filename
    # Adjust the formatting of J if needed (e.g., 0.1 instead of 0.10)
    filename = f"Bo0.0-We10-J{J}-Oh0.01-MAXlevel10-epsilon0.01.csv"
    filepath = os.path.join(results_folder, filename)
    
    # Read the data
    df = pd.read_csv(filepath)
    
    # Assume the CSV has columns "t" and "vcm" (adjust if different)
    t = df["t"]
    vcm = df["vcm"]
    
    # Plot the curve
    plt.plot(t, vcm, label=f"J = {J}")

# Add legend, axis labels, and title
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel("vcm")
plt.title("vcm vs. t for Different J Values")

# Save the figure as a PNG file
plt.savefig("figure_vcm_vs_t_for_different_J.png", dpi=300)



# Create a new figure
plt.figure(figsize=(8, 6))

for J in J_values:
    # Construct the CSV filename
    # Adjust the formatting of J if needed (e.g., 0.1 instead of 0.10)
    filename = f"Bo0.0-We10-J{J}-Oh0.01-MAXlevel10-epsilon0.01.csv"
    filepath = os.path.join(results_folder, filename)
    
    # Read the data
    df = pd.read_csv(filepath)
    
    # Assume the CSV has columns "t" and "vcm" (adjust if different)
    t = df["t"]
    vcm = df["vcm"]
    dvcm_dt = np.gradient(vcm, t)
    dvcm_dt[0]=0
    
    # Plot the curve
    plt.plot(t, dvcm_dt, label=f"J = {J}")

# Add legend, axis labels, and title
plt.legend()
plt.xlabel("Time (t)")
plt.ylabel("dvcm_dt")
plt.title("F_vcm vs. t for Different J Values")

# Save the figure as a PNG file
plt.savefig("figure_F_vcm_vs_t_for_different_J.png", dpi=300)