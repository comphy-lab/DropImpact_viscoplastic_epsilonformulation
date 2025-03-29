import pandas as pd
import numpy as np
import os
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
# Define parameter lists
Js = [0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10]
a0s = [0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8]

# Fixed parameters
Bo = 0.0
Wes = [1,5,10,50]
Oh = 0.0035
MAXlevel = 9


currentDir=os.path.dirname(os.path.abspath(__file__))
results_folder = os.path.join(currentDir, "Results_Running")
def extract(We):
    results = []
    for J in Js:
        for a0 in a0s:
            # Construct the filename using the parameters
            filename = f"Bo{Bo}-We{We}-J{J}-Oh{Oh}-a{a0}-MAXlevel{MAXlevel}.csv"
            filepath = os.path.join(results_folder, filename)
            if not os.path.exists(filepath):
                print(f"{filename} does not exist")        
                continue
            # Read the data
            df = pd.read_csv(filepath)
            
            # Extract time and vcm1
            t = df["t"].values
            vcm1 = df["vcm1"].values
            window_length = 7  # 必须为奇数
            polyorder = 3       # 必须小于 window_length
            vcm1_smooth = savgol_filter(vcm1, window_length=window_length, polyorder=polyorder)
            dvcm1_dt = np.gradient(vcm1_smooth, t)
            dvcm1_dt=dvcm1_dt/4.0
            # Compute dvcm1_dt
            dvcm1_dt = np.gradient(vcm1, t)
            
            # Find maximum of dvcm1_dt
            dvcm1_dt_max = dvcm1_dt.max()
            
            # Find maximum of Rmax
            Rmax_max = df["Rmax"].max()
            
            # Append (J, a0, dvcm1_dt_max, Rmax_max, We) to results
            results.append([J, a0, dvcm1_dt_max, Rmax_max])

    # Create a DataFrame from results
    df_summary = pd.DataFrame(results, columns=["J", "a0", "dvcm1_dt_max", "Rmax_max"])

    # Save to CSV
    df_summary.to_csv(os.path.join(currentDir, f"summary_{We}.csv"), index=False)
    print("Summary CSV file has been created: summary.csv")

        # -------------------------
    # Plot dvcm1_dt_max vs. a0
    # -------------------------
    plt.figure()  # new figure for dvcm1_dt_max vs a0
    for J in sorted(df_summary["J"].unique()):
        sub_df = df_summary[df_summary["J"] == J].sort_values("a0")
        plt.plot(
            sub_df["a0"], 
            sub_df["dvcm1_dt_max"], 
            marker='o', 
            label=f"J={J}"
        )
    plt.xlabel("a0")
    plt.xscale('log')
    plt.ylabel("dvcm1_dt_max")
    plt.title(f"dvcm1_dt_max vs a0 (We={We})")
    plt.legend()
    plt.savefig(os.path.join(currentDir, f"dvcm1_dt_max_We{We}.png"), dpi=300)
    plt.close()

    # ---------------------
    # Plot Rmax_max vs. a0
    # ---------------------
    plt.figure()  # new figure for Rmax_max vs a0
    for J in sorted(df_summary["J"].unique()):
        sub_df = df_summary[df_summary["J"] == J].sort_values("a0")
        plt.plot(
            sub_df["a0"], 
            sub_df["Rmax_max"], 
            marker='s', 
            label=f"J={J}"
        )
    plt.xlabel("a0")
    plt.xscale('log')
    plt.ylabel("Rmax_max")
    plt.title(f"Rmax_max vs a0 (We={We})")
    plt.legend()
    plt.savefig(os.path.join(currentDir, f"Rmax_max_We{We}.png"), dpi=300)
    plt.close()


for We in Wes:
    extract(We)