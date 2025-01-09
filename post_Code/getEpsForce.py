import os
import subprocess
import multiprocessing as mp
import csv
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import glob

def compile_executable():
    """
    Compiles the C program getResults.c.
    Raises RuntimeError if the compilation fails.
    """
    if os.path.exists("getResults"):
        os.remove("getResults")
    cmd = "qcc -w -Wall -O2 -disable-dimensions getResults.c -o getResults -lm"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError("Compilation of getResults.c failed!")

def run_executable(executable, args):
    args_str = [str(a) for a in args]
    try:
        completed = subprocess.run(
            [executable, *args_str],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Split the stderr output into lines
        lines = completed.stderr.strip().split("\n")

        if lines:
            # Split the first line by whitespace
            parts = lines[0].split()
            # Convert each part to float
            return tuple(float(x) for x in parts)
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"Error running {executable} with args {args}: {e}")

    # If there's an error or no lines, return an empty tuple
    return ()

def process_EpsForce(ti, tsnap,Oh,We,J):
    t = tsnap * ti
    filepath = f"intermediate/snapshot-{t:.4f}"  
    # print(f"[PID {mp.current_process().pid}] Processing: {filepath}", flush=True)  
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None    
    data = run_executable("./getEpsForce",[filepath,Oh,We,J])
    # If any value in 'data' is None, we consider it a parsing failure.
    if any(v is None for v in data):
        return None 
    # print(f"[PID {mp.current_process().pid}] Successfully processed: {filepath}", flush=True)   
    return data

def getEpsForces(tsnap, tmax, epsForce_csv, Oh, We, J, CPUStoUse):
    if os.path.exists(epsForce_csv):
        os.remove(epsForce_csv)
    file_list = sorted(glob.glob("intermediate/snapshot-*"))
    nsteps = len(file_list)  # This is now the actual number of snapshot files
    # nsteps = int(tmax / tsnap)    
    if os.path.exists(epsForce_csv):
        os.remove(epsForce_csv)
    process_func = partial(process_EpsForce, tsnap=tsnap,Oh=Oh,We=We,J=J)
    # Parallel processing using a process pool
    # with mp.Pool(processes=5) as pool:
    #     results = pool.map(process_func, range(0, nsteps + 1))
    results = []
    with mp.Pool(processes=CPUStoUse) as pool:
        # Use imap instead of map so we can iterate results as they come in
        # Wrap the iterator with tqdm for a progress bar
        for res in tqdm(
            pool.imap(process_func, range(nsteps + 1)),
            total=nsteps + 1,
            desc="Processing snapshots"
        ):
            results.append(res)

    # Write results to a CSV file
    header = ["t","mv","vcm","mv1","vcm1","epsOh", "epsJ","pforce","betaMax","vBeta","Hmax","vH"]    
    with open(epsForce_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in results:
            if row is not None:
                writer.writerow(row)
    print(f"Parallel processing finished. Results have been written to {epsForce_csv}.")
    return results

def plot_forces(epsForce_csv, output_force_fig):
    if os.path.exists(output_force_fig):
        os.remove(output_force_fig)
    # 1. Read the CSV file using pandas
    data = pd.read_csv(epsForce_csv)
    # 2. Extract specific columns using column names
    tp = data["t"].to_numpy()  # Column "t" for time points
    mv = data["mv"].to_numpy()  # Column "mv"
    mv1 = data["mv1"].to_numpy()  # Column "mv1"
    pforce = data["pforce"].to_numpy()  # Column "pforce"

    F_vb = np.gradient(mv, tp)
    F_vb1= np.gradient(mv1, tp)
    data["F_vb"] = F_vb
    data["F_vb1"] = F_vb1
    data.to_csv(epsForce_csv, index=False)
    print(f"Updated CSV with F_vb and F_vb1 saved as {epsForce_csv}.")


    # Create a figure and axis
    plt.figure(figsize=(10, 6))
    # Plot all three forces on the same axes
    # plt.plot(tp, pforce, color="green", label="F_p", linestyle="--")
    plt.plot(tp, F_vb, color="blue", label="F_vb", linestyle=":")
    plt.plot(tp, F_vb1, color="red", label="F_vb_main", linestyle="-")    

    # Add labels, legend, and title
    plt.xlabel("Time")
    plt.ylabel("Forces")
    plt.title(f"{output_force_fig}")
    plt.legend(loc="best")
    plt.grid(True)

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_force_fig, dpi=300)
    plt.close()
    print(f"Combined plot saved as {output_force_fig}.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--tMAX', type=float, default=10.0, help='tMAX')
    parser.add_argument('--tSNAP', type=float, default=0.01, help='tSNAP')
    parser.add_argument('--We', type=float, default=10.0, help='We')
    parser.add_argument('--Oh', type=float, default=0.01, help='Oh')
    parser.add_argument('--J', type=float, default=0.0, help='J')
    args = parser.parse_args()
    CPUStoUse = args.CPUs
    tsnap = args.tSNAP
    tmax = args.tMAX
    We = args.We
    Oh = args.Oh
    J = args.J
    #  file dirs
    current_folder = os.path.basename(os.path.dirname(__file__))
    output_results_csv = f"results_{current_folder}.csv"
    output_force_fig=f"figure_force_{current_folder}.png"
    output_top_fig=f"figure_tip_{current_folder}.png"
    epsForce_csv=f"epsForce_{current_folder}.csv"
    # Step 1: Process and write results
    # results = process_and_write_results(tsnap, tmax, output_results_csv)
    results=getEpsForces(tsnap, tmax, epsForce_csv, Oh, We, J,CPUStoUse)
    plot_forces(epsForce_csv, output_force_fig)
    # Step 2: Compute derivatives and write to a separate CSV
    # tp, F_vb, F_xb, pforce = compute_and_write_derivatives(results, output_force_csv)

    # Step 3: Plot force results
    # plot_forces(tp, F_vb, F_xb, pforce, output_force_fig)

    # Step 4: Plot vTP and kappaTip
    # plot_Top(output_results_csv,output_top_fig)
