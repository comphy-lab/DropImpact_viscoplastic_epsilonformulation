import os
import subprocess
import multiprocessing as mp
import csv
import numpy as np
from functools import partial
import matplotlib.pyplot as plt

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
    """
    Runs the specified executable with the given arguments and attempts to parse
    all floating-point values from the first line of the stderr output.

    If parsing or execution fails, this function returns an empty tuple.

    Parameters
    ----------
    executable : str
        The name or path of the executable, e.g. "./getResults".
    args : list
        A list of arguments passed to the executable, e.g. ["filename"].

    Returns
    -------
    tuple of float
        A tuple containing all parsed floats from the first line of stderr.
        Returns an empty tuple if parsing or execution fails.
    """
    try:
        completed = subprocess.run(
            [executable, *args],
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


def process_time_step(ti, tsnap):
    """
    Processes a single time step 'ti':
      1. Compute the actual time t = ti * tsnap
      2. Build the path to the intermediate file
      3. If the file exists, parse it with get_results
      4. Return the parsed data if successful, otherwise None
    """
    t = tsnap * ti
    filepath = f"intermediate/snapshot-{t:.4f}"    
    if not os.path.exists(filepath):
        # print(f"File not found: {filepath}")
        return None    
    data = run_executable("./getResults",filepath)
    # If any value in 'data' is None, we consider it a parsing failure.
    if any(v is None for v in data):
        return None
    
    return data

def process_and_write_results(tsnap, tmax, output_results_csv):
    """
    Compiles the C program, processes the results in parallel, and writes the results to a CSV file.
    """
    nsteps = int(tmax / tsnap)
    
    # Remove existing files if present
    for f in [output_results_csv, "getResults"]:
        if os.path.exists(f):
            os.remove(f)

    # Compile the C program
    compile_executable()

    # Use partial to fix the 'tsnap' argument, so Pool.map only needs 'ti'
    process_func = partial(process_time_step, tsnap=tsnap)

    # Parallel processing using a process pool
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(process_func, range(1, nsteps + 1))

    # Write results to a CSV file
    header = [
        "tp", "n", "xTP", "vTP", "kappaTip", "xToP", "vToP", "kappamax", "xcmax", "ycmax",
        "thetap", "R_max", "R_max_bottom", "x_top_droplet", "r_top_droplet", "u_top_droplet", "vb", "xb", "ke", "pforce"
    ]
    
    with open(output_results_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in results:
            if row is not None:
                writer.writerow(row)

    print(f"Parallel processing finished. Results have been written to {output_results_csv}.")
    return results


def compute_and_write_derivatives(results, output_force_csv):
    """
    Computes derivatives (F_vb and F_xb) and writes them along with tp and pforce to a new CSV file.
    """
    # Extract data
    tp = np.array([row[0] for row in results if row is not None])  # Column 1: tp
    vb = np.array([row[16] for row in results if row is not None])  # Column 17: vb
    xb = np.array([row[17] for row in results if row is not None])  # Column 18: xb
    pforce = np.array([row[19] for row in results if row is not None])  # Column 20: pforce

    # Compute derivatives
    F_vb = np.gradient(vb, tp)  # First derivative of vb
    F_xb = np.gradient(np.gradient(xb, tp), tp)  # Second derivative of xb

    # Prepare output data
    output_data = zip(tp, F_vb, F_xb, pforce)

    # Write to CSV file
    
    output_header = ["tp", "F_vb", "F_xb", "pforce"]

    with open(output_force_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(output_header)  # Write the header
        writer.writerows(output_data)  # Write the data rows

    print(f"Results have been written to {output_force_csv}.")
    return tp, F_vb, F_xb, pforce


def plot_forces(tp, F_vb, F_xb, pforce, output_force_fig):
    """
    Generates and saves the plot for F_vb, F_xb, and pforce against tp in a single coordinate system with different colors.
    """
    # Create a figure and axis
    plt.figure(figsize=(10, 6))

    # Plot all three forces on the same axes
    plt.plot(tp, F_vb, color="blue", label="F_vb", linestyle="-")
    plt.plot(tp, F_xb, color="green", label="F_xb", linestyle="--")
    plt.plot(tp, pforce, color="red", label="pforce", linestyle=":")

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

def plot_Top(output_results_csv,output_top_fig):
    """
    Reads data from csv_filename and plots:
      1. tp vs. vTP
      2. tp vs. kappaTip
    Then saves them as PNG files.
    """
    # Lists to store the data
    tp_data = []
    vTP_data = []
    kappaTip_data = []
    
    # Read CSV data
    with open(output_results_csv, "r", newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # You may add error handling for missing or invalid data if necessary
            tp_data.append(float(row["tp"]))
            vTP_data.append(float(row["vTP"]))
            kappaTip_data.append(float(row["kappaTip"]))

    # Create a figure and a set of subplots
    fig, ax1 = plt.subplots(figsize=(6, 4))
    
    # Plot vTP on the left y-axis
    color_1 = 'blue'
    ax1.set_xlabel("Time")
    ax1.set_ylabel("velocity_Tip", color=color_1)
    ax1.plot(tp_data, vTP_data, marker='o', linestyle='-', color=color_1, label='vTP')
    ax1.tick_params(axis='y', labelcolor=color_1)

    # Create a twin of the first axes that shares the same x-axis
    ax2 = ax1.twinx()
    
    # Plot kappaTip on the right y-axis
    color_2 = 'red'
    ax2.set_ylabel("kappa_Tip", color=color_2)
    ax2.plot(tp_data, kappaTip_data, marker='s', linestyle='-', color=color_2, label='kappaTip')
    ax2.tick_params(axis='y', labelcolor=color_2)

    # Optional: set a single title for both axes
    plt.title(f"{output_top_fig}")

    # Improve spacing
    fig.tight_layout()

    # Save the figure
    plt.savefig(output_top_fig, dpi=300)
    plt.close()

    print(f"Plot created and saved as {output_top_fig}.")

def process_EpsForce(ti, tsnap,Oh,We,J):
    t = tsnap * ti
    filepath = f"intermediate/snapshot-{t:.4f}"    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None    
    data = run_executable("./getEpsForce",[filepath,Oh,We,J])
    # If any value in 'data' is None, we consider it a parsing failure.
    if any(v is None for v in data):
        return None    
    return data

def getEpsForces(tsnap, tmax, epsForce_csv, Oh, We, J):
    nsteps = int(tmax / tsnap)    
    if os.path.exists(epsForce_csv):
        os.remove(epsForce_csv)
    process_func = partial(process_EpsForce, tsnap=tsnap,Oh=Oh,We=We,J=J)
    # Parallel processing using a process pool
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(process_func, range(1, nsteps + 1))
    for i in range(nsteps + 1):
        results=process_func(i)
    # Write results to a CSV file
    header = ["t","epsOh","epsJ","pforce","betaMax","vBeta","Hmax","vH"]    
    with open(epsForce_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in results:
            if row is not None:
                writer.writerow(row)
    print(f"Parallel processing finished. Results have been written to {epsForce_csv}.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--tMAX', type=float, default=10.0, help='tMAX')
    parser.add_argument('--tSNAP', type=float, default=0.01, help='tSNAP')
    parser.add_argument('--We', type=float, default=10.0, help='We')
    parser.add_argument('--Oh', type=float, default=0.01, help='Oh')
    parser.add_argument('--J', type=float, default=0.0, help='J')
    args = parser.parse_args()
    tsnap = args.tSNAP
    tmax = args.tMAX
    We = args.We
    Oh = args.Oh
    J = args.J
    #  file dirs
    current_folder = os.path.basename(os.path.dirname(__file__))
    output_results_csv = f"results_{current_folder}.csv"
    output_force_csv = f"force_{current_folder}.csv"
    output_force_fig=f"figure_force_{current_folder}.png"
    output_top_fig=f"figure_tip_{current_folder}.png"
    epsForce_csv=f"epsForce_{current_folder}.csv"
    # Step 1: Process and write results
    # results = process_and_write_results(tsnap, tmax, output_results_csv)
    getEpsForces(tsnap, tmax, epsForce_csv, Oh, We, J)
    # Step 2: Compute derivatives and write to a separate CSV
    # tp, F_vb, F_xb, pforce = compute_and_write_derivatives(results, output_force_csv)

    # Step 3: Plot force results
    # plot_forces(tp, F_vb, F_xb, pforce, output_force_fig)

    # Step 4: Plot vTP and kappaTip
    # plot_Top(output_results_csv,output_top_fig)
