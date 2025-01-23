import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Configuration parameters
J_values_Pre = [0,0.1]
maxLevel = 10
currentDir=os.path.dirname(os.path.abspath(__file__))
results_folder = os.path.join(currentDir, "Results_Running")
J_values=[]
# Data storage dictionary
data = {}
plt.rc('font', size=14)          # Default text size
plt.rc('axes', titlesize=18)     # Title font size
plt.rc('axes', labelsize=16)     # X and Y label size
plt.rc('xtick', labelsize=14)    # X tick label size
plt.rc('ytick', labelsize=14)    # Y tick label size
plt.rc('legend', fontsize=14)    # Legend font size
plt.rc('figure', titlesize=20)   # Figure title font size
# Load data
for J in J_values_Pre:
    filename = f"Bo0.0-We10-J{J}-Oh0.01-MAXlevel{maxLevel}-epsilon0.01.csv"
    filepath = os.path.join(results_folder, filename)
    if not os.path.exists(filepath):
        print(f"{filename} does not exist")
        continue
    if J==0:
        J_values.append(0.0009)
    else:
        J_values.append(J)
    df = pd.read_csv(filepath)
    t = df["t"]
    vcm1 = df["vcm1"]
    dvcm1_dt = np.gradient(vcm1, t)
    dvcm1_dt[0] = 0
    uH = df["uH"]
    vR = df["vR"]
    ke = df["ke"]

    mask1 = (t >= 0) & (t <= 2.5)
    dvcm1_dt_sub1 = dvcm1_dt[mask1]
    t_sub1 = t[mask1]
    if len(dvcm1_dt_sub1) > 0:
        F1 = np.max(dvcm1_dt_sub1)
        F1_time = t_sub1[np.argmax(dvcm1_dt_sub1)]
        mask2 = (t >= 2.5) & (t <= 5)

    dvcm1_dt_sub2 = dvcm1_dt[mask2]
    t_sub2 = t[mask2].values
    if len(dvcm1_dt_sub2) > 0:
        F2 = np.max(dvcm1_dt_sub2)
        F2_time = t_sub2[np.argmax(dvcm1_dt_sub2)]

    Zmax = df["Zmax"].values
    if len(uH) > 0:
        uHMAX = np.max(uH)
        uHMAX_t = t[np.argmax(uH)]

    Rmax = df["Rmax"].values
    if len(Rmax) > 0:
        RmaxMAX = np.max(Rmax)
        RmaxMAX_t = t[np.argmax(Rmax)]

    if len(vR) > 0:
        vRMAX = np.max(vR)
        vRMAX_t = t[np.argmax(vR)]

    data[J] = {"t": t, "vcm1": vcm1, "dvcm1_dt": dvcm1_dt, "uH": uH, "Zmax":Zmax, "vR": vR, "Rmax":Rmax, "F1_t":[F1,F1_time], "F2_t":[F2,F2_time], "uHMAX_t":[uHMAX,uHMAX_t], "vRMAX_t":[vRMAX,vRMAX_t], "RmaxMAX_t":[RmaxMAX,RmaxMAX_t], "ke":ke}

# Plotting function
def plot_figure(x, y, labels, xlabel, ylabel, title, save_as, xlim=None, ylim=None, xticks=None, yticks=None):
    plt.figure(figsize=(8, 6))
    for label, (x_vals, y_vals) in labels.items():
        plt.plot(x_vals, y_vals, label=f"J = {label}")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if xticks:
        plt.xticks(np.arange(*xticks))
    if yticks:
        plt.yticks(np.arange(*yticks))
    plt.savefig(save_as, dpi=300)
    plt.close()

# Plot mv-t
labels_mv = {J: (data[J]["t"], data[J]["vcm1"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_mv,
    xlabel="t", ylabel="mv",
    title="mv-t for Different J Values",
    save_as=os.path.join(results_folder, "figure_mv_vs_t_for_different_J.png")
)

# Plot F-t
labels_F = {J: (data[J]["t"], data[J]["dvcm1_dt"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_F,
    xlabel="t", ylabel="d(mv)/dt",
    title="F-t for Different J Values",
    ylim=(-1, 5),
    save_as=os.path.join(results_folder, f"figure_F_vs_t_for_different_J-{maxLevel}.png")
)

# Plot inset F-t
plot_figure(
    x=None, y=None, labels=labels_F,
    xlabel="t", ylabel="d(mv)/dt",
    title="Inset F-t for Different J Values",
    save_as=os.path.join(results_folder, f"figure_inset_F_vs_t_for_different_J-{maxLevel}.png"),
    xlim=(3.0, 3.8), ylim=(-1, 5),
    xticks=(3.0, 3.8, 0.2), yticks=(-1, 5, 1)
)

# Plot uH-t
labels_uH = {J: (data[J]["t"], data[J]["uH"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_uH,
    xlabel="t", ylabel="uH",
    title="uH-t for Different J Values",
    # ylim=(-2, 20),
    save_as=os.path.join(results_folder, f"figure_uH_vs_t_for_different_J-{maxLevel}.png")
)

# Plot vR-t
labels_vR = {J: (data[J]["t"], data[J]["vR"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_vR,
    xlabel="t", ylabel="vR",
    title="vR-t for Different J Values",
    # ylim=(-2, 20),
    save_as=os.path.join(results_folder, f"figure_vR_vs_t_for_different_J-{maxLevel}.png")
)
# Plot ke-t
labels_ke = {J: (data[J]["t"], data[J]["ke"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_ke,
    xlabel="t", ylabel="ke",
    title="ke-t for Different J Values",
    save_as=os.path.join(results_folder, f"figure_ke_vs_t_for_different_J-{maxLevel}.png")
)
# Plot Rmax-t
labels_Rmax = {J: (data[J]["t"], data[J]["Rmax"]) for J in data}
plot_figure(
    x=None, y=None, labels=labels_Rmax,
    xlabel="t", ylabel="Rmax",
    title="Rmax-t for Different J Values",
    # ylim=(-2, 20),
    save_as=os.path.join(results_folder, f"figure_Rmax_vs_t_for_different_J-{maxLevel}.png")
)
####################################################################################################
# Plot F1, F2 vs J
F1_values = [data[J]["F1_t"][0] for J in data]
F2_values = [data[J]["F2_t"][0] for J in data]
plt.figure(figsize=(8, 6))
plt.plot(J_values, F1_values, 'o-', label='F1', color='blue')
plt.plot(J_values, F2_values, 's-', label='F2', color='red')
plt.xlabel('J')
plt.ylabel('F')
plt.xscale("log")
plt.title('F1 and F2 vs J')
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(results_folder, f"figure_All_F1_F2_vs_J-{maxLevel}.png"), dpi=300)
plt.close()

# Plot uHMAX vs J
uHMAX_values = [data[J]["uHMAX_t"][0] for J in data]
plt.figure(figsize=(8, 6))
plt.plot(J_values, uHMAX_values, '^-', label='uHMAX', color='green')
plt.xlabel('J')
plt.ylabel('uHMAX')
plt.title('uHMAX vs J')
plt.xscale("log")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(results_folder, f"figure_All_uHMAX_vs_J-{maxLevel}.png"), dpi=300)
plt.close()

# Plot uHMAX vs J
vRMAX_values = [data[J]["vRMAX_t"][0] for J in data]
plt.figure(figsize=(8, 6))
plt.plot(J_values, vRMAX_values, '^-', label='vRMAX', color='green')
plt.xlabel('J')
plt.ylabel('vRMAX')
plt.title('vRMAX vs J')
plt.xscale("log")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(results_folder, f"figure_All_vRMAX_vs_J-{maxLevel}.png"), dpi=300)
plt.close()

RmaxMAX_values = [data[J]["RmaxMAX_t"][0] for J in data]
plt.figure(figsize=(8, 6))
plt.plot(J_values, RmaxMAX_values, '^-', label='Rmax', color='green')
plt.xlabel('J')
plt.ylabel('Rmax')
plt.title('RmaxMAX vs J')
plt.xscale("log")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(results_folder, f"figure_All_RmaxMAX_vs_J-{maxLevel}.png"), dpi=300)
plt.close()

# Plot F1_time, F2_time, uHMAX_t vs J
F1_times = [data[J]["F1_t"][1] for J in data]
F2_times = [data[J]["F2_t"][1] for J in data]
uHMAX_times = [data[J]["uHMAX_t"][1] for J in data]
vRMAX_times = [data[J]["vRMAX_t"][1] for J in data]
RmaxMAX_times = [data[J]["RmaxMAX_t"][1] for J in data]

fig, ax = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.05})

# Top part of the plot
ax[0].plot(J_values, F2_times, 's-', label='F2 time', color='red')
ax[0].plot(J_values, uHMAX_times, '^-', label='uHMAX time', color='green')
ax[0].set_ylim(3.2, 3.7)
ax[0].set_xticklabels([])
# ax[0].legend()
ax[0].grid(True)

# Bottom part of the plot
ax[1].plot(J_values, F1_times, 'o-', label='F1 time', color='blue')
ax[1].plot(J_values, vRMAX_times, 'o-', label='vRMAX time', color='yellow')
ax[1].plot(J_values, RmaxMAX_times, 'o-', label='RmaxMAX time', color='yellow')
ax[1].set_ylim(0.3, 0.5)
# ax[1].legend()
ax[1].grid(True)
ax[1].set_xlabel('J')

# Add diagonal lines to show break in y-axis
d = 0.015  # Size of diagonal line
kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
ax[0].plot((-d, +d), (-d, +d), **kwargs)  # Top-left diagonal
ax[0].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # Top-right diagonal

kwargs.update(transform=ax[1].transAxes)
ax[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # Bottom-left diagonal
ax[1].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # Bottom-right diagonal

fig.legend(
    loc='upper center', 
    bbox_to_anchor=(0.5, 0.98), 
    ncol=3, 
    frameon=False
)
# Save the plot
fig.supylabel('Time')
# fig.suptitle('F1_time, F2_time, and uHMAX_time vs J')
plt.savefig(os.path.join(results_folder, f"figure_All_times_vs_J-{maxLevel}.png"), dpi=300)
plt.close()