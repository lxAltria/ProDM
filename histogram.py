import numpy as np
from matplotlib import pyplot as plt

def to_weights(weights, max_weight=4):
    max_val = np.max(np.abs(weights))

    abs_weights = np.abs(weights)
    log_weights = np.log2(max_val / abs_weights, where=abs_weights != 0)
    log_weights = np.floor(log_weights).astype(np.int32)

    int_weights = np.where(abs_weights == 0, 0, max_weight - log_weights)
    int_weights = np.where(log_weights > max_weight, 0, int_weights)
            
    return int_weights

def weight_hist(weights):
    arr = [0, 0, 0, 0, 0]
    for weight in weights:
        arr[weight] += 1
    for i in range(5):
        print(f"weight = {i}: {arr[i]}")


# import your data
vx = np.fromfile("./Hurricane_f32/data/VelocityX.dat", dtype=np.float32)
vy = np.fromfile("./Hurricane_f32/data/VelocityY.dat", dtype=np.float32)
vz = np.fromfile("./Hurricane_f32/data/VelocityZ.dat", dtype=np.float32)
vtot = (vx**2 + vy**2 + vz**2)

vtot_5 = np.where(vtot > 0, 1 / vtot**0.5, 0)
vtot_4 = np.where(vtot > 0, 1 / vtot**0.4, 0)
vtot_3 = np.where(vtot > 0, 1 / vtot**0.3, 0)
vtot_2 = np.where(vtot > 0, 1 / vtot**0.2, 0)
vtot_1 = np.where(vtot > 0, 1 / vtot**0.1, 0)

vtots = [vtot_5, vtot_4, vtot_3, vtot_2, vtot_1]
labels = ["$v_{tot}^{0.5}$", "$v_{tot}^{0.4}$", "$v_{tot}^{0.3}$", "$v_{tot}^{0.2}$", "$v_{tot}^{0.1}$"]
colors = ['yellow', 'blue', 'green', 'purple', 'orange']

fig, axes = plt.subplots(5, 2, figsize=(15, 20), sharex='col', sharey=False)

for i, (vtot, label, color) in enumerate(zip(vtots, labels, colors)):
    # histogram of vtot
    axes[i, 0].hist(vtot, bins=200, range=(0,2), color=color, alpha=0.7, edgecolor='black')
    axes[i, 0].set_title(f"Histogram of {label} (Vtot)", fontsize=14)
    axes[i, 0].set_ylabel("Frequency", fontsize=12)
    axes[i, 0].grid(True, linestyle='--', alpha=0.6)

    # histogram of weight
    weights = to_weights(vtot)
    print(f"weights of vtot_{5 - i}")
    weight_hist(weights)
    axes[i, 1].hist(weights, bins=5, color=color, alpha=0.7, edgecolor='black')
    axes[i, 1].set_title(f"Histogram of {label} (Weights)", fontsize=14)
    axes[i, 1].grid(True, linestyle='--', alpha=0.6)

for ax in axes[-1, :]:
    ax.set_xlabel("$v_{tot}$ or Weights", fontsize=12)

plt.tight_layout()
plt.show()