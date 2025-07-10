import pandas as pd
import matplotlib.pyplot as plt

# Custom colors for each phase
phase_colors = {
    "Fe2O3_0": "#FFADAF",      # Pink
    "MoO3_0": "#89CEDF",       # Light blue
    "Fe2MoO43_0": "#F99527",   # Orange
    "liquid_0": "#B7A2CE"      # Purple
}

# Load the CSV data
df = pd.read_csv("phase_fractions.csv")

# Plot each phase with its color
for phase in df.columns[1:]:
    plt.plot(df["Time"], df[phase], label=phase, color=phase_colors.get(phase, None))

plt.xlabel("Time")
plt.ylabel("Phase Fraction")
plt.title("Phase Fractions vs Time")
plt.legend()
plt.tight_layout()
plt.savefig("phase_fractions_plot.png", dpi=150)
plt.show()