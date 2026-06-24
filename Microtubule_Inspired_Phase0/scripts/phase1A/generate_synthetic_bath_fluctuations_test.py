from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("runs/phase1A/day015_bath_parameter_mapping/synthetic_test")
OUT.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(123)

time_ps = np.linspace(0, 20, 2001)
dt = time_ps[1] - time_ps[0]

def ou_process(mean, sigma, tau_ps):
    x = np.zeros_like(time_ps)
    x[0] = mean
    a = np.exp(-dt / tau_ps)
    noise = sigma * np.sqrt(1 - a*a)
    for i in range(1, len(time_ps)):
        x[i] = mean + a * (x[i-1] - mean) + noise * rng.normal()
    return x

df = pd.DataFrame({
    "time_ps": time_ps,
    "PYR2": ou_process(3.7550, 0.005, 0.20),
    "PYR3": ou_process(3.7775, 0.005, 0.15),
    "PYR4": ou_process(3.7775, 0.006, 0.10),
    "PYR5": ou_process(3.7775, 0.006, 0.10),
    "J23":  ou_process(0.00237, 0.0005, 0.20),
    "J34":  ou_process(0.01700, 0.0015, 0.15),
    "J45":  ou_process(0.03300, 0.0020, 0.10),
})

path = OUT / "synthetic_md_fluctuations.csv"
df.to_csv(path, index=False)
print(path)
