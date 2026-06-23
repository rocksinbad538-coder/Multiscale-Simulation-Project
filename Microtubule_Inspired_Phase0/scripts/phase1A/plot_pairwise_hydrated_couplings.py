from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "runs/phase1A/day014_pairwise_summary/hydrated_pairwise_results.csv"
)

plt.figure(figsize=(6,4))
plt.bar(df["pair"], df["Jeff_meV"])
plt.ylabel("Effective hydrated coupling (meV)")
plt.tight_layout()

out = Path(
    "runs/phase1A/day014_pairwise_summary/hydrated_couplings_summary.png"
)
plt.savefig(out, dpi=300)

print(out)
