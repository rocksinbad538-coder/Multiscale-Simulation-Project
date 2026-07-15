from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_exciton_model")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514  # eV*fs

H = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
])

labels = ["PYR2","PYR3","PYR4","PYR5"]

gammas_ps = [1,5,10,50]

def rhs(t,y,H,gamma_fs):
    n=4
    rho=y.reshape((n,n)).astype(complex)

    comm = H@rho - rho@H

    drho = -1j/HBAR*comm

    for i in range(n):
        for j in range(n):
            if i!=j:
                drho[i,j] -= gamma_fs*rho[i,j]

    return drho.reshape(-1)

for gamma_ps in gammas_ps:

    gamma_fs = gamma_ps/1000.0

    rho0 = np.zeros((4,4),dtype=complex)

    # initial excitation on PYR5
    rho0[3,3]=1.0

    t_eval=np.linspace(0,10000,1001) # fs = 10 ps

    sol=solve_ivp(
        rhs,
        [0,10000],
        rho0.reshape(-1),
        t_eval=t_eval,
        args=(H,gamma_fs),
        rtol=1e-8,
        atol=1e-10
    )

    pops=[]

    for k in range(len(sol.t)):
        rho=sol.y[:,k].reshape((4,4))
        pops.append(np.real(np.diag(rho)))

    pops=np.array(pops)

    df=pd.DataFrame({
        "time_fs":sol.t,
        "PYR2":pops[:,0],
        "PYR3":pops[:,1],
        "PYR4":pops[:,2],
        "PYR5":pops[:,3],
    })

    csvfile=OUT/f"haken_strobl_gamma{gamma_ps}ps.csv"
    df.to_csv(csvfile,index=False)

    plt.figure(figsize=(6,4))
    for i,lbl in enumerate(labels):
        plt.plot(sol.t/1000,pops[:,i],label=lbl)

    plt.xlabel("time (ps)")
    plt.ylabel("population")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        OUT/f"haken_strobl_gamma{gamma_ps}ps.png",
        dpi=300
    )
    plt.close()

print("Finished Haken-Strobl runs.")
