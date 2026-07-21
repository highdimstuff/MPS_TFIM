# Tensor Network: A Study of Matrix Product States (MPS)

This repository contains the internship report and numerical implementation code developed during my summer research internship at **IIT (ISM) Dhanbad** under the guidance of **Prof. Sudipto Singha Roy**[cite: 1]. 

The project investigates the theoretical foundations of Tensor Networks and applies **Matrix Product States (MPS)** to study quantum phase transitions in the **Transverse Field Ising Model (TFIM)** using the **ITensor** library in **Julia**[cite: 1].

---

## 📂 Repository Contents

* `denison_report.pdf`: The complete formal internship report covering theoretical derivations (area law, canonical forms, DMRG, imaginary time evolution) and numerical findings[cite: 1].
* `tfim_benchmark.jl`: A clean, commented Julia script utilizing ITensor to implement DMRG and benchmark the 1D TFIM[cite: 1].

---

## 🚀 Project Overview & Key Findings

1. **Theoretical Formalism:** Reviewed the exponential growth of many-body Hilbert spaces and how the 1D area law of entanglement allows efficient simulation via MPS[cite: 1].
2. **Numerical Implementation:** Benchmarked a finite spin-$1/2$ chain ($N = 30$) across the paramagnetic and ferromagnetic regimes of the TFIM[cite: 1].
3. **Addressing Degeneracy (GHZ States):** Observed that sweeping into the ferromagnetic region from the paramagnetic side via adiabatic methods yields a zero-magnetization, high-entanglement state corresponding to a $\mathbb{Z}_2$ symmetry-preserving GHZ-like superposition[cite: 1]. Successfully resolved this by introducing a weak longitudinal symmetry-breaking perturbation ($\epsilon = 10^{-6}$) to recover the expected classical phase transition graph[cite: 1].

---

## 🛠️ Requirements & Tools

* **Language:** Julia[cite: 1]
* **Libraries:** `ITensors.jl`, `ITensorMPS.jl`[cite: 1]

---

## 👤 Author

* **Moirangthem Denison Meitei**  
  M.S. Quantum Technology, IISER Pune[cite: 1]
