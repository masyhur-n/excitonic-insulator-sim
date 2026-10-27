# excitonic-insulator-sim

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23684703.svg)](https://doi.org/10.5281/zenodo.23684703)

Code and data for the article

> **N. Masyhur, "Can excitonic-insulator switches beat CMOS? A benchmarked feasibility screen across the BEC–BCS crossover"** (manuscript in preparation; source in [`paper/`](paper/)).

Excitonic insulators and exciton condensates have been proposed as switches that beat the
60 mV/decade thermionic limit of CMOS. This repository tests that proposal with a chain of
models. Each ingredient first reproduces a published number (`run_benchmarks.py`). The models
are then used to evaluate subthreshold swing, condensation temperature and switching energy
at 300–350 K.

**Main results**
- Screening-limited ceiling: T_c ≤ Ry*/(4 r_s,c²). No hBN-encapsulated double layer reaches 300 K.
- Graphene double layers at the BiSFET operating point do not condense under normal-state or superfluid-state screening.
- Ta₂NiSe₅: the order is about 70% excitonic but needs the lattice. Only discrete symmetries are broken, so there is no superfluid.
- The only sub-CMOS switching (an abrupt mean-field onset) disappears with Gaussian fluctuations: S = 100–470 mV/dec, and ≥ 285 mV/dec with gate electrostatics.

## Reproducing the paper

```bash
pip install -r requirements.txt
python run_benchmarks.py          # Table I: 15 pass, 2 documented deviations (~5 min; --fast skips 2 scans)
python paper/make_figures.py      # Figs. 1–4 -> paper/figs/
python screening.py               # superfluid-onset line r_s,c(d/a_B), Fig. 2a (~8 min)
python tc_ceiling.py              # screening ceiling, Eq. (3)
python dirac.py                   # graphene / BiSFET operating point
python ta2nise5_chain.py          # Ta2NiSe5 three-chain model (~3 min)
python device_scoring.py 350      # mean-field device screen
python fluctuations.py 350        # Gaussian-fluctuation switching test (also: 300)
python real_gate.py               # gate electrostatics, 300 K
python disorder.py; python energy.py
python ml_v2.py 3000              # surrogate with uncertainty + inverse design (~2 min)
```

| manuscript item | script | notes |
|---|---|---|
| Table I | `run_benchmarks.py` | each row lists the literature value and source |
| Fig. 1, Sec. III A | `paper/make_figures.py` → `ei_model.py` | 2D crossover |
| Fig. 2, Eq. (3), Sec. III B | `screening.py`, `tc_ceiling.py` | [docs/SCREENING.md](docs/SCREENING.md) |
| Sec. III C | `dirac.py` | [docs/DIRAC.md](docs/DIRAC.md) |
| Fig. 4, Sec. III D | `ta2nise5_chain.py` | [docs/TA2NISE5.md](docs/TA2NISE5.md) |
| Fig. 3, Sec. III E | `device_scoring.py`, `fluctuations.py`, `real_gate.py`, `energy.py`, `disorder.py` | [docs/FLUCTUATIONS.md](docs/FLUCTUATIONS.md), [docs/GATE_MASS_DISORDER.md](docs/GATE_MASS_DISORDER.md), [docs/DEVICE_SCORING.md](docs/DEVICE_SCORING.md) |
| Sec. III F | `ml_v2.py` | `ei_dataset_v2.csv` |

## Repository layout

| path | content |
|---|---|
| `ei_model.py` | 2D e–h gap equation (Rytova–Keldysh, layer separation), T*, BKT bound |
| `benchmark_meanfield.py`, `test_model.py` | crossover mean-field benchmarks (Randeria & Taylor 2014), sanity checks |
| `screening.py`, `tc_ceiling.py` | fixed-density bilayer, superfluid-state screening (Neilson et al. 2014), T_c ceiling |
| `dirac.py` | Dirac double layer with US / NS / SS screening |
| `ta2nise5_chain.py` | Ta₂NiSe₅ three-chain model with exciton–phonon coupling (Kaneko et al. 2013) |
| `fluctuations.py` | Gaussian (self-consistent Hartree) order-parameter fluctuations |
| `real_gate.py`, `disorder.py`, `energy.py`, `device_scoring.py` | device metrics |
| `ml_v2.py`, `ei_dataset_v2.csv` | ML surrogate (bootstrap ensemble) and its training data |
| `run_benchmarks.py` | regenerates Table I |
| `paper/` | manuscript (REVTeX), bibliography, figures and figure script |
| `docs/` | method notes, results summary ([RESULTS_SUMMARY.md](docs/RESULTS_SUMMARY.md)) and the reference verification log ([REFERENCES.md](docs/REFERENCES.md)) |
| `results/` | figures written by the scripts |
| `materials.py`, `run_crossover.py`, `make_dataset.py`, `train_surrogate.py`, `ei_dataset.csv` | first-generation screen with illustrative parameters (superseded by `ml_v2.py`; kept for transparency) |

## Model summary

```
Δ(k) = Σ_k' V(|k−k'|) Δ(k') / 2E(k') · tanh(E(k')/2k_BT),   E = sqrt(ξ² + Δ²),   ξ = Eg/2 + ħ²k²/4m_r
V(q) = 2πe² e^(−qd) / [ε q (1 + r0 q)]        (unscreened; superfluid-state screening in screening.py)
```

| regime | condition | signature |
|---|---|---|
| BCS | Eg < 0 (semimetal) | Δ(k) peaks at k_F; v_k² has a Fermi edge |
| BEC | 0 < Eg < E_b | Δ(k) peaks at k = 0; tightly bound pairs |
| band insulator | Eg > E_b | Δ = 0 |

Numerics: radial Nyström discretisation, an analytic angular Coulomb integral
(elliptic K, plus E for the Dirac chirality factor), and cell-averaging of the
integrable logarithmic singularity.

## Limitations
- Static mean field plus Gaussian fluctuations. No vortex dynamics and no vertex corrections
  beyond those tested by Neilson et al. (2014).
- No intralayer exchange in the continuum model. It is bracketed using Zeng & MacDonald (2020).
- Two Ta₂NiSe₅ numbers deviate by 15–28% from Kaneko et al. (2013). Their supplementary
  equations were not used.
- No transport, contact or circuit modelling.

The verification status of every reference is in [docs/REFERENCES.md](docs/REFERENCES.md).

## AI assistance
Code development, literature retrieval and drafting were assisted by an AI coding assistant
(Claude, Anthropic). This is also disclosed in the manuscript. All results are checked
against the published benchmarks in `run_benchmarks.py`.

## License
MIT (see [LICENSE](LICENSE)).

## Citation
Software: N. Masyhur, *excitonic-insulator-sim*, Zenodo, doi:[10.5281/zenodo.23684703](https://doi.org/10.5281/zenodo.23684703) (all versions; v1.0.0: doi:10.5281/zenodo.23684704). See [CITATION.cff](CITATION.cff).
