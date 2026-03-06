#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot tau from tau_band_xxx.dat produced by your Fortran tau_calc output:

Expected format:
  Comment lines begin with '#'
  One line:
    # global_iband = <int>
  One line:
    # T_list(K): 30.000 40.000 ... 300.000
  Data lines:
    ikx iky ikz  E(eV)  DOS(1/eV)  tau_ps(T1) tau_ps(T2) ... tau_ps(TnT)

Features:
1) tau(E) at fixed temperature T_fix (choose nearest T in T_list)
2) tau(T) at fixed k-point (or fixed E by nearest row), plotted on log-log axes
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot tau(E) at fixed T and tau(T) at fixed k from tau_band_xxx.dat"
    )
    p.add_argument("--input", type=Path, required=True, help="tau_band_xxx.dat file")
    p.add_argument("--T-fix", type=float, default=60.0, help="Temperature (K) for tau(E) plot")
    p.add_argument("--E-fix", type=float, default=None, help="Energy (eV) to choose nearest k-row for tau(T) plot")

    # Choose a specific k row by integer indices
    p.add_argument("--k", type=int, nargs=3, default=None, metavar=("IKX", "IKY", "IKZ"),
                   help="Choose exact (ikx iky ikz) row for tau(T) plot")

    p.add_argument("--xlim-E", type=float, nargs=2, default=None,
                   help="x-range for tau(E) plot in eV, e.g. --xlim-E -0.1 0.1")
    p.add_argument("--ylim-tauE", type=float, nargs=2, default=None,
                   help="y-range for tau(E) plot (ps), e.g. --ylim-tauE 0 500")

    p.add_argument("--tauE-unit", choices=["ps", "1e-11s"], default="ps",
                   help="y unit for tau(E): ps or tau/(1e-11 s) (dimensionless)")

    p.add_argument("--out-prefix", type=str, default=None,
               help="Output prefix. If omitted, use band<iband> (e.g. band5).")
    p.add_argument("--show", action="store_true", help="Show figures interactively")

    return p.parse_args()


def _try_parse_global_iband(line: str) -> Optional[int]:
    # e.g. "# global_iband =          5"
    if "global_iband" not in line:
        return None
    # take last token that can be int
    toks = line.replace("=", " ").split()
    for t in toks[::-1]:
        try:
            return int(t)
        except ValueError:
            pass
    return None


def _parse_T_list(line: str) -> Optional[np.ndarray]:
    # e.g. "# T_list(K): 30.000 40.000 ..."
    if "T_list" not in line:
        return None
    if ":" not in line:
        return None
    tail = line.split(":", 1)[1].strip()
    arr = np.fromstring(tail, sep=" ", dtype=float)
    if arr.size == 0:
        return None
    return arr


def load_tau_band_table(path: Path) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns:
      iband_global: int (or -1 if missing)
      Tlist: (nT,)
      ik_int: (N,3) int
      E_eV: (N,)
      DOS_1eV: (N,)
      tau_ps: (N,nT)
    """
    iband_global = -1
    Tlist: Optional[np.ndarray] = None

    rows = []
    with path.open("r", encoding="ascii", errors="ignore") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("#"):
                if iband_global < 0:
                    ib = _try_parse_global_iband(s)
                    if ib is not None:
                        iband_global = ib
                if Tlist is None:
                    tarr = _parse_T_list(s)
                    if tarr is not None:
                        Tlist = tarr
                continue

            # data line
            vals = np.fromstring(s, sep=" ", dtype=float)
            if vals.size < 6:
                continue
            rows.append(vals)

    if Tlist is None:
        raise ValueError(f"Failed to find T_list(K) in header of {path}")

    nT = Tlist.size
    data = np.vstack(rows) if rows else np.empty((0, 6 + nT), dtype=float)
    if data.shape[0] == 0:
        raise ValueError(f"No data rows found in {path}")

    # Expect columns: ikx iky ikz E DOS tau(T1..Tn)
    if data.shape[1] < 5 + nT:
        raise ValueError(
            f"Data columns ({data.shape[1]}) < expected (5 + nT = {5+nT}). "
            f"Check file format or T_list length."
        )

    ik_int = data[:, 0:3].astype(int)
    E_eV = data[:, 3].astype(float)
    DOS_1eV = data[:, 4].astype(float)
    tau_ps = data[:, 5:5 + nT].astype(float)

    return iband_global, Tlist, ik_int, E_eV, DOS_1eV, tau_ps


def find_closest_T_index(Tlist: np.ndarray, T_fix: float) -> int:
    return int(np.argmin(np.abs(Tlist - T_fix)))


def find_row_by_k(ik_int: np.ndarray, k_tuple: Tuple[int, int, int]) -> int:
    m = np.all(ik_int == np.array(k_tuple, dtype=int)[None, :], axis=1)
    idxs = np.where(m)[0]
    if idxs.size == 0:
        raise ValueError(f"k-point {k_tuple} not found in file.")
    return int(idxs[0])


def find_row_by_E(E_eV: np.ndarray, E_fix: float) -> int:
    return int(np.argmin(np.abs(E_eV - E_fix)))


def main() -> None:
    args = parse_args()

    iband, Tlist, ik_int, E_eV, DOS_1eV, tau_ps = load_tau_band_table(args.input)
    prefix = args.out_prefix if args.out_prefix is not None else f"band{iband}"
    nT = Tlist.size
    N = E_eV.size

    print(f"Loaded: {args.input}")
    print(f"  global_iband(n) = {iband}")
    print(f"  N rows = {N}, nT = {nT}, T range = [{Tlist.min():.1f}, {Tlist.max():.1f}] K")
    print(f"  E range = [{E_eV.min():.6f}, {E_eV.max():.6f}] eV")

    # ---------------------------
    # Plot 1: tau(E) at fixed T
    # ---------------------------
    tidx = find_closest_T_index(Tlist, args.T_fix)
    T_used = Tlist[tidx]
    tauE_ps = tau_ps[:, tidx]

    # Optionally convert y unit
    if args.tauE_unit == "ps":
        y_tauE = tauE_ps
        ylab = r'$\tau$ (ps)'
    else:
        # tau/(1e-11 s) dimensionless; tau_ps * 1e-12 / 1e-11 = tau_ps * 0.1
        y_tauE = tauE_ps * 0.1
        ylab = r'$\tau$ / $(10^{-11}\,\mathrm{s})$'

    # Sort by energy for a clean curve
    order = np.argsort(E_eV)
    E_sorted = E_eV[order]
    y_sorted = y_tauE[order]

    plt.figure(figsize=(6.2, 4.2))
    plt.plot(E_sorted * 1e3, y_sorted)  # meV on x
    plt.xlabel(r'Energy $E$ (meV)')
    plt.ylabel(ylab)
    plt.title(f"n={iband}  tau(E) at T={T_used:.1f} K")

    if args.xlim_E is not None:
        plt.xlim(args.xlim_E[0] * 1e3, args.xlim_E[1] * 1e3)
    if args.ylim_tauE is not None:
        plt.ylim(args.ylim_tauE[0], args.ylim_tauE[1])

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out1 = f"{prefix}_tau_vs_E_T{T_used:.0f}K.png"
    plt.savefig(out1, dpi=300)
    print(f"Saved {out1}")

    # ---------------------------------------
    # Plot 2: tau(T) at fixed k (log-log)
    # ---------------------------------------
    if args.k is not None:
        row_idx = find_row_by_k(ik_int, tuple(args.k))
        choose_desc = f"k={tuple(args.k)}"
    elif args.E_fix is not None:
        row_idx = find_row_by_E(E_eV, float(args.E_fix))
        choose_desc = f"E≈{E_eV[row_idx]:.6f} eV (nearest)"
    else:
        # default: choose the row with energy closest to 0 eV
        row_idx = find_row_by_E(E_eV, 0.0)
        choose_desc = f"E≈{E_eV[row_idx]:.6f} eV (nearest to 0)"

    tauT_ps = tau_ps[row_idx, :]           # (nT,)
    tauT_s = tauT_ps * 1e-12               # seconds

    plt.figure(figsize=(6.2, 4.2))
    plt.plot(Tlist, tauT_s, marker="o", markersize=3, linewidth=1.2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Temperature (K)")
    plt.ylabel(r'$\tau$ (s)')
    plt.title(f"n={iband}  tau(T) at {choose_desc}")

    # y-range a bit padded but safe
    pos = tauT_s[tauT_s > 0]
    if pos.size > 0:
        ymin = max(pos.min() * 0.6, 1e-18)
        ymax = pos.max() * 1.6
        plt.ylim(ymin, ymax)

    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    out2 = f"{prefix}_tau_vs_T_loglog.png"
    plt.savefig(out2, dpi=300)
    print(f"Saved {out2} (row {row_idx}, {choose_desc})")

    if args.show:
        plt.show()


if __name__ == "__main__":
    # 用法示例
    # 固定温度60K画 tau-E， 固定 k 点 1 1 82 画 tau-T（loglog），选择的band编号为5
    # python plot_tau_band.py --input tau_band_5.dat --T-fix 60  --k 1 1 82 

    # 固定温度60K画 tau-E， 选定能量0.092eV,不知道 k 点，就用能量选最近的一行，画 tau-T（loglog），选择的band编号为5
    # python plot_tau_band.py --input tau_band_5.dat --T-fix 60 --E-fix 0.092 

    main()