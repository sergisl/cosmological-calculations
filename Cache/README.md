# xAlpha Cache Directory

This directory contains pre-computed `.mx` binary files for expensive xAct/xPand
calculations (second-order metric and scalar EOMs in Newtonian gauge).

Pre-shipping these files means external users load the results instantly rather than
waiting for the full tensor computation (which can take several minutes).

## Files

| File | Symbol | Description |
|---|---|---|
| `metricEOMsQuad.mx` | `metricEOMsQuad` | Second-order metric EOM, Newtonian gauge |
| `scalarEOMQuad.mx` | `scalarEOMQuad` | Second-order scalar EOM, Newtonian gauge |

## Usage

Use the package utilities in your notebook:

```mathematica
(* Load from cache, or compute and save if not present *)
metricEOMsQuad = xAlphaLoadOrCompute["metricEOMsQuad",
  -2 * MyToxPand[metricEOMs /. Xto\[CurlyPhi], "NewtonGauge", 2] /. HorndeskiDerivativeReplacements2
];

scalarEOMQuad = xAlphaLoadOrCompute["scalarEOMQuad",
  -2 * MyToxPand[scalarEOM /. Xto\[CurlyPhi], "NewtonGauge", 2] /. HorndeskiDerivativeReplacements2
];
```

## Notes

- `.mx` files are Mathematica binary format — they are version and platform independent
  for *expressions*, but require the same xAct symbol contexts to be loaded first.
- Always load Setup.wl and call `SetupFluidModel` before loading these files.
- To regenerate: delete the `.mx` file and re-run the notebook cell.
