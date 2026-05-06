# xAlpha

**xAlpha** is an [xAct](https://www.xact.es/) package for computing EFT alpha parameters ($\alpha_M$, $\alpha_K$, $\alpha_B$ and their derivatives) in luminal Horndeski gravity, and for systematically extracting and converting cosmological perturbation coefficients into EFT language.

## What it does

Given a Horndeski theory defined by functions $K(\phi, X)$, $G_3(\phi, X)$, $G_4(\phi)$, xAlpha:

- Computes the EFT alpha parameters $(\alpha_M, \alpha_K, \alpha_B)$ and all their derivatives with respect to time, $X$, and $\phi$ up to third order
- Replaces Horndeski $G$-function derivatives in perturbation equations with alpha parameters (`GToAlphas`)
- Extracts and stores the linear and quadratic perturbation coefficients from the $(0,0)$, $(0,i)$, $(i,j)$, and scalar field equations
- Converts coefficients between raw Horndeski, alpha-language, and gamma-language EFT bases
- Validates that extracted coefficients fully reconstruct the original equations

## Dependencies

- [Mathematica](https://www.wolfram.com/mathematica/) ≥ 13.0
- [xAct](https://www.xact.es/) ≥ 1.3.0: `xCore`, `xTensor`, `xCoba`, `xTras`, `xPand`

## Package structure

```
xAlpha/
├── Kernel/
│   └── init.m                  ← enables << xAct`xAlpha` syntax
├── xAlpha.m                    ← main loader: BeginPackage, loads sub-files
├── Setup.wl                    ← tensor/function definitions (load AFTER manifold)
├── HorndeskiDefinitions.wl     ← Horndeski derivative replacement rules
├── Background.wl               ← background EOM rules, cached rule accessor
├── Perturbations.wl            ← GetRemoveScalarWrapper, perturbation helpers
├── AlphaFunctions.wl           ← core EFT alpha computation and extraction
└── README.md
```

> **Note on `Setup.wl`:** It contains `DefTensor`, `DefScalarFunction`, etc. calls that require the xAct manifold, metric, and xPand slicing to already be defined. It is therefore **not** loaded automatically by `<< xAct\`xAlpha\``. You must load it manually in your notebook after setting up the geometry (see Quick Start below).

## Installation

Place the `xAlpha/` directory inside your xAct installation:

```
$UserBaseDirectory/Applications/xAct/xAlpha/
```

After that, `<< xAct\`xAlpha\`` will work from any notebook.

Alternatively, load directly by path:
```mathematica
Get["/path/to/xAlpha/xAlpha.m"];
```

## Quick start

```mathematica
(* 1. Load xAct suite *)
<< xAct`xCoba`;
<< xAct`xTras`;
<< xAct`xPand`;

(* 2. Define manifold and xPand slicing *)
DefManifold[Md, 4, {b, c, d, e, f}];
DefMetric[-1, g[-b, -c], CD, {";", "\[Del]"}];
$ConformalTime = False;
$CommuteCovDsOnScalars = True;
SetSlicing[g, n, h, cd, {"|", "\!\(\*OverscriptBox[\(D\), \(_\)]\)"}, "FLFlat"];
DefMetricFields[g, dg, h];
DefMatterFields[u, du, h];

(* 3. Load xAlpha core package *)
<< xAct`xAlpha`;          (* or: Get["/path/to/xAlpha/xAlpha.m"] *)

(* 4. Load Setup.wl AFTER manifold is defined *)
Get[FileNameJoin[{$xAlphaDir, "Setup.wl"}]];

(* 5. Initialise fluid model and helpers *)
SetupFluidModel["PerfectFluid"];
HorndeskiDerivativeReplacements1 = GetHorndeskiRules[];
XtoPhi = MakeXtoPhiRule[];

(* 6. Compute all EFT alpha parameters and their derivatives *)
ComputeAllAlphas[\[CurlyPhi], X, Hh, HorndeskiDerivativeReplacements1, phidotToX]
(* Populates: Mstarexpr, aMexpr, aKexpr, aBexpr,
              aBdotexpr, aKdotexpr, aMdotexpr, ...,
              aBXexpr, aKXexpr, aBXdotexpr, ..., etc. *)

(* 7. Transform a perturbation equation from G-language to alpha language *)
GToAlphas[myPertExpr, phidotToX]

(* 8. Extract coefficients from a collected equation *)
ExtractAndTransformCoefficients[collectedEq00, pertHeads, phidotToX, removeScalarWrapper]
(* Stores: A[i,h], Aalpha[i,h], Agamma[i,h]  for i=1,2,3 *)

(* 9. Display results as a formatted table *)
ShowPaperTables[Aalpha, Balpha, CCalpha, DDalpha,
                Aalpha, Balpha, CCalpha, DDalpha]
```

## Key user-facing functions

| Function | Description |
|---|---|
| `ComputeAllAlphas[phi, X, Hh, rules, phidotToX]` | Computes and stores all alpha parameters and derivatives |
| `GToAlphas[expr, phidotToX]` | Replaces all Horndeski derivatives with alpha parameters |
| `RemoveK[expr, phidotToX]` | Substitutes K-function derivatives only |
| `RemoveG3[expr, phidotToX]` | Substitutes G3-function derivatives only |
| `RemoveG4[expr, phidotToX]` | Substitutes G4-function derivatives only |
| `ExtractAndTransformCoefficients[...]` | Extracts linear A (00) coefficients |
| `ExtractAndTransformBCoefficients[...]` | Extracts linear B (0i) coefficients |
| `ExtractAndTransformCCoefficients[...]` | Extracts linear C (ij) coefficients |
| `ExtractAndTransformDCoefficients[...]` | Extracts linear D (scalar) coefficients |
| `ExtractAndTransformQuadratic*[...]` | Quadratic variants of the above |
| `ConvertAllCoefficientsToGammaLanguage[pertHeads]` | Converts all alpha→gamma basis |
| `ApplyConservationEOMRules[]` | Rules from background conservation equations |
| `ApplyEnergyRescalingRules[]` | Rules rescaling varE/varP by Mstar² |
| `ShowPaperTables[...]` | Formatted table of all extracted coefficients |
| `ValidateAllCoefficients[...]` | Checks coefficients reconstruct original equations |
| `MyToxPand[expr, gauge, order]` | Wrapper around xPand's ToxPand |
| `SetupFluidModel["PerfectFluid"]` | Defines stress-energy tensor for chosen fluid |
| `MakeXtoPhiRule[]` | Rule X[] → scalar gradient squared |
| `MakePhidotToXRules[]` | Rules φ'ⁿ → X[]^(n/2) |

Full documentation is available via `?FunctionName` in Mathematica after loading the package.

## EFT bases

**Alpha language** ($\alpha_K, \alpha_B, \alpha_M$ and X/φ/time derivatives): the standard Bellini–Sawicki parametrisation of Horndeski gravity on a cosmological background.

**Gamma language** ($\gamma_K, \gamma_B, \gamma_D, \gamma_M, \ldots$): an alternative algebraic basis related to alpha parameters by the equations in `GenerateAlphaToGammaRules[]`.

## Notes on `KK`

The Horndeski kinetic function $K(\phi, X)$ is defined as `KK` internally (and printed as $K$) to avoid a clash with Mathematica's built-in elliptic integral `K[m]`.

## Authors

Sergi Sirera

## Citation

If you use xAlpha in your research, please cite this repository and the xAct suite:
> J. M. Martín-García et al., *xAct: Efficient tensor computer algebra for the Wolfram Language*, https://www.xact.es/

## License

GPL-2.0 — see [LICENSE](LICENSE).

## What it does

Given a Horndeski theory defined by functions $K(\phi, X)$, $G_3(\phi, X)$, $G_4(\phi)$, xAlpha:

- Computes background equations of motion on a FLRW cosmology
- Derives the EFT alpha parameters and all their derivatives (w.r.t. time, $X$, $\phi$)
- Transforms perturbation equations from $G$-function language into alpha/EFT language
- Extracts and validates linear and quadratic perturbation coefficients in multiple gauge choices

## Dependencies

- [Mathematica](https://www.wolfram.com/mathematica/) ≥ 13.0
- [xAct](https://www.xact.es/) ≥ 1.3.0, specifically: `xCore`, `xTensor`, `xCoba`, `xTras`, `xPand`

## Installation

Place the `xAlpha/` directory inside your xAct installation folder, at the same level as `xTensor`, `xPand`, etc.:

```
xAct/
├── xAlpha/
│   ├── Kernel/
│   │   └── init.m
│   ├── xAlpha.m
│   ├── Setup.wl
│   ├── HorndeskiDefinitions.wl
│   ├── Background.wl
│   ├── AlphaFunctions.wl
│   ├── Perturbations.wl
│   └── Examples/
├── xTensor/
├── xPand/
└── ...
```

## Quick start

```mathematica
(* 1. Load xAct packages *)
<< xAct`xCoba`;
<< xAct`xTras`;
<< xAct`xPand`;

(* 2. Suppress xAct definition messages (optional but recommended) *)
$DefInfoQ = False;
Off[General::shdw];

(* 3. Define manifold and xPand slicing *)
DefManifold[Md, 4, {b, c, d, e, f}];
DefMetric[-1, g[-b, -c], CD, {";", "\[Del]"}];
$ConformalTime = False;
$CommuteCovDsOnScalars = True;
SetSlicing[g, n, h, cd, {"|", "\!\(\*OverscriptBox[\(D\), \(_\)]\)"}, "FLFlat"];
DefMetricFields[g, dg, h];
DefMatterFields[u, du, h];

(* 4. Load xAlpha — lazy packages safe to load before manifold is fully set *)
Get["/path/to/xAlpha/xAlpha.m"];  (* or << xAct`xAlpha` if installed *)

(* 5. Load Setup.wl AFTER manifold is defined *)
Get[FileNameJoin[{$xAlphaDir, "Setup.wl"}]];

(* 6. Initialise *)
SetupFluidModel["PerfectFluid"];
HorndeskiDerivativeReplacements1 = GetHorndeskiRules[];
Xto\[CurlyPhi] = MakeXtoPhiRule[];

(* 7. Compute EFT alpha parameters *)
ComputeAllAlphas[phi, X, Hh, HorndeskiDerivativeReplacements1, phidotToX]

(* 8. Transform a perturbation equation to alpha language *)
GToAlphas[myExpr, phidotToX]

(* 9. Load pre-computed second-order EOMs (cached .mx files ship with the package) *)
metricEOMsQuad = xAlphaLoadOrCompute["metricEOMsQuad",
  -2 * MyToxPand[metricEOMs /. Xto\[CurlyPhi], "NewtonGauge", 2] /. HorndeskiDerivativeReplacements2
];
scalarEOMQuad = xAlphaLoadOrCompute["scalarEOMQuad",
  -2 * MyToxPand[scalarEOM /. Xto\[CurlyPhi], "NewtonGauge", 2] /. HorndeskiDerivativeReplacements2
];
```

See `Examples/` for a complete worked notebook.
## Authors

Sergi Sirera <!-- add co-authors if any -->

## Citation

If you use xAlpha in your research, please cite this repository and the xAct suite:
> J. M. Martín-García et al., *xAct: Efficient tensor computer algebra for the Wolfram Language*, https://www.xact.es/

## License

GPL-2.0 — see [LICENSE](LICENSE).
