# xAlpha

<p align="left">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="Assets/xAlpha-dark.png?v=1">
    <source media="(prefers-color-scheme: light)" srcset="Assets/xAlpha-light.png?v=1">
    <img alt="xAlpha logo" src="Assets/xAlpha-light.png" width="170">
  </picture>
</p>

A [Mathematica](https://www.wolfram.com/mathematica/) package built on [xAct](https://www.xact.es/)/[xPand](https://www2.iap.fr/users/pitrou/xpand.htm) for computing cosmological perturbation equations in luminal Horndeski gravity and expressing them in terms of the EFT alpha parameters $\alpha_M$, $\alpha_K$, $\alpha_B$ and their derivatives.

Given a luminal Horndeski theory defined by functions $K(\phi, X)$, $G_3(\phi, X)$, $G_4(\phi)$, xAlpha:

- Computes equations of motion for a cosmological FLRW background to quadratic order in perturbations using xPand
- Computes the EFT alpha parameters $(\alpha_M, \alpha_K, \alpha_B)$ and all their derivatives with respect to time, $X$, and $\phi$ up to third order
- Replaces Horndeski $G$-function derivatives in perturbation equations with alpha parameters
- Extracts and stores the linear and quadratic perturbation coefficients from the $(0,0)$, $(0,i)$, $(i,j)$, and scalar field equations
- Converts coefficients between raw Horndeski, alpha-basis, and gamma-basis EFT bases
- Validates that extracted coefficients fully reconstruct the original equations

---

## Dependencies

- [Mathematica](https://www.wolfram.com/mathematica/) ≥ 13.0
- [xAct](https://www.xact.es/) ≥ 1.3.0: `xCore`, `xTensor`, `xPand`

---

## Installation

Run this once in any Mathematica notebook to install xAlpha system-wide alongside other xAct packages:

```mathematica
src  = "/path/to/xAlpha";   (* folder containing xAlpha.m *)
dest = FileNameJoin[{$UserBaseDirectory, "Applications", "xAct", "xAlpha"}];
RunProcess[{"ln", "-s", src, dest}];
```

After this, `<< xAct`xAlpha`` works from any notebook. Alternatively, load directly:

```mathematica
Get["/path/to/xAlpha/xAlpha.m"];
```

> **Note on `Setup.wl`:** It must be loaded manually *after* the xAct manifold and xPand slicing are defined, since it contains `DefTensor` calls that depend on them. See the worked example notebook in `Examples/`.

---

## Examples

See [`Examples/luminalHorndeski-cosmoperts.nb`](Examples/luminalHorndeski-cosmoperts.nb) for a complete worked example.

---

## Notes on `KK`

The Horndeski kinetic function $K(\phi, X)$ is defined as `KK` internally (printed as $K$) to avoid a clash with Mathematica's built-in `K[m]`.

---

## Citation

If you use this code, please cite:

```bibtex
@article{Sirera2026master,
  author  = {Sirera, Sergi and Baker, Tessa and Hallam, James and Naidoo, Krishna},
  title   = {{A Master Equation for Screening in Luminal Horndeski Gravity}},
  journal = {arXiv},
  year    = {2026},
  note    = {arXiv:XXXX.XXXXX},
  url     = {https://arxiv.org/abs/XXXX.XXXXX},
}
```

Please also cite the xAct suite:
> J. M. Martín-García et al., *xAct: Efficient tensor computer algebra for the Wolfram Language*, https://www.xact.es/

---

## AI assistance disclosure

Parts of this codebase were developed with the assistance of GitHub Copilot. All physics formulations and scientific results were designed and verified by the authors.

---

## Contact

For questions or comments, please contact sergi.sirera@port.ac.uk.
