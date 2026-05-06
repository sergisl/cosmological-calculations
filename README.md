# xAlpha

<p align="left">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="Assets/xAlpha-dark.png?v=1">
    <source media="(prefers-color-scheme: light)" srcset="Assets/xAlpha-light.png?v=1">
    <img alt="xAlpha logo" src="Assets/xAlpha-light.png" width="170">
  </picture>
</p>

A Mathematica package built on [xAct](https://www.xact.es/)/[xPand](https://www2.iap.fr/users/pitrou/xpand.htm) for computing cosmological perturbation equations in luminal Horndeski gravity and expressing them in terms of the EFT alpha parameters $\alpha_M$, $\alpha_K$, $\alpha_B$ and their derivatives.

Given a luminal Horndeski theory defined by functions $K(\phi, X)$, $G_3(\phi, X)$, $G_4(\phi)$, xAlpha:

- Computes equations of motion for a FLRW background to quadratic order in perturbations using xPand
- Computes the alpha parameters $(\alpha_M, \alpha_K, \alpha_B)$ and all their derivatives with respect to time, $X$, and $\phi$ up to third order
- Extracts and stores the linear and quadratic perturbation coefficients from the metric $(0,0)$, $(0,i)$, $(i,j)$, and scalar field equations
- Systematically replaces Horndeski $G$-function derivatives with alpha parameters
- Validates that extracted coefficients fully reconstruct the original equations

---

## Dependencies

- [Mathematica](https://www.wolfram.com/mathematica/) ≥ 13.0
- [xAct](https://www.xact.es/) ≥ 1.3.0: `xCore`, `xTensor`, `xPand`

---

## Installation

Clone the repo directly into your xAct applications folder:

```bash
cd "$HOME/Library/Mathematica/Applications/xAct"
git clone https://github.com/sergisl/xAlpha.git xAlpha
```

After this, `<< xAct``xAlpha`` ` works from any Mathematica notebook, just like `<< xAct``xPand`` `. Alternatively, load directly by path without installing:

```mathematica
Get["/path/to/xAlpha/xAlpha.m"];
```

> **Note on `Setup.wl`:** It must be loaded manually *after* the xAct manifold and xPand slicing are defined, since it contains `DefTensor` calls that depend on them. See the worked example notebook in `Examples/`.

---

## Examples

See [`luminalHorndeski-cosmoperts`](Examples/luminalHorndeski-cosmoperts.nb) for a complete worked notebook reproducing the results in [2605.XXXXX](https://arxiv.org/abs/2605.XXXXX).

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

Please also cite xAct and xPand.

---

## Contact

For questions or comments, please contact sergi.sirera@port.ac.uk.

> Parts of this codebase were developed with the assistance of GitHub Copilot. All physics formulations and scientific results were designed and verified by the authors.
