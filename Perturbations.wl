(* ::Package:: *)

(* ========================================================================== *)
(* Horndeski Perturbations Package                                           *)
(* ========================================================================== *)
(* This package provides:                                                     *)
(*   - Definitions for scalar perturbation variables                          *)
(*   - Rules for variable choices (Q, vx, dphi)                               *)
(*   - Second-order perturbation EOM computation                              *)
(* ========================================================================== *)

(* Sub-package loaded by xAct`xAlpha` — no nested BeginPackage *)

DefinePerturbationVariables::usage =
  "DefinePerturbationVariables[] defines projected tensors for scalar perturbations.";

GetGaugeTransformationRules::usage = 
  "GetGaugeTransformationRules[var] returns gauge transformation rules for the given scalar variable choice (\"Q\", \"vx\", or \"dphi\").";

pertHeads::usage = 
  "pertHeads is the list of perturbation variable heads for Q gauge: {Q, Φh, Ψh}.";

pertHeadsvx::usage = 
  "pertHeadsvx is the list of perturbation variable heads for vx gauge: {vx, Φh, Ψh}.";

pertHeadsdphi::usage = 
  "pertHeadsdphi is the list of perturbation variable heads for dphi gauge: {φ, Φh, Ψh}.";

ah0::usage = 
  "ah0 is the background scale factor slot: ah[LI[0], LI[0]].";

base::usage = 
  "base[var] returns the first-order base term: var[LI[1], LI[0]].";

tder::usage = 
  "tder[var] returns the first time derivative: var[LI[1], LI[1]].";

ttder::usage = 
  "ttder[var] returns the second time derivative: var[LI[1], LI[2]].";

lap::usage = 
  "lap[var, cc] returns the Laplacian: cd[-cc][cd[cc][base[var]]].";

gradUp::usage = 
  "gradUp[var, cc] returns the upper-index gradient: cd[cc][base[var]].";

gradDown::usage = 
  "gradDown[var, cc] returns the lower-index gradient: cd[-cc][base[var]].";

gradtUp::usage = 
  "gradtUp[var, cc] returns the upper-index time-gradient: cd[cc][tder[var]].";

gradtDown::usage = 
  "gradtDown[var, cc] returns the lower-index time-gradient: cd[-cc][tder[var]].";

gpair::usage = 
  "gpair[a, b, cc] returns the gradient pair: gradDown[a, cc]*gradUp[b, cc].";

deltaTerm::usage = 
  "deltaTerm[var] returns the metric-projected base term: h[-b, -c] base[var].";

deltaT::usage = 
  "deltaT[var] returns the metric-projected first time derivative: h[-b, -c] tder[var].";

deltaTT::usage = 
  "deltaTT[var] returns the metric-projected second time derivative: h[-b, -c] ttder[var].";

deltaLap::usage = 
  "deltaLap[var, cc] returns the metric-projected Laplacian: h[-b, -c] lap[var, cc].";

gradgrad::usage = 
  "gradgrad[var, b, c] returns the double gradient: cd[-b][cd[-c][base[var]]].";

gradgradUp::usage = 
  "gradgradUp[var, b, c] returns the double upper gradient: cd[b][cd[c][base[var]]].";

deltagpair::usage = 
  "deltagpair[a, bb, cc] returns the metric-projected gradient pair: h[-b, -c] gradDown[a, cc]*gradUp[bb, cc].";

deltagtpair::usage = 
  "deltagtpair[a, bb, cc] returns the metric-projected time-gradient pair: h[-b, -c] gradtDown[a, cc]*gradUp[bb, cc].";

gradagradb::usage = 
  "gradagradb[a, bb, cc, dd] returns the double gradient product: gradDown[a, cc]*gradDown[bb, dd].";

gradagradtb::usage = 
  "gradagradtb[a, bb, cc, dd] returns the mixed gradient-time gradient product: gradDown[a, cc]*gradtDown[bb, dd].";

Dij::usage = 
  "Dij[a, bb] returns the anisotropic stress term: base[a]*(cd[-c][cd[-b][base[bb]]] - deltaLap[bb, d]).";

HorndeskiDerivativeReplacements2::usage = 
  "HorndeskiDerivativeReplacements2 provides derivative replacement rules for second-order perturbations.";

CanonicaliseTermOrientations::usage =
  "CanonicaliseTermOrientations[expr] applies canonical index orientation rules to quadratic gradient \
structures in perturbation expressions, ensuring consistent form for Collect[].";

GetRemoveScalarWrapper::usage =
  "GetRemoveScalarWrapper[] returns the list of replacement rules that strip Scalar[] wrappers \
from Horndeski function arguments (e.g. KK[Scalar[\[CurlyPhi]], X] :> KK[\[CurlyPhi], X]), \
needed before coefficient extraction.";

GetGaugeTransformationRules::badvar =
  "Unknown scalar variable choice: `1`. Supported: \"Q\", \"vx\", \"dphi\".";

Begin["xAct`xAlpha`Private`"]

(* ========================================================================== *)
(* Helper: HorndeskiDerivativeReplacements2                                  *)
(* ========================================================================== *)

HorndeskiDerivativeReplacements2 :=
  Global`HorndeskiDerivativeReplacements1 /. (lhs_ -> rhs_) :>
    (lhs -> Global`MyToxPand[rhs /. Global`Xto\[CurlyPhi], "NewtonGauge", 0]);


(* ========================================================================== *)
(* Scalar Perturbation Variable Definitions                                   *)
(* ========================================================================== *)

DefinePerturbationVariables[] := With[{
  DefProjectedTensor = Symbol["xAct`xPand`DefProjectedTensor"],
  PrintAs = Symbol["xAct`xTensor`PrintAs"],
  Q = Symbol["Global`Q"],
  vx = Symbol["Global`vx"],
  h = Symbol["Global`h"],
  cdcd\[CurlyPhi] = Symbol["Global`cdcd\[CurlyPhi]"],
  cdcdQ = Symbol["Global`cdcdQ"],
  cdcdvx = Symbol["Global`cdcdvx"],
  cdcd\[Phi] = Symbol["Global`cdcd\[Phi]"],
  cdcd\[Psi] = Symbol["Global`cdcd\[Psi]"],
  cdcd2\[CurlyPhi] = Symbol["Global`cdcd2\[CurlyPhi]"],
  cdcd2Q = Symbol["Global`cdcd2Q"],
  cdcd2vx = Symbol["Global`cdcd2vx"],
  cdcd2\[Phi] = Symbol["Global`cdcd2\[Phi]"],
  cdcd2\[Psi] = Symbol["Global`cdcd2\[Psi]"],
  cd\[CurlyPhi] = Symbol["Global`cd\[CurlyPhi]"],
  cdQ = Symbol["Global`cdQ"],
  cdvx = Symbol["Global`cdvx"],
  cd\[Phi] = Symbol["Global`cd\[Phi]"],
  cd\[Psi] = Symbol["Global`cd\[Psi]"]
  },
  
  DefProjectedTensor[Q[], h];
  DefProjectedTensor[vx[], h, PrintAs -> "\!\(\*SubscriptBox[\(v\), \(X\)]\)"];
  DefProjectedTensor[cdcd\[CurlyPhi][], h, PrintAs -> "\!\(\*SuperscriptBox[\(\[Del]\), \(2\)]\)\[CurlyPhi]"];
  DefProjectedTensor[cdcdQ[], h, PrintAs -> "\!\(\*SuperscriptBox[\(\[Del]\), \(2\)]\)Q"];
  DefProjectedTensor[cdcdvx[], h, PrintAs -> "\!\(\*SuperscriptBox[\(\[Del]\), \(2\)]\)\!\(\*SubscriptBox[\(v\), \(X\)]\)"];
  DefProjectedTensor[cdcd\[Phi][], h, PrintAs -> "\!\(\*SuperscriptBox[\(\[Del]\), \(2\)]\)\[Phi]"];
  DefProjectedTensor[cdcd\[Psi][], h, PrintAs -> "\!\(\*SuperscriptBox[\(\[Del]\), \(2\)]\)\[Psi]"];
  DefProjectedTensor[cdcd2\[CurlyPhi][], h, PrintAs -> "\[Del]\[Del]\[CurlyPhi]"];
  DefProjectedTensor[cdcd2Q[], h, PrintAs -> "\[Del]\[Del]Q"];
  DefProjectedTensor[cdcd2vx[], h, PrintAs -> "\[Del]\[Del]\!\(\*SubscriptBox[\(v\), \(X\)]\)"];
  DefProjectedTensor[cdcd2\[Phi][], h, PrintAs -> "\[Del]\[Del]\[Phi]"];
  DefProjectedTensor[cdcd2\[Psi][], h, PrintAs -> "\[Del]\[Del]\[Psi]"];
  DefProjectedTensor[cd\[CurlyPhi][], h, PrintAs -> "\[Del]\[CurlyPhi]"];
  DefProjectedTensor[cdQ[], h, PrintAs -> "\[Del]Q"];
  DefProjectedTensor[cdvx[], h, PrintAs -> "\[Del]\!\(\*SubscriptBox[\(v\), \(X\)]\)"];
  DefProjectedTensor[cd\[Phi][], h, PrintAs -> "\[Del]\[Phi]"];
  DefProjectedTensor[cd\[Psi][], h, PrintAs -> "\[Del]\[Psi]"];
];

(* Replacement rules for scalar variables *)
(* Replacement rules for scalar variables *)
phiToQ = With[{
  LI = Symbol["xAct`xTensor`LI"],
  Hh = Symbol["Global`Hh"],
  Q = Symbol["Global`Q"],
  \[CurlyPhi] = Symbol["xAct`xPand`\[CurlyPhi]"]
  },
  
  {\[CurlyPhi][LI[1], LI[0]] -> \[CurlyPhi][LI[0], LI[1]]/
     Hh[LI[0], LI[0]] Q[LI[1], LI[0]], \[CurlyPhi][LI[1],
     LI[1]] -> \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[1]] + \[CurlyPhi][LI[0], LI[2]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[0]] - \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[1]] Q[LI[1], LI[0]], \[CurlyPhi][
     LI[1], LI[2]] -> \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[2]] + \[CurlyPhi][LI[0], LI[2]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[1]] - \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[1]] Q[LI[1],
       LI[1]] + \[CurlyPhi][LI[0], LI[2]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[1]] + \[CurlyPhi][LI[0], LI[3]]/
      Hh[LI[0], LI[0]] Q[LI[1], LI[0]] - \[CurlyPhi][LI[0], LI[2]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[1]] Q[LI[1],
       LI[0]] - \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[1]] Q[LI[1],
       LI[1]] - \[CurlyPhi][LI[0], LI[2]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[1]] Q[LI[1],
       LI[0]] - \[CurlyPhi][LI[0], LI[1]]/
      Hh[LI[0], LI[0]]^2 Hh[LI[0], LI[2]] Q[LI[1], LI[0]] + (
      2 \[CurlyPhi][LI[0], LI[1]])/
      Hh[LI[0], LI[0]]^3 Hh[LI[0], LI[1]]^2 Q[LI[1], LI[0]]}
];

phiTovx = With[{
  LI = Symbol["xAct`xTensor`LI"],
  vx = Symbol["Global`vx"],
  \[CurlyPhi] = Symbol["xAct`xPand`\[CurlyPhi]"]
  },
  \[CurlyPhi] = Symbol["xAct`xPand`\[CurlyPhi]"];
  
  {\[CurlyPhi][LI[1],
     LI[0]] -> -\[CurlyPhi][LI[0], LI[1]] vx[LI[1],
      LI[0]], \[CurlyPhi][LI[1],
     LI[1]] -> -\[CurlyPhi][LI[0], LI[1]] vx[LI[1],
       LI[1]] - \[CurlyPhi][LI[0], LI[2]] vx[LI[1],
       LI[0]], \[CurlyPhi][LI[1],
     LI[2]] -> -\[CurlyPhi][LI[0], LI[1]] vx[LI[1],
       LI[2]] - \[CurlyPhi][LI[0], LI[2]] vx[LI[1],
       LI[1]] - \[CurlyPhi][LI[0], LI[2]] vx[LI[1],
       LI[1]] - \[CurlyPhi][LI[0], LI[3]] vx[LI[1], LI[0]]}
];

(* ========================================================================== *)
(* Gauge Transformation Helper                                               *)
(* ========================================================================== *)

GetGaugeTransformationRules[var_] := Switch[var,
  "Q", phiToQ,
  "vx", phiTovx,
  "dphi", {},
  _,
  Message[GetGaugeTransformationRules::badvar, var];
  $Failed
];

(* ========================================================================== *)
(* Perturbation Variable Lists                                               *)
(* ========================================================================== *)

With[{
  Q = Symbol["Global`Q"],
  vx = Symbol["Global`vx"],
  φ = Symbol["xAct`xPand`φ"],
  Φh = Symbol["Global`\[Phi]h"],
  Ψh = Symbol["Global`\[Psi]h"]
},

pertHeads = {Q, Φh, Ψh};
pertHeadsvx = {vx, Φh, Ψh};
pertHeadsdphi = {φ, Φh, Ψh};

]; (* End With *)

(* ========================================================================== *)
(* Perturbation Structure Helper Functions                                   *)
(* ========================================================================== *)

With[{
  LI = Symbol["xAct`xTensor`LI"],
  cd = Symbol["xAct`xTensor`cd"],
  h = Symbol["Global`h"],
  ah = Symbol["Global`ah"]
},

(* Background scale factor *)
ah0 := ah[LI[0], LI[0]];

(* Basic perturbation terms *)
base[var_] := var[LI[1], LI[0]];
tder[var_] := var[LI[1], LI[1]];
ttder[var_] := var[LI[1], LI[2]];

(* Spatial derivatives *)
lap[var_, cc_] := cd[-cc][cd[cc][base[var]]];
gradUp[var_, cc_] := cd[cc][base[var]];
gradDown[var_, cc_] := cd[-cc][base[var]];
gradtUp[var_, cc_] := cd[cc][tder[var]];
gradtDown[var_, cc_] := cd[-cc][tder[var]];
gpair[a_, b_, cc_] := gradDown[a, cc]*gradUp[b, cc];

(* Metric projection helpers *)
deltaTerm[var_] := h[-b, -c] base[var];
deltaT[var_] := h[-b, -c] tder[var];
deltaTT[var_] := h[-b, -c] ttder[var];
deltaLap[var_, cc_] := h[-b, -c] lap[var, cc];
gradgrad[var_, b_, c_] := cd[-b][cd[-c][base[var]]];
gradgradUp[var_, b_, c_] := cd[b][cd[c][base[var]]];
deltagpair[a_, bb_, cc_] := h[-b, -c] gradDown[a, cc]*gradUp[bb, cc];
deltagtpair[a_, bb_, cc_] := h[-b, -c] gradtDown[a, cc]*gradUp[bb, cc];
gradagradb[a_, bb_, cc_, dd_] := gradDown[a, cc]*gradDown[bb, dd];
gradagradtb[a_, bb_, cc_, dd_] := gradDown[a, cc]*gradtDown[bb, dd];
Dij[a_, bb_] := base[a]*(cd[-c][cd[-b][base[bb]]] - deltaLap[bb, d]);

]; (* End With *)

(* ========================================================================== *)
(* Canonical index orientation for quadratic gradient structures              *)
(* ========================================================================== *)

CanonicaliseTermOrientations[expr_] :=
  Module[{cdSym, LISym},
    cdSym  = Symbol["Global`cd"];
    LISym  = Symbol["xAct`xTensor`LI"];
    expr //. {
      (* 1. D9/B2 type: cd[-b][tder] cd[b][base]
         Forces the time derivative to take the lower index. *)
      cdSym[b_Symbol][tder : (YY_[LISym[1], LISym[1]] | _[YY_[LISym[1], LISym[1]]])]*
         cdSym[m_][base : (XX_[LISym[1], LISym[0]] | _[XX_[LISym[1], LISym[0]]])] /;
        m === -b :> cdSym[-b][tder]*cdSym[b][base],
      (* 2. D7/B1 type: cd[-b][base1] cd[b][base2]
         Alphabetically first field gets the lower index. *)
      cdSym[b_Symbol][t1 : (XX_[LISym[1], LISym[0]] | _[XX_[LISym[1], LISym[0]]])]*
         cdSym[m_][t2 : (YY_[LISym[1], LISym[0]] | _[YY_[LISym[1], LISym[0]]])] /;
        m === -b && OrderedQ[{XX, YY}] && XX =!= YY :> cdSym[-b][t1]*cdSym[b][t2],
      (* 3. D11 type: cd[-b][tder1] cd[b][tder2]
         Alphabetically first field gets the lower index. *)
      cdSym[b_Symbol][t1 : (XX_[LISym[1], LISym[1]] | _[XX_[LISym[1], LISym[1]]])]*
         cdSym[m_][t2 : (YY_[LISym[1], LISym[1]] | _[YY_[LISym[1], LISym[1]]])] /;
        m === -b && OrderedQ[{XX, YY}] && XX =!= YY :> cdSym[-b][t1]*cdSym[b][t2],
      (* 4. CC8 type: cd[-b][base1] cd[-c][base2]
         Forces idx1 to alphabetically precede idx2 before checking field order. *)
      cdSym[idx1_][t1 : (XX_[LISym[1], LISym[0]] | _[XX_[LISym[1], LISym[0]]])]*
         cdSym[idx2_][t2 : (YY_[LISym[1], LISym[0]] | _[YY_[LISym[1], LISym[0]]])] /;
        OrderedQ[{idx1, idx2}] && !OrderedQ[{XX, YY}] :> cdSym[idx1][t2]*cdSym[idx2][t1]
    }
  ];

(* ========================================================================== *)
(* removeScalarWrapper rules                                                  *)
(* ========================================================================== *)

GetRemoveScalarWrapper[] :=
  With[{Sc = Symbol["xAct`xTensor`Scalar"]},
    Flatten[{
      Table[
        With[{fn = Symbol["Global`" <> name]},
          HoldPattern[fn[Sc[phi_], XX_]] :> fn[phi, XX]
        ],
        {name, {"KK", "KX", "KXX", "KXXX", "Kphi", "KphiX", "KphiXX",
                "Kphiphi", "KphiphiX", "G3", "G3X", "G3XX", "G3XXX",
                "G3phi", "G3phiX", "G3phiXX", "G3phiXXX",
                "G3phiphi", "G3phiphiX", "G3phiphiphi"}}
      ],
      Table[
        With[{fn = Symbol["Global`" <> name]},
          HoldPattern[fn[Sc[phi_]]] :> fn[phi]
        ],
        {name, {"G4", "G4phi", "G4phiphi"}}
      ]
    }]
  ];

End[]
