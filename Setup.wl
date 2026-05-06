(* ::Package:: *)

(* ========================================================================== *)
(* xAlpha: Setup Script                                                       *)
(* ========================================================================== *)
(* Defines all xAct tensors and functions for Horndeski gravity.              *)
(* Must be loaded AFTER DefManifold / SetSlicing / DefMetricFields.           *)
(* ========================================================================== *)

If[$xAlphaVerbose, Print["[xAlpha] Loading Setup.wl..."]]; 

(* Usage strings for Setup.wl user-facing helpers *)
MyToxPand::usage = "MyToxPand[expr, gauge, order] expands expr using xPand's ToxPand with the metric perturbation dg, 4-velocity u and its perturbation du, and spatial metric h, for the specified gauge (e.g. \"NewtonGauge\", \"SynchronousGauge\") and perturbation order (1 or 2). This is the standard xAct community wrapper around ToxPand.";
org::usage = "org[expr] applies NoScalar, ContractMetric, and collects in powers of $PerturbationParameter via ToCanonical. Shorthand for organising xPand perturbation output into canonical form.";
collect::usage = "collect[expr] applies NoScalar and collects in powers of $PerturbationParameter with Identity (no canonicalisation). Lighter-weight alternative to org[] when index canonicalisation is not needed.";
MakeXtoPhiRule::usage = "MakeXtoPhiRule[] constructs and returns the replacement rule X[] -> -(1/2) Scalar[g^{ab} \[Del]_a\[CurlyPhi] \[Del]_b\[CurlyPhi]], expressing the kinetic variable X in terms of the scalar field gradient. Used to convert X-dependent expressions back to metric/field language before applying ToxPand.";
MakePhidotToXRules::usage = "MakePhidotToXRules[] constructs and returns a list of replacement rules mapping powers of the background scalar field time derivative \[CurlyPhi][LI[0],LI[1]]^n to powers of X[], using X = \[CurlyPhi]'^2/2. Also handles Scalar[\[CurlyPhi]']^n and inverse powers. Returns a flat rule list suitable for use with //.";
SetupFluidModel::usage = "SetupFluidModel[model] defines the stress-energy tensor T[-b,-c] for the specified fluid model. Supported values: \"Dust\" (density only), \"PerfectFluid\" (density + isotropic pressure), \"ImperfectFluid\" (density + pressure + anisotropic stress). Sets $FluidModel, $Density, $Pressure, $AnisotropicStress flags. Must be called after the manifold and metric are defined."; 

Block[{$DefInfoQ = False},

(* ========================================================================== *)
(* xPand Helper Functions                                                     *)
(* ========================================================================== *)

MyToxPand[expr_, gauge_, order_] := ToxPand[expr, dg, u, du, h, gauge, order];
org[expr_] := NoScalar@Collect[ContractMetric[expr], $PerturbationParameter, ToCanonical];
collect[expr_] := NoScalar@Collect[expr, $PerturbationParameter, Identity];

(* ========================================================================== *)
(* Scalar Field Definitions                                                   *)
(* ========================================================================== *)

DefTensor[X[], Md];
DefTensor[tX[], Md, PrintAs -> "\!\(\*OverscriptBox[\(X\), \(~\)]\)"];
DefProjectedTensor[tphi[], h, PrintAs -> "\!\(\*OverscriptBox[\(\[CurlyPhi]\), \(~\)]\)"];

(* ========================================================================== *)
(* Horndeski Functions                                                        *)
(* ========================================================================== *)

(* Base functions *)
DefScalarFunction[KK, PrintAs -> "K"];
DefScalarFunction[G3, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3\)]\)"];
DefScalarFunction[G4, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4\)]\)"];

(* K derivatives *)
DefScalarFunction[Kphi, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\)]\)"];
DefScalarFunction[KX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(X\)]\)"];
DefScalarFunction[Kphiphi, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[KXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(XX\)]\)"];
DefScalarFunction[KphiX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]X\)]\)"];
DefScalarFunction[Kphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[KphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]X\)]\)"];
DefScalarFunction[KphiXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]XX\)]\)"];
DefScalarFunction[KXXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(XXX\)]\)"];
DefScalarFunction[KXXXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(XXXX\)]\)"];
DefScalarFunction[Kphiphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[KphiphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]X\)]\)"];
DefScalarFunction[KphiphiXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]XX\)]\)"];
DefScalarFunction[KphiXXX, PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]XXX\)]\)"];

(* G3 derivatives *)
DefScalarFunction[G3phi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\)]\)"];
DefScalarFunction[G3X, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  X\)]\)"];
DefScalarFunction[G3phiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[G3XX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  XX\)]\)"];
DefScalarFunction[G3phiX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]X\)]\)"];
DefScalarFunction[G3phiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[G3phiphiX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]X\)]\)"];
DefScalarFunction[G3phiXX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]XX\)]\)"];
DefScalarFunction[G3XXX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  XXX\)]\)"];
DefScalarFunction[G3XXXX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  XXXX\)]\)"];
DefScalarFunction[G3phiphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[G3phiphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]X\)]\)"];
DefScalarFunction[G3phiphiXX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]XX\)]\)"];
DefScalarFunction[G3phiXXX, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]XXX\)]\)"];

(* G4 derivatives *)
DefScalarFunction[G4phi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\)]\)"];
DefScalarFunction[G4phiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[G4phiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
DefScalarFunction[G4phiphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];

(* ========================================================================== *)
(* Dimensionless Horndeski Functions                                         *)
(* ========================================================================== *)

(* tK functions *)
DefScalarFunction[tK, PrintAs -> "\!\(\*OverscriptBox[\(K\), \(~\)]\)"];
DefScalarFunction[tKphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tKphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tKphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tKX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(X\)], \(~\)]\)"];
DefScalarFunction[tKXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(XX\)], \(~\)]\)"];
DefScalarFunction[tKXXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(XXX\)], \(~\)]\)"];
DefScalarFunction[tKphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tKphiphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tKphiXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]XX\)], \(~\)]\)"];
DefScalarFunction[tKXXXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(XXXX\)], \(~\)]\)"];
DefScalarFunction[tKphiphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tKphiphiphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tKphiphiXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]XX\)], \(~\)]\)"];
DefScalarFunction[tKphiXXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(K\), \(\[CurlyPhi]XXX\)], \(~\)]\)"];

(* tG3 functions *)
DefScalarFunction[tG3, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3\)], \(~\)]\)"];
DefScalarFunction[tG3phi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG3phiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG3phiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG3X, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  X\)], \(~\)]\)"];
DefScalarFunction[tG3XX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  XX\)], \(~\)]\)"];
DefScalarFunction[tG3XXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  XXX\)], \(~\)]\)"];
DefScalarFunction[tG3phiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tG3phiphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tG3phiXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]XX\)], \(~\)]\)"];
DefScalarFunction[tG3XXXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  XXXX\)], \(~\)]\)"];
DefScalarFunction[tG3phiphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG3phiphiphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]X\)], \(~\)]\)"];
DefScalarFunction[tG3phiphiXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]XX\)], \(~\)]\)"];
DefScalarFunction[tG3phiXXX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(3  \[CurlyPhi]XXX\)], \(~\)]\)"];

(* tG4 functions *)
DefScalarFunction[tG4, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(4\)], \(~\)]\)"];
DefScalarFunction[tG4phi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(4  \[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG4phiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];
DefScalarFunction[tG4phiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)], \(~\)]\)"];

(* Dimensionless fluid quantities *)
DefScalarFunction[trho, PrintAs -> "\!\(\*OverscriptBox[\(\[Rho]\), \(~\)]\)"];
DefScalarFunction[tp, PrintAs -> "\!\(\*OverscriptBox[\(p\), \(~\)]\)"];

(* ========================================================================== *)
(* Physical Scales                                                            *)
(* ========================================================================== *)

DefConstantSymbol[Mpl, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(Pl\)]\)"];
DefConstantSymbol[H0, PrintAs -> "\!\(\*SubscriptBox[\(H\), \(0\)]\)"];
DefConstantSymbol[Ms, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(s\)]\)"];
DefConstantSymbol[MK, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(K\)]\)"];
DefConstantSymbol[MG3, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(G3\)]\)"];
DefConstantSymbol[MG4, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(G4\)]\)"];

(* ========================================================================== *)
(* EFT Alpha Parameters                                                       *)
(* ========================================================================== *)

(* Basic alpha parameters *)
DefScalarFunction[Mstar, PrintAs -> "\!\(\*SubscriptBox[\(M\), \(*\)]\)"];
DefScalarFunction[aM, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\)]\)"];
DefScalarFunction[aK, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\)]\)"];
DefScalarFunction[aB, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\)]\)"];

DefScalarFunction[aMdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(M\)], \(.\)]\)"];
DefScalarFunction[aKdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(K\)], \(.\)]\)"];
DefScalarFunction[aBdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(B\)], \(.\)]\)"];
DefScalarFunction[aMddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(M\)], \(..\)]\)"];
DefScalarFunction[aKddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(K\)], \(..\)]\)"];
DefScalarFunction[aBddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(B\)], \(..\)]\)"];
DefScalarFunction[aMdddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(M\)], \(...\)]\)"];
DefScalarFunction[aKdddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(K\)], \(...\)]\)"];
DefScalarFunction[aBdddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(B\)], \(...\)]\)"];

(* X derivatives *)
DefScalarFunction[aMX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(MX\)]\)"];
DefScalarFunction[aBX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(BX\)]\)"];
DefScalarFunction[aKX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(KX\)]\)"];
DefScalarFunction[aMXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(MXX\)]\)"];
DefScalarFunction[aBXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(BXX\)]\)"];
DefScalarFunction[aKXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(KXX\)]\)"];
DefScalarFunction[aMXXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(MXXX\)]\)"];
DefScalarFunction[aBXXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(BXXX\)]\)"];
DefScalarFunction[aKXXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(KXXX\)]\)"];

(* Phi derivatives *)
DefScalarFunction[aMphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]\)]\)"];
DefScalarFunction[aBphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]\)]\)"];
DefScalarFunction[aKphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]\)]\)"];
DefScalarFunction[aMphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]\[Phi]\)]\)"];
DefScalarFunction[aBphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]\[Phi]\)]\)"];
DefScalarFunction[aKphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]\[Phi]\)]\)"];
DefScalarFunction[aMphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]\[Phi]\[Phi]\)]\)"];
DefScalarFunction[aBphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]\[Phi]\[Phi]\)]\)"];
DefScalarFunction[aKphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]\[Phi]\[Phi]\)]\)"];

(* Mixed derivatives *)
DefScalarFunction[aMphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]X\)]\)"];
DefScalarFunction[aBphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]X\)]\)"];
DefScalarFunction[aKphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]X\)]\)"];
DefScalarFunction[aMphiXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]XX\)]\)"];
DefScalarFunction[aBphiXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]XX\)]\)"];
DefScalarFunction[aKphiXX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]XX\)]\)"];
DefScalarFunction[aMphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(M\[Phi]\[Phi]X\)]\)"];
DefScalarFunction[aBphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\[Phi]\[Phi]X\)]\)"];
DefScalarFunction[aKphiphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\[Phi]\[Phi]X\)]\)"];

(* Time derivatives of alpha derivatives *)
DefScalarFunction[aKXdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(KX\)], \(.\)]\)"];
DefScalarFunction[aBXXdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(BXX\)], \(.\)]\)"];
DefScalarFunction[aBXdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(BX\)], \(.\)]\)"];
DefScalarFunction[aMXdot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(MX\)], \(.\)]\)"];
DefScalarFunction[aKXddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(KX\)], \(..\)]\)"];
DefScalarFunction[aBXddot, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Alpha]\), \(BX\)], \(..\)]\)"];

(* Gamma functions (used in alpha-to-gamma conversion) *)
DefScalarFunction[gamK, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(K\)]\)"];
DefScalarFunction[gamM, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(M\)]\)"];
DefScalarFunction[gamB, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(B\)]\)"];
DefScalarFunction[gamE, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(E\)]\)"];
DefScalarFunction[gamX, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(X\)]\)"];
DefScalarFunction[gamA, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(A\)]\)"];
DefScalarFunction[gamD, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(D\)]\)"];
DefScalarFunction[gamC, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(C\)]\)"];
DefScalarFunction[gamF, PrintAs -> "\!\(\*SubscriptBox[\(\[Gamma]\), \(F\)]\)"];
DefScalarFunction[dotgamK, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(K\)], \(.\)]\)"];
DefScalarFunction[dotgamM, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(M\)], \(.\)]\)"];
DefScalarFunction[dotgamB, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(B\)], \(.\)]\)"];
DefScalarFunction[dotgamX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(X\)], \(.\)]\)"];
DefScalarFunction[dotgamD, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(D\)], \(.\)]\)"];
DefScalarFunction[dotgamF, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[Gamma]\), \(F\)], \(.\)]\)"];
DefScalarFunction[Dc, PrintAs -> "\!\(\*SubscriptBox[\(D\), \(c\)]\)"];

(* ========================================================================== *)
(* Background Metric EOM Combinations                                         *)
(* ========================================================================== *)

(* Dimensionful combinations *)
DefScalarFunction[varE, PrintAs -> "\[ScriptCapitalE]"];
DefScalarFunction[varP, PrintAs -> "\[ScriptCapitalP]"];
DefScalarFunction[varS, PrintAs -> "\[ScriptCapitalS]"];
DefScalarFunction[varEdot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(.\)]\)"];
DefScalarFunction[varPdot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalP]\), \(.\)]\)"];
DefScalarFunction[varSdot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalS]\), \(.\)]\)"];
DefScalarFunction[varEddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(..\)]\)"];
DefScalarFunction[varPddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalP]\), \(..\)]\)"];
DefScalarFunction[varSddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalS]\), \(..\)]\)"];
DefScalarFunction[varEdddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(...\)]\)"];
DefScalarFunction[varPdddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalP]\), \(...\)]\)"];
DefScalarFunction[varSdddot, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalS]\), \(...\)]\)"];
DefScalarFunction[varEphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\)]\)"];
DefScalarFunction[varPphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\)]\)"];
DefScalarFunction[varSphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\)]\)"];
DefScalarFunction[varEphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varPphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varSphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varEphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varPphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varSphiphiphi, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\[Phi]\[Phi]\)]\)"];
DefScalarFunction[varEphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]X\)]\)"];
DefScalarFunction[varPphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]X\)]\)"];
DefScalarFunction[varSphiX, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]X\)]\)"];

(* Dimensionless combinations *)
DefScalarFunction[tvarE, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(~\)]\)"];
DefScalarFunction[tvarP, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalP]\), \(~\)]\)"];
DefScalarFunction[tvarS, PrintAs -> "\!\(\*OverscriptBox[\(\[ScriptCapitalS]\), \(~\)]\)"];
DefScalarFunction[tvarEdot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalE]\), \(~\)], \(.\)]\)"];
DefScalarFunction[tvarPdot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalP]\), \(~\)], \(.\)]\)"];
DefScalarFunction[tvarSdot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalS]\), \(~\)], \(.\)]\)"];
DefScalarFunction[tvarEddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalE]\), \(~\)], \(..\)]\)"];
DefScalarFunction[tvarPddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalP]\), \(~\)], \(..\)]\)"];
DefScalarFunction[tvarSddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalS]\), \(~\)], \(..\)]\)"];
DefScalarFunction[tvarEdddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalE]\), \(~\)], \(...\)]\)"];
DefScalarFunction[tvarPdddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalP]\), \(~\)], \(...\)]\)"];
DefScalarFunction[tvarSdddot, PrintAs -> "\!\(\*OverscriptBox[OverscriptBox[\(\[ScriptCapitalS]\), \(~\)], \(...\)]\)"];
DefScalarFunction[tvarEphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarPphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarSphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarEphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarPphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarSphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarEphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarPphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarSphiphiphi, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]\[Phi]\)], \(~\)]\)"];
DefScalarFunction[tvarEphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalE]\), \(\[Phi]X\)], \(~\)]\)"];
DefScalarFunction[tvarPphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalP]\), \(\[Phi]X\)], \(~\)]\)"];
DefScalarFunction[tvarSphiX, PrintAs -> "\!\(\*OverscriptBox[SubscriptBox[\(\[ScriptCapitalS]\), \(\[Phi]X\)], \(~\)]\)"]

]; (* End Block[$DefInfoQ = False] *)

If[$xAlphaVerbose, Print["[xAlpha] Setup complete: Horndeski tensors and functions defined."]];

(* ========================================================================== *)
(* Stress-Energy Tensor Setup Function                                        *)
(* ========================================================================== *)

SetupFluidModel[model_String] := Module[{},
  (* Define tensor components *)
  DefTensor[T[-b, -c], Md];
  DefTensor[Ts[], Md, PrintAs -> "T"];
  DefTensor[PiRaw[-b, -c], Md, Symmetric[{-b, -c}], PrintAs -> "\[CapitalPi]^{raw}"];
  DefTensor[PiST[-b, -c], Md, Symmetric[{-b, -c}], PrintAs -> "\[CapitalPi]"];
  IndexSet[PiST[b_, c_], PiRaw[b, c] - (1/4) g[b, c] PiRaw[d, -d]];
  
  (* Set fluid model *)
  $FluidModel = model;
  
  (* Set flags based on model *)
  If[model === "Dust",
    $Density = True; $Pressure = False; $AnisotropicStress = False;
  ];
  If[model === "PerfectFluid",
    $Density = True; $Pressure = True; $AnisotropicStress = False;
  ];
  If[model === "ImperfectFluid",
    $Density = True; $Pressure = True; $AnisotropicStress = True;
  ];
  
  (* Define T_ab based on fluid model *)
  IndexSet[T[b_, c_], (
    If[$Density, \[Rho][u][] u[b] u[c], 0] +
    If[$Pressure, P[u][] (g[b, c] + u[b] u[c]), 0] +
    If[$AnisotropicStress, PiST[b, c], 0]
  )];
  
  If[$xAlphaVerbose, Print["[xAlpha] Fluid model: ", model]];
];

(* ========================================================================== *)
(* Conversion Rule Helpers                                                    *)
(* ========================================================================== *)

MakeXtoPhiRule[] := MakeRule[{X[], -(1/2) Scalar[g[d, e] CD[-d][\[CurlyPhi][]] CD[-e][\[CurlyPhi][]]]}];

(* phidotToX is built lazily on first call to avoid ToCanonical errors at load time *)
MakePhidotToXRules[] := Flatten[{
  (* Bare patterns *)
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^2, 2 X[]}], 
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^3, 2 \[CurlyPhi][LI[0], LI[1]] X[]}], 
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^4, 4 X[]^2}], 
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^5, 4 \[CurlyPhi][LI[0], LI[1]] X[]^2}], 
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^6, 8 X[]^3}], 
  MakeRule[{\[CurlyPhi][LI[0], LI[1]]^7, 6 \[CurlyPhi][LI[0], LI[1]] X[]^3}], 
  MakeRule[{1/(\[CurlyPhi][LI[0], LI[1]]^2), 1/(2 X[])}],
  {c_.*\[CurlyPhi][LI[0], LI[1]]^-2 :> c/(2 X[]), 
   c_.*\[CurlyPhi][LI[0], LI[1]]^-3 :> c/(2 X[]*\[CurlyPhi][LI[0], LI[1]]^-1)},
  (* Patterns inside Scalar wrapper *)
  {Scalar[\[CurlyPhi][LI[0], LI[1]]]^2/2 :> X[],
   Scalar[\[CurlyPhi][LI[0], LI[1]]]^2 :> 2 X[],
   Scalar[\[CurlyPhi][LI[0], LI[1]]]^3 :> 2 Scalar[\[CurlyPhi][LI[0], LI[1]]] X[],
   Scalar[\[CurlyPhi][LI[0], LI[1]]]^4 :> 4 X[]^2,
   1/(Scalar[\[CurlyPhi][LI[0], LI[1]]]^2) :> 1/(2 X[])}
}];

(* Initialise phidotToX on load (manifold is available at this point) *)
phidotToX := phidotToX = MakePhidotToXRules[];

(* ========================================================================== *)
(* Dimensionless Replacement Rules                                            *)
(* ========================================================================== *)

HorndeskiDimensionlessReplacement = {
  (* Base Functions *)
  \[CurlyPhi][LI[0], LI[0]] :> Ms tphi[LI[0], LI[0]], 
  \[CurlyPhi][LI[0], LI[1]] :> Ms tphi[LI[0], LI[1]], 
  \[CurlyPhi][LI[0], LI[2]] :> Ms tphi[LI[0], LI[2]], 
  \[CurlyPhi][LI[0], LI[3]] :> Ms tphi[LI[0], LI[3]], 
  X[] :> Ms^2 H0^2 tX[], 
  KK[a_, b_] :> MK^4 tK[a, b], 
  G3[a_, b_] :> MG3 tG3[a, b], 
  G4[a_] :> MG4^2 tG4[a],
  
  (* First Derivatives *)
  Kphi[a_, b_] :> 1/Ms MK^4 tKphi[a, b], 
  KX[a_, b_] :> 1/(Ms^2 H0^2) MK^4 tKX[a, b], 
  G3phi[a_, b_] :> 1/Ms MG3 tG3phi[a, b], 
  G3X[a_, b_] :> 1/(Ms^2 H0^2) MG3 tG3X[a, b], 
  G4phi[a_] :> 1/Ms MG4^2 tG4phi[a],
  
  (* Second Derivatives *)
  Kphiphi[a_, b_] :> (1/Ms)^2 MK^4 tKphiphi[a, b], 
  KphiX[a_, b_] :> 1/Ms 1/(Ms^2 H0^2) MK^4 tKphiX[a, b], 
  KXX[a_, b_] :> (1/(Ms^2 H0^2))^2 MK^4 tKXX[a, b],
  G3phiphi[a_, b_] :> (1/Ms)^2 MG3 tG3phiphi[a, b], 
  G3phiX[a_, b_] :> 1/Ms 1/(Ms^2 H0^2) MG3 tG3phiX[a, b], 
  G3XX[a_, b_] :> (1/(Ms^2 H0^2))^2 MG3 tG3XX[a, b],
  G4phiphi[a_] :> (1/Ms)^2 MG4^2 tG4phiphi[a],
  
  (* Third Derivatives *)
  Kphiphiphi[a_, b_] :> (1/Ms)^3 MK^4 tKphiphiphi[a, b], 
  KphiphiX[a_, b_] :> (1/Ms)^2 1/(Ms^2 H0^2) MK^4 tKphiphiX[a, b], 
  KphiXX[a_, b_] :> 1/Ms (1/(Ms^2 H0^2))^2 MK^4 tKphiXX[a, b], 
  KXXX[a_, b_] :> (1/(Ms^2 H0^2))^3 MK^4 tKXXX[a, b],
  G3phiphiphi[a_, b_] :> (1/Ms)^3 MG3 tG3phiphiphi[a, b], 
  G3phiphiX[a_, b_] :> (1/Ms)^2 1/(Ms^2 H0^2) MG3 tG3phiphiX[a, b], 
  G3phiXX[a_, b_] :> 1/Ms (1/(Ms^2 H0^2))^2 MG3 tG3phiXX[a, b], 
  G3XXX[a_, b_] :> (1/(Ms^2 H0^2))^3 MG3 tG3XXX[a, b],
  G4phiphiphi[a_] :> (1/Ms)^3 MG4^2 tG4phiphiphi[a],
  
  (* Fourth Derivatives *)
  Kphiphiphiphi[a_, b_] :> (1/Ms)^4 MK^4 tKphiphiphiphi[a, b],
  KphiphiphiX[a_, b_] :> (1/Ms)^3 1/(Ms^2 H0^2) MK^4 tKphiphiphiX[a, b],
  KphiphiXX[a_, b_] :> (1/Ms)^2 (1/(Ms^2 H0^2))^2 MK^4 tKphiphiXX[a, b],
  KphiXXX[a_, b_] :> 1/Ms (1/(Ms^2 H0^2))^3 MK^4 tKphiXXX[a, b],
  G3phiphiphiphi[a_, b_] :> (1/Ms)^4 MG3 tG3phiphiphiphi[a, b],
  G3phiphiphiX[a_, b_] :> (1/Ms)^3 1/(Ms^2 H0^2) MG3 tG3phiphiphiX[a, b],
  G3phiphiXX[a_, b_] :> (1/Ms)^2 (1/(Ms^2 H0^2))^2 MG3 tG3phiphiXX[a, b],
  G3phiXXX[a_, b_] :> 1/Ms (1/(Ms^2 H0^2))^3 MG3 tG3phiXXX[a, b]
};

CosmologicalMassScaleMatching = {
  Ms -> Mpl, 
  MK -> (Mpl H0)^(1/2), 
  MG3 -> Mpl, 
  MG4 -> Mpl
};

NormaliseMstar = {
  varE -> Mstar[]^2 tvarE, varEdot -> Mstar[]^2 tvarEdot, varEddot -> Mstar[]^2 tvarEddot, varEdddot -> Mstar[]^2 tvarEdddot,
  varEphi -> Mstar[]^2 tvarEphi, varEphiphi -> Mstar[]^2 tvarEphiphi,
  varEphiphiphi -> Mstar[]^2 tvarEphiphiphi, varEphiX -> Mstar[]^2 tvarEphiX,
  varP -> Mstar[]^2 tvarP, varPdot -> Mstar[]^2 tvarPdot, varPddot -> Mstar[]^2 tvarPddot, varPdddot -> Mstar[]^2 tvarPdddot,
  varPphi -> Mstar[]^2 tvarPphi, varPphiphi -> Mstar[]^2 tvarPphiphi,
  varPphiphiphi -> Mstar[]^2 tvarPphiphiphi, varPphiX -> Mstar[]^2 tvarPphiX,
  varS -> Mstar[]^2 tvarS, varSdot -> Mstar[]^2 tvarSdot, varSddot -> Mstar[]^2 tvarSddot, varSdddot -> Mstar[]^2 tvarSdddot,
  varSphi -> Mstar[]^2 tvarSphi, varSphiphi -> Mstar[]^2 tvarSphiphi,
  varSphiphiphi -> Mstar[]^2 tvarSphiphiphi, varSphiX -> Mstar[]^2 tvarSphiX,
  \[Rho]u[LI[0], LI[0]] -> Mstar[]^2 trho[],
  Pu[LI[0], LI[0]] -> Mstar[]^2 tp[]
};

If[$xAlphaVerbose, Print["[xAlpha] Setup.wl loaded successfully"]];
