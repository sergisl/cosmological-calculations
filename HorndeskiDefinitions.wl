(* ::Package:: *)

(* ========================================================================== *)
(* Horndeski Theory Definitions Package                                       *)
(* ========================================================================== *)
(* This package contains xAct definitions for Horndeski gravity computations *)
(* To use: Load after defining manifold, metric, and xPand setup             *)
(* ========================================================================== *)

(* Sub-package loaded by xAct`xAlpha` — no nested BeginPackage *)

(* ========================================================================== *)
(* Public API                                                                 *)
(* ========================================================================== *)

SetupFluidModel::usage = "SetupFluidModel[model] sets up stress-energy tensor. Options: \"Dust\", \"PerfectFluid\", \"ImperfectFluid\"";
SetupHorndeskiFunctions::usage = "SetupHorndeskiFunctions[] defines all Horndeski functions and their derivatives";
SetupScalarField::usage = "SetupScalarField[] defines kinetic term X and tilded quantities";

Begin["xAct`xAlpha`Private`"]

(* ========================================================================== *)
(* Stress-Energy Tensor Setup                                                *)
(* ========================================================================== *)

SetupFluidModel[model_String] := Module[{density, pressure, anisoStress},
  (* Define tensor components - Md, g, b, c, d, u must be defined in Global context *)
  DefTensor[T[-Global`b, -Global`c], Global`Md];
  DefTensor[Ts[], Global`Md, PrintAs -> "T"];
  DefTensor[PiRaw[-Global`b, -Global`c], Global`Md, Symmetric[{-Global`b, -Global`c}], 
    PrintAs -> "\[CapitalPi]^{raw}"];
  DefTensor[PiST[-Global`b, -Global`c], Global`Md, Symmetric[{-Global`b, -Global`c}], 
    PrintAs -> "\[CapitalPi]"];
  IndexSet[PiST[Global`b_, Global`c_], 
    PiRaw[Global`b, Global`c] - (1/4) Global`g[Global`b, Global`c] PiRaw[Global`d, -Global`d]];
  
  (* Set flags based on model *)
  {density, pressure, anisoStress} = Which[
    model === "Dust", {True, False, False},
    model === "PerfectFluid", {True, True, False},
    model === "ImperfectFluid", {True, True, True},
    True, (Message[SetupFluidModel::badmodel, model]; {True, True, False})
  ];
  
  (* Store in global variables *)
  Global`$FluidModel = model;
  Global`$Density = density;
  Global`$Pressure = pressure;
  Global`$AnisotropicStress = anisoStress;
  
  (* Define T_ab based on fluid model *)
  IndexSet[T[Global`b_, Global`c_], (
    If[density, \[Rho][Global`u][] Global`u[Global`b] Global`u[Global`c], 0] +
    If[pressure, P[Global`u][] (Global`g[Global`b, Global`c] + Global`u[Global`b] Global`u[Global`c]), 0] +
    If[anisoStress, PiST[Global`b, Global`c], 0]
  )];
  
  If[$xAlphaVerbose, Print["Fluid model set to: ", model]];
];

SetupFluidModel::badmodel = "Unknown fluid model `1`. Using PerfectFluid.";


(* ========================================================================== *)
(* Scalar Field and Kinetic Term Setup                                       *)
(* ========================================================================== *)

SetupScalarField[] := (
  (* Define tensors - Md and h must already be defined in Global context *)
  DefTensor[X[], Global`Md];
  DefTensor[tX[], Global`Md, PrintAs -> "\!\(\*OverscriptBox[\(X\), \(~\)]\)"];
  DefProjectedTensor[tphi[], Global`h, 
    PrintAs -> "\!\(\*OverscriptBox[\(\[CurlyPhi]\), \(~\)]\)"];
  
  If[$xAlphaVerbose, Print["Scalar field tensors defined"]];
);


(* ========================================================================== *)
(* Horndeski Functions Setup                                                 *)
(* ========================================================================== *)

SetupHorndeskiFunctions[] := Module[{},
  (* Base functions - let xAct create them in appropriate context *)
  DefScalarFunction[KK, PrintAs -> "K"];
  DefScalarFunction[G3, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3\)]\)"];
  DefScalarFunction[G4, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4\)]\)"];
  
  (* K derivatives *)
  DefScalarFunction[Kphi, 
    PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\)]\)"];
  DefScalarFunction[KX, 
    PrintAs -> "\!\(\*SubscriptBox[\(K\), \(X\)]\)"];
  DefScalarFunction[Kphiphi, 
    PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]\[CurlyPhi]\)]\)"];
  DefScalarFunction[KXX, 
    PrintAs -> "\!\(\*SubscriptBox[\(K\), \(XX\)]\)"];
  DefScalarFunction[KphiX, 
    PrintAs -> "\!\(\*SubscriptBox[\(K\), \(\[CurlyPhi]X\)]\)"];
  
  (* G3 derivatives *)
  DefScalarFunction[G3phi, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\)]\)"];
  DefScalarFunction[G3X, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  X\)]\)"];
  DefScalarFunction[G3phiphi, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]\[CurlyPhi]\)]\)"];
  DefScalarFunction[G3XX, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  XX\)]\)"];
  DefScalarFunction[G3phiX, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3  \[CurlyPhi]X\)]\)"];
  
  (* G4 derivatives *)
  DefScalarFunction[G4phi, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\)]\)"];
  DefScalarFunction[G4phiphi, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\)]\)"];
  DefScalarFunction[G4phiphiphi, 
    PrintAs -> "\!\(\*SubscriptBox[\(G\), \(4  \[CurlyPhi]\[CurlyPhi]\[CurlyPhi]\)]\)"];
  
  If[$xAlphaVerbose, Print["Horndeski functions defined"]];
];

End[]
