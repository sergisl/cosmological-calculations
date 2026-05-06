(* ::Package:: *)

(* ========================================================================== *)
(* xAlpha: EFT alpha parameters for Horndeski gravity                        *)
(* An xAct package for cosmological perturbation theory                      *)
(* ========================================================================== *)
(* Authors:                                                                   *)
(* License: GPL-2.0                                                           *)
(* ========================================================================== *)

BeginPackage["xAct`xAlpha`",
  {"xAct`xCore`", "xAct`xTensor`", "xAct`xCoba`", "xAct`xPand`",
   "xAct`xTras`"}
]

(* Version info *)
xAct`xAlpha`$Version     = {"0.1.0", {2026, 4, 28}};
xAct`xAlpha`$xActVersion = {"1.3.0", {2025, 12, 29}};

If[Unevaluated[$xAlphaVerbose] === $xAlphaVerbose, $xAlphaVerbose = False];

(* ========================================================================== *)
(* Load sub-packages                                                          *)
(* ========================================================================== *)
(* NOTE: Setup.wl is NOT loaded here because it contains top-level DefTensor *)
(* and DefProjectedTensor calls that require the user's manifold, metric,    *)
(* and xPand slicing to already be defined. Load it manually in your         *)
(* notebook AFTER DefManifold / SetSlicing / DefMetricFields, via:           *)
(*   Get[FileNameJoin[{$xAlphaDir, "Setup.wl"}]]                             *)
(* ========================================================================== *)

$xAlphaDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{$xAlphaDir, "HorndeskiDefinitions.wl"}]];
Get[FileNameJoin[{$xAlphaDir, "Background.wl"}]];
Get[FileNameJoin[{$xAlphaDir, "Perturbations.wl"}]];
Get[FileNameJoin[{$xAlphaDir, "AlphaFunctions.wl"}]];

EndPackage[]
