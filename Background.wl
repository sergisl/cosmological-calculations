(* ::Package:: *)

(* ========================================================================== *)
(* Horndeski Background EOMs Package                                         *)
(* ========================================================================== *)
(* Provides tools for computing background Equations of Motion (EOMs),       *)
(* their derivatives, and simplification rules (Scalar EOM, rho/p subs).     *)
(* ========================================================================== *)

(* Sub-package loaded by xAct`xAlpha` — no nested BeginPackage *)

GetHorndeskiRules::usage = 
  "GetHorndeskiRules[] returns the list of replacement rules to convert abstract Derivatives of Horndeski functions (like Derivative[1][G4][phi]) into named symbols (like G4phi).";

ComputeBackgroundEOMs::usage =
  "ComputeBackgroundEOMs[lagrangian] is a placeholder for the background EOM computation logic usually performed in the main notebook using xPand.";

ComputeBackgroundQuantityDerivatives::usage =
  "ComputeBackgroundQuantityDerivatives[varEexpr, varPexpr, varSexpr, rules] computes all required derivatives of the background quantities (varE, varP, varS) with respect to \[CurlyPhi], X, and time. Returns an Association.";

(* NEW EXPORTS *)
ApplyBackgroundEOMSimplifications::usage = 
  "ApplyBackgroundEOMSimplifications[expr, phidotToX] simplifies an expression by applying the Scalar EOM (eliminating phiddot/phidddot) and Background EOMs (varE -> -rho, varP -> -p).";

GetSimpleBackgroundEOMRules::usage = 
  "GetSimpleBackgroundEOMRules[] returns simple substitution rules for energy and pressure: {varE -> -trho, varP -> -tp}. Used by ApplyBackgroundEOMSimplifications. Note: the xAct-tensor form (matching \[Rho]u[LI[0],LI[0]] etc.) is GetBackgroundEOMRules[] from AlphaFunctions.wl.";

GetScalarEOMRules::usage = 
  "GetScalarEOMRules[phidotToX] solves the global scalar equation (varS=0) to return replacement rules for phiddot and phidddot.";

xAlphaLoadOrCompute::usage =
  "xAlphaLoadOrCompute[\"name\", expr] loads the pre-computed result for \"name\" from the xAlpha Cache directory \
if an .mx file exists there; otherwise evaluates expr (held until needed), saves the result to Cache, and returns it. \
Requires $xAlphaDir to be set (done automatically when xAlpha.m is loaded).";

xAlphaSaveCache::usage =
  "xAlphaSaveCache[\"name\", expr] saves expr to the xAlpha Cache directory as \"name\".mx \
and returns expr. Use this to manually refresh a cached result.";

Begin["xAct`xAlpha`Private`"]

(* ========================================================================== *)
(* 1. Force Global Context for Horndeski Symbols                              *)
(* ========================================================================== *)

KK = Symbol["Global`KK"];
G3 = Symbol["Global`G3"];
G4 = Symbol["Global`G4"];

Kphi = Symbol["Global`Kphi"]; KX = Symbol["Global`KX"]; Kphiphi = Symbol["Global`Kphiphi"];
KXX = Symbol["Global`KXX"]; KXXX = Symbol["Global`KXXX"]; KphiX = Symbol["Global`KphiX"];
KphiXX = Symbol["Global`KphiXX"]; KphiphiX = Symbol["Global`KphiphiX"];
KphiphiphiX = Symbol["Global`KphiphiphiX"]; KphiphiXX = Symbol["Global`KphiphiXX"];
KphiXXX = Symbol["Global`KphiXXX"]; KXXXX = Symbol["Global`KXXXX"];

G3phi = Symbol["Global`G3phi"]; G3X = Symbol["Global`G3X"]; G3phiphi = Symbol["Global`G3phiphi"];
G3XX = Symbol["Global`G3XX"]; G3XXX = Symbol["Global`G3XXX"]; G3phiX = Symbol["Global`G3phiX"];
G3phiXX = Symbol["Global`G3phiXX"]; G3phiphiphi = Symbol["Global`G3phiphiphi"];
G3phiphiX = Symbol["Global`G3phiphiX"]; G3phiphiphiphi = Symbol["Global`G3phiphiphiphi"];
G3phiphiphiX = Symbol["Global`G3phiphiphiX"]; G3phiphiXX = Symbol["Global`G3phiphiXX"];
G3phiXXX = Symbol["Global`G3phiXXX"]; G3XXXX = Symbol["Global`G3XXXX"];

G4phi = Symbol["Global`G4phi"]; G4phiphi = Symbol["Global`G4phiphi"];
G4phiphiphi = Symbol["Global`G4phiphiphi"]; G4phiphiphiphi = Symbol["Global`G4phiphiphiphi"];

(* ========================================================================== *)
(* 2. Horndeski Derivative Replacement Rules                                  *)
(* ========================================================================== *)

ClearAll[horndeskiSymbolFromOrder, horndeskiOrderMap, horndeskiRulesK, 
  horndeskiRulesG3, horndeskiRulesG4, HorndeskiDerivativeReplacements1];

horndeskiSymbolFromOrder[base_, m_, n_] := Module[{name, sym},
  Which[
    base === "K" && m == 0 && n == 0, KK,
    base === "G3" && m == 0 && n == 0, G3,
    base === "G4" && m == 0 && n == 0, G4,
    True,
    name = base <> StringRepeat["phi", m] <> StringRepeat["X", n];
    sym = Symbol["Global`" <> name];
    sym
  ]
];

horndeskiOrderMap[base_, symbols_List] := Module[{baseLen},
  baseLen = StringLength[base];
  Association@
   Map[
    Function[sym,
      sym -> If[sym === KK || sym === G3 || sym === G4,
        {0, 0},
        With[{suffix = StringDrop[SymbolName[sym], baseLen]},
          {StringCount[suffix, "phi"], StringCount[suffix, "X"]}
        ]
      ]
    ],
    symbols
   ]
];

$KSymbols = {
   KK, Kphi, KX, Kphiphi, KphiX, KXX, Symbol["Global`Kphiphiphi"], KphiphiX, KphiXX, KXXX,
   KphiphiphiX, KphiphiXX, KphiXXX, KXXXX
};
$G3Symbols = {
   G3, G3phi, G3X, G3phiphi, G3phiX, G3XX, G3phiphiphi, G3phiphiX, G3phiXX,
   G3XXX, G3phiphiphiphi, G3phiphiphiX, G3phiphiXX, G3phiXXX, G3XXXX
};
$G4Symbols = {G4, G4phi, G4phiphi, G4phiphiphi, G4phiphiphiphi};

$KOrder = horndeskiOrderMap["K", $KSymbols];
$G3Order = horndeskiOrderMap["G3", $G3Symbols];
$G4Order = horndeskiOrderMap["G4", $G4Symbols];

horndeskiRulesK := Module[{phiSym, XSym},
  phiSym = Symbol["xAct`xPand`\[CurlyPhi]"];
  XSym = Symbol["Global`X"];
  {
   KK[pphi_, XX_] :> KK[phiSym[], XSym[]],
   Derivative[dp_, dx_][sym_][pphi_, XX_] /; KeyExistsQ[$KOrder, sym] :>
    With[{ord = $KOrder[sym] + {dp, dx},
      tgt = horndeskiSymbolFromOrder["K", Sequence @@ ($KOrder[sym] + {dp, dx})]},
     If[KeyExistsQ[$KOrder, tgt], tgt[phiSym[], XSym[]], Derivative[dp, dx][sym][pphi, XX]]
    ]
   }
];

horndeskiRulesG3 := Module[{phiSym, XSym},
  phiSym = Symbol["xAct`xPand`\[CurlyPhi]"];
  XSym = Symbol["Global`X"];
  {
   G3[pphi_, XX_] :> G3[phiSym[], XSym[]],
   Derivative[dp_, dx_][sym_][pphi_, XX_] /; KeyExistsQ[$G3Order, sym] :>
    With[{ord = $G3Order[sym] + {dp, dx},
      tgt = horndeskiSymbolFromOrder["G3", Sequence @@ ($G3Order[sym] + {dp, dx})]},
     If[KeyExistsQ[$G3Order, tgt], tgt[phiSym[], XSym[]], Derivative[dp, dx][sym][pphi, XX]]
    ]
   }
];

horndeskiRulesG4 := Module[{phiSym},
  phiSym = Symbol["xAct`xPand`\[CurlyPhi]"];
  {
   G4[pphi_] :> G4[phiSym[]],
   Derivative[dp_][sym_][pphi_] /; KeyExistsQ[$G4Order, sym] :>
    With[{ord = $G4Order[sym] + {dp, 0},
      tgt = horndeskiSymbolFromOrder["G4", Sequence @@ ($G4Order[sym] + {dp, 0})]},
     If[KeyExistsQ[$G4Order, tgt], tgt[phiSym[]], Derivative[dp][sym][pphi]]
    ]
   }
];

HorndeskiDerivativeReplacements1 := Flatten[{horndeskiRulesK, horndeskiRulesG3, horndeskiRulesG4}];
GetHorndeskiRules[] := HorndeskiDerivativeReplacements1;

ComputeBackgroundEOMs[lagrangian_] := Module[{}, <|"Note" -> "Use xPand in notebook."|>];


(* ========================================================================== *)
(* 3. Compute Derivatives of Background Quantities                            *)
(* ========================================================================== *)

ComputeBackgroundQuantityDerivatives[varEexpr_, varPexpr_, varSexpr_, rules_] := Module[
  {phi, X, Hh, LI, Scalar, 
   rhsE, rhsP, rhsS,
   varEphi, varPphi, varSphi, 
   varEphiphi, varEphiphiphi, varPphiphi,
   varEphiX, varPphiX, 
   varEdot, varPdot, varSdot, 
   varEddot, varPddot, varSddot,
   varEdddot, varPdddot, varSdddot,
   scalarDerivRule, phidot, phiddot, phidddot, TimeDeriv},
  
  (* FIX: Use the xPand symbol for phi to match G-function arguments *)
  phi = Symbol["xAct`xPand`\[CurlyPhi]"]; 
  
  X = Symbol["Global`X"];
  Hh = Symbol["Global`Hh"];
  LI = Symbol["xAct`xTensor`LI"];
  Scalar = Symbol["xAct`xTensor`Scalar"];
  
  (* Rule to clean up Derivative[1][Scalar][...] -> 1 *)
  scalarDerivRule = {Derivative[1][Scalar][_] -> 1, Scalar[x_] :> x};
  
  (* Helper to extract RHS of equation *)
  GetRHS[eq_] := If[Head[eq] === Equal, eq[[2]], eq];

  rhsE = GetRHS[varEexpr] /. rules;
  rhsP = GetRHS[varPexpr] /. rules;
  rhsS = GetRHS[varSexpr] /. rules;
  
  (* --- 1. Phi Derivatives --- *)
  varEphi = System`D[rhsE, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  varPphi = System`D[rhsP, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  varSphi = System`D[rhsS, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  
  varEphiphi = System`D[varEphi, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  varPphiphi = System`D[varPphi, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  varEphiphiphi = System`D[varEphiphi, phi[LI[0], LI[0]]] /. scalarDerivRule /. rules;
  
  (* --- 2. Mixed Phi-X Derivatives --- *)
  phidot = phi[LI[0], LI[1]];
  varEphiX = (System`D[varEphi, X[]] + System`D[varEphi, phidot] * (1/phidot)) /. scalarDerivRule /. rules;
  varPphiX = (System`D[varPphi, X[]] + System`D[varPphi, phidot] * (1/phidot)) /. scalarDerivRule /. rules;
  
  (* --- 3. Time Derivatives --- *)
  phiddot = phi[LI[0], LI[2]];
  phidddot = phi[LI[0], LI[3]];
  
  TimeDeriv[expr_] := (
      System`D[expr, phi[LI[0], LI[0]]] * phidot +
      System`D[expr, phi[LI[0], LI[1]]] * phiddot +
      System`D[expr, phi[LI[0], LI[2]]] * phidddot +
      System`D[expr, X[]] * (phidot * phiddot) + 
      System`D[expr, Hh[LI[0], LI[0]]] * Hh[LI[0], LI[1]] +
      System`D[expr, Hh[LI[0], LI[1]]] * Hh[LI[0], LI[2]]
  ) /. scalarDerivRule /. rules;
             
  varEdot = TimeDeriv[rhsE];
  varPdot = TimeDeriv[rhsP];
  varSdot = TimeDeriv[rhsS];
  
  (* --- 4. Second Time Derivatives --- *)
  varEddot = TimeDeriv[varEdot];
  varPddot = TimeDeriv[varPdot];
  varSddot = TimeDeriv[varSdot];
  
  (* --- 5. Third Time Derivatives --- *)
  varEdddot = TimeDeriv[varEddot];
  varPdddot = TimeDeriv[varPddot];
  varSdddot = TimeDeriv[varSddot];
  
  (* Return as association *)
  <|
    "varEphi" -> (Symbol["Global`varEphi"] == varEphi),
    "varPphi" -> (Symbol["Global`varPphi"] == varPphi),
    "varSphi" -> (Symbol["Global`varSphi"] == varSphi),
    "varEphiphi" -> (Symbol["Global`varEphiphi"] == varEphiphi),
    "varEphiphiphi" -> (Symbol["Global`varEphiphiphi"] == varEphiphiphi),
    "varPphiphi" -> (Symbol["Global`varPphiphi"] == varPphiphi),
    "varEphiX" -> (Symbol["Global`varEphiX"] == varEphiX),
    "varPphiX" -> (Symbol["Global`varPphiX"] == varPphiX),
    "varEdot" -> (Symbol["Global`varEdot"] == varEdot),
    "varPdot" -> (Symbol["Global`varPdot"] == varPdot),
    "varSdot" -> (Symbol["Global`varSdot"] == varSdot),
    "varEddot" -> (Symbol["Global`varEddot"] == varEddot),
    "varPddot" -> (Symbol["Global`varPddot"] == varPddot),
    "varSddot" -> (Symbol["Global`varSddot"] == varSddot),
    "varEdddot" -> (Symbol["Global`varEdddot"] == varEdddot),
    "varPdddot" -> (Symbol["Global`varPdddot"] == varPdddot),
    "varSdddot" -> (Symbol["Global`varSdddot"] == varSdddot)
  |>
];

(* ========================================================================== *)
(* 4. Background Simplifications (Scalar EOM & Substitution)                  *)
(* ========================================================================== *)

GetSimpleBackgroundEOMRules[] := GetSimpleBackgroundEOMRules[] = 
  Module[{trho, tp, varE, varP},
    trho = Symbol["Global`trho"][];
    tp = Symbol["Global`tp"][];
    varE = Symbol["Global`varE"]; 
    varP = Symbol["Global`varP"]; 
    {varE -> -trho, varP -> -tp}
  ];

GetScalarEOMRules[phidotToX_] := 
  Module[{varSexpr, varSdotexpr, phiddot, phidddot, phi, LI, solvePhiddot, solvePhidddot, rules},
    
    varSexpr = Symbol["Global`varSexpr"];
    varSdotexpr = Symbol["Global`varSdotexpr"];
    
    (* Use xPand phi to match the derivatives *)
    phi = Symbol["xAct`xPand`\[CurlyPhi]"];
    LI = Symbol["xAct`xTensor`LI"];
    
    phiddot = phi[LI[0], LI[2]];
    phidddot = phi[LI[0], LI[3]];
    
    (* 1. Solve varS = 0 for phiddot (Eliminate 2nd derivative) *)
    (* We assume varSexpr is "varS == ...". We take the RHS. *)
    solvePhiddot = Solve[(If[Head[varSexpr]===Equal, varSexpr[[2]], varSexpr] //. phidotToX) == 0, phiddot];
    
    If[Length[solvePhiddot] == 0, 
       If[$xAlphaVerbose, Print["Warning: Could not solve varS for phiddot (scalar EOM)."]];
       solvePhiddot = {}
    ];
    
    (* 2. Solve varSdot = 0 for phidddot (Eliminate 3rd derivative) *)
    solvePhidddot = {};
    If[ValueQ[varSdotexpr],
       (* varSdot depends on phiddot. Substitute phiddot solution first. *)
       solvePhidddot = Solve[(If[Head[varSdotexpr]===Equal, varSdotexpr[[2]], varSdotexpr] //. phidotToX //. solvePhiddot) == 0, phidddot]
    ];

    (* Combine and simplify *)
    rules = Flatten[{solvePhiddot, solvePhidddot}];
    rules = rules //. phidotToX // Simplify;
    
    rules
  ];

(* Master function to apply Background EOMs *)
ApplyBackgroundEOMSimplifications[expr_, phidotToX_] :=
  Module[{res, scalarRules, backgroundRules, scalarDerivRule},
    
    (* 1. Clean artifacts *)
    scalarDerivRule = {Derivative[1][Symbol["xAct`xTensor`Scalar"]][_] -> 1, 
                       Symbol["xAct`xTensor`Scalar"][x_] :> x};
    res = expr /. scalarDerivRule;
    
    (* 2. Apply Scalar EOM (eliminate phiddot, phidddot) *)
    scalarRules = GetScalarEOMRules[phidotToX];
    res = res //. scalarRules;
    
    (* 3. Apply Background EOM (varE -> -rho, varP -> -p) *)
    backgroundRules = GetSimpleBackgroundEOMRules[];
    res = res //. backgroundRules;
    
    Simplify[res]
  ];

(* ========================================================================== *)
(* Cache Utilities                                                            *)
(* ========================================================================== *)

(* Internal helper: full path to a cache file *)
xAlphaCacheFile[name_String] :=
  FileNameJoin[{Symbol["xAct`xAlpha`$xAlphaDir"], "Cache", name <> ".mx"}];

(* Load from Cache/name.mx if it exists; otherwise evaluate computeExpr,
   save to cache, and return the result. *)
SetAttributes[xAlphaLoadOrCompute, HoldRest];
xAlphaLoadOrCompute[name_String, computeExpr_] :=
  Module[{cacheFile, result},
    cacheFile = xAlphaCacheFile[name];
    If[FileExistsQ[cacheFile],
      If[$xAlphaVerbose, Print["[xAlpha] Loading cached ", name, " from ", cacheFile]];
      Get[cacheFile];
      Symbol["Global`" <> name]
      ,
      If[$xAlphaVerbose, Print["[xAlpha] Computing ", name, " (no cache found)..."]];
      result = computeExpr;
      xAlphaSaveCache[name, result];
      result
    ]
  ];

(* Save expr to Cache/name.mx and return expr. *)
xAlphaSaveCache[name_String, expr_] :=
  Module[{cacheFile, tmpSym},
    cacheFile = xAlphaCacheFile[name];
    (* Assign to a Global symbol matching the name so DumpSave works *)
    tmpSym = Symbol["Global`" <> name];
    tmpSym = expr;
    DumpSave[cacheFile, tmpSym];
    If[$xAlphaVerbose, Print["[xAlpha] Saved cache: ", cacheFile]];
    expr
  ];

End[]