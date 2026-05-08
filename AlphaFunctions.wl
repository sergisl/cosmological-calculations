(* ::Package:: *)

(* ========================================================================== *)
(* Horndeski Alpha Functions Package                                          *)
(* ========================================================================== *)

(* Sub-package loaded by xAct`xAlpha` — no nested BeginPackage *)

(* --- Derivatives & Base Alphas --- *)
ComputeAlphaM::usage = "ComputeAlphaM[G4, G4phi, Hh, phidot] computes the kinetic braiding parameter \[Alpha]_M = 2 \[CurlyPhi]' G_{4\[CurlyPhi]} / (2 H G_4), where primes denote cosmic-time derivatives. Returns a symbolic expression in the Horndeski functions.";
ComputeAlphaK::usage = "ComputeAlphaK[G4, G3phi, G3phiX, G3X, G3XX, KX, KXX, Hh, phidot, X] computes the kineticity parameter \[Alpha]_K from the Horndeski functions K, G3, G4 and their derivatives evaluated at the background. Returns a symbolic expression.";
ComputeAlphaB::usage = "ComputeAlphaB[G4, G4phi, G3X, Hh, phidot, X] computes the braiding parameter \[Alpha]_B = \[CurlyPhi]' (X G_{3X} - G_{4\[CurlyPhi]}) / (H G_4). Returns a symbolic expression.";

ComputeAlphaTimeDerivative::usage = "ComputeAlphaTimeDerivative[alphaExpr, phi, phidot, phiddot, Hh, Hdot, X, replacementRules] computes the cosmic-time derivative of an alpha parameter expression using the chain rule: d\[Alpha]/dt = (\[PartialD]\[Alpha]/\[PartialD]\[CurlyPhi])\[CurlyPhi]' + (\[PartialD]\[Alpha]/\[PartialD]\[CurlyPhi]')\[CurlyPhi]'' + (\[PartialD]\[Alpha]/\[PartialD]X)X' + (\[PartialD]\[Alpha]/\[PartialD]H)H'. alphaExpr may be a bare expression or an Equal equation (rhs is used). replacementRules should be HorndeskiDerivativeReplacements1. Returns a symbolic expression.";
ComputeAlphaXDerivative::usage = "ComputeAlphaXDerivative[alphaExpr, X, phidot, replacementRules] computes the X-derivative of an alpha parameter expression, accounting for the dependence of X = \[CurlyPhi]'^2/2 on \[CurlyPhi]'. alphaExpr may be a bare expression or an Equal equation (rhs is used). replacementRules should be HorndeskiDerivativeReplacements1. Returns a symbolic expression.";
ComputeAlphaPhiDerivative::usage = "ComputeAlphaPhiDerivative[alphaExpr, phi, replacementRules] computes the \[CurlyPhi]-derivative of an alpha parameter expression. alphaExpr may be a bare expression or an Equal equation (rhs is used). replacementRules should be HorndeskiDerivativeReplacements1. Returns a symbolic expression.";
ComputeAllAlphas::usage = "ComputeAllAlphas[phi, X, Hh, HorndeskiDerivativeReplacements1, phidotToX] computes and stores all EFT alpha parameters and their derivatives as Global` equations. Populates: Mstarexpr, aMexpr, aKexpr, aBexpr; first/second/third time derivatives (aBdotexpr, etc.); X-derivatives (aBXexpr, aKXexpr, aMXexpr, ...); \[CurlyPhi]-derivatives (aBphiexpr, ...); and mixed time+X derivatives (aBXdotexpr, etc.).";

(* --- Transformations --- *)
RemoveK::usage = "RemoveK[expr, phidotToX] substitutes all Horndeski K-function derivatives (KX, KXX, Kphi, KphiX, ...) in expr with the corresponding combinations of EFT alpha parameters (\[Alpha]_K, \[Alpha]_B and their derivatives), using equations previously computed by ComputeAllAlphas[]. Returns {transformedExpr, substitutionRules}.";
RemoveG3::usage = "RemoveG3[expr, phidotToX] substitutes all G3 derivatives (G3X, G3XX, G3phi, G3phiX, ...) in expr with combinations of EFT alpha parameters (\[Alpha]_B and its X and \[CurlyPhi] derivatives), using equations from ComputeAllAlphas[]. Returns {transformedExpr, substitutionRules}.";
RemoveG4::usage = "RemoveG4[expr, phidotToX] substitutes all G4 derivatives (G4, G4phi, G4phiphi, ...) in expr with the Planck mass Mstar[] and its time derivatives (\[Alpha]_M and its derivatives), using equations from ComputeAllAlphas[]. Returns {transformedExpr, substitutionRules}.";
GToAlphas::usage = "GToAlphas[expr, phidotToX] applies the full transformation chain RemoveK \[Rule] RemoveG3 \[Rule] RemoveG4 to expr, replacing all Horndeski G-function derivatives with EFT alpha parameters. Returns the transformed expression. Convenience wrapper around RemoveK, RemoveG3, RemoveG4.";

ApplyConservationEOMRules::usage = "ApplyConservationEOMRules[] returns a list of replacement rules derived from the background energy and momentum conservation equations, expressing \!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(.\)]\), \!\(\*OverscriptBox[\(\[ScriptCapitalP]\), \(.\)]\) and higher derivatives in terms of \[ScriptCapitalE], \[ScriptCapitalP], H and the alpha parameters. Apply with //. to an expression to eliminate time derivatives of the fluid variables.";
SimplifyWithConservation::usage = "SimplifyWithConservation[expr] applies ApplyConservationEOMRules[] to expr and calls Simplify[]. Convenience wrapper.";
ApplyEnergyRescalingRules::usage = "ApplyEnergyRescalingRules[] returns replacement rules that rescale the dimensionful background EOM combinations varE, varP, varS (and their derivatives) by Mstar[]^2, rewriting them as Mstar[]^2 * tvarE, etc. Apply after GToAlphas[] to put coefficients in standard dimensionless form.";

(* --- Extraction --- *)
ExtractAndTransformCoefficients::usage = "ExtractAndTransformCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts and transforms the three linear coefficients A[1,h], A[2,h], A[3,h] from the (0,0) metric perturbation equation collectedEqn for each perturbation head h in pertHeads. The three coefficients correspond to: A[1,h] (prefactor 1/H^2, multiplies h[LI[1],LI[0]]), A[2,h] (prefactor 1/H, multiplies h'), A[3,h] (prefactor -a^2, multiplies \[Del]^2 h). Stores raw (Horndeski) results as A[i,h], alpha-language results as Aalpha[i,h], and gamma-language results as Agamma[i,h] in Global` context.";
ExtractAndTransformQuadraticCoefficients::usage = "ExtractAndTransformQuadraticCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts the six quadratic A coefficients A[1..6, h1, h2] from the (0,0) equation for all pairs {h1,h2} in pertHeads. Stores raw/alpha/gamma results in A[i,h1,h2], Aalpha[i,h1,h2], Agamma[i,h1,h2].";
ExtractAndTransformBCoefficients::usage = "ExtractAndTransformBCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts the two linear B coefficients B[1,h], B[2,h] from the (0,i) momentum constraint equation. B[1,h]: prefactor 1/H, multiplies \[Del]_c h; B[2,h]: prefactor -1, multiplies \[Del]_c h'. Stores results in B[i,h], Balpha[i,h], Bgamma[i,h].";
ExtractAndTransformQuadraticBCoefficients::usage = "ExtractAndTransformQuadraticBCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts quadratic B coefficients from the (0,i) equation for all pairs {h1,h2}. Stores results in B[i,h1,h2], Balpha[i,h1,h2], Bgamma[i,h1,h2].";
ExtractAndTransformCCoefficients::usage = "ExtractAndTransformCCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts linear C coefficients from the (i,j) traceless spatial equation. Stores results in CC[i,h], CCalpha[i,h], Cgamma[i,h].";
ExtractAndTransformQuadraticCCoefficients::usage = "ExtractAndTransformQuadraticCCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts quadratic C coefficients from the (i,j) equation. Stores results in CC[i,h1,h2], CCalpha[i,h1,h2], Cgamma[i,h1,h2].";
ExtractAndTransformDCoefficients::usage = "ExtractAndTransformDCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts linear D coefficients from the scalar field perturbation equation. Stores results in DD[i,h], DDalpha[i,h], Dgamma[i,h].";
ExtractAndTransformQuadraticDCoefficients::usage = "ExtractAndTransformQuadraticDCoefficients[collectedEqn, pertHeads, phidotToX, removeScalarWrapper] extracts quadratic D coefficients from the scalar field equation. Stores results in DD[i,h1,h2], DDalpha[i,h1,h2], Dgamma[i,h1,h2].";

CanonicalizeGradientProducts::usage = "CanonicalizeGradientProducts[expr] canonicalizes index placement in gradient products to ensure unique representation: cd[upper][tder] cd[lower][base] -> cd[lower][tder] cd[upper][base].";

DefineCollectionPatterns::usage = "DefineCollectionPatterns[hh] defines patterns for term collection like lap[hh], gradDown[hh], etc.";
DefineCollectionPatternsQuadratic::usage = "DefineCollectionPatternsQuadratic[h1, h2] defines quadratic patterns like gpair[h1, h2], deltagpair[h1, h2], etc.";
DefineAllCollectionPatterns::usage = "DefineAllCollectionPatterns[pertHeadsList] defines all linear and quadratic patterns for the given perturbation heads.";
GetCollectionRules::usage = "GetCollectionRules[pertHeadsList] returns replacement rules to convert expanded terms to pattern notation for display.";

ConvertAllCoefficientsToGammaLanguage::usage = "ConvertAllCoefficientsToGammaLanguage[pertHeadsList] converts extracted Alpha coefficients into Gamma basis.";

(* --- Validation --- *)
ReconstructEquationFromCoefficients::usage = "ReconstructEquationFromCoefficients[equationType, pertHead, removeScalarWrapper] reconstructs equation from extracted raw coefficients. equationType: \"linear-A\", \"quadratic-A\", \"linear-B\", \"quadratic-B\", \"linear-C\", \"quadratic-C\", \"linear-D\", \"quadratic-D\". removeScalarWrapper: rules to remove Scalar[] wrappers.";
ValidateCoefficients::usage = "ValidateCoefficients[originalEqn, equationType, pertHead(s), coeffType] validates that extracted coefficients fully reconstruct the original equation. Returns {isComplete, residual}.";
ValidateAllCoefficients::usage = "ValidateAllCoefficients[equations, pertHeads, coeffType] systematically validates all coefficients. equations is an Association with keys: \"eq00Linear\", \"eq00Quadratic\", \"eq0iLinear\", \"eq0iQuadratic\", \"eqijLinear\", \"eqijQuadratic\", \"eqScalarLinear\", \"eqScalarQuadratic\". Returns {results, failedValidations}.";
ValidateFullLinearEquation::usage = "ValidateFullLinearEquation[equationType, originalEqn, pertHeads, removeScalarWrapper] validates full linear equation reconstruction and returns {isComplete, residual}.";
ValidateFullQuadraticEquation::usage = "ValidateFullQuadraticEquation[equationType, originalEqn, pertHeads, removeScalarWrapper] validates full quadratic equation reconstruction and returns {isComplete, residual}.";
ValidateAllFullEquations::usage = "ValidateAllFullEquations[equations, pertHeads, removeScalarWrapper] validates full reconstructions for linear/quadratic A,B,C,D and returns an Association with results and summary.";

ShowPaperTables::usage =
  "ShowPaperTables[linA, linB, linC, linD, quadA, quadB, quadC, quadD] prints two \n\
formatted Grid tables summarising the linear (first-order) and quadratic \n\
(second-order) EFT coefficients extracted from the A/B/C/D equation arrays. \n\
Green entries indicate non-zero coefficients; grey entries are zero. \n\
Typically called as ShowPaperTables[Aalpha, Balpha, CCalpha, DDalpha, \n\
  Aalpha, Balpha, CCalpha, DDalpha].";

GetBackgroundEOMRules::usage = "GetBackgroundEOMRules[] returns replacement rules derived from the Friedmann and scalar field background equations of motion, expressing combinations of H, H', H'', \[CurlyPhi]', \[CurlyPhi]'' in terms of the fluid energy/pressure varE, varP and alpha parameters. Apply with /. to simplify background-dependent expressions.";
GenerateAlphaToGammaRules::usage = "GenerateAlphaToGammaRules[] constructs and returns the full list of replacement rules that convert alpha-language parameters (\[Alpha]_K, \[Alpha]_B, \[Alpha]_M, \[Alpha]_{BX}, ...) into gamma-language parameters (\[Gamma]_K, \[Gamma]_B, \[Gamma]_D, ...) via the algebraic relations between the two EFT bases. Returns a flat list of Rule objects.";
ApplyGammaSubstitution::usage = "ApplyGammaSubstitution[expr, gammaRules, eomRules] applies gammaRules (from GenerateAlphaToGammaRules[]) to expr, optionally also applying eomRules for background EOM simplification. Returns the expression in the gamma basis.";
ConvertAllCoefficientsToGammaLanguage::usage = "ConvertAllCoefficientsToGammaLanguage[pertHeadsList] applies GenerateAlphaToGammaRules[] to all stored Aalpha[i,h], Balpha[i,h], CCalpha[i,h], DDalpha[i,h] coefficients (linear and quadratic) for each head in pertHeadsList, updating the corresponding Agamma, Bgamma, Cgamma, Dgamma arrays in Global` context.";
TransformCoeffToAlphaLanguage::usage = "TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper] applies the full pipeline: RemoveK \[Rule] RemoveG3 \[Rule] RemoveG4 \[Rule] background EOM rules \[Rule] conservation \[Rule] energy rescaling, with intermediate Scalar[] wrapper removal and phidotToX normalisation at each stage. Returns the final expression in alpha-parameter language, ready to store as an Aalpha/Balpha/etc. coefficient. This is the internal workhorse called by all ExtractAndTransform* functions.";
ToArgAgnosticRules::usage = "ToArgAgnosticRules[rules] converts a list of rules of the form f[a,b]->rhs to HoldPattern[f[_,_]]->rhs, so that each rule fires regardless of the concrete argument form (bare symbol, indexed tensor, function call, etc.). Used internally by RemoveK, RemoveG3, RemoveG4 to ensure Horndeski function substitutions match in all contexts.";
SafeCoefficient::usage = "SafeCoefficient[expr, patt] returns Coefficient[expr, patt] unless patt===0, in which case it returns 0 (avoiding the Mathematica gotcha where Coefficient[expr,0] returns the full expression). Used internally by the ExtractAndTransform* coefficient extraction functions.";
SafeRuleList::usage = "SafeRuleList[rules] validates a rule list: flattens it and returns it only if every element is a Rule or RuleDelayed; returns {} otherwise. Used as a guard when passing rule sets between functions.";

(* Forward declarations: anchor these cross-file symbols in the public
   xAct`xAlpha` context before Begin[Private], so that references inside
   the Private block always resolve to the same symbol as the definition. *)
GetRemoveScalarWrapper;
GetHorndeskiRules;
GetBackgroundEOMRules;

Begin["xAct`xAlpha`Private`"];

(* ========================================================================== *)
(* Helper: solve for the first concrete instance of a given head in eq.       *)
(* SolveForHead finds e.g. KXX[phi0, X0] from the equation and solves for it.*)
(* ========================================================================== *)
SolveForHead[eq_, head_] :=
  Module[{instances},
    instances = DeleteDuplicates @ With[{h = head}, Cases[eq, _h, Infinity]];
    If[Length[instances] === 0, {}, Solve[eq, First @ instances]]
  ];

(* ========================================================================== *)
(* Helper: full K→G3→G4→conservation→rescaling pipeline.                      *)
(* ========================================================================== *)
TransformCoeffToAlphaLanguage[coeff_, phidotToX_, removeScalarWrapper_] :=
  Module[{c = coeff, rules},
    {c, rules} = RemoveK[c, phidotToX];
    c = c //. removeScalarWrapper // Simplify;
    c = c //. phidotToX // Simplify;
    {c, rules} = RemoveG3[c, phidotToX];
    c = c //. removeScalarWrapper // Simplify;
    c = c //. phidotToX // Simplify;
    {c, rules} = RemoveG4[c, phidotToX];
    c = c //. removeScalarWrapper // Simplify;
    c = c //. phidotToX // Simplify;
    c = c /. GetBackgroundEOMRules[] //. ApplyConservationEOMRules[] //. ApplyEnergyRescalingRules[] // Simplify;
    ApplyReverseXRuleToCoeff[Simplify[c]]
  ];

(* ========================================================================== *)
(* 0. Background EOM & Alpha EOM Caching Logic                                *)
(* ========================================================================== *)

$AlphaEOMRulesCache = Null;

GetBackgroundEOMRules[] := GetBackgroundEOMRules[] = {
  Symbol["Global`\[Rho]u"][Symbol["xAct`xTensor`LI"][0], Symbol["xAct`xTensor`LI"][0]] -> -Symbol["Global`varE"],
  Symbol["Global`Pu"][Symbol["xAct`xTensor`LI"][0], Symbol["xAct`xTensor`LI"][0]] -> -Symbol["Global`varP"]
};

(* ========================================================================== *)
(* Conservation Equation Rules for Energy Density                             *)
(* ========================================================================== *)

ApplyConservationEOMRules[] := ApplyConservationEOMRules[] = Module[{Hh, varE, varP, varEdot, varPdot, varEddot, varPddot, varEdddot, varPdddot, H, Hdot, Hddot},
  Hh = Symbol["Global`Hh"]; varE = Symbol["Global`varE"]; varP = Symbol["Global`varP"];
  varEdot = Symbol["Global`varEdot"]; varPdot = Symbol["Global`varPdot"];
  varEddot = Symbol["Global`varEddot"]; varPddot = Symbol["Global`varPddot"];
  varEdddot = Symbol["Global`varEdddot"]; varPdddot = Symbol["Global`varPdddot"];
  H = Hh[Symbol["xAct`xTensor`LI"][0], Symbol["xAct`xTensor`LI"][0]];
  Hdot = Hh[Symbol["xAct`xTensor`LI"][0], Symbol["xAct`xTensor`LI"][1]];
  Hddot = Hh[Symbol["xAct`xTensor`LI"][0], Symbol["xAct`xTensor`LI"][2]];
  
  {
    (* Conservation Equation (continuity equation): Ė = -3H(E + P) *)
    varEdot -> -3 H (varE + varP),
    
    (* First time derivative of conservation equation: Ë = 3(E+P)(3H²-Ḣ) - 3H Ṗ *)
    varEddot -> 3 (varE + varP) (3 H^2 - Hdot) - 3 H varPdot,
    
    (* Second time derivative of conservation equation: E... = -3(E+P)(Ḧ-9HḢ+9H³) + 3Ṗ(3H²-2Ḣ) - 3HṖ *)
    varEdddot -> -3 (varE + varP) (Hddot - 9 H Hdot + 9 H^3) + 3 varPdot (3 H^2 - 2 Hdot) - 3 H varPdddot
  }
];

(* Convenience function to apply conservation rules to a coefficient *)
SimplifyWithConservation[coeffAlpha_] := 
  Simplify[coeffAlpha //. ApplyConservationEOMRules[]];

SafeRuleList[rules_] := Module[{r = Flatten[{rules}]},
  If[VectorQ[r, MatchQ[#, _Rule | _RuleDelayed] &], r, {}]
];

(* Shared private helper: avoids Coefficient[expr, 0] which returns the whole expr *)
SafeCoefficient[expr_, patt_] := If[patt === 0, 0, Coefficient[expr, patt]];

(* Shared private helper: convert f[a,b]->rhs rules to HoldPattern[f[_,_]]->rhs so *)
(* they fire regardless of the concrete argument form (bare X vs X[], indexed phi, etc.) *)
ToArgAgnosticRules[rules_List] :=
  rules /. Rule[f_[_, _], rhs_] :> Rule[HoldPattern[f[_, _]], rhs];

(* ========================================================================== *)
(* Rescaling Rules for Energy Density Variables                               *)
(* ========================================================================== *)

ApplyEnergyRescalingRules[] := ApplyEnergyRescalingRules[] = Module[{Mstar, varE, varP, varEdot, varPdot, varEddot, varPddot, varEdddot, varPdddot},
  Mstar = Symbol["Global`Mstar"]; 
  varE = Symbol["Global`varE"]; varP = Symbol["Global`varP"];
  varEdot = Symbol["Global`varEdot"]; varPdot = Symbol["Global`varPdot"];
  varEddot = Symbol["Global`varEddot"]; varPddot = Symbol["Global`varPddot"];
  varEdddot = Symbol["Global`varEdddot"]; varPdddot = Symbol["Global`varPdddot"];
  
  {
    (* Rescale energy density and pressure: E = Mstar^2 * tE *)
    varE -> Mstar[]^2 * Symbol["Global`tvarE"],
    varP -> Mstar[]^2 * Symbol["Global`tvarP"],
    
    (* Rescale time derivatives *)
    varEdot -> Mstar[]^2 * Symbol["Global`tvarEdot"],
    varPdot -> Mstar[]^2 * Symbol["Global`tvarPdot"],
    varEddot -> Mstar[]^2 * Symbol["Global`tvarEddot"],
    varPddot -> Mstar[]^2 * Symbol["Global`tvarPddot"],
    varEdddot -> Mstar[]^2 * Symbol["Global`tvarEdddot"],
    varPdddot -> Mstar[]^2 * Symbol["Global`tvarPdddot"]
  }
];

(* ========================================================================== *)
(* Canonicalization for gradient products                                    *)
(* ========================================================================== *)

CanonicalizeGradientProducts[expr_] := Module[{cd, LI, b, XX, YY, result},
  cd = Symbol["Global`cd"];
  LI = Symbol["xAct`xTensor`LI"];
  
  result = expr //. {
    (* D9-type: cd[upper][tder] cd[lower][base] -> cd[lower][tder] cd[upper][base] *)
    cd[b_][XX_[LI[1], LI[1]]] * cd[-b_][YY_[LI[1], LI[0]]] :> 
      cd[-b][XX[LI[1], LI[1]]] * cd[b][YY[LI[1], LI[0]]],
    (* Also catch opposite multiplicative ordering after canonical term sorting *)
    cd[-b_][YY_[LI[1], LI[0]]] * cd[b_][XX_[LI[1], LI[1]]] :> 
      cd[-b][XX[LI[1], LI[1]]] * cd[b][YY[LI[1], LI[0]]],
    (* Catch remaining order: cd[upper][base] cd[lower][tder] *)
    cd[b_][YY_[LI[1], LI[0]]] * cd[-b_][XX_[LI[1], LI[1]]] :>
      cd[-b][XX[LI[1], LI[1]]] * cd[b][YY[LI[1], LI[0]]],
    
    (* D11-type: cd[upper][tder1] cd[lower][tder2] -> cd[lower][tder1] cd[upper][tder2] *)
    cd[b_][XX_[LI[1], LI[1]]] * cd[-b_][YY_[LI[1], LI[1]]] :> 
      cd[-b][XX[LI[1], LI[1]]] * cd[b][YY[LI[1], LI[1]]],

    (* C8-type and similar: canonicalize field ordering for two lower-index gradients *)
    cd[idx1_][t1 : (XX_[LI[1], LI[0]] | _[XX_[LI[1], LI[0]]])] *
      cd[idx2_][t2 : (YY_[LI[1], LI[0]] | _[YY_[LI[1], LI[0]]])] /;
      OrderedQ[{idx1, idx2}] && ! OrderedQ[{XX, YY}] :>
      cd[idx1][t2] * cd[idx2][t1]
  };
  
  result
];

(* ========================================================================== *)

DefineCollectionPatterns[hh_] := Quiet[Block[{},
  (* Define simple evaluation rules with explicit Global context *)
  ToExpression["base[" <> SymbolName[Unevaluated[hh]] <> "] := " <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]"];
  ToExpression["tder[" <> SymbolName[Unevaluated[hh]] <> "] := " <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]"];
  ToExpression["ttder[" <> SymbolName[Unevaluated[hh]] <> "] := " <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[2]]"];
  ToExpression["lap[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`b][Global`cd[Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]]"];
  ToExpression["gradDown[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradUp[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradtDown[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]]"];
  ToExpression["gradgrad[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`c][Global`cd[-Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]]"];
  ToExpression["gradgradUp[" <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[Global`b][Global`cd[Global`c][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]]"];
  ToExpression["deltaTerm[" <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]"];
  ToExpression["deltaT[" <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]"];
  ToExpression["deltaTT[" <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[2]]"];
  ToExpression["deltaLap[" <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*Global`cd[-Global`d][Global`cd[Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]]"];
  
  (* Quadratic patterns *)
  ToExpression["gpair[" <> SymbolName[Unevaluated[hh]] <> ", " <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["deltagpair[" <> SymbolName[Unevaluated[hh]] <> ", " <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*Global`cd[-Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["deltagtpair[" <> SymbolName[Unevaluated[hh]] <> "] := Global`h[-Global`b, -Global`c]*Global`cd[-Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradagradb[" <> SymbolName[Unevaluated[hh]] <> ", " <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[-Global`c][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradagradtb[" <> SymbolName[Unevaluated[hh]] <> ", " <> SymbolName[Unevaluated[hh]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[-Global`c][" <> SymbolName[Unevaluated[hh]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]]"];
  
], {Set::write}];

DefineCollectionPatternsQuadratic[h1_, h2_] := Quiet[Block[{},
  (* Define simple evaluation rules with explicit Global context *)
  ToExpression["gpair[" <> SymbolName[Unevaluated[h1]] <> ", " <> SymbolName[Unevaluated[h2]] <> "] := Global`cd[-Global`d][" <> SymbolName[Unevaluated[h1]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[h2]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["deltagpair[" <> SymbolName[Unevaluated[h1]] <> ", " <> SymbolName[Unevaluated[h2]] <> "] := Global`h[-Global`b, -Global`c]*Global`cd[-Global`d][" <> SymbolName[Unevaluated[h1]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[h2]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["deltagtpair[" <> SymbolName[Unevaluated[h1]] <> ", " <> SymbolName[Unevaluated[h2]] <> "] := Global`h[-Global`b, -Global`c]*Global`cd[-Global`d][" <> SymbolName[Unevaluated[h1]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]]*Global`cd[Global`d][" <> SymbolName[Unevaluated[h2]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradagradb[" <> SymbolName[Unevaluated[h1]] <> ", " <> SymbolName[Unevaluated[h2]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[h1]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[-Global`c][" <> SymbolName[Unevaluated[h2]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]"];
  ToExpression["gradagradtb[" <> SymbolName[Unevaluated[h1]] <> ", " <> SymbolName[Unevaluated[h2]] <> "] := Global`cd[-Global`b][" <> SymbolName[Unevaluated[h1]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[0]]]*Global`cd[-Global`c][" <> SymbolName[Unevaluated[h2]] <> "[xAct`xTensor`LI[1], xAct`xTensor`LI[1]]]"];
  
], {Set::write}];

DefineAllCollectionPatterns[pertHeadsList_] := Module[{},
  (* Define linear patterns for each perturbation head *)
  Do[DefineCollectionPatterns[hh], {hh, pertHeadsList}];
  
  (* Define quadratic patterns for all pairs *)
  Do[
    DefineCollectionPatternsQuadratic[h1, h2],
    {h1, pertHeadsList}, {h2, pertHeadsList}
  ];
  
  Null
];

GetCollectionRules[pertHeadsList_] := Module[{LI, cd, h, rules},
  LI = Symbol["xAct`xTensor`LI"];
  cd = Symbol["Global`cd"];
  h = Symbol["Global`h"];
  
  rules = {};
  
  (* Linear patterns for each perturbation head *)
  Do[
    With[{hhVal = hh},
      rules = Join[rules, {
        (* Laplacian - match any dummy index *)
        cd[Pattern[idx1, Blank[]]][cd[-Pattern[idx1, Blank[]]][hhVal[LI[1], LI[0]]]] :> lap[hhVal],
        cd[-Pattern[idx1, Blank[]]][cd[Pattern[idx1, Blank[]]][hhVal[LI[1], LI[0]]]] :> lap[hhVal],
        (* Gradients *)
        cd[-Pattern[idx1, Blank[]]][hhVal[LI[1], LI[0]]] :> gradDown[hhVal],
        cd[Pattern[idx1, Blank[]]][hhVal[LI[1], LI[0]]] :> gradUp[hhVal],
        cd[-Pattern[idx1, Blank[]]][hhVal[LI[1], LI[1]]] :> gradtDown[hhVal],
        (* Double gradients *)
        cd[-Pattern[idx2, Blank[]]][cd[-Pattern[idx1, Blank[]]][hhVal[LI[1], LI[0]]]] :> gradgrad[hhVal],
        cd[Pattern[idx1, Blank[]]][cd[Pattern[idx2, Blank[]]][hhVal[LI[1], LI[0]]]] :> gradgradUp[hhVal],
        (* Delta terms *)
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*hhVal[LI[1], LI[0]] :> deltaTerm[hhVal],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*hhVal[LI[1], LI[1]] :> deltaT[hhVal],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*hhVal[LI[1], LI[2]] :> deltaTT[hhVal],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[-Pattern[idx3, Blank[]]][cd[Pattern[idx3, Blank[]]][hhVal[LI[1], LI[0]]]] :> deltaLap[hhVal],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[Pattern[idx3, Blank[]]][cd[-Pattern[idx3, Blank[]]][hhVal[LI[1], LI[0]]]] :> deltaLap[hhVal]
      }]
    ],
    {hh, pertHeadsList}
  ];
  
  (* Quadratic patterns for all pairs *)
  Do[
    With[{h1Val = h1, h2Val = h2},
      rules = Join[rules, {
        cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]]*cd[Pattern[idx1, Blank[]]][h2Val[LI[1], LI[0]]] :> gpair[h1Val, h2Val],
        cd[Pattern[idx1, Blank[]]][h2Val[LI[1], LI[0]]]*cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]] :> gpair[h1Val, h2Val],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[-Pattern[idx3, Blank[]]][h1Val[LI[1], LI[0]]]*cd[Pattern[idx3, Blank[]]][h2Val[LI[1], LI[0]]] :> deltagpair[h1Val, h2Val],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[Pattern[idx3, Blank[]]][h2Val[LI[1], LI[0]]]*cd[-Pattern[idx3, Blank[]]][h1Val[LI[1], LI[0]]] :> deltagpair[h1Val, h2Val],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[-Pattern[idx3, Blank[]]][h1Val[LI[1], LI[1]]]*cd[Pattern[idx3, Blank[]]][h2Val[LI[1], LI[0]]] :> deltagtpair[h1Val, h2Val],
        h[-Pattern[idx1, Blank[]], -Pattern[idx2, Blank[]]]*cd[Pattern[idx3, Blank[]]][h2Val[LI[1], LI[0]]]*cd[-Pattern[idx3, Blank[]]][h1Val[LI[1], LI[1]]] :> deltagtpair[h1Val, h2Val],
        cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]]*cd[-Pattern[idx2, Blank[]]][h2Val[LI[1], LI[0]]] :> gradagradb[h1Val, h2Val],
        cd[-Pattern[idx2, Blank[]]][h2Val[LI[1], LI[0]]]*cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]] :> gradagradb[h1Val, h2Val],
        cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]]*cd[-Pattern[idx2, Blank[]]][h2Val[LI[1], LI[1]]] :> gradagradtb[h1Val, h2Val],
        cd[-Pattern[idx2, Blank[]]][h2Val[LI[1], LI[1]]]*cd[-Pattern[idx1, Blank[]]][h1Val[LI[1], LI[0]]] :> gradagradtb[h1Val, h2Val]
      }]
    ],
    {h1, pertHeadsList}, {h2, pertHeadsList}
  ];
  
  rules
];

(* ========================================================================== *)
(* Reverse X Rule - Convert φ̇^n to X^(n/2) notation                          *)
(* ========================================================================== *)

(* Rule to convert powers of phidot back to X *)
(* φ[0,1]^2 -> 2X[], φ[0,1]^3 -> 2X[] φ[0,1], φ[0,1]^4 -> 4X[]^2, etc. *)
(* Also handles negative powers: φ[0,1]^(-2) -> 1/(2X[]), φ[0,1]^(-3) -> φ[0,1]/(4X[]^2) *)
GetReverseXRule[] := Module[{phi, LI, Xvar},
  phi = Symbol["xAct`xPand`\[CurlyPhi]"];
  LI = Symbol["xAct`xTensor`LI"];
  Xvar = Symbol["Global`X"][];
  {
    phi[LI[0], LI[1]]^n_Integer /; (n >= 2 || n <= -2) :> (2 * Xvar)^Quotient[n, 2] * phi[LI[0], LI[1]]^Mod[n, 2]
  }
];

(* Apply reverse X rule as the final step of coefficient transformation *)
ApplyReverseXRuleToCoeff[coeff_] := Module[{result, phi, LI},
  phi = Symbol["xAct`xPand`\[CurlyPhi]"];
  LI = Symbol["xAct`xTensor`LI"];
  
  result = coeff;
  (* Expand first to allow cancellations of common φ[0,1] factors *)
  result = Expand[result];
  (* Simplify to cancel common factors in numerator/denominator *)
  result = Simplify[result];
  (* Together to combine fractions and enable more cancellations *)
  result = Together[result];
  (* Apply the rule multiple times to ensure all nested powers are converted *)
  Do[result = result //. GetReverseXRule[], {5}];
  (* Final aggressive simplification *)
  result = Simplify[result, TimeConstraint -> 10];
  (* If single φ[0,1] terms remain, try to factor them out *)
  result = Collect[result, phi[LI[0], LI[1]], Simplify];
  result
];

(* Solves the Scalar EOMs *IN THE ALPHA BASIS* for extreme speed *)
GetAlphaEOMRules[phidotToX_, removeScalarWrapper_] := Module[
  {varSexpr, varSdotexpr, phi, LI, phiddot, phidddot, sol2, sol3, alphaVarS, alphaVarSdot},
  
  If[$AlphaEOMRulesCache =!= Null, Return[$AlphaEOMRulesCache]];
  
  phi = Symbol["xAct`xPand`\[CurlyPhi]"];
  LI = Symbol["xAct`xTensor`LI"];
  phiddot = phi[LI[0], LI[2]];
  phidddot = phi[LI[0], LI[3]];
  
  sol2 = {};
  If[ValueQ[Symbol["Global`varSexpr"]],
    varSexpr = Symbol["Global`varSexpr"];
    alphaVarS = If[Head[varSexpr]===Equal, varSexpr[[2]], varSexpr] /. removeScalarWrapper // Simplify;
    (* Translate equation to Alphas BEFORE solving *)
    alphaVarS = GToAlphas[alphaVarS, phidotToX];
    sol2 = Quiet@Solve[(alphaVarS //. phidotToX) == 0, phiddot];
    If[Length[sol2]==0, sol2 = {}];
  ];
  
  sol3 = {};
  If[ValueQ[Symbol["Global`varSdotexpr"]],
    varSdotexpr = Symbol["Global`varSdotexpr"];
    alphaVarSdot = If[Head[varSdotexpr]===Equal, varSdotexpr[[2]], varSdotexpr] /. removeScalarWrapper // Simplify;
    alphaVarSdot = GToAlphas[alphaVarSdot, phidotToX];
    sol3 = Quiet@Solve[(alphaVarSdot //. phidotToX //. sol2) == 0, phidddot];
    If[Length[sol3]==0, sol3 = {}];
  ];
  
  $AlphaEOMRulesCache = Flatten[{sol2, sol3}] //. phidotToX;
  $AlphaEOMRulesCache
];

(* ========================================================================== *)
(* 1. Base Alpha Parameters & Derivatives                                     *)
(* ========================================================================== *)

ComputeAlphaM[G4_, G4phi_, Hh_, phidot_] := (2 phidot*G4phi)/(2*Hh*G4);
ComputeAlphaK[G4_, G3phi_, G3phiX_, G3X_, G3XX_, KX_, KXX_, Hh_, phidot_, X_] := (2 X/(2*Hh^2*G4))*(KX + 2*X*KXX - 2*G3phi - 2*X*G3phiX + 6*phidot*Hh*(G3X + X*G3XX));
ComputeAlphaB[G4_, G4phi_, G3X_, Hh_, phidot_, X_] := (2 phidot/(2*Hh*G4))*(X*G3X - G4phi);

ComputeAlphaTimeDerivative[alphaExpr_, phi_, phidot_, phiddot_, Hh_, Hdot_, X_, replacementRules_] :=
  Module[{rhs, result},
    rhs = If[Head[alphaExpr] === Equal, alphaExpr[[2]], alphaExpr];
    result = D[rhs, phi]*phidot + D[rhs, phidot]*phiddot + D[rhs, X]*phidot*phiddot + D[rhs, Hh]*Hdot;
    result /. replacementRules // Simplify
  ];

ComputeAlphaXDerivative[alphaExpr_, X_, phidot_, replacementRules_] :=
  Module[{rhs, result},
    rhs = If[Head[alphaExpr] === Equal, alphaExpr[[2]], alphaExpr];
    result = D[rhs, X] + D[rhs, phidot]*phidot/(2 X);
    result /. replacementRules
  ];

ComputeAlphaPhiDerivative[alphaExpr_, phi_, replacementRules_] :=
  Module[{rhs},
    rhs = If[Head[alphaExpr] === Equal, alphaExpr[[2]], alphaExpr];
    D[rhs, phi] /. replacementRules
  ];

(* ========================================================================== *)
(* 2. Compute All Alphas                                                      *)
(* ========================================================================== *)

ComputeAllAlphas[phi_, X_, Hh_, HorndeskiDerivativeReplacements1_, phidotToX_] :=
  Module[{LI, phi0, phidot, phiddot, Hh0, Hdot, X0, G4val, G4phival, G3phival, G3phiXval, G3Xval, G3XXval, KXval, KXXval,
          MstarexprVal, aMexprVal, aKexprVal, aBexprVal, aBdotexprVal, aBddotexprVal, aKdotexprVal, aKddotexprVal, aMdotexprVal, aMddotexprVal, aMdddotexprVal,
          aBXexprVal, aKXexprVal, aMXexprVal, aKXXexprVal, aMXXexprVal, aBXXexprVal, aKXXXexprVal, aMXXXexprVal, aBXXXexprVal,
          aBphiexprVal, aMphiexprVal, aKphiexprVal, aBphiphiexprVal, aMphiphiexprVal, aKphiphiexprVal, aBphiphiphiexprVal, aMphiphiphiexprVal, aKphiphiphiexprVal,
          aBphiXexprVal, aKphiXexprVal, aBphiphiXexprVal, aBphiXXexprVal,
          Mstar, aM, aB, aK, aBdot, aBddot, aKdot, aKddot, aMdot, aMddot, aMdddot,
          aBX, aKX, aMX, aKXX, aMXX, aBXX, aKXXX, aMXXX, aBXXX, aBphi, aMphi, aKphi, aBphiphi, aMphiphi, aKphiphi,
          aBphiphiphi, aMphiphiphi, aKphiphiphi, aBphiX, aKphiX, aBphiphiX, aBphiXX, G4, G4phi, G3phi, G3phiX, G3X, G3XX, KX, KXX},
    
    $AlphaEOMRulesCache = Null; (* Reset Cache *)

    LI = Symbol["xAct`xTensor`LI"]; Mstar = Symbol["Global`Mstar"];
    aM = Symbol["Global`aM"]; aB = Symbol["Global`aB"]; aK = Symbol["Global`aK"];
    aBdot = Symbol["Global`aBdot"]; aBddot = Symbol["Global`aBddot"]; aKdot = Symbol["Global`aKdot"]; aKddot = Symbol["Global`aKddot"];
    aMdot = Symbol["Global`aMdot"]; aMddot = Symbol["Global`aMddot"]; aMdddot = Symbol["Global`aMdddot"];
    
    aBX = Symbol["Global`aBX"]; aKX = Symbol["Global`aKX"]; aMX = Symbol["Global`aMX"];
    aKXX = Symbol["Global`aKXX"]; aMXX = Symbol["Global`aMXX"]; aBXX = Symbol["Global`aBXX"];
    aKXXX = Symbol["Global`aKXXX"]; aMXXX = Symbol["Global`aMXXX"]; aBXXX = Symbol["Global`aBXXX"];
    
    aBphi = Symbol["Global`aBphi"]; aMphi = Symbol["Global`aMphi"]; aKphi = Symbol["Global`aKphi"];
    aBphiphi = Symbol["Global`aBphiphi"]; aMphiphi = Symbol["Global`aMphiphi"]; aKphiphi = Symbol["Global`aKphiphi"];
    aBphiphiphi = Symbol["Global`aBphiphiphi"]; aMphiphiphi = Symbol["Global`aMphiphiphi"]; aKphiphiphi = Symbol["Global`aKphiphiphi"];
    aBphiX = Symbol["Global`aBphiX"]; aKphiX = Symbol["Global`aKphiX"]; aBphiphiX = Symbol["Global`aBphiphiX"]; aBphiXX = Symbol["Global`aBphiXX"];
    G4 = Symbol["Global`G4"]; G4phi = Symbol["Global`G4phi"]; G3phi = Symbol["Global`G3phi"]; G3phiX = Symbol["Global`G3phiX"];
    G3X = Symbol["Global`G3X"]; G3XX = Symbol["Global`G3XX"]; KX = Symbol["Global`KX"]; KXX = Symbol["Global`KXX"];
    
    phi0 = phi[LI[0], LI[0]]; phidot = phi[LI[0], LI[1]]; phiddot = phi[LI[0], LI[2]]; Hh0 = Hh[LI[0], LI[0]]; Hdot = Hh[LI[0], LI[1]]; X0 = X[];
    G4val = G4[phi0]; G4phival = G4phi[phi0]; G3phival = G3phi[phi0, X0]; G3phiXval = G3phiX[phi0, X0]; G3Xval = G3X[phi0, X0]; G3XXval = G3XX[phi0, X0]; KXval = KX[phi0, X0]; KXXval = KXX[phi0, X0];
    
    ToExpression["Global`Mstarexpr", StandardForm, Function[sym, sym = (Mstar[]^2 == 2 G4val), HoldFirst]];
    ToExpression["Global`aMexpr", StandardForm, Function[sym, sym = (aM[] == ComputeAlphaM[G4val, G4phival, Hh0, phidot]), HoldFirst]];
    ToExpression["Global`aKexpr", StandardForm, Function[sym, sym = (aK[] == ComputeAlphaK[G4val, G3phival, G3phiXval, G3Xval, G3XXval, KXval, KXXval, Hh0, phidot, X0]), HoldFirst]];
    ToExpression["Global`aBexpr", StandardForm, Function[sym, sym = (aB[] == ComputeAlphaB[G4val, G4phival, G3Xval, Hh0, phidot, X0]), HoldFirst]];
    
    ToExpression["Global`aBdotexpr", StandardForm, Function[sym, sym = (aBdot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aBexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aBddotexpr", StandardForm, Function[sym, sym = (aBddot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aBdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aKdotexpr", StandardForm, Function[sym, sym = (aKdot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aKexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aKddotexpr", StandardForm, Function[sym, sym = (aKddot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aKdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aMdotexpr", StandardForm, Function[sym, sym = (aMdot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aMexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aMddotexpr", StandardForm, Function[sym, sym = (aMddot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aMdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aMdddotexpr", StandardForm, Function[sym, sym = (aMdddot[] == (ComputeAlphaTimeDerivative[Symbol["Global`aMddotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    
    ToExpression["Global`aBXexpr", StandardForm, Function[sym, sym = (aBX[] == ComputeAlphaXDerivative[Symbol["Global`aBexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKXexpr", StandardForm, Function[sym, sym = (aKX[] == ComputeAlphaXDerivative[Symbol["Global`aKexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMXexpr", StandardForm, Function[sym, sym = (aMX[] == ComputeAlphaXDerivative[Symbol["Global`aMexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKXXexpr", StandardForm, Function[sym, sym = (aKXX[] == ComputeAlphaXDerivative[Symbol["Global`aKXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMXXexpr", StandardForm, Function[sym, sym = (aMXX[] == ComputeAlphaXDerivative[Symbol["Global`aMXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBXXexpr", StandardForm, Function[sym, sym = (aBXX[] == ComputeAlphaXDerivative[Symbol["Global`aBXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKXXXexpr", StandardForm, Function[sym, sym = (aKXXX[] == ComputeAlphaXDerivative[Symbol["Global`aKXXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMXXXexpr", StandardForm, Function[sym, sym = (aMXXX[] == ComputeAlphaXDerivative[Symbol["Global`aMXXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBXXXexpr", StandardForm, Function[sym, sym = (aBXXX[] == ComputeAlphaXDerivative[Symbol["Global`aBXXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    
    (* Time derivatives of alpha X-derivatives *)
    ToExpression["Global`aKXdotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aKXdot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aKXexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aBXdotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aBXdot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aBXexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aMXdotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aMXdot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aMXexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aKXddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aKXddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aKXdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aBXddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aBXddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aBXdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aMXddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aMXddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aMXdotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    
    (* Third time derivatives of alphas and alpha X-derivatives *)
    ToExpression["Global`aKdddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aKdddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aKddotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aBdddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aBdddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aBddotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    ToExpression["Global`aKXdddotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aKXdddot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aKXddotexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    
    (* X-derivatives of K functions needed for higher order rules *)
    ToExpression["Global`aKXXXXexpr", StandardForm, Function[sym, sym = (Symbol["Global`aKXXXX"][] == ComputeAlphaXDerivative[Symbol["Global`aKXXXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBXXXXexpr", StandardForm, Function[sym, sym = (Symbol["Global`aBXXXX"][] == ComputeAlphaXDerivative[Symbol["Global`aBXXXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    
    (* Time derivatives of higher-order X-derivatives *)
    ToExpression["Global`aBXXdotexpr", StandardForm, Function[sym, sym = (Symbol["Global`aBXXdot"][] == (ComputeAlphaTimeDerivative[Symbol["Global`aBXXexpr"][[2]], phi0, phidot, phiddot, Hh0, Hdot, X0, HorndeskiDerivativeReplacements1] //. phidotToX)), HoldFirst]];
    
    ToExpression["Global`aBphiexpr", StandardForm, Function[sym, sym = (aBphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aBexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMphiexpr", StandardForm, Function[sym, sym = (aMphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aMexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKphiexpr", StandardForm, Function[sym, sym = (aKphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aKexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBphiphiexpr", StandardForm, Function[sym, sym = (aBphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aBphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMphiphiexpr", StandardForm, Function[sym, sym = (aMphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aMphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKphiphiexpr", StandardForm, Function[sym, sym = (aKphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aKphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBphiphiphiexpr", StandardForm, Function[sym, sym = (aBphiphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aBphiphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aMphiphiphiexpr", StandardForm, Function[sym, sym = (aMphiphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aMphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKphiphiphiexpr", StandardForm, Function[sym, sym = (aKphiphiphi[] == ComputeAlphaPhiDerivative[Symbol["Global`aKphiphiexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    
    ToExpression["Global`aBphiXexpr", StandardForm, Function[sym, sym = (aBphiX[] == ComputeAlphaPhiDerivative[Symbol["Global`aBXexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aKphiXexpr", StandardForm, Function[sym, sym = (aKphiX[] == ComputeAlphaPhiDerivative[Symbol["Global`aKXexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBphiphiXexpr", StandardForm, Function[sym, sym = (aBphiphiX[] == ComputeAlphaPhiDerivative[Symbol["Global`aBphiXexpr"], phi0, HorndeskiDerivativeReplacements1]), HoldFirst]];
    ToExpression["Global`aBphiXXexpr", StandardForm, Function[sym, sym = (aBphiXX[] == ComputeAlphaXDerivative[Symbol["Global`aBphiXexpr"], X0, phidot, HorndeskiDerivativeReplacements1]), HoldFirst]];
    
    If[$xAlphaVerbose, Print["Alpha functions computed successfully"]];
  ];

(* ========================================================================== *)
(* 3. Transformations (K, G3, G4)                                             *)
(* ========================================================================== *)

PrepEq[eq_, derivRules_, phidotToX_] :=
  Module[{result},
    result = eq /. derivRules //. phidotToX;
    (* If Equal head survived substitution, return the equation.
       If it collapsed (e.g. both sides identical → True), return
       $Failed so callers can detect and skip this equation safely.
       Never fabricate "result == 0" — that would set the wrong thing to zero. *)
    If[Head[result] === Equal, result, $Failed]
  ];

RemoveK[expr_, phidotToX_] := Module[{SolveKXX, SolveKX, SolveK, SolveKphiXX, SolveKphiX, SolveKphi, SolveKXXX, SolveKphiphiX, SolveKXXXX, SolveKphiphiphiX,
         allRules, result, solveResult, 
         aKexpr, aKdotexpr, aKddotexpr, aKXexpr, aKXdotexpr, aKXXexpr,
         varPexpr, varPdotexpr, varPddotexpr, varPdddotexpr, varPphiphiexpr,
         varEexpr, varEdotexpr, varEddotexpr, varEdddotexpr,
         KK, KX, KXX, Kphi, KphiX, KphiXX, KXXX, KXXXX, Kphiphi, Kphiphiphi, KphiphiX, KphiphiphiX, KphiphiXX, KphiXXX, derivRules},
    
    (* Correctly map all global symbols to local variables *)
    aKexpr = Symbol["Global`aKexpr"]; aKdotexpr = Symbol["Global`aKdotexpr"]; aKddotexpr = Symbol["Global`aKddotexpr"];
    aKXexpr = Symbol["Global`aKXexpr"]; aKXdotexpr = Symbol["Global`aKXdotexpr"]; 
    aKXXexpr = Symbol["Global`aKXXexpr"];
    
    varPexpr = Symbol["Global`varPexpr"]; varPdotexpr = Symbol["Global`varPdotexpr"]; 
    varPddotexpr = Symbol["Global`varPddotexpr"]; varPdddotexpr = Symbol["Global`varPdddotexpr"]; varPphiphiexpr = Symbol["Global`varPphiphiexpr"];
    
    varEexpr = Symbol["Global`varEexpr"]; varEdotexpr = Symbol["Global`varEdotexpr"]; 
    varEddotexpr = Symbol["Global`varEddotexpr"]; varEdddotexpr = Symbol["Global`varEdddotexpr"];
    
    KK = Symbol["Global`KK"]; KX = Symbol["Global`KX"]; KXX = Symbol["Global`KXX"]; Kphi = Symbol["Global`Kphi"]; 
    KphiX = Symbol["Global`KphiX"]; KphiXX = Symbol["Global`KphiXX"]; KXXX = Symbol["Global`KXXX"]; KXXXX = Symbol["Global`KXXXX"];
    Kphiphi = Symbol["Global`Kphiphi"]; Kphiphiphi = Symbol["Global`Kphiphiphi"]; KphiphiX = Symbol["Global`KphiphiX"]; KphiphiphiX = Symbol["Global`KphiphiphiX"];
    KphiphiXX = Symbol["Global`KphiphiXX"]; KphiXXX = Symbol["Global`KphiXXX"];
    derivRules = GetHorndeskiRules[];
    
    If[!ValueQ[aKexpr], Return[{expr, {}}]]; allRules = {};

    Quiet[
      (* K base from varPexpr (pressure constraint) *)
      If[ValueQ[varPexpr], solveResult = With[{$pe = PrepEq[varPexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KK]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KX from varEexpr (energy constraint) *)
      If[ValueQ[varEexpr], solveResult = With[{$pe = PrepEq[varEexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KXX from aKexpr (alpha function) *)
      solveResult = With[{$pe = PrepEq[aKexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]];
      
      (* KXXX from aKXexpr (alpha X-derivative) *)
      If[ValueQ[aKXexpr], solveResult = With[{$pe = PrepEq[aKXexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KXXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KXXXX from aKXXexpr (alpha XX-derivative) *)
      If[ValueQ[aKXXexpr], solveResult = With[{$pe = PrepEq[aKXXexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KXXXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* Kphi from varPdotexpr (time-derivative of pressure constraint) *)
      If[ValueQ[varPdotexpr], solveResult = With[{$pe = PrepEq[varPdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, Kphi]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KphiX from varEdotexpr (time-derivative of energy constraint) *)
      If[ValueQ[varEdotexpr], solveResult = With[{$pe = PrepEq[varEdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KphiXX from aKdotexpr (alpha time-derivative) *)
      solveResult = With[{$pe = PrepEq[aKdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]];
      
      (* KphiXXX from aKXdotexpr (time-derivative of alpha X-derivative) *)
      If[ValueQ[aKXdotexpr], solveResult = With[{$pe = PrepEq[aKXdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiXXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* Kphiphi from varPddotexpr (second time-derivative of pressure constraint) *)
      If[ValueQ[varPddotexpr], solveResult = With[{$pe = PrepEq[varPddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, Kphiphi]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];

      (* Kphiphiphi from varPdddotexpr (third time-derivative of pressure constraint) *)
      If[ValueQ[varPdddotexpr], solveResult = With[{$pe = PrepEq[varPdddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, Kphiphiphi]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KphiphiX from varEddotexpr (second time-derivative of energy constraint) *)
      If[ValueQ[varEddotexpr], solveResult = With[{$pe = PrepEq[varEddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiphiX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KphiphiphiX from varEdddotexpr (third time-derivative of energy constraint) *)
      If[ValueQ[varEdddotexpr], solveResult = With[{$pe = PrepEq[varEdddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiphiphiX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
      (* KphiphiXX from aKddotexpr (second time-derivative of alpha function) *)
      If[ValueQ[aKddotexpr], solveResult = With[{$pe = PrepEq[aKddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, KphiphiXX]]]; If[Length[solveResult]>0, AppendTo[allRules, solveResult[[1]]]]];
      
    , {Set::nosym}];
    
    allRules = Flatten[allRules]; allRules = Select[allRules, MatchQ[#, _Rule | _RuleDelayed] &];
    allRules = ToArgAgnosticRules[allRules];
    result = expr;
    result = result /. derivRules;
    result = result //. phidotToX;
    result = result //. allRules;
    result = result //. phidotToX // Simplify;
    With[{rsw = GetRemoveScalarWrapper[]}, If[rsw =!= {}, result = result /. rsw]];
    {result, allRules}
];

RemoveG3[expr_, phidotToX_] := Module[{SolveG3X, SolveG3XX, SolveG3XXX, SolveG3XXXX, SolveG3phiphi, SolveG3phiphiphi, SolveG3phiphiphiphi, SolveG3phiX, SolveG3phiphiX, SolveG3phiXX,
         SolveG3phiphiXX, SolveG3phiXXX, SolveG3phiphiphiX, allRules, result, aBexpr, aBXexpr, aBXXexpr, aBXXXexpr, aBdotexpr, aBddotexpr, aBdddotexpr, aBphiexpr, aBphiphiexpr, aBphiXexpr, aBXdotexpr, aBXXdotexpr, aBXddotexpr,
         varSexpr, varSdotexpr, varSddotexpr,
         G3X, G3XX, G3XXX, G3XXXX, G3phi, G3phiphi, G3phiX, G3phiphiX, G3phiXX, G3phiphiphi,  G3phiphiphiphi, G3phiphiXX, G3phiXXX, G3phiphiphiX, derivRules},
    
    aBexpr = Symbol["Global`aBexpr"]; aBXexpr = Symbol["Global`aBXexpr"]; aBXXexpr = Symbol["Global`aBXXexpr"]; aBXXXexpr = Symbol["Global`aBXXXexpr"];
    aBdotexpr = Symbol["Global`aBdotexpr"]; aBddotexpr = Symbol["Global`aBddotexpr"]; aBdddotexpr = Symbol["Global`aBdddotexpr"]; aBphiexpr = Symbol["Global`aBphiexpr"]; aBphiphiexpr = Symbol["Global`aBphiphiexpr"]; 
    aBphiXexpr = Symbol["Global`aBphiXexpr"]; aBXdotexpr = Symbol["Global`aBXdotexpr"]; aBXXdotexpr = Symbol["Global`aBXXdotexpr"]; aBXddotexpr = Symbol["Global`aBXddotexpr"];
    varSexpr = Symbol["Global`varSexpr"]; varSdotexpr = Symbol["Global`varSdotexpr"]; varSddotexpr = Symbol["Global`varSddotexpr"];
    
    G3X = Symbol["Global`G3X"]; G3XX = Symbol["Global`G3XX"]; G3XXX = Symbol["Global`G3XXX"]; G3XXXX = Symbol["Global`G3XXXX"];
    G3phi = Symbol["Global`G3phi"]; G3phiphi = Symbol["Global`G3phiphi"]; G3phiX = Symbol["Global`G3phiX"]; 
    G3phiphiX = Symbol["Global`G3phiphiX"]; G3phiXX = Symbol["Global`G3phiXX"]; G3phiphiphi = Symbol["Global`G3phiphiphi"]; G3phiphiphiphi = Symbol["Global`G3phiphiphiphi"]; G3phiphiphiX = Symbol["Global`G3phiphiphiX"];
    G3phiphiXX = Symbol["Global`G3phiphiXX"]; G3phiXXX = Symbol["Global`G3phiXXX"]; G3phiphiphiX = Symbol["Global`G3phiphiphiX"];
    derivRules = GetHorndeskiRules[];
    
    Quiet[
      (* Solve for each G3 function independently *)
      SolveG3X = With[{$pe = PrepEq[aBexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3X]]];
      SolveG3XX = With[{$pe = PrepEq[aBXexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3XX]]];
      SolveG3XXX = With[{$pe = PrepEq[aBXXexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3XXX]]];
      SolveG3XXXX = If[ValueQ[aBXXXexpr],
        With[{$pe = PrepEq[aBXXXexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3XXXX]]],
        {}
      ];
      SolveG3phiX = With[{$pe = PrepEq[aBdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiX]]];
      SolveG3phiphiX = With[{$pe = PrepEq[aBddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiX]]];
      SolveG3phiXX = If[ValueQ[aBXdotexpr],
        With[{$pe = PrepEq[aBXdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiXX]]],
        {}
      ];
      
      (* G3phiphi from varSexpr (scalar field equation) *)
      SolveG3phiphi = If[ValueQ[varSexpr],
        With[{$pe = PrepEq[varSexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphi]]],
        With[{$pe = PrepEq[aBphiexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphi]]]
      ];
      
      (* G3phiphiphi from varSdotexpr (time-derivative of scalar equation) *)
      SolveG3phiphiphi = If[ValueQ[varSdotexpr],
        With[{$pe = PrepEq[varSdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiphi]]],
        With[{$pe = PrepEq[aBphiphiexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiphi]]]
      ];

      (* G3phiphiphiphi from varSddotexpr (time-derivative of scalar equation) *)
      SolveG3phiphiphiphi = If[ValueQ[varSddotexpr],
        With[{$pe = PrepEq[varSddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiphiphi]]],
        With[{$pe = PrepEq[aBphiphiphiexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiphiphi]]]
      ];
      
      (* G3phiphiXX from aBXddotexpr (alpha X second time-derivative) *)
      SolveG3phiphiXX = If[ValueQ[aBXddotexpr],
        With[{$pe = PrepEq[aBXddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiXX]]],
        {}
      ];
      
      (* G3phiXXX from aBXXdotexpr (alpha XX time-derivative) *)
      SolveG3phiXXX = If[ValueQ[aBXXdotexpr],
        With[{$pe = PrepEq[aBXXdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiXXX]]],
        {}
      ];
      
      (* G3phiphiphiX from aBdddotexpr (alpha third time-derivative) *)
      SolveG3phiphiphiX = If[ValueQ[aBdddotexpr],
        With[{$pe = PrepEq[aBdddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G3phiphiphiX]]],
        {}
      ];
    , {Set::nosym}];
    
    allRules = Flatten[{SolveG3XXXX, SolveG3XXX, SolveG3XX, SolveG3X, SolveG3phiphiX, SolveG3phiX, SolveG3phiXXX, SolveG3phiXX, SolveG3phiphiXX, SolveG3phiphiphiX, SolveG3phiphi, SolveG3phiphiphi, SolveG3phiphiphiphi}];
    allRules = Select[allRules, MatchQ[#, _Rule | _RuleDelayed] &];
    allRules = ToArgAgnosticRules[allRules];
    result = expr;
    result = result /. derivRules;
    result = result //. phidotToX;
    result = result //. allRules;
    result = result //. phidotToX // Simplify;
    With[{rsw = GetRemoveScalarWrapper[]}, If[rsw =!= {}, result = result /. rsw]];
    {result, allRules}
];

RemoveG4[expr_, phidotToX_] := Module[{SolveG4, SolveG4phi, SolveG4phiphi, SolveG4phiphiphi, SolveG4phiphiphiphi, SolveG4phiphiphiphiphi,
         allRules, result, Mstarexpr, aMexpr, aMdotexpr, aMddotexpr, aMdddotexpr, aMphiphiexpr, aMphiphiphiexpr,
         G4, G4phi, G4phiphi, G4phiphiphi, G4phiphiphiphi, derivRules},
    
    Mstarexpr = Symbol["Global`Mstarexpr"]; aMexpr = Symbol["Global`aMexpr"]; aMdotexpr = Symbol["Global`aMdotexpr"];
    aMddotexpr = Symbol["Global`aMddotexpr"]; aMdddotexpr = Symbol["Global`aMdddotexpr"]; aMphiphiexpr = Symbol["Global`aMphiphiexpr"]; aMphiphiphiexpr = Symbol["Global`aMphiphiphiexpr"];
    G4 = Symbol["Global`G4"]; G4phi = Symbol["Global`G4phi"]; G4phiphi = Symbol["Global`G4phiphi"];
    G4phiphiphi = Symbol["Global`G4phiphiphi"]; G4phiphiphiphi = Symbol["Global`G4phiphiphiphi"]; G4phiphiphiphiphi = Symbol["Global`G4phiphiphiphiphi"];
    derivRules = GetHorndeskiRules[];
    
    Quiet[
      (* Solve for highest order first *)
      SolveG4phiphiphiphi = If[ValueQ[aMdddotexpr],
        With[{$pe = PrepEq[aMdddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G4phiphiphiphi]]],
        {}
      ];
      SolveG4phiphiphi = With[{$pe = PrepEq[aMddotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G4phiphiphi]]];
      SolveG4phiphi = With[{$pe = PrepEq[aMdotexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G4phiphi]]];
      SolveG4phi = With[{$pe = PrepEq[aMexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G4phi]]];
      SolveG4 = With[{$pe = PrepEq[Mstarexpr, derivRules, phidotToX]}, If[$pe === $Failed, {}, SolveForHead[$pe, G4]]];
    , {Set::nosym}];
    
    allRules = Flatten[{SolveG4phiphiphiphi, SolveG4phiphiphi, SolveG4phiphi, SolveG4phi, SolveG4}];
    allRules = Select[allRules, MatchQ[#, _Rule | _RuleDelayed] &];
    allRules = ToArgAgnosticRules[allRules];
    result = expr;
    result = result /. derivRules;
    result = result //. phidotToX;
    result = result //. allRules;
    result = result //. phidotToX // Simplify;
    With[{rsw = GetRemoveScalarWrapper[]}, If[rsw =!= {}, result = result /. rsw]];
    {result, allRules}
];

GToAlphas[expr_, phidotToX_] := 
  Module[{result, kRules, g3Rules, g4Rules},
    {result, kRules} = RemoveK[expr, phidotToX];
    {result, g3Rules} = RemoveG3[result, phidotToX];
    {result, g4Rules} = RemoveG4[result, phidotToX];
    result
  ];

(* ========================================================================== *)
(* 5. Extraction Functions                                                    *)
(* ========================================================================== *)

ExtractAndTransformCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, A, Aalpha, Agamma, gammaRules, coeff, coeffAlpha, basePatt, tderPatt, lapPatt, d, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"];
    cd = Symbol["Global`cd"]; A = Symbol["Global`A"]; Aalpha = Symbol["Global`Aalpha"]; Agamma = Symbol["Global`Agamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]];
    d = Symbol["Global`d"]; horndeskiRules = GetHorndeskiRules[];
    
    Block[{d = Symbol["Global`d"]}, Do[
      basePatt = hh[LI[1], LI[0]]; tderPatt = hh[LI[1], LI[1]]; lapPatt = cd[-d][cd[d][hh[LI[1], LI[0]]]];
      
      coeff = (1/Hh[LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, basePatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules];
      A[1, hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      
      Aalpha[1, hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = (1/Hh[LI[0], LI[0]]) Simplify@Coefficient[collectedEqn, tderPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules];
      A[2, hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      
      Aalpha[2, hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(ah[LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, lapPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules];
      A[3, hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      
      Aalpha[3, hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
    , {hh, pertHeads}];];
    Do[Agamma[i, hh] = ApplyGammaSubstitution[Aalpha[i, hh], gammaRules], {i, 1, 3}, {hh, pertHeads}];
    Print["A coefficients extracted and transformed."];
];

ExtractAndTransformQuadraticCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, A, Aalpha, Agamma, gammaRules, coeff, coeffAlpha, base1, base2, tder1, tder2, lap1, lap2, gpair12, gpair21, d, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"];
    A = Symbol["Global`A"]; Aalpha = Symbol["Global`Aalpha"]; Agamma = Symbol["Global`Agamma"];
    gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]];
    d = Symbol["Global`d"]; horndeskiRules = GetHorndeskiRules[];
    
    Block[{d = Symbol["Global`d"]}, Do[
      base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]]; tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
      lap1 = cd[-d][cd[d][h1[LI[1], LI[0]]]]; lap2 = cd[-d][cd[d][h2[LI[1], LI[0]]]]; gpair12 = cd[-d][base1] * cd[d][base2]; gpair21 = cd[-d][base2] * cd[d][base1];
      
      coeff = (2/Hh[LI[0], LI[0]]^2) Simplify@Mean@{Coefficient[collectedEqn, base1 * base2], Coefficient[collectedEqn, base2 * base1]};
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[1, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      
      (* Transform with cleaning between steps, not cascading all rules *)
      Aalpha[1, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = (1/Hh[LI[0], LI[0]]) Simplify@Coefficient[collectedEqn, tder1 * base2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[2, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Aalpha[2, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(1/3) Simplify@Mean@{Coefficient[collectedEqn, tder1 * tder2], Coefficient[collectedEqn, tder2 * tder1]};
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[3, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Aalpha[3, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(ah[LI[0], LI[0]]^2)/2 Simplify@Coefficient[collectedEqn, base1 * lap2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[4, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Aalpha[4, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      (* With all-pairs reconstruction, keep A5 directional (h1,h2) and let (h2,h1)
        be reconstructed by the separate ordered pair in the outer sum. *)
      coeff = -(2 ah[LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, gpair12];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[5, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Aalpha[5, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = (ah[LI[0], LI[0]]^2)*Hh[LI[0], LI[0]] Simplify@Coefficient[collectedEqn, tder1 * lap2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; A[6, h1, h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Aalpha[6, h1, h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
    , {h1, pertHeads}, {h2, pertHeads}];];
    Do[Agamma[i, h1, h2] = ApplyGammaSubstitution[Aalpha[i, h1, h2], gammaRules], {i, 1, 6}, {h1, pertHeads}, {h2, pertHeads}];
    Print["A quadratic extracted and transformed."];
];

ExtractAndTransformBCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, cd, c, B, Balpha, Bgamma, gammaRules, coeff, coeffAlpha, gradDownPatt, gradtDownPatt, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; cd = Symbol["Global`cd"]; c = Symbol["Global`c"]; B = Symbol["Global`B"]; Balpha = Symbol["Global`Balpha"]; Bgamma = Symbol["Global`Bgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    Block[{c = Symbol["Global`c"]}, Do[
      gradDownPatt = cd[-c][hh[LI[1], LI[0]]]; gradtDownPatt = cd[-c][hh[LI[1], LI[1]]];
      
      coeff = (1/Hh[LI[0], LI[0]]) Simplify@Coefficient[collectedEqn, gradDownPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[1, hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[1, hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -Simplify@Coefficient[collectedEqn, gradtDownPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[2, hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[2, hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
    , {hh, pertHeads}];];
    Do[Bgamma[i, hh] = ApplyGammaSubstitution[Balpha[i, hh], gammaRules], {i, 1, 2}, {hh, pertHeads}];
    Print["B coefficients extracted and transformed."];
];

ExtractAndTransformQuadraticBCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, c, d, B, Balpha, Bgamma, gammaRules, coeff, coeffAlpha, base1, tder1, lap2, gradDown1, gradDown2, gradtDown2, gradUp1, gradgrad2, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"]; c = Symbol["Global`c"]; d = Symbol["Global`d"]; B = Symbol["Global`B"]; Balpha = Symbol["Global`Balpha"]; Bgamma = Symbol["Global`Bgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    Block[{c = Symbol["Global`c"], d = Symbol["Global`d"]}, Do[
      base1=h1[LI[1],LI[0]]; tder1=h1[LI[1],LI[1]]; gradDown1=cd[-c][h1[LI[1],LI[0]]]; gradUp1=cd[d][h1[LI[1],LI[0]]]; gradDown2=cd[-c][h2[LI[1],LI[0]]]; gradtDown2=cd[-c][h2[LI[1],LI[1]]]; lap2=cd[-d][cd[d][h2[LI[1],LI[0]]]]; gradgrad2=cd[-d][cd[-c][h2[LI[1],LI[0]]]];
      
      coeff = (1/Hh[LI[0],LI[0]]) Simplify@Coefficient[collectedEqn, base1 * gradDown2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[1,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[1,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -Simplify@Coefficient[collectedEqn, tder1 * gradDown2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[2,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[2,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -Simplify@Coefficient[collectedEqn, base1 * gradtDown2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[3,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[3,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = (Hh[LI[0],LI[0]]) Simplify@Coefficient[collectedEqn, tder1 * gradtDown2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[4,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[4,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(Hh[LI[0],LI[0]] ah[LI[0],LI[0]]^2)/2 Simplify@(Coefficient[collectedEqn, gradDown1*lap2] - Coefficient[collectedEqn, gradUp1*gradgrad2]);
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; B[5,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      Balpha[5,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
    , {h1, pertHeads}, {h2, pertHeads}];];
    Do[Bgamma[i, h1, h2] = ApplyGammaSubstitution[Balpha[i, h1, h2], gammaRules], {i, 1, 5}, {h1, pertHeads}, {h2, pertHeads}];
    Print["B quadratic extracted and transformed."];
];

ExtractAndTransformCCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, h, b, c, d, CC, CCalpha, Cgamma, gammaRules, coeff, coeffAlpha, deltaTermPatt, deltaTPatt, deltaTTPatt, gradgradPatt, lapPatt, tracelessPatt, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"]; h = Symbol["Global`h"]; b = Symbol["Global`b"]; c = Symbol["Global`c"]; d = Symbol["Global`d"]; CC = Symbol["Global`CC"]; CCalpha = Symbol["Global`CCalpha"]; Cgamma = Symbol["Global`Cgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    
    Block[{h = Symbol["Global`h"], b = Symbol["Global`b"], c = Symbol["Global`c"], d = Symbol["Global`d"]}, Do[
      deltaTermPatt = h[-b,-c]*hh[LI[1],LI[0]]; deltaTPatt = h[-b,-c]*hh[LI[1],LI[1]]; deltaTTPatt = h[-b,-c]*hh[LI[1],LI[2]]; 
      gradgradPatt = cd[-c][cd[-b][hh[LI[1],LI[0]]]];
      lapPatt = cd[-d][cd[d][hh[LI[1],LI[0]]]];
      
      coeff = 1/(Hh[LI[0],LI[0]]^2 ah[LI[0],LI[0]]^2) Simplify@Coefficient[collectedEqn, deltaTermPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[1,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      CCalpha[1,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = 1/(Hh[LI[0],LI[0]] ah[LI[0],LI[0]]^2) Simplify@Coefficient[collectedEqn, deltaTPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[2,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      CCalpha[2,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -1/(ah[LI[0],LI[0]]^2) Simplify@Coefficient[collectedEqn, deltaTTPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[3,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      CCalpha[3,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      (* Extract traceless part: h[-b,-c]*cd[-d][cd[d][hh]] - cd[-c][cd[-b][hh]] *)
      coeff = -Simplify@(SafeCoefficient[collectedEqn, deltaTermPatt*lapPatt] - SafeCoefficient[collectedEqn, gradgradPatt]);
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[4,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      CCalpha[4,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
    , {hh, pertHeads}];];
    Do[Cgamma[i, hh] = ApplyGammaSubstitution[CCalpha[i, hh], gammaRules], {i, 1, 4}, {hh, pertHeads}];
    Print["C coefficients extracted and transformed."];
];

ExtractAndTransformQuadraticCCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, h, b, c, d, CC, CCalpha, Cgamma, gammaRules, coeff, coeffAlpha, base1, base2, tder1, tder2, ttder2, deltaTerm1, deltaTerm2, deltaT1, deltaT2, gradgrad2, deltagpair12, deltagpair21, deltagtpair12, gradagradb12, gradagradb21, gradagradtb12, normalizedEqn, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"]; h = Symbol["Global`h"]; b = Symbol["Global`b"]; c = Symbol["Global`c"]; d = Symbol["Global`d"]; CC = Symbol["Global`CC"]; CCalpha = Symbol["Global`CCalpha"]; Cgamma = Symbol["Global`Cgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    
    Block[{h = Symbol["Global`h"], b = Symbol["Global`b"], c = Symbol["Global`c"], d = Symbol["Global`d"]},
      normalizedEqn = CanonicalizeGradientProducts[collectedEqn /. {Symbol["Global`b"] -> b, Symbol["Global`c"] -> c}];
      Do[
        base1=h1[LI[1],LI[0]]; base2=h2[LI[1],LI[0]]; tder1=h1[LI[1],LI[1]]; tder2=h2[LI[1],LI[1]]; ttder2=h2[LI[1],LI[2]];
        deltaTerm1=h[-b,-c]*base1; deltaTerm2=h[-b,-c]*base2; deltaT1=h[-b,-c]*tder1; deltaT2=h[-b,-c]*tder2;
        gradgrad2=cd[-c][cd[-b][base2]]; deltagpair12=h[-b,-c]*cd[-d][base1]*cd[d][base2]; deltagpair21=h[-b,-c]*cd[-d][base2]*cd[d][base1];
        deltagtpair12=h[-b,-c]*cd[-d][tder1]*cd[d][base2]; gradagradb12=cd[-b][base1]*cd[-c][base2]; gradagradb21=cd[-b][base2]*cd[-c][base1];
        gradagradtb12=cd[-b][base1]*cd[-c][tder2];
        
        coeff = -1/(Hh[LI[0],LI[0]]^2 ah[LI[0],LI[0]]^2) Simplify@Mean@{SafeCoefficient[normalizedEqn, deltaTerm1*base2], SafeCoefficient[normalizedEqn, deltaTerm2*base1]};
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[1,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[1,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = 1/(Hh[LI[0],LI[0]] ah[LI[0],LI[0]]^2) Simplify@SafeCoefficient[normalizedEqn, deltaTerm1*tder2];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[2,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[2,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = -1/(ah[LI[0],LI[0]]^2) Simplify@Mean@{SafeCoefficient[normalizedEqn, deltaT1*tder2], SafeCoefficient[normalizedEqn, deltaT2*tder1]};
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[3,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[3,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = 1/(ah[LI[0],LI[0]]^2) Simplify@SafeCoefficient[normalizedEqn, deltaTerm1*ttder2];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[4,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[4,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = Hh[LI[0],LI[0]]/(ah[LI[0],LI[0]]^2) Simplify@SafeCoefficient[normalizedEqn, deltaT1*ttder2];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[5,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[5,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = Simplify@Mean@{SafeCoefficient[normalizedEqn, deltagpair12], SafeCoefficient[normalizedEqn, deltagpair21]};
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[6,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[6,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = (Hh[LI[0],LI[0]])/2 Simplify@SafeCoefficient[normalizedEqn, deltagtpair12];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[7,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[7,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = Simplify@SafeCoefficient[normalizedEqn, gradagradb12];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[8,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[8,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = -Hh[LI[0],LI[0]] Simplify@SafeCoefficient[normalizedEqn, gradagradtb12];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[9,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[9,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
        
        coeff = -Simplify@SafeCoefficient[normalizedEqn, base1*gradgrad2];
        If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; CC[10,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
        CCalpha[10,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      , {h1, pertHeads}, {h2, pertHeads}];];
      Do[Cgamma[i, h1, h2] = ApplyGammaSubstitution[CCalpha[i, h1, h2], gammaRules], {i, 1, 10}, {h1, pertHeads}, {h2, pertHeads}];
      Print["C quadratic extracted and transformed."];
];

ExtractAndTransformDCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, b, phi, delta, DD, DDalpha, Dgamma, gammaRules, coeff, coeffAlpha, basePatt, tderPatt, ttderPatt, lapPatt, horndeskiRules},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"]; b = Symbol["Global`b"]; phi = Symbol["xAct`xPand`\[CurlyPhi]"]; delta = Symbol["xAct`xTensor`delta"]; DD = Symbol["Global`DD"]; DDalpha = Symbol["Global`DDalpha"]; Dgamma = Symbol["Global`Dgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    
    Block[{b = Symbol["Global`b"]}, Do[
      basePatt = hh[LI[1],LI[0]]; tderPatt = hh[LI[1],LI[1]]; ttderPatt = hh[LI[1],LI[2]]; lapPatt = cd[-b][cd[b][hh[LI[1],LI[0]]]];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*1/(Hh[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, basePatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[1,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[1,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*1/(Hh[LI[0],LI[0]])*1/(delta[-LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, tderPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[2,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[2,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*1/(delta[-LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, ttderPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[3,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[3,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(ah[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) Simplify@Coefficient[collectedEqn, lapPatt];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[4,hh] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[4,hh] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
    , {hh, pertHeads}];];
    Do[Dgamma[i, hh] = ApplyGammaSubstitution[DDalpha[i, hh], gammaRules], {i, 1, 4}, {hh, pertHeads}];
    Print["D coefficients extracted and transformed."];
];

ExtractAndTransformQuadraticDCoefficients[collectedEqn_, pertHeads_, phidotToX_, removeScalarWrapper_] :=
  Module[{LI, Hh, ah, cd, b, c, d, phi, delta, DD, DDalpha, Dgamma, gammaRules, coeff, coeffAlpha, base1, base2, tder1, tder2, ttder1, ttder2, ttderPair12, lap1, lap2, gradPair12, gradPair21, gradtPair12, gradtPair21, gradtPairSwap12, gradtderPair12, gradgrad1, gradgrad2, gradgradUp1a, gradgradUp1b, gradgradUp2a, gradgradUp2b, horndeskiRules, canonicalEqn, d9TaggedEqn, d9Tag, d11TaggedEqn, d11Tag},
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"]; cd = Symbol["Global`cd"]; b = Symbol["Global`b"]; c = Symbol["Global`c"]; d = Symbol["Global`d"]; phi = Symbol["xAct`xPand`\[CurlyPhi]"]; delta = Symbol["xAct`xTensor`delta"]; DD = Symbol["Global`DD"]; DDalpha = Symbol["Global`DDalpha"]; Dgamma = Symbol["Global`Dgamma"]; gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]]; horndeskiRules = GetHorndeskiRules[];
    
    (* Canonicalize gradient products to ensure unique representation *)
    canonicalEqn = CanonicalizeGradientProducts[collectedEqn];
    
    Block[{b = Symbol["Global`b"], c = Symbol["Global`c"], d = Symbol["Global`d"]}, Do[
      base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]];
      tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
      ttder1 = h1[LI[1], LI[2]]; ttder2 = h2[LI[1], LI[2]];
      ttderPair12 = tder1 * ttder2;
      lap1 = cd[-b][cd[b][base1]]; lap2 = cd[-b][cd[b][base2]];
      gradPair12 = cd[-b][base1] * cd[b][base2]; gradPair21 = cd[-b][base2] * cd[b][base1];
      gradtPair12 = cd[-b][tder1] * cd[b][base2];
      gradtPair21 = cd[-b][tder2] * cd[b][base1];
      gradtPairSwap12 = cd[-b][base2] * cd[b][tder1];
      gradgrad1 = cd[-c][cd[-b][base1]]; gradgrad2 = cd[-c][cd[-b][base2]];
      gradgradUp1a = cd[b][cd[c][base1]]; gradgradUp1b = cd[c][cd[b][base1]];
      gradgradUp2a = cd[b][cd[c][base2]]; gradgradUp2b = cd[c][cd[b][base2]];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*1/(Hh[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) Simplify@Mean@{SafeCoefficient[collectedEqn, base1*base2], SafeCoefficient[collectedEqn, base2*base1]};
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[1,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[1,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*1/(3*Hh[LI[0],LI[0]])*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, base1*tder2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[2,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[2,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(1/3)*1/(delta[-LI[0], LI[0]]^2) Simplify@Mean@{SafeCoefficient[collectedEqn, tder1*tder2], SafeCoefficient[collectedEqn, tder2*tder1]};
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[3,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[3,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(1/3)*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, base1*ttder2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[4,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[4,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(-1/3)*(Hh[LI[0], LI[0]])*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, ttderPair12];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[5,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[5,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
      
      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(ah[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, base1*lap2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[6,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[6,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

            (* D7 kept directional to match collected scalar equation ordering:
         cd[-b][base1] cd[b][base2] for ordered pair {h1,h2}. *)
            coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(ah[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) *
              SafeCoefficient[collectedEqn, gradPair12] // Simplify;
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[7,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[7,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(-ah[LI[0],LI[0]]^2*Hh[LI[0],LI[0]])*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, tder1*lap2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[8,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[8,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

            (* D9: extract from canonicalized expression using a single index-agnostic pattern *)
            (* This avoids double counting equivalent orderings in an Orderless product. *)
            d9Tag = Unique["d9Tag$"];
            d9TaggedEqn = canonicalEqn /. cd[-Pattern[idx, Blank[]]][tder1] * cd[Pattern[idx, Blank[]]][base2] :> d9Tag;
            coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(-1/2)*(ah[LI[0],LI[0]]^2*Hh[LI[0],LI[0]])*1/(delta[-LI[0], LI[0]]^2) *
              SafeCoefficient[d9TaggedEqn, d9Tag] // Simplify;
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[9,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[9,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

      coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(-ah[LI[0],LI[0]]^2*Hh[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) Simplify@SafeCoefficient[collectedEqn, ttder1*lap2];
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[10,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[10,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

            (* D11: extract from canonicalized expression using one orientation to avoid QQ overcounting *)
            d11Tag = Unique["d11Tag$"];
            d11TaggedEqn = canonicalEqn /. cd[-Pattern[idx, Blank[]]][tder1] * cd[Pattern[idx, Blank[]]][tder2] :> d11Tag;
            coeff = -(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]]))*(ah[LI[0],LI[0]]^2*Hh[LI[0],LI[0]]^2)*1/(delta[-LI[0], LI[0]]^2) *
              SafeCoefficient[d11TaggedEqn, d11Tag] // Simplify;
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[11,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[11,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];

      coeff = 2(phi[LI[0],LI[1]]/(Hh[LI[0],LI[0]])) (Hh[LI[0],LI[0]]^2 ah[LI[0],LI[0]]^4)*1/(delta[-LI[0], LI[0]]^2) *
              If[h1 === h2,
                SafeCoefficient[collectedEqn, cd[-b][cd[b][base1]] * cd[-c][cd[c][base2]]],
                SafeCoefficient[collectedEqn, cd[-b][cd[b][base1]] * cd[-c][cd[c][base2]]] +
                SafeCoefficient[collectedEqn, cd[-b][cd[b][base2]] * cd[-c][cd[c][base1]]]
              ] // Simplify;
      If[ValueQ[horndeskiRules], coeff = coeff /. horndeskiRules]; DD[12,h1,h2] = coeff /. GetBackgroundEOMRules[] //. GetReverseXRule[] //. removeScalarWrapper // Simplify;
      DDalpha[12,h1,h2] = TransformCoeffToAlphaLanguage[coeff, phidotToX, removeScalarWrapper];
    , {h1, pertHeads}, {h2, pertHeads}];];
    Do[Dgamma[i, h1, h2] = ApplyGammaSubstitution[DDalpha[i, h1, h2], gammaRules], {i, 1, 12}, {h1, pertHeads}, {h2, pertHeads}];
    Print["D quadratic extracted and transformed."];
];

(* ========================================================================== *)
(* Gamma Conversion                                                           *)
(* ========================================================================== *)

ApplyGammaSubstitution[expr_, gammaRules_List, eomRules_List:{}] :=
  Module[{aK, aB, aM, gamA, allRules, out},
    aK = Symbol["Global`aK"]; aB = Symbol["Global`aB"]; aM = Symbol["Global`aM"]; gamA = Symbol["Global`gamA"];
    allRules = Join[eomRules, gammaRules];

    (* Apply all rules to a fixed point to catch nested dependencies *)
    out = FixedPoint[Simplify[# /. allRules] &, expr, 8];

    (* Enforce the primary gammaA definition last for robust cancellation *)
    out = Simplify[out /. (aK[] -> gamA - 6 aB[] - 2 aM[])];

    (* Remove leftover linear identity combinations, if any *)
    out = Simplify[Expand[out] /. {
      c_. * (gamA - 6 aB[] - aK[] - 2 aM[]) :> 0,
      c_. * (-gamA + 6 aB[] + aK[] + 2 aM[]) :> 0
    }];

    out
  ];

GenerateAlphaToGammaRules[] := 
  Module[{LI, Hh, X, phi, 
          aK, aKdot, aKddot, aM, aMdot, aMddot, aB, aBdot, aBddot,
          aBX, aBXdot, aBXX, aBXXdot, aKX, aKXdot,
          tvarE, tvarP,
          gamK, dotgamK, gamM, dotgamM, gamB, dotgamB, gamE, gamX, dotgamX, gamA, gamD, dotgamD, gamF, dotgamF, gamC, Dc},
    
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; X = Symbol["Global`X"]; phi = Symbol["xAct`xPand`\[CurlyPhi]"];
    aK = Symbol["Global`aK"]; aKdot = Symbol["Global`aKdot"]; aKddot = Symbol["Global`aKddot"];
    aM = Symbol["Global`aM"]; aMdot = Symbol["Global`aMdot"]; aMddot = Symbol["Global`aMddot"];
    aB = Symbol["Global`aB"]; aBdot = Symbol["Global`aBdot"]; aBddot = Symbol["Global`aBddot"];
    aBX = Symbol["Global`aBX"]; aBXdot = Symbol["Global`aBXdot"]; aBXX = Symbol["Global`aBXX"]; aBXXdot = Symbol["Global`aBXXdot"];
    aKX = Symbol["Global`aKX"]; aKXdot = Symbol["Global`aKXdot"];
    tvarE = Symbol["Global`tvarE"]; tvarP = Symbol["Global`tvarP"];
    gamK = Symbol["Global`gamK"]; dotgamK = Symbol["Global`dotgamK"];
    gamM = Symbol["Global`gamM"]; dotgamM = Symbol["Global`dotgamM"];
    gamB = Symbol["Global`gamB"]; dotgamB = Symbol["Global`dotgamB"];
    gamE = Symbol["Global`gamE"]; gamX = Symbol["Global`gamX"]; dotgamX = Symbol["Global`dotgamX"];
    gamA = Symbol["Global`gamA"]; gamD = Symbol["Global`gamD"]; dotgamD = Symbol["Global`dotgamD"];
    gamF = Symbol["Global`gamF"]; dotgamF = Symbol["Global`dotgamF"]; gamC = Symbol["Global`gamC"]; Dc = Symbol["Global`Dc"];

    Flatten[{
      Solve[gamK == aKdot[]/Hh[LI[0], LI[0]] + aK[]*(3 + aM[]), aKdot[]][[1]],
      Solve[dotgamK == aKddot[]/Hh[LI[0], LI[0]] - aKdot[] Hh[LI[0], LI[1]]/Hh[LI[0], LI[0]]^2 + aKdot[]*(3 + aM[]) + aK[]*aMdot[], aKddot[]][[1]],
      Solve[gamM == aMdot[]/Hh[LI[0], LI[0]] + aM[]*(3 + aM[]), aMdot[]][[1]],
      Solve[dotgamM == aMddot[]/Hh[LI[0], LI[0]] + Hh[LI[0], LI[0]] aB[] (gamM - aM[] (3 + aM[])) (3 + 2 aM[] - (Hh[LI[0], LI[1]]/Hh[LI[0], LI[0]]^2)), aMddot[]][[1]],
      Solve[gamB == aBdot[]/Hh[LI[0], LI[0]] + aB[]*(3 + aM[]) + (Hh[LI[0], LI[1]]/Hh[LI[0], LI[0]]^2) aB[], aBdot[]][[1]],
      Solve[dotgamB == aBddot[]/Hh[LI[0], LI[0]] + Hh[LI[0], LI[0]] aB[] ((gamB/aB[] - 3 - 2 aM[] - (Hh[LI[0], LI[1]]/Hh[LI[0], LI[0]]^2)) (3 + aM[]) + gamM + (Hh[LI[0], LI[2]]/Hh[LI[0], LI[0]]^3) - (2 Hh[LI[0], LI[1]]^2/Hh[LI[0], LI[0]]^4)), aBddot[]][[1]],
      (tvarE + tvarP) -> (2 Hh[LI[0], LI[1]] - gamE Hh[LI[0], LI[0]]^2),
      (tvarE[] + tvarP[]) -> (2 Hh[LI[0], LI[1]] - gamE Hh[LI[0], LI[0]]^2),
      Solve[gamE == (2 Hh[LI[0], LI[1]] - tvarE - tvarP)/Hh[LI[0], LI[0]]^2, tvarE][[1]],
      Solve[gamE == (2 Hh[LI[0], LI[1]] - tvarE[] - tvarP[])/Hh[LI[0], LI[0]]^2, tvarE[]][[1]],
      Solve[4 gamX == 81*aB[] + 4*(2*aK[] + 6*aM[] - 27*aBX[]*X[] - 2*aKX[]*X[] + 3*aBXX[]*X[]^2), aBXX[]][[1]],
      Solve[4 dotgamX == 81*aBdot[] + 4*(2*aKdot[] + 6*aMdot[] - 27*aBXdot[]*X[] - 27*aBX[]*phi[LI[0], LI[2]]*phi[LI[0], LI[1]] - 2*aKXdot[]*X[] - 2*aKX[]*phi[LI[0], LI[2]]*phi[LI[0], LI[1]] + 3*aBXXdot[]*X[]^2 + 6*aBXX[]*X[]*phi[LI[0], LI[2]]*phi[LI[0], LI[1]]), aBXXdot[]][[1]],
      Solve[gamD == aB[] - 2 X[] aBX[], aBX[]][[1]],
      Solve[dotgamD == Hh[LI[0], LI[0]] (gamB - 3 aB[] - aB[] aM[]) - (Hh[LI[0], LI[1]]/Hh[LI[0], LI[0]]) aB[] - (2 phi[LI[0], LI[2]]/phi[LI[0], LI[1]]) (-gamD + aB[]) - 2 X[] aBXdot[], aBXdot[]][[1]],
      Solve[gamF == aK[] - 2 X[] aKX[], aKX[]][[1]],
      Solve[dotgamF == aKdot[] - 2 X[] aKXdot[] - 2 phi[LI[0], LI[2]] phi[LI[0], LI[1]] aKX[], aKXdot[]][[1]],
      Solve[gamC == 3 aB[] + aK[] + 3 aM[] - X[] aKX[], aKX[]][[1]],
      Solve[gamA == 6*aB[] + aK[] + 2*aM[], aK[]][[1]],
      Solve[2*Dc == 3*aB[]^2 + 2*aK[], aK[]][[1]]
    }]
  ];

ConvertAllCoefficientsToGammaLanguage[pertHeadsList_] :=
  Module[{LI, hh, h1, h2, Aalpha, Agamma, Balpha, Bgamma, Calpha, Cgamma, Dalpha, Dgamma, gammaRules, eomRules, allRules},
    
    Aalpha = Symbol["Global`Aalpha"]; Balpha = Symbol["Global`Balpha"]; Calpha = Symbol["Global`Calpha"]; Dalpha = Symbol["Global`Dalpha"];
    Agamma = Symbol["Global`Agamma"]; Bgamma = Symbol["Global`Bgamma"]; Cgamma = Symbol["Global`Cgamma"]; Dgamma = Symbol["Global`Dgamma"];
    gammaRules = SafeRuleList[GenerateAlphaToGammaRules[]];
    eomRules = SafeRuleList[GetBackgroundEOMRules[]];
    allRules = Join[eomRules, gammaRules];
    
    (* Apply both EOM rules and gamma rules, then simplify more aggressively *)
    Do[Agamma[i, hh] = ApplyGammaSubstitution[Aalpha[i, hh], gammaRules, eomRules]; Symbol["Global`Agammaalphas"][i, hh] = Agamma[i, hh];, {i, 1, 3}, {hh, pertHeadsList}];
    Do[Agamma[i, h1, h2] = ApplyGammaSubstitution[Aalpha[i, h1, h2], gammaRules, eomRules]; Symbol["Global`Agammaalphas"][i, h1, h2] = Agamma[i, h1, h2];, {i, 1, 6}, {h1, pertHeadsList}, {h2, pertHeadsList}];
    Do[Bgamma[i, hh] = ApplyGammaSubstitution[Balpha[i, hh], gammaRules, eomRules]; Symbol["Global`Bgammaalphas"][i, hh] = Bgamma[i, hh];, {i, 1, 2}, {hh, pertHeadsList}];
    Do[Bgamma[i, h1, h2] = ApplyGammaSubstitution[Balpha[i, h1, h2], gammaRules, eomRules]; Symbol["Global`Bgammaalphas"][i, h1, h2] = Bgamma[i, h1, h2];, {i, 1, 5}, {h1, pertHeadsList}, {h2, pertHeadsList}];
    Do[Cgamma[i, hh] = ApplyGammaSubstitution[CCalpha[i, hh], gammaRules, eomRules]; Symbol["Global`Cgammaalphas"][i, hh] = Cgamma[i, hh];, {i, 1, 4}, {hh, pertHeadsList}];
    Do[Cgamma[i, h1, h2] = ApplyGammaSubstitution[CCalpha[i, h1, h2], gammaRules, eomRules]; Symbol["Global`Cgammaalphas"][i, h1, h2] = Cgamma[i, h1, h2];, {i, 1, 10}, {h1, pertHeadsList}, {h2, pertHeadsList}];
    Do[Dgamma[i, hh] = ApplyGammaSubstitution[DDalpha[i, hh], gammaRules, eomRules]; Symbol["Global`Dgammaalphas"][i, hh] = Dgamma[i, hh];, {i, 1, 4}, {hh, pertHeadsList}];
    Do[Dgamma[i, h1, h2] = ApplyGammaSubstitution[DDalpha[i, h1, h2], gammaRules, eomRules]; Symbol["Global`Dgammaalphas"][i, h1, h2] = Dgamma[i, h1, h2];, {i, 1, 12}, {h1, pertHeadsList}, {h2, pertHeadsList}];
    
    If[$xAlphaVerbose, Print["Success! Coefficients converted to Gamma language."]];
  ];

(* ========================================================================== *)
(* Validation: Completeness Check                                            *)
(* ========================================================================== *)
(* RECONSTRUCTION AND VALIDATION FUNCTIONS *)
(* ========================================================================== *)

(* Reconstruct complete linear equation by summing over all perturbations *)
ReconstructFullLinearEquation[equationType_String, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{reconstructed, hh},
    reconstructed = 0;
    
    Do[
      reconstructed += ReconstructEquationFromCoefficients[equationType, hh, removeScalarWrapper];
    , {hh, pertHeads}];
    
    Simplify[reconstructed]
  ];

(* Reconstruct complete quadratic equation by summing over all perturbation pairs *)
ReconstructFullQuadraticEquation[equationType_String, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{reconstructed, h1, h2},
    reconstructed = 0;

    Do[
      reconstructed += ReconstructEquationFromCoefficients[equationType, {h1, h2}, removeScalarWrapper];
    , {h1, pertHeads}, {h2, pertHeads}];
    
    Simplify[reconstructed]
  ];

(* Validate full linear equation reconstruction *)
ValidateFullLinearEquation[equationType_String, originalEqn_, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{reconstructed, residual, isComplete, phidotToX, originalCanon, reconstructedCanon},
    
    reconstructed = ReconstructFullLinearEquation[equationType, pertHeads, removeScalarWrapper];
    originalCanon = CanonicalizeGradientProducts[originalEqn];
    reconstructedCanon = CanonicalizeGradientProducts[reconstructed];
    
    phidotToX = GetReverseXRule[];
    residual1 = Simplify[(originalCanon - reconstructedCanon) //. phidotToX];
    residual = Simplify[CanonicalizeGradientProducts[(residual1) //. phidotToX]];
    isComplete = (residual === 0 || Simplify[residual] === 0);
    
    {isComplete, residual}
  ];

(* Validate full quadratic equation reconstruction *)
ValidateFullQuadraticEquation[equationType_String, originalEqn_, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{reconstructed, residual, isComplete, phidotToX, originalCanon, reconstructedCanon},
    
    reconstructed = ReconstructFullQuadraticEquation[equationType, pertHeads, removeScalarWrapper];
    originalCanon = CanonicalizeGradientProducts[originalEqn];
    reconstructedCanon = CanonicalizeGradientProducts[reconstructed];
    
    phidotToX = GetReverseXRule[];
    residual1 = Simplify[(originalCanon - reconstructedCanon) //. phidotToX];
    residual = Simplify[CanonicalizeGradientProducts[(residual1) //. phidotToX]];
    isComplete = (residual === 0 || Simplify[residual] === 0);
    
    {isComplete, residual}
  ];

(* Validate all full equations at once *)
ValidateAllFullEquations[equations_Association, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{results = <||>, isComplete, residual},
    
    (* Validate 00 equation - linear *)
    If[KeyExistsQ[equations, "eq00Linear"],
      {isComplete, residual} = ValidateFullLinearEquation["linear-A", equations["eq00Linear"], pertHeads, removeScalarWrapper];
      results["00-linear"] = {isComplete, residual};
      Print["A linear equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate 00 equation - quadratic *)
    If[KeyExistsQ[equations, "eq00Quadratic"],
      {isComplete, residual} = ValidateFullQuadraticEquation["quadratic-A", equations["eq00Quadratic"], pertHeads, removeScalarWrapper];
      results["00-quadratic"] = {isComplete, residual};
      Print["A quadratic equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate 0i equation - linear *)
    If[KeyExistsQ[equations, "eq0iLinear"],
      {isComplete, residual} = ValidateFullLinearEquation["linear-B", equations["eq0iLinear"], pertHeads, removeScalarWrapper];
      results["0i-linear"] = {isComplete, residual};
      Print["B linear equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate 0i equation - quadratic *)
    If[KeyExistsQ[equations, "eq0iQuadratic"],
      {isComplete, residual} = ValidateFullQuadraticEquation["quadratic-B", equations["eq0iQuadratic"], pertHeads, removeScalarWrapper];
      results["0i-quadratic"] = {isComplete, residual};
      Print["B quadratic equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate ij equation - linear *)
    If[KeyExistsQ[equations, "eqijLinear"],
      {isComplete, residual} = ValidateFullLinearEquation["linear-C", equations["eqijLinear"], pertHeads, removeScalarWrapper];
      results["ij-linear"] = {isComplete, residual};
      Print["C linear equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate ij equation - quadratic *)
    If[KeyExistsQ[equations, "eqijQuadratic"],
      {isComplete, residual} = ValidateFullQuadraticEquation["quadratic-C", equations["eqijQuadratic"], pertHeads, removeScalarWrapper];
      results["ij-quadratic"] = {isComplete, residual};
      Print["C quadratic equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate scalar equation - linear *)
    If[KeyExistsQ[equations, "eqScalarLinear"],
      {isComplete, residual} = ValidateFullLinearEquation["linear-D", equations["eqScalarLinear"], pertHeads, removeScalarWrapper];
      results["scalar-linear"] = {isComplete, residual};
      Print["D linear equation reconstructed and validated: ", isComplete];
    ];
    
    (* Validate scalar equation - quadratic *)
    If[KeyExistsQ[equations, "eqScalarQuadratic"],
      {isComplete, residual} = ValidateFullQuadraticEquation["quadratic-D", equations["eqScalarQuadratic"], pertHeads, removeScalarWrapper];
      results["scalar-quadratic"] = {isComplete, residual};
      Print["D quadratic equation reconstructed and validated: ", isComplete];
    ];

    (* Count successes and failures *)
    Module[{passed = 0, failed = 0},
      Do[
        If[results[key][[1]], passed++, failed++];
      , {key, Keys[results]}];

      Print["Validation summary -> Passed: ", passed, ", Failed: ", failed];
    ];
  ];

(* ========================================================================== *)

(* Helper function to extract only terms containing a specific perturbation *)
ExtractTermsForPerturbation[equation_, pertHead_] :=
  Module[{LI, terms, filteredTerms},
    LI = Symbol["xAct`xTensor`LI"];
    
    (* Split into sum of terms *)
    terms = If[Head[equation] === Plus, List @@ equation, {equation}];
    
    (* Keep only terms that contain the specified perturbation *)
    filteredTerms = Select[terms, !FreeQ[#, pertHead[LI[_], LI[_]]] &];
    
    (* Return sum of filtered terms *)
    If[Length[filteredTerms] == 0, 0, Total[filteredTerms]]
  ];

ReconstructEquationFromCoefficients[equationType_String, pertHead_, removeScalarWrapper_:{}] :=
  Module[{LI, Hh, ah, cd, phi, h, delta, d, b, c, reconstructed, basePatt, tderPatt, lapPatt, ttderPatt,
          A, B, CC, DD},
    
    LI = Symbol["xAct`xTensor`LI"]; Hh = Symbol["Global`Hh"]; ah = Symbol["Global`ah"];
    cd = Symbol["Global`cd"]; phi = Symbol["xAct`xPand`\[CurlyPhi]"]; h = Symbol["Global`h"]; delta = Symbol["xAct`xTensor`delta"]; 
    d = Symbol["Global`d"]; b = Symbol["Global`b"]; c = Symbol["Global`c"];
    A = Symbol["Global`A"]; B = Symbol["Global`B"]; CC = Symbol["Global`CC"]; DD = Symbol["Global`DD"];
    
    reconstructed = 0;
    
    Which[
      (* Linear A coefficients (00 equation) *)
      (* A[1,h] = (1/H^2) * coeff, A[2,h] = (1/H) * coeff, A[3,h] = (-a^2) * coeff *)
      equationType == "linear-A",
      Block[{d = Symbol["Global`d"]},
        basePatt = pertHead[LI[1], LI[0]];
        tderPatt = pertHead[LI[1], LI[1]];
        lapPatt = cd[-d][cd[d][pertHead[LI[1], LI[0]]]];
        
        reconstructed = Hh[LI[0], LI[0]]^2 * A[1, pertHead] * basePatt +
                       Hh[LI[0], LI[0]] * A[2, pertHead] * tderPatt -
                       (1/ah[LI[0], LI[0]]^2) * A[3, pertHead] * lapPatt;
      ],
      
      (* Quadratic A coefficients *)
      equationType == "quadratic-A",
      Block[{d = Symbol["Global`d"], h1, h2, base1, base2, tder1, tder2, lap1, lap2, gpair12, symWeight},
        {h1, h2} = pertHead;
        base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]];
        tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
        lap1 = cd[-d][cd[d][h1[LI[1], LI[0]]]]; lap2 = cd[-d][cd[d][h2[LI[1], LI[0]]]];
        gpair12 = cd[-d][base1] * cd[d][base2];
        gpair21 = cd[-d][base2] * cd[d][base1];
        
        (* A1 and A3 basis monomials are symmetric under h1<->h2, so with all-pairs summation
           use half-weight on off-diagonal pairs and full weight on diagonal pairs. *)
        symWeight = If[h1 === h2, 1, 1/2];
        
        (* A[i,h1,h2] coefficients have various H and a prefactors *)
        reconstructed = symWeight * (1/2) * Hh[LI[0], LI[0]]^2 * A[1, h1, h2] * base1 * base2 +
                       Hh[LI[0], LI[0]] * A[2, h1, h2] * tder1 * base2 -
                       symWeight * 3 * A[3, h1, h2] * tder1 * tder2 -
                       (2/ah[LI[0], LI[0]]^2) * A[4, h1, h2] * base1 * lap2 -
                       (1/(2*ah[LI[0], LI[0]]^2)) * A[5, h1, h2] * gpair12 +
                       (1/(Hh[LI[0], LI[0]]*ah[LI[0], LI[0]]^2)) * A[6, h1, h2] * tder1 * lap2;
      ],
      
      (* Linear B coefficients (0i equation) *)
      (* B[1,h] = (1/H^3) * coeff, B[2,h] = (1/H) * coeff *)
      equationType == "linear-B",
      Block[{c = Symbol["Global`c"], gradPatt, gradTderPatt},
        gradPatt = cd[-c][pertHead[LI[1], LI[0]]];
        gradTderPatt = cd[-c][pertHead[LI[1], LI[1]]];
        
        reconstructed = Hh[LI[0], LI[0]] * B[1, pertHead] * gradPatt - 
                       B[2, pertHead] * gradTderPatt;
      ],
      
      (* Quadratic B coefficients *)
      (* B[1,h1,h2] = (1/H) * coeff{base1*gradDown2}, B[2,h1,h2] = -coeff{tder1*gradDown2} *)
      (* B[3,h1,h2] = -coeff{base1*gradtDown2}, B[4,h1,h2] = H * coeff{tder1*gradtDown2} *)
      (* B[5,h1,h2] = -(H*a²/2) * (coeff{gradDown1*lap2} - coeff{gradUp1*gradgrad2}) *)
      equationType == "quadratic-B",
        Block[{c = Symbol["Global`c"], h1, h2, base1, base2, tder1, tder2, lap1, lap2,
          gradDown1, gradDown2, gradTder1, gradTder2, gradUp1, gradUp2, 
          gradgradDown1, gradgradDown2},
        {h1, h2} = pertHead;
        base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]];
        tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
        lap1 = cd[-d][cd[d][h1[LI[1], LI[0]]]]; lap2 = cd[-d][cd[d][h2[LI[1], LI[0]]]];
        gradDown1 = cd[-c][base1]; gradDown2 = cd[-c][base2];
        gradUp1 = cd[d][base1]; gradUp2 = cd[d][base2];
        gradTder1 = cd[-c][tder1]; gradTder2 = cd[-c][tder2];
        gradgradDown1 = cd[-d][cd[-c][base1]]; gradgradDown2 = cd[-d][cd[-c][base2]];
        
          (* Multiply back the extraction prefactors: H for B[1], identity for B[2], identity for B[3], 1/H for B[4], -1/(H*a²) for B[5] *)
          (* Note: quadratic-B loops over ALL pairs, so no need for swapped-index terms here *)
          reconstructed = Hh[LI[0], LI[0]] * B[1, h1, h2] * base1 * gradDown2 -
                   B[2, h1, h2] * tder1 * gradDown2 -
                   B[3, h1, h2] * base1 * gradTder2 +
                   (1/Hh[LI[0], LI[0]]) * B[4, h1, h2] * tder1 * gradTder2 -
                   (1/(Hh[LI[0], LI[0]]*ah[LI[0], LI[0]]^2)) * B[5, h1, h2] * ( gradDown1*lap2 - gradUp1*gradgradDown2 )
      ],
      
      (* Linear C coefficients (ij equation) *)
      (* CC[i,h] = (1/H^3, 1/H, 1/H, 1/H) * coeff *)
      equationType == "linear-C",
      Block[{basePatt, tderPatt, lapPatt, hbc},
        basePatt = pertHead[LI[1], LI[0]];
        lap = cd[-d][cd[d][pertHead[LI[1], LI[0]]]];
        tderPatt = pertHead[LI[1], LI[1]];
        ttderPatt = pertHead[LI[1], LI[2]];
        lapPatt = cd[-d][cd[d][pertHead[LI[1], LI[0]]]];
        gradDown = cd[-c][basePatt];
        gradUp = cd[d][basePatt];
        gradgradDown = cd[-c][cd[-b][basePatt]];
        hbc = h[-b, -c];
        
        reconstructed = Hh[LI[0], LI[0]]^2 * ah[LI[0], LI[0]]^2 * CC[1, pertHead] * hbc * basePatt +
                       Hh[LI[0], LI[0]] * ah[LI[0], LI[0]]^2 * CC[2, pertHead] * hbc * tderPatt -
                       ah[LI[0], LI[0]]^2 * CC[3, pertHead] * hbc * ttderPatt -
                       CC[4, pertHead] * (hbc * lap - gradgradDown );
      ],
      
      (* Quadratic C coefficients *)
      equationType == "quadratic-C",
         Block[{h1, h2, base1, base2, tder1, tder2, lap1, lap2,
           gradgrad1, gradgrad2, gpair12, gpair21, tgpair12, tgpair21, tgpairSwap12, hbc, symWeight},
        {h1, h2} = pertHead;
        base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]];
        tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
        ttder1 = h1[LI[1], LI[2]]; ttder2 = h2[LI[1], LI[2]];
        lap1 = cd[-d][cd[d][h1[LI[1], LI[0]]]]; lap2 = cd[-d][cd[d][h2[LI[1], LI[0]]]];
        gradgrad1 = cd[-c][cd[-b][base1]]; gradgrad2 = cd[-c][cd[-b][base2]];
        gpair12 = cd[-b][base1] * cd[-c][base2]; gpair21 = cd[-b][base2] * cd[-c][base1];
        tgpair12 = cd[-b][base1] * cd[-c][tder2]; tgpair21 = cd[-b][base2] * cd[-c][tder1];
        tgpairSwap12 = cd[-b][tder1] * cd[-c][base2];
        gpaircontr12 = cd[-d][base1] * cd[d][base2]; gpaircontr21 = cd[-d][base2] * cd[d][base1];
        tgpaircontr12 = cd[-d][base1] * cd[d][tder2]; tgpaircontr21 = cd[-d][base2] * cd[d][tder1];
         hbc = h[-b, -c];
        symWeight = If[h1 === h2, 1, 1/2];
        
         reconstructed = - symWeight * Hh[LI[0], LI[0]]^2 * ah[LI[0], LI[0]]^2 * CC[1, h1, h2] * hbc * base1 * base2 +
               Hh[LI[0], LI[0]] * ah[LI[0], LI[0]]^2 * CC[2, h1, h2] * hbc * base1 * tder2 -
               symWeight * ah[LI[0], LI[0]]^2 * CC[3, h1, h2] * hbc * tder1 * tder2 +
               ah[LI[0], LI[0]]^2 * CC[4, h1, h2] * hbc * base1 * ttder2 +
                     (ah[LI[0], LI[0]]^2/Hh[LI[0], LI[0]]) * CC[5, h1, h2] * hbc * tder1 * ttder2 +
                     CC[6, h1, h2] * hbc * gpaircontr12 +
                       (2/Hh[LI[0], LI[0]]) * CC[7, h1, h2] * hbc * tgpaircontr12 +
                     CC[8, h1, h2] * gpair12 -
                       (1/Hh[LI[0], LI[0]]) * CC[9, h1, h2] * (tgpair12 + tgpairSwap12) +
                       CC[10, h1, h2] * ( hbc * base1 * lap2 - base1 * gradgrad2 );
      ],
      
      (* Linear D coefficients (scalar equation) *)
      (* DD[1,h] = (1/H^2) * coeff, DD[2,h] = (1/H) * coeff, DD[3,h] = 1 * coeff, DD[4,h] = 1 * coeff *)
      equationType == "linear-D",
      Block[{b = Symbol["Global`b"], basePatt, tderPatt, ttderPatt, lapPatt},
        basePatt = pertHead[LI[1], LI[0]];
        tderPatt = pertHead[LI[1], LI[1]];
        ttderPatt = pertHead[LI[1], LI[2]];
        lapPatt = cd[-b][cd[b][pertHead[LI[1], LI[0]]]];
        
        reconstructed = - delta[-LI[0], LI[0]]^2 * (Hh[LI[0], LI[0]]/phi[LI[0], LI[1]]) * (Hh[LI[0], LI[0]]^2 * DD[1, pertHead] * basePatt +
                       Hh[LI[0], LI[0]] * DD[2, pertHead] * tderPatt +
                       DD[3, pertHead] * ttderPatt +
                       (1/ah[LI[0], LI[0]]^2) * DD[4, pertHead] * lapPatt );
      ],
      
      (* Quadratic D coefficients *)
      equationType == "quadratic-D",
            Block[{b, c, d, h1, h2, base1, base2, tder1, tder2, ttder1, ttder2, ttderPair12, lap1, lap2, gradPair12, gradPair21, gradtPair12, gradtPairSwap12, gradtPairSym12, gradtderPair12,
              gradgrad1, gradgrad2, gradgradUp1a, gradgradUp1b, gradgradUp2a, gradgradUp2b, lap2c, symWeight, gradtNorm},
        {h1, h2} = pertHead;
        b = Symbol["Global`b"]; c = Symbol["Global`c"]; d = Symbol["Global`d"];
        base1 = h1[LI[1], LI[0]]; base2 = h2[LI[1], LI[0]];
        tder1 = h1[LI[1], LI[1]]; tder2 = h2[LI[1], LI[1]];
        ttder1 = h1[LI[1], LI[2]]; ttder2 = h2[LI[1], LI[2]];
        ttderPair12 = tder1 * ttder2;
        gradtderPair12 = cd[-b][tder1] * cd[b][tder2];
        gradPair12 = cd[-b][base1] * cd[b][base2]; gradPair21 = cd[-b][base2] * cd[b][base1];
        gradtPair12 = cd[-b][tder1] * cd[b][base2];
        symWeight = If[h1 === h2, 1, 1/2];
        
        (* Define laplacian patterns for base fields only *)
        lap1 = cd[-b][cd[b][base1]]; lap2 = cd[-b][cd[b][base2]];
        lap2c = cd[-c][cd[c][base2]];
        gradgrad1 = cd[-c][cd[-b][base1]]; gradgrad2 = cd[-c][cd[-b][base2]];
        gradgradUp2a = cd[b][cd[c][base2]]; gradgradUp2b = cd[c][cd[b][base2]];
        
        reconstructed = - delta[-LI[0], LI[0]]^2 * (Hh[LI[0], LI[0]]/phi[LI[0], LI[1]]) * ( 
                       symWeight * Hh[LI[0], LI[0]]^2 * DD[1, h1, h2] * base1 * base2 +
                       3 * Hh[LI[0], LI[0]] * DD[2, h1, h2] * base1 * tder2 +
                       3 * symWeight * DD[3, h1, h2] * tder1 * tder2 +
                       3 * DD[4, h1, h2] * base1 * ttder2 -
                       (3/Hh[LI[0], LI[0]]) * DD[5, h1, h2] * ttderPair12 +
                       (1/ah[LI[0], LI[0]]^2) * DD[6, h1, h2] * base1 * lap2 +
                       (1/ah[LI[0], LI[0]]^2) * DD[7, h1, h2] * gradPair12 -
                       (1/(Hh[LI[0], LI[0]]*ah[LI[0], LI[0]]^2)) * DD[8, h1, h2] * tder1 * lap2 -
                       (2/(Hh[LI[0], LI[0]]*ah[LI[0], LI[0]]^2)) * DD[9, h1, h2] * gradtPair12 -
                       (1/(Hh[LI[0], LI[0]]^2*ah[LI[0], LI[0]]^2)) * DD[10, h1, h2] * ttder2 * lap2 +
                       symWeight * (1/(Hh[LI[0], LI[0]]^2*ah[LI[0], LI[0]]^2)) * DD[11, h1, h2] * gradtderPair12 -
                       symWeight * (1/(2 * Hh[LI[0], LI[0]]^2 * ah[LI[0], LI[0]]^4)) * DD[12, h1, h2] * (lap1 * lap2c - gradgrad1 * gradgradUp2b) );
      ],
      
      True,
      If[$xAlphaVerbose, Print["Error: Unknown equation type: ", equationType]];
      Return[$Failed];
    ];
    
    reconstructed
  ];

ValidateCoefficients[originalEqn_, equationType_String, pertHead_, removeScalarWrapper_:{}] :=
  Module[{reconstructed, residual, isComplete, filteredEqn, filteredCanon, reconstructedCanon},
    
    (* Filter original equation to only include terms with this perturbation *)
    filteredEqn = ExtractTermsForPerturbation[originalEqn, pertHead];
    filteredCanon = CanonicalizeGradientProducts[filteredEqn];
    
    reconstructed = ReconstructEquationFromCoefficients[equationType, pertHead, removeScalarWrapper];
    reconstructedCanon = CanonicalizeGradientProducts[reconstructed];
    
    If[reconstructed === $Failed,
      Return[{False, "Reconstruction failed"}]
    ];
    
    residual = Simplify[CanonicalizeGradientProducts[filteredCanon - reconstructedCanon]];
    isComplete = (residual === 0 || Simplify[residual] === 0);
    
    {isComplete, residual}
  ];

ValidateAllCoefficients[equations_Association, pertHeads_List, removeScalarWrapper_:{}] :=
  Module[{results, eq00Lin, eq00Quad, eq0iLin, eq0iQuad, eqijLin, eqijQuad, eqScalarLin, eqScalarQuad,
          isComplete, residual, failedValidations = {}},
    
    (* Extract equations from association *)
    eq00Lin = Lookup[equations, "eq00Linear", None];
    eq00Quad = Lookup[equations, "eq00Quadratic", None];
    eq0iLin = Lookup[equations, "eq0iLinear", None];
    eq0iQuad = Lookup[equations, "eq0iQuadratic", None];
    eqijLin = Lookup[equations, "eqijLinear", None];
    eqijQuad = Lookup[equations, "eqijQuadratic", None];
    eqScalarLin = Lookup[equations, "eqScalarLinear", None];
    eqScalarQuad = Lookup[equations, "eqScalarQuadratic", None];
    
    results = <||>;
    
    (* ===== Linear A Coefficients (00 equation) ===== *)
    If[eq00Lin =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eq00Lin, "linear-A", hh, removeScalarWrapper];
        results["A-linear-" <> ToString[hh]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"A-linear", hh, residual}];
        ];
      , {hh, pertHeads}];
    ];
    
    (* ===== Quadratic A Coefficients (00 equation) ===== *)
    If[eq00Quad =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eq00Quad, "quadratic-A", {h1, h2}, removeScalarWrapper];
        results["A-quadratic-" <> ToString[h1] <> "-" <> ToString[h2]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"A-quadratic", {h1, h2}, residual}];
        ];
      , {h1, pertHeads}, {h2, pertHeads}];
    ];
    
    (* ===== Linear B Coefficients (0i equation) ===== *)
    If[eq0iLin =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eq0iLin, "linear-B", hh, removeScalarWrapper];
        results["B-linear-" <> ToString[hh]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"B-linear", hh, residual}];
        ];
      , {hh, pertHeads}];
    ];
    
    (* ===== Quadratic B Coefficients (0i equation) ===== *)
    If[eq0iQuad =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eq0iQuad, "quadratic-B", {h1, h2}, removeScalarWrapper];
        results["B-quadratic-" <> ToString[h1] <> "-" <> ToString[h2]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"B-quadratic", {h1, h2}, residual}];
        ];
      , {h1, pertHeads}, {h2, pertHeads}];
    ];
    
    (* ===== Linear C Coefficients (ij equation) ===== *)
    If[eqijLin =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eqijLin, "linear-C", hh, removeScalarWrapper];
        results["C-linear-" <> ToString[hh]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"C-linear", hh, residual}];
        ];
      , {hh, pertHeads}];
    ];
    
    (* ===== Quadratic C Coefficients (ij equation) ===== *)
    If[eqijQuad =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eqijQuad, "quadratic-C", {h1, h2}, removeScalarWrapper];
        results["C-quadratic-" <> ToString[h1] <> "-" <> ToString[h2]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"C-quadratic", {h1, h2}, residual}];
        ];
      , {h1, pertHeads}, {h2, pertHeads}];
    ];
    
    (* ===== Linear D Coefficients (scalar equation) ===== *)
    If[eqScalarLin =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eqScalarLin, "linear-D", hh, removeScalarWrapper];
        results["D-linear-" <> ToString[hh]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"D-linear", hh, residual}];
        ];
      , {hh, pertHeads}];
    ];
    
    (* ===== Quadratic D Coefficients (scalar equation) ===== *)
    If[eqScalarQuad =!= None,
      Do[
        {isComplete, residual} = ValidateCoefficients[eqScalarQuad, "quadratic-D", {h1, h2}, removeScalarWrapper];
        results["D-quadratic-" <> ToString[h1] <> "-" <> ToString[h2]] = {isComplete, residual};
        If[!isComplete, 
          AppendTo[failedValidations, {"D-quadratic", {h1, h2}, residual}];
        ];
      , {h1, pertHeads}, {h2, pertHeads}];
    ];
    
    (* ===== Summary ===== *)
    If[$xAlphaVerbose, Print["Coefficient validation summary -> Passed: ", Count[Values[results], {True, _}], ", Failed: ", Count[Values[results], {False, _}]]];
    
    {results, failedValidations}
  ];


(* Note: display formatting for alpha/gamma symbols is handled by PrintAs
   in Setup.wl via DefScalarFunction. Do not add Format[] rules here for
   those symbols — it causes Symbol::symname errors in xAct internals. *)

(* ========================================================================== *)
(* ShowPaperTables: display linear and quadratic coefficient tables           *)
(* ========================================================================== *)

ShowPaperTables[linA_, linB_, linC_, linD_, quadA_, quadB_, quadC_, quadD_] :=
  Module[{phih, psih, Qvar, qData, lData, checkLin, checkQuad,
    linDefs, quadDefs, linRows, quadRows, gridLin, gridQuad,
    eqnE1, eqnA1, eqnP1, eqnS1, eqnE2, eqnA2, eqnP2, eqnS2,
    Y1, Y1d, Y1dd, Y2, Y2d, Y2dd, Y0, Y0d, Y0dd,
    del2, di, dj, dij, DijSym},

    phih  = Symbol["Global`\[Phi]h"];
    psih  = Symbol["Global`\[Psi]h"];
    Qvar  = Symbol["Global`Q"];

    lData = {{phih, "\[Phi]"}, {psih, "\[Psi]"}, {Qvar, "Q"}};
    qData = {
      {phih, phih, "\[Phi]\[Phi]"}, {psih, psih, "\[Psi]\[Psi]"},
      {Qvar, Qvar, "QQ"}, {phih, psih, "\[Phi]\[Psi]"},
      {phih, Qvar, "\[Phi]Q"}, {psih, Qvar, "\[Psi]Q"},
      {psih, phih, "\[Psi]\[Phi]"}, {Qvar, phih, "Q\[Phi]"},
      {Qvar, psih, "Q\[Psi]"}};

    (* Column labels *)
    eqnE1 = Style[Superscript["\[ScriptE]", "(1)"], Bold, 14];
    eqnA1 = Style[Subsuperscript["\[ScriptA]", "i", "(1)"], Bold, 14];
    eqnP1 = Style[Subsuperscript["\[ScriptP]", "ij", "(1)"], Bold, 14];
    eqnS1 = Style[Superscript["\[ScriptS]", "(1)"], Bold, 14];
    eqnE2 = Style[Superscript["\[ScriptE]", "(2)"], Bold, 14];
    eqnA2 = Style[Subsuperscript["\[ScriptA]", "i", "(2)"], Bold, 14];
    eqnP2 = Style[Subsuperscript["\[ScriptP]", "ij", "(2)"], Bold, 14];
    eqnS2 = Style[Superscript["\[ScriptS]", "(2)"], Bold, 14];

    (* Formatting shorthands *)
    Y1   = Superscript["Y", "a"];
    Y1d  = Superscript[Overscript["Y", "."], "a"];
    Y1dd = Superscript[Overscript["Y", ".."], "a"];
    Y2   = Superscript["Y", "b"];
    Y2d  = Superscript[Overscript["Y", "."], "b"];
    Y2dd = Superscript[Overscript["Y", ".."], "b"];
    Y0   = "Y";
    Y0d  = Superscript[Overscript["Y", "."], ""];
    Y0dd = Superscript[Overscript["Y", ".."], ""];
    del2 = Superscript["\[Del]", "2"];
    di   = Subscript["\[PartialD]", "i"];
    dj   = Subscript["\[PartialD]", "j"];
    dij  = Subscript["\[Delta]", "ij"];
    DijSym = Subscript["\[ScriptD]", "ij"];

    (* Evaluators *)
    checkLin[arr_, letter_, idx_, i_] :=
      Module[{val = arr[idx, lData[[i, 1]]]},
        If[val === 0 || PossibleZeroQ[val],
          Style["0", Gray, 14],
          Style[Subsuperscript[letter, idx, lData[[i, 2]]], Bold, Darker[Green], 14]]];

    checkQuad[arr_, letter_, isSym_, idx_, i_] :=
      Module[{val},
        If[isSym && i >= 7,
          Return[If[i == 7,
            Item[Style["(sym)", Italic, Gray], Alignment -> Center],
            SpanFromLeft]]];
        val = arr[idx, qData[[i, 1]], qData[[i, 2]]];
        If[val === 0 || PossibleZeroQ[val],
          Style["0", Gray, 14],
          Style[Subsuperscript[letter, idx, qData[[i, 3]]], Bold, Darker[Green], 14]]];

    (* Linear row definitions: {EqnLabel, Term, Array, Letter, Index} *)
    linDefs = {
      {eqnE1, Y1,                    linA, "A", 1},
      {SpanFromAbove, Y1d,            linA, "A", 2},
      {SpanFromAbove, Row[{del2,Y1}], linA, "A", 3},
      {eqnA1, Row[{di,Y1}],          linB, "B", 1},
      {SpanFromAbove, Row[{di,Y1d}],  linB, "B", 2},
      {eqnP1, Row[{dij,Y1}],         linC, "C", 1},
      {SpanFromAbove, Row[{dij,Y1d}], linC, "C", 2},
      {SpanFromAbove, Row[{dij,Y1dd}],linC, "C", 3},
      {SpanFromAbove, Row[{DijSym,Y1}],linC,"C", 4},
      {eqnS1, Y1,                    linD, "D", 1},
      {SpanFromAbove, Y1d,            linD, "D", 2},
      {SpanFromAbove, Y1dd,           linD, "D", 3},
      {SpanFromAbove, Row[{del2,Y1}], linD, "D", 4}};

    (* Quadratic row definitions: {EqnLabel, Term, Array, Letter, isSym, Index} *)
    quadDefs = {
      {eqnE2, Row[{Y1," ",Y2}],                               quadA,"A",True, 1},
      {SpanFromAbove, Row[{Y1d," ",Y2}],                      quadA,"A",False,2},
      {SpanFromAbove, Row[{Y1d," ",Y2d}],                     quadA,"A",True, 3},
      {SpanFromAbove, Row[{Y1,del2,Y2}],                      quadA,"A",False,4},
      {SpanFromAbove, Row[{"\[PartialD]Y"," ","\[PartialD]Y"}],quadA,"A",True, 5},
      {SpanFromAbove, Row[{Y1d,del2,Y2}],                     quadA,"A",False,6},
      {eqnA2, Row[{Y1,di,Y2}],                                quadB,"B",False,1},
      {SpanFromAbove, Row[{Y1d,di,Y2}],                       quadB,"B",False,2},
      {SpanFromAbove, Row[{Y1,di,Y2d}],                       quadB,"B",False,3},
      {SpanFromAbove, Row[{Y1d,di,Y2d}],                      quadB,"B",False,4},
      {SpanFromAbove, Row[{"\[PartialD]Y"," ","\[ScriptD]"," Y"}],quadB,"B",False,5},
      {eqnP2, Row[{dij,Y0," ",Y0}],                          quadC,"C",True, 1},
      {SpanFromAbove, Row[{dij,Y0," ",Y0d}],                  quadC,"C",False,2},
      {SpanFromAbove, Row[{dij,Y0d," ",Y0d}],                 quadC,"C",True, 3},
      {SpanFromAbove, Row[{dij,Y0," ",Y0dd}],                 quadC,"C",False,4},
      {SpanFromAbove, Row[{dij,Y0d," ",Y0dd}],                quadC,"C",False,5},
      {SpanFromAbove, Row[{dij,"\[PartialD]Y"," ","\[PartialD]Y"}],quadC,"C",True, 6},
      {SpanFromAbove, Row[{dij,"\[PartialD]Y"," ","\[PartialD]",Y0d}],quadC,"C",False,7},
      {SpanFromAbove, Row[{di,Y0," ",dj,Y0}],                 quadC,"C",True, 8},
      {SpanFromAbove, Row[{di,Y0," ",dj,Y0d}],                quadC,"C",False,9},
      {SpanFromAbove, Row[{Y0," ",DijSym,Y0}],                quadC,"C",False,10},
      {eqnS2, Row[{Y1," ",Y2}],                               quadD,"D",True, 1},
      {SpanFromAbove, Row[{Y1," ",Y2d}],                      quadD,"D",False,2},
      {SpanFromAbove, Row[{Y1d," ",Y2d}],                     quadD,"D",True, 3},
      {SpanFromAbove, Row[{Y1," ",Y2dd}],                     quadD,"D",False,4},
      {SpanFromAbove, Row[{Y1d," ",Y2dd}],                    quadD,"D",False,5},
      {SpanFromAbove, Row[{Y1,del2,Y2}],                      quadD,"D",False,6},
      {SpanFromAbove, Row[{di,Y1," ",di,Y2}],                 quadD,"D",True, 7},
      {SpanFromAbove, Row[{Y1d,del2,Y2}],                     quadD,"D",False,8},
      {SpanFromAbove, Row[{di,Y1d," ",di,Y2}],                quadD,"D",False,9},
      {SpanFromAbove, Row[{Y1dd,del2,Y2}],                    quadD,"D",False,10},
      {SpanFromAbove, Row[{di,Y1d," ",di,Y2d}],               quadD,"D",True, 11},
      {SpanFromAbove, Superscript["\[ScriptY]","a(2)"],        quadD,"D",True, 12}};

    (* Build rows *)
    linRows = Join[
      {{Style["Eqn",Bold], Style["Term",Bold],
        Style["\[Phi]",Bold,14], Style["\[Psi]",Bold,14], Style["Q",Bold,14]}},
      Map[{#[[1]], #[[2]],
          Sequence @@ Table[checkLin[#[[3]], #[[4]], #[[5]], i], {i, 3}]} &,
        linDefs]];

    quadRows = Join[
      {Prepend[
          Map[Style[#, Bold, 14] &, qData[[All, 3]]],
          Style["Term", Bold]] ~Prepend~ Style["Eqn", Bold]},
      Map[{#[[1]], #[[2]],
          Sequence @@ Table[checkQuad[#[[3]], #[[4]], #[[5]], #[[6]], i], {i, 9}]} &,
        quadDefs]];

    (* Render *)
    gridLin = Grid[linRows,
      Frame -> All,
      Dividers -> {All, {1->True, 2->True, 5->True, 7->True, 11->True, -1->True}},
      Spacings -> {2, 1.5},
      Background -> {{LightGray, LightGray, None}, {LightGray, None}}];

    gridQuad = Grid[quadRows,
      Frame -> All,
      Dividers -> {All, {1->True, 2->True, 8->True, 13->True, 23->True, -1->True}},
      Spacings -> {2, 1.5},
      Background -> {{LightGray, LightGray, None}, {LightGray, None}}];

    Print[Panel@Column[
      {Style["TABLE 1: Linear Coefficients", Bold, 16], gridLin},
      Alignment -> Center]];
    Print["\n"];
    Print[Panel@Column[
      {Style["TABLE 2: Second-Order Coefficients", Bold, 16], gridQuad},
      Alignment -> Center]];
  ];

End[];