(* ::Package:: *)

(* ::Package:: *)
(**)


BeginPackage["NLSMAlgebra`"];

(* --- 1. FUNCTION EXPORTS --- *)
FinalSimplify::usage = "FinalSimplify[expr] fully simplifies a tensor expression.";
Para::usage = "Para[expr] projects onto the longitudinal component.";
Perp::usage = "Perp[expr] projects onto the transverse component.";
MySimplify::usage = "MySimplify[expr] applies basic algebraic simplification.";
TraceSimplify::usage = "TraceSimplify[expr] handles trace identities and constants.";
Exchange::usage = "Exchange[expr, x, y] swaps variables.";

(* --- 2. SYMBOL EXPORTS --- *)
W::usage = "W field symbol.";
Phi::usage = "Phi field symbol (English P).";
\[Phi]::usage = "Phi field symbol (Small Greek).";
\[CapitalPhi]::usage = "Phi field symbol (Capital Greek).";
\[CapitalEpsilon]::usage = "Epsilon/Ubareu symbol.";
\[CapitalPi]::usage = "Pi field symbol.";
\[Sigma]::usage = "Sigma field symbol.";
\[CapitalEth]::usage = "Diffusion Model.";

Der::usage = "Derivative symbol wrapper.";
Id::usage = "Identity operator.";
t::usage = "Coupling t.";
\[Lambda]::usage = "Coupling Lambda.";
\[Nu]0::usage = "Parameter nu0.";
Dp::usage = "Parameter Dp.";
D0::usage = "Parameter D0.";
\[CapitalEth]0::usage = "Parameter Eth0.";

Begin["`Private`"];

(* --- CONFIGURATION: Rules & Commutations --- *)

(* 1. Final Substitution Rules (Run at the very end) (**)
rule = {
  Tr[Subscript[\[Sigma], 3]] -> 0,
  Tr[Subscript[\[Sigma], 1]] -> 0, 
  Tr[\[CapitalPi][0]] -> 0, 
  Tr[\[CapitalPi][0] . Subscript[\[Sigma], 3]] -> 0, 
  Tr[Subscript[\[Sigma], 3] . \[CapitalPi][0]] -> 0, 
  Tr[Id] -> 0, Tr[1] -> 0, Tr[0] -> 0, 
  Tr[Der[x_]] -> 0, 
  Tr[Der[x_] . Der[y_]] -> 0, 
  Tr[Der[x_] . Subscript[\[Sigma], 3]] -> 0, 
  Tr[Der[x_] . Der[y_] . Subscript[\[Sigma], 3]] -> 0
};*)


rule = {
  Tr[Subscript[\[Sigma], 3]] -> 0,
  Tr[Subscript[\[Sigma], 1]] -> 0, 
  Tr[Id] -> 0, Tr[1] -> 0, Tr[0] -> 0, 
  Tr[Der[x_]] -> 0, 
  Tr[Der[x_] . Der[y_]] -> 0, 
  Tr[Der[x_] . Subscript[\[Sigma], 3]] -> 0, 
  Tr[Der[x_] . Der[y_] . Subscript[\[Sigma], 3]] -> 0,
  Tr[y___ . \[CapitalPi][x_] . z___]:>Tr[y . z](DR[x]+DA[x])+Tr[y . Subscript[\[Sigma], 3] . z](DR[x]-DA[x])
};

commutationRules = {
  (* --- 1. EXPLICIT 2-TERM CASES (The Fix) --- *)
  (* Catches "S3 . S3" when they are the ONLY things in the Dot *)
  Subscript[\[Sigma], 3] . Subscript[\[Sigma], 3] :> 1,
  1 . x_ :> x, 
  x_ . 1 :> x,

  (* --- 2. SEQUENCE CASES (Chains) --- *)
  (* Catches "A . 1 . B" -> "A . B" *)
  Dot[1, rest__] :> Dot[rest],         
  Dot[front__, 1] :> Dot[front],       
  Dot[front__, 1, rest__] :> Dot[front, rest],
  
  (* Catches "A . S3 . S3 . B" -> "A . B" *)
  Dot[front___, Subscript[\[Sigma], 3], Subscript[\[Sigma], 3], rest___] :> Dot[front, rest],

  (* --- 3. STANDARD LOGIC --- *)
  Tr[x___ . minus . y___] :> -Tr[x . y],
  
  Subscript[\[Sigma], 3] . \[CapitalPi][arg_] :> \[CapitalPi][arg] . Subscript[\[Sigma], 3],
  Subscript[\[Sigma], 3] . Der[arg_] :> Der[arg] . Subscript[\[Sigma], 3],
  Subscript[\[Sigma], 3] . Para[any_] :> Para[any] . Subscript[\[Sigma], 3],
  Subscript[\[Sigma], 3] . Perp[any_] :> minus . Perp[any] . Subscript[\[Sigma], 3],
  
  Tr[Subscript[\[Sigma], 3] . rest__] :> Tr[rest . Subscript[\[Sigma], 3]]
};

(* --- PROJECTION LOGIC --- *)

Para[Perp[x_]] := 0;
Perp[Para[x_]] := 0;
Para[Para[x_]] := Para[x];
Perp[Perp[x_]] := Perp[x];

Para[Der[x_]] := Der[x];
Perp[Der[x_]] := 0;

Para[\[CapitalPi][x_]] := \[CapitalPi][x];
Perp[\[CapitalPi][x_]] := 0;

Para[Subscript[\[Sigma], 3]] := Subscript[\[Sigma], 3];
Perp[Subscript[\[Sigma], 3]] := 0;

Para[Subscript[\[Sigma], 1]] := 0;
Perp[Subscript[\[Sigma], 1]] := Subscript[\[Sigma], 1];

Para[Id] := 1; Perp[Id] := 0;
Para[1] := 1; Perp[1] := 0;
Para[0] := 0; Perp[0] := 0;

Para[x_ . y_] := Para[x] . Para[y] + Perp[x] . Perp[y];
Perp[x_ . y_] := Para[x] . Perp[y] + Perp[x] . Para[y];

(* --- SIMPLIFICATION --- *)

MySimplify[A_ + B_] := MySimplify[A] + MySimplify[B];
MySimplify[A_ B_] := MySimplify[A] MySimplify[B];
MySimplify[A_] := A;
MySimplify[Tr[x_]] := \[Sigma]3Simplify[Tr[x]];

\[Sigma]3Simplify[expr_] := expr //. commutationRules;

TraceSimplify[x_ y_] := TraceSimplify[x] TraceSimplify[y];
TraceSimplify[x_^n_] := (TraceSimplify[x])^n;
TraceSimplify[x_ + y_] := TraceSimplify[x] + TraceSimplify[y];
TraceSimplify[Tr[x_ + y_]] := TraceSimplify[Tr[x]] + TraceSimplify[Tr[y]];

(* CRITICAL FIX: Changed = to := so it evaluates at runtime, not definition time *)
TraceSimplify[Tr[x___ . \[CapitalPi][0] . y___]] := TraceSimplify[Tr[x . y]] \[CapitalEth][0];

TraceSimplify[Tr[x_]] := MySimplify[Tr[x]];
TraceSimplify[x_?NumericQ] := x;
TraceSimplify[s_Symbol] := s;
TraceSimplify[s_Symbol[x_]] := s[x];

(* --- CYCLIC OPERATIONS (Minimal & Fast) --- *)

cycle\[CapitalPhi]toFirst[x_ y_] := cycle\[CapitalPhi]toFirst[x] cycle\[CapitalPhi]toFirst[y];
cycle\[CapitalPhi]toFirst[x_^n_] := (cycle\[CapitalPhi]toFirst[x])^n;
cycle\[CapitalPhi]toFirst[x_ + y_] := cycle\[CapitalPhi]toFirst[x] + cycle\[CapitalPhi]toFirst[y];

(* Target CapitalPhi to normalize cyclic order *)
cycle\[CapitalPhi]toFirst[Tr[x___ . target_ . y___]] /; (!FreeQ[target, \[CapitalPhi]]) := 
   Tr[target . y . x];

cycle\[CapitalPhi]toFirst[Tr[x_]] := Tr[x];
cycle\[CapitalPhi]toFirst[x_?NumericQ] := x;
cycle\[CapitalPhi]toFirst[s_Symbol] := s;
cycle\[CapitalPhi]toFirst[s_Symbol[x_]] := s[x];

(* --- UTILITIES --- *)
Exchange[expr_, x_, y_] := (expr /. {x -> y1, y -> x}) /. {y1 -> y};
FinalSimplify[x_] := ((cycle\[CapitalPhi]toFirst[TraceSimplify[TensorExpand[x]]]) //. rule) // FullSimplify;

End[];
EndPackage[];
