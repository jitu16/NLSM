(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMAlgebra`"];

(* --- 1. FUNCTION EXPORTS --- *)
FinalSimplify::usage = "FinalSimplify[expr] fully simplifies a tensor expression.";
Para::usage = "Para[expr] projects onto the longitudinal component.";
Perp::usage = "Perp[expr] projects onto the transverse component.";
MySimplify::usage = "MySimplify[expr] applies basic algebraic simplification.";
TraceSimplify::usage = "TraceSimplify[expr] handles trace identities and constants.";
Exchange::usage = "Exchange[expr, x, y] swaps variables.";

(* --- 2. SYMBOL EXPORTS (THE SAFETY LOCK) --- *)
(* By defining usage here, we force these to be the SAME symbol 
   in your notebook and the package. *)

(* Fields *)
W::usage = "W field symbol.";
Phi::usage = "Phi field symbol (English P).";
\[Phi]::usage = "Phi field symbol (Small Greek).";
\[CapitalPhi]::usage = "Phi field symbol (Capital Greek).";
\[CapitalEpsilon]::usage = "Epsilon/Ubareu symbol.";
\[CapitalPi]::usage = "Pi field symbol.";
\[Sigma]::usage = "Sigma field symbol.";

(* Operators & Params *)
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

rule = {
  Tr[Subscript[\[Sigma], 3]] -> 0, 
  Tr[\[CapitalPi][0]] -> 0, 
  Tr[\[CapitalPi][0] . Subscript[\[Sigma], 3]] -> 0, 
  Tr[Subscript[\[Sigma], 3] . \[CapitalPi][0]] -> 0, 
  Tr[Id] -> 0, Tr[1] -> 0, Tr[0] -> 0, 
  Tr[Der[x_]] -> 0, 
  Tr[Der[x_] . Der[y_]] -> 0, 
  Tr[Der[x_] . Subscript[\[Sigma], 3]] -> 0, 
  Tr[Der[x_] . Der[y_] . Subscript[\[Sigma], 3]] -> 0
};

commutationRules = {
  Subscript[\[Sigma], 3] . Subscript[\[Sigma], 3] :> 1, 
  1 . x_ :> x, x_ . 1 :> x, 
  Tr[x___ . minus . y___] :> -Tr[x . y], 
  Tr[Subscript[\[Sigma], 3] . (rest___)] /; ! MatchQ[{rest}, {Subscript[\[Sigma], 3], ___}] :> Tr[rest . Subscript[\[Sigma], 3]], 
  Subscript[\[Sigma], 3] . \[CapitalPi][arg_] :> \[CapitalPi][arg] . Subscript[\[Sigma], 3], 
  Subscript[\[Sigma], 3] . Der[arg_] :> Der[arg] . Subscript[\[Sigma], 3], 
  Subscript[\[Sigma], 3] . Para[any_] :> Para[any] . Subscript[\[Sigma], 3], 
  Subscript[\[Sigma], 3] . Perp[any_] :> minus . Perp[any] . Subscript[\[Sigma], 3]
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
TraceSimplify[Tr[x___ . \[CapitalPi][0] . y___]] = TraceSimplify[Tr[x . y]] \[CapitalPi][0];
TraceSimplify[Tr[x_]] := MySimplify[Tr[x]];
TraceSimplify[x_?NumericQ] := x;

(* Identity rules for constants *)
TraceSimplify[\[Lambda]] = \[Lambda]; TraceSimplify[\[Lambda]^n_] := \[Lambda]^n;
TraceSimplify[t] := t; TraceSimplify[t^n_] := t^n;
TraceSimplify[\[Nu]0] := \[Nu]0; TraceSimplify[\[Nu]0^n_] := \[Nu]0^n;
TraceSimplify[Dp] := Dp; TraceSimplify[Dp^n_] := Dp^n;
TraceSimplify[D0] := D0; TraceSimplify[D0^n_] := D0^n;
TraceSimplify[\[CapitalEth]0] := \[CapitalEth]0; TraceSimplify[\[CapitalEth]0^n_] := \[CapitalEth]0^n;
TraceSimplify[\[CapitalPi][0]] := \[CapitalPi][0]; TraceSimplify[\[CapitalPi][0]^n_] := \[CapitalPi][0]^n;

(* --- CYCLIC OPERATIONS --- *)
(* Now robust because \[CapitalPhi] and Phi are explicitly exported above *)

cycle\[CapitalPhi]toFirst[x_ y_] := cycle\[CapitalPhi]toFirst[x] cycle\[CapitalPhi]toFirst[y];
cycle\[CapitalPhi]toFirst[x_^n_] := (cycle\[CapitalPhi]toFirst[x])^n;
cycle\[CapitalPhi]toFirst[x_ + y_] := cycle\[CapitalPhi]toFirst[x] + cycle\[CapitalPhi]toFirst[y];

(* This pattern matches properly now *)
cycle\[CapitalPhi]toFirst[Tr[x___ . \[CapitalPhi][a_] . y___]] := Tr[\[CapitalPhi][a] . y . x];
cycle\[CapitalPhi]toFirst[Tr[x___ . Phi[a_] . y___]] := Tr[Phi[a] . y . x];

cycle\[CapitalPhi]toFirst[Tr[x_]] := Tr[x];

cycle\[CapitalPhi]toFirst[x_?NumericQ] := x;
cycle\[CapitalPhi]toFirst[t] := t; cycle\[CapitalPhi]toFirst[t^n_] := t^n;
cycle\[CapitalPhi]toFirst[\[Lambda]] := \[Lambda]; cycle\[CapitalPhi]toFirst[\[Lambda]^n_] := \[Lambda]^n;
cycle\[CapitalPhi]toFirst[\[Nu]0] := \[Nu]0; cycle\[CapitalPhi]toFirst[\[Nu]0^n_] := \[Nu]0^n;
cycle\[CapitalPhi]toFirst[Dp] := Dp; cycle\[CapitalPhi]toFirst[Dp^n_] := Dp^n;
cycle\[CapitalPhi]toFirst[D0] := D0; cycle\[CapitalPhi]toFirst[D0^n_] := D0^n;
cycle\[CapitalPhi]toFirst[\[CapitalEth]0] := \[CapitalEth]0; cycle\[CapitalPhi]toFirst[\[CapitalEth]0^n_] := \[CapitalEth]0^n;
cycle\[CapitalPhi]toFirst[\[CapitalPi][0]] := \[CapitalPi][0]; cycle\[CapitalPhi]toFirst[\[CapitalPi][0]^n_] := \[CapitalPi][0]^n;

(* --- UTILITIES --- *)
Exchange[expr_, x_, y_] := (expr /. {x -> y1, y -> x}) /. {y1 -> y};
FinalSimplify[x_] := ((cycle\[CapitalPhi]toFirst[TraceSimplify[TensorExpand[x]]]) //. rule) // FullSimplify;

End[];
EndPackage[];
