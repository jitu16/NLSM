(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMLisp`"];

(* Public Exported Symbols *)
ToLispString::usage = "ToLispString[expr] converts a Mathematica expression to Lisp S-expression.";
EvaluateModel::usage = "EvaluateModel[lispCommand] runs the external Lisp script.";

Begin["`Private`"];

(* --- 1. ROBUST TRANSLATOR LOGIC --- *)
ClearAll[recToLisp];

(* HELPER: safeName
   Safely gets the name of a symbol, handling special characters. *)
safeName[s_Symbol] := SymbolName[s];
safeName[other_] := "";

(* --- SPECIFIC FIELDS (Robust Matching) --- *)

(* 1. STANDARD FIELDS (No Derivatives) *)
(* Matches Phi[x], W[x], etc. by name *)
recToLisp[head_[x_]] /; MemberQ[{"Phi", "CapitalPhi", "\[CapitalPhi]"}, safeName[head]] := 
  "(PHI " <> recToLisp[x] <> ")";

recToLisp[head_[x_]] /; MemberQ[{"Phi", "\[Phi]"}, safeName[head]] := 
  "(VARPHI " <> recToLisp[x] <> ")";

recToLisp[head_[x_]] /; MemberQ[{"W"}, safeName[head]] := 
  "(W " <> recToLisp[x] <> ")";

recToLisp[head_[x_]] /; MemberQ[{"CapitalEpsilon", "\[CapitalEpsilon]"}, safeName[head]] := 
  "(UBAREU " <> recToLisp[x] <> ")";

(* 2. FIELDS WITH DERIVATIVES (The Fix) *)
(* We match d_[y_] and check if d is named "Der" *)

recToLisp[head_[x_, d_[y_]]] /; (MemberQ[{"Phi", "CapitalPhi", "\[CapitalPhi]"}, safeName[head]] && safeName[d] === "Der") := 
  "(PHI " <> recToLisp[x] <> " (d " <> recToLisp[y] <> "))";

recToLisp[head_[x_, d_[y_]]] /; (MemberQ[{"Phi", "\[Phi]"}, safeName[head]] && safeName[d] === "Der") := 
  "(VARPHI " <> recToLisp[x] <> " (d " <> recToLisp[y] <> "))";

recToLisp[head_[x_, d_[y_]]] /; (MemberQ[{"W"}, safeName[head]] && safeName[d] === "Der") := 
  "(W " <> recToLisp[x] <> " (d " <> recToLisp[y] <> "))";

recToLisp[head_[x_, d_[y_]]] /; (MemberQ[{"CapitalEpsilon", "\[CapitalEpsilon]"}, safeName[head]] && safeName[d] === "Der") := 
  "(UBAREU " <> recToLisp[x] <> " (d " <> recToLisp[y] <> "))";


(* --- SIGMA HANDLING --- *)
(* Robustly catches Sigma in subscript, regardless of context *)
recToLisp[Subscript[s_Symbol, 3]] /; MemberQ[{"Sigma", "\[Sigma]"}, safeName[s]] := "SIGMA3";
recToLisp[s_Symbol] /; safeName[s] == "\[Sigma]3" := "SIGMA3";


(* --- OPERATIONS --- *)
recToLisp[d_Dot] := StringRiffle[recToLisp /@ (List @@ d), " "];
recToLisp[Tr[content_]] := "(tr " <> recToLisp[content] <> ")";
recToLisp[t_Times] := "(* " <> StringRiffle[recToLisp /@ (List @@ t), " "] <> ")";
recToLisp[p_Plus] := "(+ " <> StringRiffle[recToLisp /@ (List @@ p), " "] <> ")";
recToLisp[Power[x_, n_]] := "(expt " <> recToLisp[x] <> " " <> recToLisp[n] <> ")";
recToLisp[Rational[n_, d_]] := "(/ " <> ToString[n] <> " " <> ToString[d] <> ")";
recToLisp[n_Integer] := ToString[n];
recToLisp[r_Real] := ToString[r];

(* --- SAFETY NET --- *)
(* 1. Catch generic symbols (x, y, r) *)
recToLisp[s_Symbol] := ToString[s];

(* 2. Catch generic Subscripts that aren't Sigma3 *)
(* Subscript[A, B] -> "A_B" style string, no commas *)
recToLisp[Subscript[a_, b_]] := recToLisp[a] <> ToString[b];

(* 3. Catch unknown functions - remove commas to prevent crashes *)
recToLisp[x_] := StringReplace[ToString[InputForm[x]], "," -> " "];

(* Wrapper *)
ToLispString[expr_] := "(main '" <> recToLisp[expr] <> ")";


(* --- EXECUTION ENGINE --- *)
EvaluateModel[lispCommand_String] := Module[{scriptPath, sbclPath, processResult, rawOutput, cleanString, data},
   (* UPDATE THESE PATHS IF NEEDED *)
   scriptPath = "/Users/zjitu/Documents/Lisp/NLSM/nlsm.lisp";
   sbclPath = "/opt/homebrew/bin/sbcl";

   processResult = RunProcess[{sbclPath, "--script", scriptPath, lispCommand}];

   If[processResult["ExitCode"] =!= 0, 
    Print["Lisp Script Failed!"];
    Print[processResult["StandardError"]];
    Return[$Failed];
   ];

   rawOutput = processResult["StandardOutput"];
   cleanString = First[StringCases[rawOutput, "<|" ~~ ___ ~~ "|>"], "" ];
   
   If[cleanString === "", 
    Print["Extraction Failed. No valid Association found."];
    Print["Raw Output: ", rawOutput];
    Return[$Failed];
   ];

   data = ToExpression[cleanString];
   Return[data];
];

End[];
EndPackage[];
