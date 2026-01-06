(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMDriver`", {"NLSMAlgebra`", "NLSMLisp`", "NLSMVisuals`"}];

(* --- EXPORTS --- *)
GenerateDiagrams::usage = "GenerateDiagrams[expr, outputDir, prefix, options...] processes a tensor expression term-by-term, generating diagrams and saving reports.
Arguments:
  expr: The sum of terms to process.
  outputDir: String path (relative or absolute) for the output folder.
  prefix: String prefix for filenames (e.g., 'Diagram').
Options:
  \"Display\" -> True (Show grids in notebook) | False (Save only, return paths)
  \"Factor\" -> (Multiplier applied to every term, default 1)";

Begin["`Private`"];

(* --- 1. THE SPLITTER: Isolates Scalar Factors from Tensor Traces --- *)
IsolateFactor[term_] := Module[{factors, traces, listContent},
  (* Case A: Product (a * b * Tr[...]) *)
  If[Head[term] === Times,
     (* Flatten product list *)
     listContent = List @@ term;
     (* Select items that DO NOT contain Tr/Trace *)
     factors = Select[listContent, FreeQ[#, Tr] && FreeQ[#, Trace] &];
     (* Select items that DO contain Tr/Trace *)
     traces  = Select[listContent, (!FreeQ[#, Tr] || !FreeQ[#, Trace]) &];
     
     (* Return {ScalarProduct, TensorProduct} *)
     {Times @@ factors, Times @@ traces},
     
  (* Case B: Single Trace or Power of Trace (Tr[...] or Tr[...]^2) *)
  If[(!FreeQ[term, Tr] || !FreeQ[term, Trace]),
     {1, term},
     
  (* Case C: Pure Scalar (No Traces) *)
     {term, 1}
  ]]
];

(* --- 2. THE WORKER: Processes a Single Term --- *)
ProcessSingleTerm[term_, termIndex_, outputDir_, prefix_, doDisplay_, globalFactor_] := 
  Module[{localFactor, traceExpr, fullFactor, lispString, data, path},
    
    (* A. Split the term *)
    {localFactor, traceExpr} = IsolateFactor[term];
    
    (* B. Combine with Global Factor *)
    fullFactor = globalFactor * localFactor;
    
    (* C. Run Lisp Engine *)
    lispString = NLSMLisp`ToLispString[traceExpr];
    data = NLSMLisp`EvaluateModel[lispString];
    
    (* D. Construct File Path (Relative if outputDir is relative) *)
    path = FileNameJoin[{outputDir, prefix <> "_" <> ToString[termIndex] <> ".pdf"}];
    
    (* E. Save Report (Always) *)
    NLSMVisuals`saveReport[data, path, fullFactor];
    
    (* F. Return Result based on Display Mode *)
    If[doDisplay,
       NLSMVisuals`generateReport[data, fullFactor], (* Return Visual Grid *)
       path (* Return File Path String *)
    ]
];

(* --- 3. THE MANAGER: Batch Processor --- *)
GenerateDiagrams[expr_, outputDir_, prefix_: "Diagram", OptionsPattern[{"Display" -> True, "Factor" -> 1}]] := 
  Module[{expandedExpr, termList, displayMode, globalFactor, results},
    
    (* Setup *)
    displayMode = OptionValue["Display"];
    globalFactor = OptionValue["Factor"];
    
    (* This creates the directory relative to current Directory[] *)
    If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
    Print["Processing terms into: ", outputDir];

    (* 1. Expand and Flatten Sums *)
    expandedExpr = Expand[expr];
    termList = If[Head[expandedExpr] === Plus, List @@ expandedExpr, {expandedExpr}];
    
    (* 2. Map Over Terms with Index *)
    results = MapIndexed[
       ProcessSingleTerm[#1, #2[[1]], outputDir, prefix, displayMode, globalFactor] &,
       termList
    ];
    
    (* 3. Output *)
    If[displayMode, 
       Column[results, Spacings -> 3, Alignment -> Center], 
       results (* Returns list of paths *)
    ]
];

End[];
EndPackage[];
