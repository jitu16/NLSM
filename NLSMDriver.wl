(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMDriver`", {"NLSMAlgebra`", "NLSMLisp`", "NLSMVisuals`"}];

(* --- EXPORTS --- *)
GenerateDiagrams::usage = "GenerateDiagrams[expr, outputDir, prefix, options...] processes a tensor expression term-by-term, generating diagrams and saving reports.";

Begin["`Private`"];

(* --- HELPER: ROBUST NAME MATCHING --- *)
safeName[s_Symbol] := SymbolName[s];
safeName[other_] := "";

(* --- 1. THE SPLITTER: Isolates Scalar Factors --- *)
IsolateFactor[term_] := Module[{factors, traces, listContent},
  If[Head[term] === Times,
     listContent = List @@ term;
     factors = Select[listContent, FreeQ[#, Tr] && FreeQ[#, Trace] &];
     traces  = Select[listContent, (!FreeQ[#, Tr] || !FreeQ[#, Trace]) &];
     {Times @@ factors, Times @@ traces},
  If[(!FreeQ[term, Tr] || !FreeQ[term, Trace]),
     {1, term},
     {term, 1}
  ]]
];

(* --- 2. THE FORMATTER: Prepares Math for Display --- *)
cleanSkeletonForTitle[expr_] := expr //. {
   tag_[_, content_] /; safeName[tag] == "Tagged" :> content, 
   d_[x_] /; safeName[d] == "Der" :> Subscript["\[Del]", x], 
   w_[pos_, extra_] /; safeName[w] == "W" :> Dot[extra, w[pos]]
};

(* --- 3. THE WORKER: Processes a Single Term --- *)
ProcessSingleTerm[term_, termIndex_, outputDir_, prefix_, doDisplay_, globalFactor_] := 
  Module[{localFactor, traceExpr, fullFactor, lispString, rawResultList, enhancedList, path},
    
    (* A. Split the term *)
    {localFactor, traceExpr} = IsolateFactor[term];
    fullFactor = globalFactor * localFactor;
    
    (* B. Run Lisp Engine -> Returns a LIST of Associations *)
    lispString = NLSMLisp`ToLispString[traceExpr];
    rawResultList = NLSMLisp`EvaluateModel[lispString];
    
    (* SANITY CHECK: Ensure we have a List. If Lisp failed, wrap in empty list. *)
    If[!ListQ[rawResultList], rawResultList = {}];

    (* C. Enhance the Data (Map over the list, preserving structure!) *)
    enhancedList = Map[
        Function[rawData,
            Module[{title, updatedDiagrams},
                
                (* 1. Pre-Calculate Title *)
                title = TraditionalForm[fullFactor * cleanSkeletonForTitle[rawData["Skeleton"]]];
                
                (* 2. Pre-Calculate Values *)
                If[KeyExistsQ[rawData, "Diagrams"] && ListQ[rawData["Diagrams"]],
                   updatedDiagrams = Map[
                       Function[d, 
                           Append[d, "DisplayValue" -> NLSMAlgebra`FinalSimplify[fullFactor * d["Value"]]]
                       ], 
                       rawData["Diagrams"]
                   ];
                   ,
                   updatedDiagrams = {};
                ];
                
                (* 3. Return Enhanced Association *)
                Append[rawData, {
                    "Diagrams" -> updatedDiagrams, 
                    "DisplayTitle" -> title
                }]
            ]
        ],
        rawResultList
    ];
    
    (* D. Save & Display (Pass the LIST to Visuals) *)
    path = FileNameJoin[{outputDir, prefix <> "_" <> ToString[termIndex] <> ".pdf"}];
    
    (* Visuals package knows how to handle a List! *)
    NLSMVisuals`saveReport[enhancedList, path];
    
    If[doDisplay,
       NLSMVisuals`generateReport[enhancedList], 
       path
    ]
];

(* --- 4. THE MANAGER: Batch Processor --- *)
GenerateDiagrams[expr_, outputDir_, prefix_: "Diagram", OptionsPattern[{"Display" -> True, "Factor" -> 1}]] := 
  Module[{expandedExpr, termList, displayMode, globalFactor, results},
    
    displayMode = OptionValue["Display"];
    globalFactor = OptionValue["Factor"];
    
    If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
    Print["Processing terms into: ", outputDir];

    expandedExpr = Expand[expr];
    termList = If[Head[expandedExpr] === Plus, List @@ expandedExpr, {expandedExpr}];
    
    results = MapIndexed[
       ProcessSingleTerm[#1, #2[[1]], outputDir, prefix, displayMode, globalFactor] &,
       termList
    ];
    
    If[displayMode, 
       Column[results, Spacings -> 3, Alignment -> Center], 
       results 
    ]
];

End[];
EndPackage[];
