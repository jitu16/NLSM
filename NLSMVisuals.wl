(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMVisuals`", {"NLSMAlgebra`"}];

(* Public Exported Symbols *)
generateReport::usage = "generateReport[data, factor:1] returns a Grid showing the skeleton equation and all diagrams with their values.";
saveReport::usage = "saveReport[data, path, factor:1] exports the report Grid to a PDF file.";
saveDiagram::usage = "saveDiagram[data, id, path] exports a specific diagram to a PDF/EPS file.";
saveAllDiagrams::usage = "saveAllDiagrams[data, pathTemplate] exports all diagrams using the template for naming.";

Begin["`Private`"];

(* --- CONFIGURATION --- *)
$HubBaseRadius = 0.01;
$HubGrowthFactor = 0.25;
$RibbonWidthOuter = 0.08;
$RibbonWidthInner = 0.06;
$HashWidth = 0.55;

(* --- HELPER: ROBUST NAME MATCHING --- *)
(* This ensures we match 'W' or 'Phi' regardless of Context *)
safeName[s_Symbol] := SymbolName[s];
safeName[other_] := "";

(* --- DRAWING PRIMITIVES --- *)
getNormal[tangent_] := {-tangent[[2]], tangent[[1]]}/Norm[tangent];

drawRibbonCurve[p1_, c1_, c2_, p2_] := {
   Black, Thickness[$RibbonWidthOuter], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}], 
   White, Thickness[$RibbonWidthInner], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}]
};

drawRedRibbonCurve[p1_, c1_, c2_, p2_] := {
   Red, Thickness[$RibbonWidthOuter], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}], 
   White, Thickness[$RibbonWidthInner], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}]
};

drawTripleHash[mid_, tangent_] := Module[{norm, spacing}, 
   norm = getNormal[tangent]; spacing = 0.09;
   {Black, Thickness[0.012], CapForm["Butt"], 
    Line[{mid - norm*$HashWidth/2, mid + norm*$HashWidth/2}], 
    Line[{mid - norm*$HashWidth/2 - tangent*spacing, mid + norm*$HashWidth/2 - tangent*spacing}], 
    Line[{mid - norm*$HashWidth/2 + tangent*spacing, mid + norm*$HashWidth/2 + tangent*spacing}]}
];

drawPhi[pos_, normal_] := Module[{splitLen = 0.7, spread = 0.6, tip1, tip2, c1, root}, 
   root = pos - normal*0.05;
   tip1 = pos + normal*splitLen + {-normal[[2]], normal[[1]]}*spread;
   tip2 = pos + normal*splitLen - {-normal[[2]], normal[[1]]}*spread;
   c1 = pos + normal*(splitLen*0.3);
   {drawRedRibbonCurve[root, c1, tip1, tip1], drawRedRibbonCurve[root, c1, tip2, tip2]}
];

drawHub[center_, radius_] := {FaceForm[White], EdgeForm[{Black, Thickness[0.01]}], Disk[center, radius]};


(* --- PARSING & TOPOLOGY (ROBUST) --- *)

getPortCoords[center_, n_, i_, radius_] := Module[{angle}, angle = 135 Degree - (360 Degree*(i - 1))/n; center + radius*{Cos[angle], Sin[angle]}];
getPortNormal[n_, i_] := Module[{angle}, angle = 135 Degree - (360 Degree*(i - 1))/n; {Cos[angle], Sin[angle]}];

(* Robust extraction of Position from W[Pos, ...] *)
extractPos[arg_] := Module[{found}, 
   If[MatchQ[arg, _Symbol | _Integer | _Plus | _Subtract], Return[arg]];
   
   (* Match W[...] by name *)
   found = FirstCase[arg, w_[p_, ___] /; safeName[w] == "W" :> p, None, Infinity];
   
   If[found =!= None, Return[found]];
   arg
];

(* Robust Skeleton Parsing *)
parseSkeleton[expr_] := Module[{wOps, phiOps, allOps}, 
   
   (* Match Tagged[id, W[pos]] by names *)
   wOps = Cases[expr, 
     tag_[id_, w_[pos_, ___]] /; (safeName[tag] == "Tagged" && safeName[w] == "W") 
     :> <|"Type" -> "W", "ID" -> id, "Pos" -> pos|>, Infinity];
   
   (* Match Phi[arg] or CapitalPhi[arg] by names *)
   phiOps = Cases[expr, 
     p_[arg_] /; MemberQ[{"Phi", "\[Phi]", "CapitalPhi", "\[CapitalPhi]"}, safeName[p]] 
     :> <|"Type" -> "Phi", "ID" -> "Phi_" <> ToString[Unique[]], "Pos" -> extractPos[arg]|>, Infinity];
   
   allOps = Join[wOps, phiOps];
   GroupBy[allOps, #Pos &]
];

drawSmartTube[p1_, n1_, p2_, n2_] := Module[{dist, c1, c2, mid, tangent, curve, sumN, detourVector}, 
   dist = EuclideanDistance[p1, p2];
   c1 = p1 + n1*2.5; c2 = p2 + n2*2.5; sumN = n1 + n2;
   If[Norm[sumN] < 0.5, 
    detourVector = Normalize[{-(p2[[2]] - p1[[2]]), (p2[[1]] - p1[[1]])}];
    If[detourVector . n1 < 0, detourVector = -detourVector];
    c1 = c1 + detourVector*3.0; c2 = c2 + detourVector*3.0;, 
    If[dist > 2.0, c1 = p1 + n1*(1.5 + dist*0.8); c2 = p2 + n2*(1.5 + dist*0.8);];
   ];
   curve = BezierFunction[{p1, c1, c2, p2}];
   mid = curve[0.5]; tangent = Normalize[curve'[0.5]];
   {drawRibbonCurve[p1, c1, c2, p2], drawTripleHash[mid, tangent]}
];

drawDiagram[skeleton_, topology_] := Module[{blocksData, distinctPos, posMap, layerTubes, layerHubs, layerHashes, globalPortMap, blockCenter, ops, nPorts, radius, coord, norm, p1, n1, p2, n2, id1, id2, tubes, hash}, 
   blocksData = parseSkeleton[skeleton];
   distinctPos = Keys[blocksData];
   posMap = AssociationThread[distinctPos -> Table[{i*6.0, 0}, {i, 0, Length[distinctPos] - 1}]];
   layerTubes = {}; layerHubs = {}; layerHashes = {}; globalPortMap = <||>;
   
   Do[
    blockCenter = posMap[pos]; ops = blocksData[pos]; nPorts = Length[ops];
    radius = $HubBaseRadius + (nPorts*$HubGrowthFactor);
    AppendTo[layerHubs, drawHub[blockCenter, radius]];
    Do[
     coord = getPortCoords[blockCenter, nPorts, i, radius];
     norm = getPortNormal[nPorts, i];
     op = ops[[i]];
     If[op["Type"] == "Phi", AppendTo[layerTubes, drawPhi[coord, norm]]];
     If[op["Type"] == "W", globalPortMap[op["ID"]] = {coord, norm}];, 
     {i, 1, nPorts}
    ];, {pos, distinctPos}
   ];
   
   Do[
    id1 = pair[[1]]; id2 = pair[[2]];
    If[KeyExistsQ[globalPortMap, id1] && KeyExistsQ[globalPortMap, id2], 
     p1 = globalPortMap[id1][[1]]; n1 = globalPortMap[id1][[2]];
     p2 = globalPortMap[id2][[1]]; n2 = globalPortMap[id2][[2]];
     {tubes, hash} = drawSmartTube[p1, n1, p2, n2];
     AppendTo[layerTubes, tubes]; AppendTo[layerHashes, hash];
    ];, {pair, topology}
   ];
   Graphics[{layerHubs, layerTubes, layerHashes}, ImageSize -> {Automatic, 200}, Background -> None, PlotRange -> All]
];


(* --- REPORTING (ROBUST) --- *)

cleanSkeletonDisplay[expr_] := expr //. {
   (* Match Tagged by name *)
   tag_[_, content_] /; safeName[tag] == "Tagged" :> content, 
   
   (* Match Der by name *)
   d_[x_] /; safeName[d] == "Der" :> Subscript["\[Del]", x], 
   
   (* Match W by name *)
   w_[pos_, extra_] /; safeName[w] == "W" :> Dot[extra, w[pos]]
};

generateReport[data_, factor_: 1] := Module[{skeleton, cleanedSkeleton, diagrams, gridItems, headerRow, titleRow}, 
   skeleton = data["Skeleton"];
   diagrams = data["Diagrams"];
   cleanedSkeleton = cleanSkeletonDisplay[skeleton];
   
   gridItems = Table[{
      Style[diag["ID"], 20, Bold], 
      drawDiagram[skeleton, diag["Topology"]], 
      Pane[FinalSimplify[factor*diag["Value"]], {350, Automatic}, ImageSizeAction -> "Scrollable"]
   }, {diag, diagrams}];
   
   headerRow = {Style["ID", Bold], Style["Diagram", Bold], Style["Value", Bold]};
   titleRow = {
      Item[Style[TraditionalForm[cleanedSkeleton], 20, Black], Background -> LightBlue, Frame -> True, Alignment -> Center], 
      SpanFromLeft, SpanFromLeft
   };
   
   Grid[Join[{titleRow, headerRow}, gridItems], Frame -> All, Spacings -> {1, 1}, ItemStyle -> {{Left, Center, Left}}, Alignment -> {Left, Center}]
];

saveAllDiagrams[data_, pathTemplate_] := Module[{skeleton, diagrams, graphic, dir, base, ext, fullPath, createdFiles}, 
   skeleton = data["Skeleton"]; diagrams = data["Diagrams"];
   dir = DirectoryName[pathTemplate]; base = FileBaseName[pathTemplate]; ext = FileExtension[pathTemplate];
   If[ext == "", ext = "pdf"]; If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   createdFiles = {};
   Do[
    graphic = drawDiagram[skeleton, diag["Topology"]];
    fullPath = FileNameJoin[{dir, base <> "_" <> ToString[diag["ID"]] <> "." <> ext}];
    Export[fullPath, graphic]; AppendTo[createdFiles, fullPath];, {diag, diagrams}
   ];
   Print["Saved ", Length[createdFiles], " diagrams to: ", dir]; createdFiles
];

saveReport[data_, path_, factor_: 1] := Module[{dir, reportGrid}, 
   dir = DirectoryName[path]; If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   reportGrid = generateReport[data, factor];
   Export[path, reportGrid]; Print["Report saved successfully to: ", path]; path
];

saveDiagram[data_, id_, path_] := Module[{skeleton, diagrams, targetDiag, graphic, dir}, 
   skeleton = data["Skeleton"]; diagrams = data["Diagrams"];
   targetDiag = FirstCase[diagrams, d_ /; d["ID"] == id, None];
   If[targetDiag === None, Print["Error: Diagram ID ", id, " not found."]; Return[$Failed]];
   dir = DirectoryName[path]; If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   graphic = drawDiagram[skeleton, targetDiag["Topology"]];
   Export[path, graphic]; Print["Saved Diagram ", id, " to: ", path]; path
];

End[];
EndPackage[];
