(* ::Package:: *)

(* ::Package:: *)

BeginPackage["NLSMVisuals`"]; (* No Algebra dependency *)

(* Public Exported Symbols *)
generateReport::usage = "generateReport[data] returns a Grid showing the skeleton equation and all diagrams with their pre-calculated values.";
saveReport::usage = "saveReport[data, path] exports the report Grid to a PDF file.";
saveDiagram::usage = "saveDiagram[data, id, path] exports a specific diagram to a PDF/EPS file.";
saveAllDiagrams::usage = "saveAllDiagrams[data, pathTemplate] exports all diagrams using the template for naming.";
drawDiagram::usage = "drawDiagram[skeleton, topology] returns the Graphics object.";

Begin["`Private`"];

(* --- CONFIGURATION --- *)
$HubBaseRadius = 0.01;
$HubGrowthFactor = 0.25;
$RibbonWidthOuter = 0.08;
$RibbonWidthInner = 0.06;
$HashWidth = 0.55;

(* --- HELPER: ROBUST NAME MATCHING --- *)
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

drawGrayRibbonCurve[p1_, c1_, c2_, p2_] := {
   Gray, Thickness[$RibbonWidthOuter], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}], 
   LightGray, Thickness[$RibbonWidthInner], CapForm["Round"], BezierCurve[{p1, c1, c2, p2}]
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

drawSmallPhi[pos_, normal_] := Module[{len = 0.7, tip, c1, curve, t, pt, tan, norm, spacing, lines},
   tip = pos + normal*len;
   c1 = pos + normal*(len*0.5);
   curve = BezierFunction[{pos, c1, tip}];
   lines = Table[
     pt = curve[t];
     tan = Normalize[curve'[t]];
     norm = {-tan[[2]], tan[[1]]}; 
     {Black, Thickness[0.005], Line[{pt - norm*0.025, pt + norm*0.025}]}
     , {t, 0.1, 0.9, 0.1}];
   {drawGrayRibbonCurve[pos, c1, tip, tip], lines}
];

drawWavyLine[p1_, p2_, amplitude_: 0.15, waveDensity_: 2.0] := Module[{dist, tan, norm, points, t, numWaves, numSegments},
   dist = EuclideanDistance[p1, p2];
   If[dist < 0.2, Return[{Gray, Thickness[0.008], CapForm["Round"], Line[{p1, p2}]}]];
   tan = Normalize[p2 - p1];
   norm = {-tan[[2]], tan[[1]]}; 
   numWaves = Max[1, Round[dist * waveDensity]];
   numSegments = numWaves * 20; 
   points = Table[
      t = i / numSegments; 
      (p1 + t*(p2 - p1)) + (amplitude * Sin[t * numWaves * 2 * Pi] * norm),
      {i, 0, numSegments}
   ];
   {Gray, Thickness[0.008], CapForm["Round"], Line[points]}
];

drawHub[center_, radius_] := {FaceForm[White], EdgeForm[{Black, Thickness[0.01]}], Disk[center, radius]};

(* --- PARSING & TOPOLOGY --- *)

getPortCoords[center_, n_, i_, radius_] := Module[{angle}, angle = 135 Degree - (360 Degree*(i - 1))/n; center + radius*{Cos[angle], Sin[angle]}];
getPortNormal[n_, i_] := Module[{angle}, angle = 135 Degree - (360 Degree*(i - 1))/n; {Cos[angle], Sin[angle]}];

extractPos[arg_] := Module[{found}, 
   If[MatchQ[arg, _Symbol | _Integer | _Plus | _Subtract], Return[arg]];
   found = FirstCase[arg, w_[p_, ___] /; safeName[w] == "W" :> p, None, Infinity];
   If[found =!= None, Return[found]];
   arg
];

parseSkeleton[expr_] := Module[{wOps, capPhiOps, smallPhiOps, allOps}, 
   wOps = Cases[expr, tag_[id_, w_[pos_, ___]] /; (safeName[tag] == "Tagged" && safeName[w] == "W") :> 
     <|"Type" -> "W", "ID" -> id, "Pos" -> pos|>, Infinity];
   capPhiOps = Cases[expr, tag_[id_, p_[arg_]] /; (safeName[tag] == "Tagged" && MemberQ[{"Phi", "CapitalPhi", "\[CapitalPhi]"}, safeName[p]]) :> 
     <|"Type" -> "CapitalPhi", "ID" -> id, "Pos" -> extractPos[arg]|>, Infinity];
   smallPhiOps = Cases[expr, tag_[id_, p_[arg_]] /; (safeName[tag] == "Tagged" && MemberQ[{"\[Phi]", "VarPhi"}, safeName[p]]) :> 
     <|"Type" -> "SmallPhi", "ID" -> id, "Pos" -> extractPos[arg]|>, Infinity];
   If[capPhiOps === {}, 
      capPhiOps = Cases[expr, p_[arg_] /; MemberQ[{"Phi", "CapitalPhi", "\[CapitalPhi]"}, safeName[p]] :> 
        <|"Type" -> "CapitalPhi", "ID" -> "Phi_" <> ToString[Unique[]], "Pos" -> extractPos[arg]|>, Infinity];
   ];
   If[smallPhiOps === {}, 
      smallPhiOps = Cases[expr, p_[arg_] /; MemberQ[{"\[Phi]", "VarPhi"}, safeName[p]] :> 
        <|"Type" -> "SmallPhi", "ID" -> "phi_" <> ToString[Unique[]], "Pos" -> extractPos[arg]|>, Infinity];
   ];
   allOps = Join[wOps, capPhiOps, smallPhiOps];
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

drawDiagram[skeleton_, topology_] := Module[{blocksData, distinctPos, posMap, layerTubes, layerHubs, layerHashes, globalPortMap, blockCenter, ops, nPorts, radius, coord, norm, p1, p2, id1, id2, type1, type2, tubes, hash, connectedIDs, idTypeMap, n1, n2}, 
   blocksData = parseSkeleton[skeleton];
   distinctPos = Keys[blocksData];
   posMap = AssociationThread[distinctPos -> Table[{i*6.0, 0}, {i, 0, Length[distinctPos] - 1}]];
   
   connectedIDs = Flatten[topology];
   idTypeMap = <||>;
   Do[Do[idTypeMap[op["ID"]] = op["Type"], {op, blocksData[pos]}], {pos, distinctPos}];

   layerTubes = {}; layerHubs = {}; layerHashes = {}; globalPortMap = <||>;
   
   Do[
    blockCenter = posMap[pos]; ops = blocksData[pos]; nPorts = Length[ops];
    radius = $HubBaseRadius + (nPorts*$HubGrowthFactor);
    AppendTo[layerHubs, drawHub[blockCenter, radius]];
    Do[
     coord = getPortCoords[blockCenter, nPorts, i, radius];
     norm = getPortNormal[nPorts, i];
     op = ops[[i]];
     globalPortMap[op["ID"]] = {coord, norm};
     If[op["Type"] == "CapitalPhi", AppendTo[layerTubes, drawPhi[coord, norm]]];
     If[op["Type"] == "SmallPhi", 
        If[!MemberQ[connectedIDs, op["ID"]], AppendTo[layerTubes, drawSmallPhi[coord, norm]]];
     ];
     ;, {i, 1, nPorts}
    ];, {pos, distinctPos}
   ];
   
   Do[
    id1 = pair[[1]]; id2 = pair[[2]];
    If[KeyExistsQ[globalPortMap, id1] && KeyExistsQ[globalPortMap, id2], 
     p1 = globalPortMap[id1][[1]]; n1 = globalPortMap[id1][[2]];
     p2 = globalPortMap[id2][[1]]; n2 = globalPortMap[id2][[2]];
     type1 = idTypeMap[id1]; type2 = idTypeMap[id2];

     If[type1 == "W" && type2 == "W",
         {tubes, hash} = drawSmartTube[p1, n1, p2, n2];
         AppendTo[layerTubes, tubes]; AppendTo[layerHashes, hash];
     ];
     If[type1 == "SmallPhi" && type2 == "SmallPhi",
         AppendTo[layerHashes, drawWavyLine[p1, p2]];
     ];
    ];, {pair, topology}
   ];
   Graphics[{layerHubs, layerTubes, layerHashes}, ImageSize -> {Automatic, 200}, Background -> None, PlotRange -> All]
];

(* --- REPORTING (SIMPLIFIED) --- *)

(* 1. MANAGER *)
generateReport[data_List] := Column[
   Table[
    Column[{
      Style["Term " <> ToString[i], "Section", Gray],
      generateReport[data[[i]]], 
      Style["____________________", Gray]
    }, Alignment -> Center]
    , {i, Length[data]}],
   Spacings -> 2, Alignment -> Center
];

(* 2. WORKER *)
generateReport[data_Association] := Module[{skeleton, diagrams, gridItems, headerRow, titleRow}, 
   skeleton = data["Skeleton"];
   diagrams = data["Diagrams"];
   
   (* Use pre-calculated Display values from Driver *)
   gridItems = Table[{
      Style[diag["ID"], 20, Bold], 
      drawDiagram[skeleton, diag["Topology"]], 
      Pane[diag["DisplayValue"], {350, Automatic}, ImageSizeAction -> "Scrollable"]
   }, {diag, diagrams}];
   
   headerRow = {Style["ID", Bold], Style["Diagram", Bold], Style["Value", Bold]};
   titleRow = {
      Item[Style[data["DisplayTitle"], 20, Black], Background -> LightBlue, Frame -> True, Alignment -> Center], 
      SpanFromLeft, SpanFromLeft
   };
   
   Grid[Join[{titleRow, headerRow}, gridItems], Frame -> All, Spacings -> {1, 1}, ItemStyle -> {{Left, Center, Left}}, Alignment -> {Left, Center}]
];

saveAllDiagrams[data_List, pathTemplate_] := Module[{ext, base, dir, termPath, allFiles = {}},
   dir = DirectoryName[pathTemplate];
   base = FileBaseName[pathTemplate];
   ext = FileExtension[pathTemplate];
   If[ext == "", ext = "pdf"];
   Do[
    termPath = FileNameJoin[{dir, base <> "_T" <> ToString[i] <> "." <> ext}];
    Join[allFiles, saveAllDiagrams[data[[i]], termPath]];
    , {i, Length[data]}];
   Flatten[allFiles]
];

saveAllDiagrams[data_Association, pathTemplate_] := Module[{skeleton, diagrams, graphic, dir, base, ext, fullPath, createdFiles}, 
   skeleton = data["Skeleton"];
   diagrams = data["Diagrams"];
   dir = DirectoryName[pathTemplate];
   base = FileBaseName[pathTemplate];
   ext = FileExtension[pathTemplate];
   If[ext == "", ext = "pdf"];
   If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   createdFiles = {};
   Do[
    graphic = drawDiagram[skeleton, diag["Topology"]];
    fullPath = FileNameJoin[{dir, base <> "_" <> ToString[diag["ID"]] <> "." <> ext}];
    Export[fullPath, graphic];
    AppendTo[createdFiles, fullPath];, 
    {diag, diagrams}
   ];
   Print["Saved ", Length[createdFiles], " diagrams to: ", dir];
   createdFiles
];

saveReport[data_List, path_] := Module[{reportColumn, dir},
   dir = DirectoryName[path];
   If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   reportColumn = generateReport[data];
   Export[path, reportColumn];
   Print["Full Report saved to: ", path];
   path
];

saveReport[data_Association, path_] := Module[{dir, reportGrid}, 
   dir = DirectoryName[path];
   If[dir =!= "" && ! DirectoryQ[dir], CreateDirectory[dir]];
   reportGrid = generateReport[data];
   Export[path, reportGrid];
   Print["Report saved successfully to: ", path];
   path
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
