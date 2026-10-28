(* ::Package:: *)
(* ------------------------------------------------------------------ *)
(*  GTPIntersections.wl                                               *)
(*  Triple intersection numbers of divisors in the Calabi-Yau three-  *)
(*  fold associated with a triangulated generalized toric polygon.    *)
(* ------------------------------------------------------------------ *)

BeginPackage["GTPIntersections`"];

ClearAll["GTPIntersections`*"];

(* ---------- combinatorics of the Symington polygon ---------------- *)
monodromyMatrix::usage      = "monodromyMatrix[{p,q}] is the SL(2,Z) monodromy of a (p,q) 7-brane.";
invariantDirection::usage   = "invariantDirection[{p,q}] is the primitive monodromy-invariant direction.";
cutImage::usage             = "cutImage[{charge,cut}] is the image of the cut polyline under the monodromy fold.";
gtpEdgeIdentifications::usage = "gtpEdgeIdentifications[brane] pairs the segments of the cut with their images.";
gtpVertexPairs::usage       = "gtpVertexPairs[brane] pairs the vertices of the cut with their images.";
gtpQuotient::usage          = "gtpQuotient[edges,branes] assembles all combinatorial data of a triangulated GTP.";
edgeKey::usage              = "edgeKey[e] is the orientation independent key of an edge.";
gtpKey::usage               = "gtpKey[x] wraps x so that it can be used as an Association key.";
triangleList::usage         = "triangleList[edges] lists the triangles of the triangulation.";
thirdVertex::usage          = "thirdVertex[t,e] is the vertex of the triangle t not on the edge e.";
trianglesContainingEdge::usage = "trianglesContainingEdge[triangles,e] selects the triangles having e as a side.";
internalEdgeSelfIntersections::usage = "internalEdgeSelfIntersections[data] gives the normal bundle degrees of the compact double curves.";

(* ---------- the algorithm ----------------------------------------- *)
gtpELines::usage            = "gtpELines[data] traces the E-lines of the exuberant divisors.";
gtpIndexSets::usage         = "gtpIndexSets[data] builds the divisor and curve index sets I and A.";
gtpToricAlpha::usage        = "gtpToricAlpha[ix] is Step 1 of the algorithm (Section 6.1.2 of the draft): the coefficients alpha^a_ij of D_i.D_j for central-fiber divisors.";
gtpToricBeta::usage         = "gtpToricBeta[ix] is Step 1 of the algorithm: the coefficients beta_ia = D_i.C_a for central-fiber divisors and double curves (Undefined outside the admissibility domain B_0).";
gtpAnsatz::usage            = "gtpAnsatz[ix] is Steps 2 and 3 of the algorithm: the admissibility domain B and the coefficients alpha and beta including the exuberant unknowns.";
gtpCompatibilityEquations::usage = "gtpCompatibilityEquations[ans] is the part of the compatibility system of Step 4 whose first divisor is compact; by default gtpIntersections imposes all admissible presentations T_{i;jk} = T_{j;ik} (see its option Presentations).";
gtpIntersections::usage     = "gtpIntersections[edges,branes] runs the whole algorithm and returns an Association with the triple intersection numbers (\"Intersections\"), the coefficients \"Alpha\" and \"Beta\", the equations, their solution and the \"Status\".  Divisors are D_1, ..., D_n in the order in which the vertices first appear in the edge list, then the exuberant divisors E_{n+1}, ...; curves are numbered as in the draft: compact double curves (in the order of the edges), exuberant curves, non-compact double curves.  Options: \"CurveNames\" -> \"Ends\" names the double curves C_{i,j} instead; \"BypassSide\" -> \"Left\" or \"Right\" fixes the side on which the deformed E-line passes the compact vertices it meets (by default the left, unless the tip of another 7-brane is there); \"Presentations\" -> \"CompactFirst\" imposes only the compatibility relations whose first divisor is compact.  When the compatibility equations leave some coefficients free without affecting D_i.D_j.D_k, the status is Solved, \"CoefficientsDetermined\" is False, \"FreeParameters\" lists them, \"SymbolicAlpha\" and \"SymbolicBeta\" keep them, and \"Alpha\", \"Beta\" use the choice \"Gauge\".";
gtpIntersectionTables::usage = "gtpIntersectionTables[res] formats the matrices D_i.D_j.D_k for compact D_i.";
gtpReport::usage            = "gtpReport[res] prints a summary of a run of the algorithm, ending with the intersection matrices (omitted with \"IntersectionTables\" -> False).";
gtpPlot::usage              = "gtpPlot[data] draws the GTP with divisor / curve labels, the normal bundle degrees, the E-lines and the shaded generalized faces.  Accepts either the output of gtpQuotient or of gtpIntersections.";
gtpPictureData::usage       = "gtpPictureData[data] collects the geometric data that gtpPlot and gtpTikZ draw.";
gtpTikZ::usage              = "gtpTikZ[data] returns, as a String, the tikzpicture reproducing gtpPlot[data].";
gtpDivisorSummary::usage    = "gtpDivisorSummary[res] tabulates the divisors, their charge and D^3.";
gtpCurveSummary::usage      = "gtpCurveSummary[res] tabulates the curves C_a.";
gtpBetaTable::usage         = "gtpBetaTable[res] displays the matrix beta_{ia} = D_i.C_a.";
gtpAlphaTable::usage        = "gtpAlphaTable[res] displays D_i.D_j = sum_a alpha^a_ij C_a.";

gtpConsistencyChecks::usage = "gtpConsistencyChecks[res] runs geometric consistency tests on the output of gtpIntersections.";

(* ---------- input ------------------------------------------------ *)
gtpDraw::usage = "gtpDraw[] asks for the size n of the grid and opens an n x n grid on which the edges of the triangulation and the 7-brane cuts are drawn with the mouse; clicking Done returns {edges, branes}, the two arguments of gtpIntersections.  A 7-brane is drawn as the two sides of its wedge, starting from the tip, and its charge is deduced from them (see gtpBraneFromSides).  gtpDraw[n] skips the question, gtpDraw[{{xmin,xmax},{ymin,ymax}}] uses that range of lattice points, and further arguments give the edges and branes to start from, e.g. {gtp, branes} = gtpDraw[5, gtp, branes] to edit them (Cancel then returns them unchanged).";
gtpBraneFromSides::usage = "gtpBraneFromSides[side1, side2] gives the 7-brane {charge, cut} whose monodromy maps one of the two polylines side1, side2 (lists of lattice points starting at the common tip of the wedge) onto the other; the cut is whichever side is mapped onto the other.  The charge is {p,q}, or {p,q,k} for k coincident branes.  Returns $Failed, with the reason, if no 7-brane does.";
gtpRun::usage = "gtpRun[] does the whole computation in one go: it opens the editor of gtpDraw, runs gtpIntersections on the GTP drawn there and prints the input, the picture of the GTP, the report of gtpReport, the triple intersection numbers and the TikZ code of the picture, with buttons that copy it.  It returns the result of gtpIntersections, with the drawn GTP under \"Input\", so write  res = gtpRun[];  and  res = gtpRun[n, Sequence @@ res[\"Input\"]];  to edit it.  gtpRun[n] and gtpRun[{{xmin,xmax},{ymin,ymax}}] open the editor on that grid, gtpRun[n, edges, branes] opens it on a GTP to be edited, and gtpRun[edges, branes] skips the editor.  Options: \"Report\" -> False prints only the picture, the intersection numbers and the TikZ code; \"PlotOptions\" and \"TikZOptions\" are lists of options for gtpPlot and gtpTikZ; \"TikZFile\" -> \"name.tex\" also writes the TikZ code to that file (next to the notebook unless the path is absolute); \"Presentations\", \"CurveNames\" and \"BypassSide\" are passed to gtpIntersections.";

Begin["`Private`"];

gtpQuotient::branes = "The second argument must be the list of 7-branes, each \
of the form {{p,q},{{x1,y1},{x2,y2},...}} (charge, then monodromy cut), or {} \
if there is none.  Received `1`.";
gtpQuotient::onebrane = "The second argument, `1`, is a single 7-brane; the \
argument must be the list of the 7-branes, so write `2` instead.";
gtpQuotient::wedge = "The generalized face swept out by the 7-brane `1` contains \
the triangle(s) `2`.  The face of a 7-brane must not be triangulated: remove the \
edges inside it.";

(* ================================================================== *)
(*  1.  Cuts, monodromy, quotient combinatorics                        *)
(* ================================================================== *)

(*  monodromyMatrix[{p,q}] is the monodromy of one (p,q) 7-brane.
    monodromyMatrix[{p,q,k}] is the monodromy of k coincident (p,q)
    7-branes, M^k = 1 + k N with N nilpotent.  A boundary segment of
    the GTP of integral length l whose 7-brane is pulled in to a vertex
    from which another such segment is also pulled in needs the second
    form: the two wedges then share a side, which is not an edge of the
    triangulation, and only the total monodromy is visible.
    Conventions relative to the draft: the charge {p,q} is the [p,q] of
    the draft, whose M_[p,q] = {{1-pq, p^2}, {-q^2, 1+pq}} acts on
    charges and M_[p,q]^T on directions in the polygon.  monodromyMatrix
    acts on directions and equals Inverse[Transpose[M_[p,q]]]: the cut is
    the side of the wedge onto which M_[p,q]^T maps the other side,
    cutImage.  E.g. for the SU(4) example of the draft, brane
    {{1,1},{{-1,4},{0,0}}}, the cut runs to v_4'' = (0,0) and cutImage
    to v_4' = (-3,3).                                                  *)
monodromyMatrix[{p_, q_}] := monodromyMatrix[{p, q, 1}];
monodromyMatrix[{p_, q_, k_}] := IdentityMatrix[2] + k {{p q, q^2}, {-p^2, -p q}};

(* M.(q,-p) = (q,-p) : the monodromy invariant line *)
invariantDirection[{p_, q_}] := With[{g = GCD[p, q]}, {q, -p}/If[g === 0, 1, g]];
invariantDirection[{p_, q_, k_}] := invariantDirection[{p, q}];

cutImage[{charge_, cut_List}] :=
  Module[{m = monodromyMatrix[charge]},
   FoldList[#1 + m . #2 &, First[cut], Differences[cut]]];

gtpEdgeIdentifications[{charge_, cut_List}] :=
  Module[{img = cutImage[{charge, cut}]},
   Transpose[{Partition[cut, 2, 1], Partition[img, 2, 1]}]];

gtpVertexPairs[{charge_, cut_List}] :=
  Module[{img = cutImage[{charge, cut}]},
   DeleteDuplicates[Select[Transpose[{cut, img}], #[[1]] =!= #[[2]] &]]];

edgeKey[e_] := Sort[e];
gtpKey[x_]  := HoldComplete[x];

(*  a triple of vertices spans a triangle of the triangulation iff its
    three sides are edges, it is non degenerate, and no further vertex
    lies in the closed triangle.  The last condition is essential: the
    outer boundary of a polygon whose interior vertices are all joined
    to it also has all three sides present (e.g. the local P^2 with an
    interior point), and must not be counted as a face.               *)
inClosedTriangle[t_List, p_] :=
 Module[{s = Sign[Det[{t[[2]] - t[[1]], t[[3]] - t[[1]]}]]},
  And @@ (# >= 0 & /@ (s*{Det[{t[[2]] - t[[1]], p - t[[1]]}],
       Det[{t[[3]] - t[[2]], p - t[[2]]}],
       Det[{t[[1]] - t[[3]], p - t[[3]]}]}))];

triangleList[edges_List] :=
 Module[{verts, edgeAssoc, candidates},
  verts = DeleteDuplicates[Flatten[edges, 1]];
  edgeAssoc = AssociationThread[gtpKey /@ edgeKey /@ edges -> True];
  candidates = Select[Subsets[verts, {3}],
    Function[t, Det[{t[[2]] - t[[1]], t[[3]] - t[[1]]}] =!= 0 &&
      And @@ Table[TrueQ[Lookup[edgeAssoc, gtpKey[edgeKey[pair]], False]],
        {pair, Subsets[t, {2}]}]]];
  Select[candidates,
   Function[t, NoneTrue[Complement[verts, t], inClosedTriangle[t, #] &]]]];

trianglesContainingEdge[triangles_List, edge_List] := Select[triangles, SubsetQ[#, edge] &];
thirdVertex[triangle_List, edge_List] := First[Complement[triangle, edge]];

makeQuotientRelation[universe_List, pairs_List] :=
 Module[{ord, g, classes, rep, classOf},
  ord = AssociationThread[gtpKey /@ universe -> Range[Length[universe]]];
  g = Graph[universe, UndirectedEdge @@@ pairs];
  classes = ConnectedComponents[g];
  rep = Association@Flatten[
     Function[class, With[{r = First[SortBy[class, Lookup[ord, gtpKey[#], Infinity] &]]},
        (gtpKey[#] -> r &) /@ class]] /@ classes, 1];
  classOf[x_] := Lookup[rep, gtpKey[x], x];
  Association["Universe" -> universe, "Pairs" -> pairs, "Classes" -> classes,
   "Representative" -> rep, "ClassOf" -> classOf,
   "Quotient" -> DeleteDuplicates[classOf /@ universe]]];

(* ------------------------------------------------------------------ *)
(*  normal bundle degrees of the double curves                         *)
(* ------------------------------------------------------------------ *)
(*  For an interior edge e = [v_i,v_j] of the triangulation, with the
    two adjacent triangles having third vertices v_k and v_l, the toric
    relation
        (C|_{D_i})^2 (v_j-v_i) + (v_k-v_i) + (v_l-v_i) = 0
    fixes (C|_{D_i})^2 = D_i.D_j.D_j.  When the edge is a cut edge the
    two adjacent triangles live in different charts and the far one has
    to be transported by the monodromy.
    selfIntersectionForRep returns {(C|_{D_{p2}})^2 , (C|_{D_{p1}})^2}
    where {p1,p2} = Ends of the representative edge, i.e. the first
    entry is beta_{p1,e} = D_{p1}.D_{p1}.D_{p2} and the second is
    beta_{p2,e}.                                                      *)

(*  The gluings of a class form a chain e_{k} ~ e_{k+1}: when two 7-branes
    sit at the same apex their wedges are adjacent and the edge they
    share is glued twice, so the chain has length > 1 and its middle
    edges carry no triangle at all.  cutStep gives the elementary
    transport across one gluing (point map on the cut, monodromy on
    vectors), transportChain composes them along the chain.            *)
cutStep[data_Association, pair_List] :=
 Module[{b = data["BraneOfCutPair"][gtpKey[pair]], img},
  img = cutImage[{braneCharge[b], braneCut[b]}];
  <|"M" -> monodromyMatrix[braneCharge[b]],
    "Map" -> AssociationThread[gtpKey /@ braneCut[b] -> img],
    "InvMap" -> AssociationThread[gtpKey /@ img -> braneCut[b]]|>];

(*  {Phi, A} : the point map and the linear map carrying the chart of
    the first edge of the chain onto the chart of the last one.        *)
transportChain[data_Association, path_List] :=
 Module[{pairs = data["CutEdgeIdPairs"], phi = (# &), lin = IdentityMatrix[2], st, fwd},
  Do[fwd = MemberQ[pairs, {path[[k]], path[[k + 1]]}];
     st = cutStep[data, If[fwd, {path[[k]], path[[k + 1]]}, {path[[k + 1]], path[[k]]}]];
     With[{f = phi, map = If[fwd, st["Map"], st["InvMap"]],
           m = If[fwd, st["M"], Inverse[st["M"]]]},
      phi = Function[x, Lookup[map, gtpKey[f[x]], Missing["NotOnCut", f[x]]]];
      lin = m . lin],
   {k, Length[path] - 1}];
  {phi, lin}];

(*  the chain of gluings joining two edges of the same class           *)
gluingPath[data_Association, from_Integer, to_Integer] :=
 Module[{eclass = data["EdgeEquivalence"]["ClassOf"], members, g, p},
  members = Select[Range[Length[data["Edges"]]], eclass[#] === eclass[from] &];
  g = Graph[members, UndirectedEdge @@@ Select[data["CutEdgeIdPairs"],
      MemberQ[members, #[[1]]] && MemberQ[members, #[[2]]] &]];
  p = Quiet@FindShortestPath[g, from, to];
  If[ListQ[p] && p =!= {}, p, None]];

selfIntersectionForRep[data_Association, repId_Integer] :=
 Module[{edges, triangles, eclass, members, nTri, ends, lEdge, rEdge, path, phi, lin,
   mono, t1, t2, p11, p12, p13, p21, p22, p23, a, b, val, vcl, vEnds},
  edges = data["Edges"][[All, "Ends"]];
  triangles = data["Triangles"];
  eclass = data["EdgeEquivalence"]["ClassOf"];
  members = Select[Range[Length[edges]], eclass[#] === repId &];
  nTri = AssociationThread[members,
     Length[trianglesContainingEdge[triangles, edges[[#]]]] & /@ members];
  ends = Select[members, nTri[#] === 1 &];
  Which[
   (* an ordinary interior edge : the two charts are the same one *)
   Length[members] === 1 && nTri[First[members]] === 2,
    lEdge = rEdge = First[members]; mono = IdentityMatrix[2];
    {t1, t2} = trianglesContainingEdge[triangles, edges[[lEdge]]];
    {p11, p12} = edges[[lEdge]]; {p21, p22} = {p11, p12},
   (* a curve glued from the two ends of a chain of cut edges *)
   Length[ends] === 2 && Total[Values[nTri]] === 2,
    {lEdge, rEdge} = ends;
    path = gluingPath[data, lEdge, rEdge];
    If[path === None, Return[Missing["NoGluingPath", <|"RepId" -> repId|>]]];
    {phi, lin} = transportChain[data, path];
    mono = Inverse[lin];
    t1 = First[trianglesContainingEdge[triangles, edges[[lEdge]]]];
    t2 = First[trianglesContainingEdge[triangles, edges[[rEdge]]]];
    {p11, p12} = edges[[lEdge]];
    {p21, p22} = phi /@ {p11, p12};
    If[Sort[{p21, p22}] =!= Sort[edges[[rEdge]]],
     Return[Missing["GluingMismatch",
       <|"RepId" -> repId, "Image" -> {p21, p22}, "Edge" -> edges[[rEdge]]|>]]],
   True,
    Return[Missing["UnexpectedIncidence", <|"RepId" -> repId, "TriangleSides" -> nTri|>]]];
  p13 = thirdVertex[t1, edges[[lEdge]]]; p23 = thirdVertex[t2, edges[[rEdge]]];
  a = p11 - p12; b = p12 - p11;
  (* val1 = (C|_{D_{p12}})^2 = D_{p11}.D_{p11}.D_{p12} = beta_{p11,a}
     val2 = (C|_{D_{p11}})^2 = beta_{p12,a}                            *)
  val = {-((p13 - p12) + mono . (p23 - p22)) . a/(a . a),
         -((p13 - p11) + mono . (p23 - p21)) . b/(b . b)};
  (* reorder so that the r-th entry belongs to the r-th end of the
     *representative* edge after the quotient (lEdge need not be the
     representative when the class comes from a cut)                   *)
  vcl = data["VertexEquivalence"]["ClassOf"];
  vEnds = vcl /@ data["Edges"][[repId, "Ends"]];
  If[vcl[p11] === vEnds[[2]] && vcl[p12] === vEnds[[1]] && vEnds[[1]] =!= vEnds[[2]],
   val = Reverse[val]];
  val];

internalEdgeSelfIntersections[data_Association] :=
  AssociationThread[data["InternalEdgeRepresentatives"],
   selfIntersectionForRep[data, #] & /@ data["InternalEdgeRepresentatives"]];

(* ================================================================== *)
(*  2.  gtpQuotient : all combinatorial data of a triangulated GTP     *)
(* ================================================================== *)
(*  edges  : list of {{x1,y1},{x2,y2}} , the 1-cells of the unimodular
             triangulation of the Symington polygon.  The interior of
             every generalized face (T-cone) is *not* triangulated: the
             face is removed and its two sides are glued by the cut.
    branes : the list of the 7-branes, each of them {charge,cut}.
             charge={p,q} is the 7-brane charge, cut is the polyline of
             lattice points along one side of the removed wedge, starting
             at its tip.  The other side is cutImage[{charge,cut}].
             Convention: the polygon lies counterclockwise from the cut.
                                                                       *)

(*  A brane is {charge, cut} , {charge, cut, u} or {charge, cut, {u1,...,uk}}
    with charge = {p,q} or {p,q,k}.  cut is the polyline along one side
    of the wedge, starting at its tip; the other side is
    cutImage[{charge,cut}].  k is the number of coincident 7-branes, so
    the number of exuberant divisors emanating from the tip; the ui are
    their E-line directions, all parallel to invariantDirection.       *)
(*  The second argument of gtpQuotient / gtpIntersections is always the
    *list* of the 7-branes, so a single brane has to be wrapped:
    {{{0,1},{{1,2},{2,0}}}} and not {{0,1},{{1,2},{2,0}}}.  An empty
    list means no 7-brane, i.e. an ordinary toric polygon.             *)
braneQ[b_] := MatchQ[b, {{_Integer, _Integer} | {_Integer, _Integer, _Integer},
     {{_Integer, _Integer} ..}, ___}] && Length[b[[2]]] >= 2;

checkBranes[bs_] := Which[
   ListQ[bs] && AllTrue[bs, braneQ], bs,
   braneQ[bs], Message[gtpQuotient::onebrane, bs, {bs}]; Abort[],
   True, Message[gtpQuotient::branes, bs]; Abort[]];

(*  The generalized face swept out by a 7-brane -- the polygon bounded by
    its cut, the image of the cut and the segment joining their far ends
    -- must not be triangulated (Conventions (i)).  braneWedgeTriangles
    lists the triangles lying in it.  For a one-segment cut the face is
    the cone on d and M.d, which are on the same side of the invariant
    line since M preserves det(u, .): it is convex.                   *)
braneWedge[b_] := Join[braneCut[b], Reverse[Rest[cutImage[{braneCharge[b], braneCut[b]}]]]];

(*  exact lattice geometry rather than the region functions: those load
    on first use, which a Dynamic of the editor must not have to do    *)
polygonArea2[w_] := Total[det2 @@@ Partition[w, 2, 1, 1]];          (* twice the signed area *)
insidePolygonQ[w_, pt_] := Module[{wn = 0, n = Length[w], a, b},     (* winding number *)
  Do[a = w[[i]]; b = w[[Mod[i, n] + 1]];
   If[a[[2]] <= pt[[2]],
    If[b[[2]] > pt[[2]] && det2[b - a, pt - a] > 0, wn++],
    If[b[[2]] <= pt[[2]] && det2[b - a, pt - a] < 0, wn--]], {i, n}];
  wn =!= 0];

braneWedgeTriangles[triangles_List, b_] := With[{w = braneWedge[b]},
  If[Length[DeleteDuplicates[w]] < 3 || polygonArea2[w] == 0, {},
   Select[triangles, insidePolygonQ[w, Mean[#]] &]]];

braneCharge[b_]       := b[[1]];
braneCut[b_]          := b[[2]];
braneMultiplicity[b_] := Which[
   Length[braneCharge[b]] === 3, braneCharge[b][[3]],
   Length[b] >= 3 && MatchQ[b[[3]], {{_, _} ..}], Length[b[[3]]],
   True, 1];
(* the canonical sign: the polygon lies counterclockwise from the cut *)
braneDirection[b_] := With[{u = invariantDirection[braneCharge[b]],
    d = braneCut[b][[2]] - braneCut[b][[1]]}, If[Det[{d, u}] > 0, u, -u]];
(* one direction per coincident 7-brane; by default they alternate, so
   that k=2 gives the two opposite rays of the invariant line *)
braneDirections[b_] := braneDirections[None, b];
(*  The counterclockwise rule assumes the triangulated part of the
    polygon sits immediately counterclockwise of the cut.  That fails
    for the second of two 7-branes sharing an apex: what lies there has
    already been removed by the first one.  Whenever the star of the tip
    settles the question -- exactly one of the two rays of the invariant
    line enters a triangle -- we let it decide, and fall back on the
    counterclockwise rule otherwise.                                   *)
braneDirectionInStar[data_, b_] :=
 Module[{u = braneDirection[b], v0 = First[braneCut[b]], star, inStar},
  If[! AssociationQ[data], Return[u]];
  star = Select[data["Triangles"], MemberQ[#, v0] &];
  inStar[w_] := AnyTrue[star, inConeQ[w, (# - v0 & /@ Complement[#, {v0}])] &];
  If[! inStar[u] && inStar[-u], -u, u]];

braneDirections[data_, b_] := Which[
   Length[b] >= 3 && MatchQ[b[[3]], {{_, _} ..}], b[[3]],
   Length[b] >= 3 && MatchQ[b[[3]], {_, _}], ConstantArray[b[[3]], braneMultiplicity[b]],
   True, With[{u = braneDirectionInStar[data, b]},
     Table[If[OddQ[j], u, -u], {j, braneMultiplicity[b]}]]];

gtpQuotient[edges_List, branesIn_] :=
 Module[{branes, vertices, vertexPairs, vertexEq, vclass, edgeObjects, edgeIdsByGeom,
   cutEdgePairs, edgePairsOnIds, braneOfPair, imageOfPoint, edgeEq, eclass,
   quotientVertices, quotientEdges, divisorIndex, triangles, repIds, sideCount,
   classIncidence, compactMask, internalEdges, externalEdges, internalReps,
   externalReps, boundaryVerts, quotientInternal, quotientExternal, selfInts, data, wedgeBad},
  branes = checkBranes[branesIn];
  vertices = DeleteDuplicates[Flatten[edges, 1]];
  triangles = triangleList[edges];
  wedgeBad = Select[Table[<|"Brane" -> b, "Triangles" -> braneWedgeTriangles[triangles, b]|>, {b, branes}],
    #["Triangles"] =!= {} &];
  Do[Message[gtpQuotient::wedge, w["Brane"], w["Triangles"]], {w, wedgeBad}];
  (* --- vertex identifications coming from the cuts --- *)
  vertexPairs = DeleteDuplicates[Join @@ (gtpVertexPairs[{braneCharge[#], braneCut[#]}] & /@ branes)];
  vertexEq = makeQuotientRelation[vertices, vertexPairs];
  vclass = vertexEq["ClassOf"];
  quotientVertices = vertexEq["Quotient"];
  divisorIndex = AssociationThread[gtpKey /@ quotientVertices -> Range[Length[quotientVertices]]];
  (* --- edge identifications coming from the cuts --- *)
  edgeObjects = MapIndexed[<|"Id" -> First[#2], "Ends" -> #1, "GeomKey" -> edgeKey[#1]|> &, edges];
  edgeIdsByGeom = GroupBy[edgeObjects, gtpKey[#["GeomKey"]] &, #[[All, "Id"]] &];
  cutEdgePairs = Join @@ (gtpEdgeIdentifications[{braneCharge[#], braneCut[#]}] & /@ branes);
  braneOfPair = Association[];
  imageOfPoint = Association[];
  Do[With[{b = branes[[bi]], img = cutImage[{braneCharge[branes[[bi]]], braneCut[branes[[bi]]]}]},
    Do[imageOfPoint[gtpKey[braneCut[b][[k]]]] = img[[k]], {k, Length[braneCut[b]]}]],
   {bi, Length[branes]}];
  edgePairsOnIds = DeleteDuplicates@Flatten[
     Table[With[{b = branes[[bi]], ids = gtpEdgeIdentifications[{braneCharge[branes[[bi]]], braneCut[branes[[bi]]]}]},
       Table[With[{lhs = Lookup[edgeIdsByGeom, gtpKey[edgeKey[ids[[k, 1]]]], {}],
          rhs = Lookup[edgeIdsByGeom, gtpKey[edgeKey[ids[[k, 2]]]], {}]},
         (braneOfPair[gtpKey[#]] = b) & /@ Tuples[{lhs, rhs}];
         Tuples[{lhs, rhs}]], {k, Length[ids]}]], {bi, Length[branes]}], 2];
  edgeEq = makeQuotientRelation[Range[Length[edges]], edgePairsOnIds];
  eclass = edgeEq["ClassOf"];
  repIds = DeleteDuplicates[eclass /@ Range[Length[edges]]];
  quotientEdges = (vclass /@ edges[[#]] &) /@ repIds;
  (* --- compact / non compact double curves ---------------------------
     an edge class is a compact P^1 iff the total number of (triangle,
     side) incidences of its members is 2 (either one interior edge, or
     two cut edges glued to each other)                                *)
  sideCount = Counts@Flatten@Table[gtpKey[edgeKey[s]], {t, triangles}, {s, Subsets[t, {2}]}];
  classIncidence = Association@Table[
     rep -> Total[Lookup[sideCount, gtpKey[edgeKey[edges[[#]]]], 0] & /@
        Select[Range[Length[edges]], eclass[#] === rep &]], {rep, repIds}];
  compactMask = (classIncidence[#] >= 2 &) /@ repIds;
  internalReps = Pick[repIds, compactMask];
  externalReps = Pick[repIds, Not /@ compactMask];
  internalEdges = Pick[quotientEdges, compactMask];
  externalEdges = Pick[quotientEdges, Not /@ compactMask];
  (* --- compact / non compact divisors ------------------------------- *)
  boundaryVerts = DeleteDuplicates[Flatten[
     (vclass /@ edges[[#]] &) /@ Select[Range[Length[edges]], MemberQ[externalReps, eclass[#]] &], 1]];
  quotientExternal = Select[quotientVertices, MemberQ[boundaryVerts, #] &];
  quotientInternal = Complement[quotientVertices, quotientExternal];
  data = <|"Branes" -> branes, "Vertices" -> vertices, "Triangles" -> triangles,
    "Edges" -> edgeObjects, "CutEdgePairs" -> cutEdgePairs, "CutEdgeIdPairs" -> edgePairsOnIds,
    "BraneOfCutPair" -> braneOfPair, "ImageOfCutPoint" -> imageOfPoint,
    "VertexEquivalence" -> vertexEq, "EdgeEquivalence" -> edgeEq,
    "VerticesQuotiented" -> quotientVertices, "DivisorIndex" -> divisorIndex,
    "QuotientEdgeRepresentatives" -> repIds, "EdgesQuotiented" -> quotientEdges,
    "ClassIncidence" -> classIncidence,
    "ExternalVertices" -> quotientExternal, "InternalVertices" -> quotientInternal,
    "ExternalDivisorIndices" -> Lookup[divisorIndex, gtpKey /@ quotientExternal],
    "InternalDivisorIndices" -> Lookup[divisorIndex, gtpKey /@ quotientInternal],
    "InternalEdges" -> internalEdges, "ExternalEdges" -> externalEdges,
    "InternalEdgeRepresentatives" -> internalReps,
    "ExternalEdgeRepresentatives" -> externalReps, "TriangulatedWedges" -> wedgeBad|>;
  selfInts = AssociationThread[internalReps, selfIntersectionForRep[data, #] & /@ internalReps];
  gtpSmoothEllipticCompletion[Append[data, "InternalEdgeSelfIntersections" -> selfInts]]];

(* Fully smoothing the three corners of a smooth P^2 fan leaves a
   smooth anticanonical elliptic curve. Triangle incidence alone cannot
   see its compactness. Recognize this geometry from the completed fan,
   not from coordinates or an example label. More general fully cut
   stars are deliberately left unsupported.

   For each E-curve, undo the OTHER corner smoothings to obtain a
   one-cut reference presentation. Its local (p,q) rule determines the
   curve class on the unchanged P^2 surface. Further boundary-corner
   smoothings preserve that numerical class. These are numerical
   records, not rays traced through the zero-area Symington polygon. *)
gtpSmoothEllipticCompletion[data_Association] :=
 Module[{bs = data["Branes"], es, tips, v, rays, faces, parentEdges, parent,
   pix, pbe, form, ps, s, d, rep, other, selfInts, records = {}, refEdges,
   ref, line, pq, degree, rootK, direction, record, compactV, externalV},
  If[Length[bs] =!= 3 || data["Triangles"] =!= {} ||
     Length[data["VerticesQuotiented"]] =!= 2 ||
     Length[data["QuotientEdgeRepresentatives"]] =!= 1 ||
     ! AllTrue[bs, Length[braneCut[#]] === 2 && braneMultiplicity[#] === 1 &], Return[data]];
  tips = DeleteDuplicates[First[braneCut[#]] & /@ bs];
  If[Length[tips] =!= 1, Return[data]]; v = First[tips];
  es = data["Edges"][[All, "Ends"]];
  rays = (# - v &) /@ DeleteCases[data["Vertices"], v];
  If[Length[es] =!= 3 || Length[rays] =!= 3 || ! AllTrue[es, MemberQ[#, v] &] ||
     ! AllTrue[rays, GCD @@ # === 1 &] || Total[rays] =!= {0, 0}, Return[data]];
  faces = braneFacePolygon /@ bs;
  If[! AllTrue[faces, Length[#] === 3 && Abs[det2[#[[2]] - v, #[[3]] - v]] === 1 &] ||
     Length[DeleteDuplicates[Sort /@ faces]] =!= 3, Return[data]];
  parentEdges = DeleteDuplicatesBy[Join[es, braneFaceSegment /@ bs], edgeKey];
  parent = gtpQuotient[parentEdges, {}]; pix = gtpIndexSets[parent];
  If[Length[pix["Ic"]] =!= 1, Return[data]]; ps = First[pix["Ic"]];
  pbe = gtpToricBeta[pix]; form = gtpCompactBoundaryForm[pix, pbe, ps];
  If[! AssociationQ[form] || ! TrueQ[form["SpansNumericalPicard"]] ||
     form["PicardRank"] =!= 1, Return[data]];
  rootK = Sqrt[form["K2"]]; If[! IntegerQ[rootK], Return[data]];
  compactV = data["VertexEquivalence"]["ClassOf"][v];
  s = data["DivisorIndex"][gtpKey[compactV]];
  other = First[DeleteCases[data["VerticesQuotiented"], compactV]];
  d = data["DivisorIndex"][gtpKey[other]];
  rep = First[data["QuotientEdgeRepresentatives"]];
  Do[
   refEdges = Select[parentEdges, edgeKey[#] =!= edgeKey[braneFaceSegment[b]] &];
   ref = gtpQuotient[refEdges, {b}]; line = First[gtpELines[ref]];
   If[! AssociationQ[line], Return[data]];
   pq = line["PQ"]; degree = Total[pq]/rootK;
   If[! AllTrue[pq, IntegerQ[#] && # > 0 &] || ! IntegerQ[degree] || degree <= 0 ||
      degree^2 =!= Times @@ pq - 1, Return[data]];
   direction = line["Direction"];
   record = <|"Brane" -> b, "Start" -> s, "StartVertex" -> v,
     "Direction" -> direction, "PQ" -> pq,
     "SpanEdges" -> ({rep, #[[2]], d} & /@ line["SpanEdges"]),
     "Crossed" -> {}, "Triangles" -> {}, "AlongEdge" -> False,
     "Path" -> {}, "EndPoint" -> v, "Termination" -> {"SmoothBoundary", rep},
     "NumericalOnly" -> True, "OtherComponent" -> d,
     "ClassDerivation" -> <|"Method" -> "Preserved class under smooth boundary-corner deformations",
       "ReferencePQ" -> pq, "AnticanonicalDegree" -> Total[pq],
       "SelfIntersection" -> Times @@ pq - 1, "DegreeInPicardGenerator" -> degree|>|>;
   AppendTo[records, record], {b, bs}];
  externalV = DeleteCases[data["ExternalVertices"], compactV];
  selfInts = data["InternalEdgeSelfIntersections"];
  selfInts[rep] = (If[# === compactV, -form["K2"], form["K2"]] & /@
     First[data["EdgesQuotiented"]]);
  <|data, "InternalVertices" -> {compactV}, "ExternalVertices" -> externalV,
    "InternalDivisorIndices" -> {s},
    "ExternalDivisorIndices" -> Lookup[data["DivisorIndex"], gtpKey /@ externalV],
    "InternalEdges" -> data["EdgesQuotiented"], "ExternalEdges" -> {},
    "InternalEdgeRepresentatives" -> {rep}, "ExternalEdgeRepresentatives" -> {},
    "InternalEdgeSelfIntersections" -> selfInts,
    "SmoothEllipticBoundaries" -> <|rep -> <|"Surface" -> s, "OtherComponent" -> d,
       "ArithmeticGenus" -> 1, "NormalizationGenus" -> 1, "NodeCount" -> 0,
       "K2" -> form["K2"], "ParentBoundaryForm" -> form|>|>,
    "NumericalELines" -> records|>];

(* ================================================================== *)
(*  3.  E-lines of the exuberant divisors                              *)
(* ================================================================== *)
(*  The tropicalisation l of the exuberant divisor E starts at the tip
    v of the cut and runs in the monodromy invariant direction u.  We
    follow it through the triangulation and record
      "Start"        the vertex v  (a divisor index),
      "Direction"    u,
      "PQ"           {p,q} with u = p u_1 + q u_2 , where u_1, u_2 are
                     the sides at v of the first triangle, ordered
                     counterclockwise (det(u_1,u_2) = 1) as in the draft,
      "SpanVectors"  {u_1, u_2},
      "SpanEdges"    the classes of the edges e_1, e_2 of the first
                     triangle, together with p and q and the divisor at
                     the far end,
      "Crossed"      the classes of the edges crossed transversally by
                     the deformed line l' (with those of the bypasses),
      "Bypassed"     the compact vertices that l' goes around,
      "Triangles"    the corridor of triangles swept by l',
      "Termination"  {"Edge",class} or {"Vertex",divisorIndex}.
      "BypassSides"  the side ("left"/"right") of each bypass.
    l' is the small deformation of l of Section 5.3.2 of the draft: at
    a compact vertex w met by l it passes on one side of w and crosses
    the edges of the star of w on that side.  data["BypassSide"] is
    "Left", "Right", or Automatic (the default): the left of the
    direction of travel, unless the tip of another 7-brane lies on the
    left of l next to w and none on the right -- passing next to another
    focus-focus singularity needs the third support rule of Section
    5.3.2, which is not implemented.  Only vertices whose class is a
    single lattice point that is not the tip of a cut are bypassed; l
    stops at any other vertex.                                         *)

(* signed area helper *)
det2[a_, b_] := a[[1]] b[[2]] - a[[2]] b[[1]];

(* ------------------------------------------------------------------ *)
(*  the generalized face (T-cone) that has been cut out.  Its boundary
    is  cut[1] -> ... -> cut[n] -> img[n] -> ... -> img[2] -> cut[1] ,
    i.e. one side of the wedge, the primitive boundary segment of the
    GTP carrying the white dots, and the other side of the wedge.     *)
braneFacePolygon[brane_] :=
 Module[{cut = braneCut[brane], img},
  img = cutImage[{braneCharge[brane], cut}];
  DeleteDuplicates[Join[cut, Reverse[Rest[img]]]]];

(*  The boundary segment of the GTP closing the generalized face has
    integral length L and carries the k 7-branes of the brane, so it is
    made of k primitive boundary segments of length L/k.  Its interior
    lattice points are white, except the k-1 junctions between
    consecutive primitive segments, which are black dots (divisors).  *)
braneBaseDots[brane_] :=
 Module[{cut = braneCut[brane], img, a, b, L, k, pts},
  img = cutImage[{braneCharge[brane], cut}];
  a = Last[cut]; b = Last[img];
  L = GCD @@ (b - a); k = braneMultiplicity[brane];
  pts = Table[{j, a + j (b - a)/L}, {j, L - 1}];
  {Cases[pts, {j_, p_} /; ! Divisible[j, L/k] :> p],
   Cases[pts, {j_, p_} /; Divisible[j, L/k] :> p]}];

braneWhiteDots[brane_] := First[braneBaseDots[brane]];
braneBlackDots[brane_] := Last[braneBaseDots[brane]];

braneFaceSegment[brane_] :=
 {Last[braneCut[brane]], Last[cutImage[{braneCharge[brane], braneCut[brane]}]]};

(* is the direction d inside the closed cone spanned by a and b (a,b
   spanning a positively oriented unimodular cone)? *)
inConeQ[d_, {a_, b_}] := Module[{s = Sign[det2[a, b]]},
   s det2[a, d] >= 0 && s det2[d, b] >= 0];

gtpELines[data_Association] :=
 Module[{branes, tri, edges, eclass, vclass, dIndex, res},
  If[KeyExistsQ[data, "NumericalELines"], Return[data["NumericalELines"]]];
  branes = data["Branes"];
  tri = data["Triangles"]; edges = data["Edges"];
  eclass = data["EdgeEquivalence"]["ClassOf"];
  vclass = data["VertexEquivalence"]["ClassOf"];
  dIndex = data["DivisorIndex"];
  res = Join @@ Table[traceELine[data, branes[[k]], #] & /@ braneDirections[data, branes[[k]]],
     {k, Length[branes]}];
  res];

(* the class of the edge {a,b} of the triangulation *)
edgeClassOfSegment[data_, seg_] := Module[{ids},
  ids = Select[Range[Length[data["Edges"]]], edgeKey[data["Edges"][[#, "Ends"]]] === edgeKey[seg] &];
  If[ids === {}, Missing["NoSuchEdge", seg], data["EdgeEquivalence"]["ClassOf"][First[ids]]]];

(*  a vertex met by an E-line is bypassed (Section 5.3.2 of the draft)
    when it is compact, is not the tip of a cut, and its class is the
    single lattice point w, so that its star lies in one chart        *)
bypassableQ[data_, w_] := With[{c = data["VertexEquivalence"]["ClassOf"][w]},
   MemberQ[data["InternalVertices"], c] && FreeQ[First[braneCut[#]] & /@ data["Branes"], w] &&
    Count[data["Vertices"], y_ /; data["VertexEquivalence"]["ClassOf"][y] === c] === 1];

traceELine[data_Association, brane_, dirIn_ : Automatic] :=
 Module[{v0, u, star, t0, w1, w2, pq, spanEdges, crossed = {}, corridor = {},
   cur, dir, tri, term = None, guard = 0, nextT, exitSeg, cand, s, alongEdge, cutIds, m,
   partner, dIndex, vcl, x, segs = {}, curPath, atVertex = None, bypassed = {}, side, sides = {},
   w, starW, nbrs, next, mode, otherTips},
  tri = data["Triangles"]; dIndex = data["DivisorIndex"];
  vcl = data["VertexEquivalence"]["ClassOf"];
  mode = Lookup[data, "BypassSide", Automatic];
  v0 = First[braneCut[brane]];
  otherTips = DeleteCases[First[braneCut[#]] & /@ data["Branes"], v0];
  u = If[dirIn === Automatic, braneDirectionInStar[data, brane], dirIn];
  (* --- the triangles of the star of v0 --- *)
  star = Select[tri, MemberQ[#, v0] &];
  (* does u run along an edge of the star ? *)
  alongEdge = SelectFirst[DeleteDuplicates[Flatten[Complement[#, {v0}] & /@ star, 1]],
    det2[# - v0, u] == 0 && (# - v0) . u > 0 &, None];
  If[alongEdge =!= None,
   (* degenerate case: l runs along the edge [v0,alongEdge] : (p,q)=(1,0) *)
   pq = {1, 0}; {w1, w2} = {alongEdge, None};
   spanEdges = {{edgeClassOfSegment[data, {v0, alongEdge}], 1, dIndex[gtpKey[vcl[alongEdge]]]}};
   corridor = Select[tri, SubsetQ[#, {v0, alongEdge}] &];
   segs = {{v0, alongEdge}}; dir = u; atVertex = alongEdge,
   (* --- generic case : find the star triangle whose cone contains u --- *)
   t0 = SelectFirst[star, inConeQ[u, (# - v0 & /@ Complement[#, {v0}])] &, None];
   If[t0 === None, Return[Missing["ELineLeavesPolygon", v0]]];
   {w1, w2} = Complement[t0, {v0}];
   If[det2[w1 - v0, w2 - v0] < 0, {w1, w2} = {w2, w1}];
   pq = LinearSolve[Transpose[{w1 - v0, w2 - v0}], u];
   spanEdges = {{edgeClassOfSegment[data, {v0, w1}], pq[[1]], dIndex[gtpKey[vcl[w1]]]},
                {edgeClassOfSegment[data, {v0, w2}], pq[[2]], dIndex[gtpKey[vcl[w2]]]}};
   cur = v0; dir = u; nextT = t0; curPath = {v0}];
  (* --- walk: through triangles, across cuts, around compact vertices --- *)
  While[term === None && guard++ < 400,
   If[atVertex =!= None,
    w = atVertex; atVertex = None;
    If[! bypassableQ[data, w], term = {"Vertex", dIndex[gtpKey[vcl[w]]]}; Break[]];
    (* l' passes w on the chosen side: it crosses the edges of the star
       of w that lie strictly on that side, and sweeps those triangles *)
    AppendTo[bypassed, w];
    starW = Select[tri, MemberQ[#, w] &];
    nbrs = DeleteDuplicates[Flatten[Complement[#, {w}] & /@ starW, 1]];
    side = Which[mode === "Left", 1, mode === "Right", -1,
      With[{t = Select[nbrs, MemberQ[otherTips, #] &]},
       AnyTrue[t, det2[dir, # - w] > 0 &] && NoneTrue[t, det2[dir, # - w] < 0 &]], -1,
      True, 1];
    AppendTo[sides, If[side === 1, "left", "right"]];
    crossed = Join[crossed, edgeClassOfSegment[data, {w, #}] & /@ Select[nbrs, side det2[dir, # - w] > 0 &]];
    corridor = Join[corridor, Select[starW, AnyTrue[Complement[#, {w}], Function[y, side det2[dir, y - w] > 0]] &]];
    next = SelectFirst[nbrs, det2[# - w, dir] == 0 && (# - w) . dir > 0 &, None];
    If[next =!= None, AppendTo[segs, {w, next}]; atVertex = next; Continue[]];
    nextT = SelectFirst[starW, inConeQ[dir, (# - w & /@ Complement[#, {w}])] &, None];
    If[nextT === None, term = {"Error", w}; Break[]];
    cur = w; curPath = {w}];
   AppendTo[corridor, nextT];
   exitSeg = eLineExit[nextT, cur, dir];
   If[exitSeg === None, term = {"Error", nextT}; AppendTo[segs, curPath]; Break[]];
   If[First[exitSeg] === "Vertex",
    AppendTo[curPath, exitSeg[[2]]]; AppendTo[segs, curPath];
    atVertex = exitSeg[[2]]; Continue[]];
   (* exitSeg = {"Edge", {a,b}, point} *)
   s = exitSeg[[2]]; x = exitSeg[[3]];
   AppendTo[curPath, x];
   AppendTo[crossed, edgeClassOfSegment[data, s]];
   (* is s a cut edge? *)
   cutIds = SelectFirst[data["CutEdgeIdPairs"],
     (edgeKey[data["Edges"][[#[[1]], "Ends"]]] === edgeKey[s] ||
      edgeKey[data["Edges"][[#[[2]], "Ends"]]] === edgeKey[s]) &, None];
   Which[
    cutIds =!= None,
    (* jump across the cut *)
    m = monodromyMatrix[braneCharge[data["BraneOfCutPair"][gtpKey[cutIds]]]];
    If[edgeKey[data["Edges"][[cutIds[[1]], "Ends"]]] === edgeKey[s],
     partner = data["Edges"][[cutIds[[2]], "Ends"]];
     x = data["ImageOfCutPoint"][gtpKey[First[s]]] + m . (x - First[s]); dir = m . dir,
     partner = data["Edges"][[cutIds[[1]], "Ends"]];
     x = First[Select[braneCut[data["BraneOfCutPair"][gtpKey[cutIds]]],
          data["ImageOfCutPoint"][gtpKey[#]] === First[s] &]] + Inverse[m] . (x - First[s]);
     dir = Inverse[m] . dir];
    AppendTo[segs, curPath]; curPath = {x};
    nextT = SelectFirst[tri, SubsetQ[#, partner] &, None];
    If[nextT === None, term = {"Edge", edgeClassOfSegment[data, partner]};
     AppendTo[segs, curPath], cur = x],
    True,
    cand = Select[tri, SubsetQ[#, s] && # =!= nextT &];
    If[cand === {},
     term = {"Edge", Last[crossed]}; AppendTo[segs, curPath],
     nextT = First[cand]; cur = x]]];
  <|"Brane" -> brane, "Start" -> dIndex[gtpKey[vcl[v0]]],
   "StartVertex" -> v0, "Direction" -> u, "PQ" -> pq,
   "SpanVectors" -> If[w2 === None, {w1 - v0}, {w1 - v0, w2 - v0}], "SpanEdges" -> spanEdges,
   "Crossed" -> crossed, "Triangles" -> DeleteDuplicates[corridor], "AlongEdge" -> (alongEdge =!= None),
   "Bypassed" -> bypassed, "BypassSides" -> sides,
   "Path" -> segs, "EndPoint" -> Last[Last[segs]], "Termination" -> term|>];

(* exit of the ray (cur,dir) from the triangle t.  cur is either a
   vertex of t or a point on one of its sides.                        *)
eLineExit[t_, cur_, dir_] :=
 Module[{sides, hits},
  sides = Subsets[t, {2}];
  hits = Table[
    Module[{a = s[[1]], b = s[[2]], den, lam, mu},
     den = det2[dir, b - a];
     If[den == 0, Nothing,
      lam = det2[a - cur, b - a]/den;   (* cur + lam dir = a + mu (b-a) *)
      mu = det2[a - cur, dir]/den;
      If[lam > 0 && 0 <= mu <= 1, {lam, mu, s, cur + lam dir}, Nothing]]],
   {s, sides}];
  hits = DeleteCases[hits, Nothing];
  If[hits === {}, Return[None]];
  hits = MinimalBy[hits, First];
  Module[{h = First[hits]},
   Which[h[[2]] == 0, {"Vertex", h[[3, 1]]},
    h[[2]] == 1, {"Vertex", h[[3, 2]]},
    True, {"Edge", h[[3]], h[[4]]}]]];

(* ================================================================== *)
(*  4.  Combinatorial input : the index sets I and A (Section 6.1.1)  *)
(* ================================================================== *)

(* Each retained triangle is a triple point. If its two sides at a
   corner become the same double curve, that curve acquires one node.
   For rational normalization this is the arithmetic genus. Recognized
   smooth elliptic boundaries instead have normalization genus one and
   no nodes. Neither genus is inferred from the normal-degree check. *)
gtpDoubleCurveGenera[data_Association] :=
 Module[{reps = data["QuotientEdgeRepresentatives"], nodes, t, a, b},
  nodes = AssociationThread[reps, ConstantArray[0, Length[reps]]];
  Do[
   a = eclassOfSeg[data, {t[[k]], t[[Mod[k, 3] + 1]]}];
   b = eclassOfSeg[data, {t[[k]], t[[Mod[k + 1, 3] + 1]]}];
   If[IntegerQ[a] && a === b, nodes[a] = nodes[a] + 1],
   {t, data["Triangles"]}, {k, 3}];
  Do[nodes[rep] = 1, {rep, Keys[Lookup[data, "SmoothEllipticBoundaries", <||>]]}];
  nodes];

(*  Curves are numbered as in the draft (Example 1 of Section 7):
    first the compact double curves A_{0,c}, in the order of the input
    edges, then the exuberant curves A_e, then the non-compact double
    curves A_{0,boundary}, and last the non-compact curves B_k of the
    recognized smooth elliptic boundaries.  The exuberant divisors are
    E_{|I_0|+k}.  data["CurveNames"] = "Ends" (option "CurveNames" of
    gtpIntersections) names the double curves C_{i,j} after the divisors
    they join instead, and the exuberant curve of E_alpha C_alpha, as in
    the SU(4) example; a prime distinguishes curves joining the same
    two divisors.  The index a of C_a in all arrays is its position in
    this list, so that the unknowns alpha^a_ij and beta_ia carry the
    indices of the draft.                                             *)
gtpIndexSets[data_Association] :=
 Module[{qv, reps, nD0, nE, nC0, nC, lines, div, recs, cur, cIdx, eIdx, genera, smooth,
   boundaryLines, nB, compactQ, toricRec, labels, dI, ends, A0, Ae, Ab},
  qv = data["VerticesQuotiented"];
  reps = data["QuotientEdgeRepresentatives"];
  nD0 = Length[qv]; nC0 = Length[reps];
  (* E-lines that cannot be traced are dropped, so that gtpPlot still
     draws such a GTP; gtpIntersections rejects it before this point *)
  lines = Select[gtpELines[data], AssociationQ]; nE = Length[lines];
  genera = gtpDoubleCurveGenera[data];
  smooth = Lookup[data, "SmoothEllipticBoundaries", <||>];
  boundaryLines = Select[Range[nE], TrueQ[Lookup[lines[[#]], "NumericalOnly", False]] &];
  nB = Length[boundaryLines]; nC = nC0 + nE + nB;
  div = Join[
    Table[<|"Index" -> i, "Type" -> "Toric", "Label" -> "D" <> ToString[i],
      "Vertex" -> qv[[i]], "Compact" -> MemberQ[data["InternalVertices"], qv[[i]]]|>, {i, nD0}],
    Table[<|"Index" -> nD0 + k, "Type" -> "Exuberant",
      "Label" -> "E" <> ToString[nD0 + k], "ELine" -> k, "Compact" -> False|>, {k, nE}]];
  compactQ[r_] := MemberQ[data["InternalEdgeRepresentatives"], r];
  toricRec[r_] := <|"Type" -> "Toric", "Rep" -> r,
     "Ends" -> data["EdgesQuotiented"][[First[FirstPosition[reps, r]]]],
     "ArithmeticGenus" -> If[compactQ[r], genera[r], Missing["NonCompact"]],
     "NormalizationGenus" -> If[KeyExistsQ[smooth, r], 1, 0],
     "NodeCount" -> If[KeyExistsQ[smooth, r], 0, genera[r]],
     "Geometry" -> If[KeyExistsQ[smooth, r], "SmoothElliptic", "RationalNormalization"],
     "Compact" -> compactQ[r]|>;
  recs = Join[
    toricRec /@ Select[reps, compactQ],
    Table[<|"Type" -> "Exuberant", "ELine" -> k,
      "Ends" -> {lines[[k, "StartVertex"]]}, "Compact" -> True,
      "ArithmeticGenus" -> If[Min[lines[[k, "PQ"]]] > 0, Times @@ (lines[[k, "PQ"]] - 1)/2, 0]|>, {k, nE}],
    toricRec /@ Select[reps, ! compactQ[#] &],
    (* B_k = D_boundary.E_k is noncompact. It cannot be replaced by
       -C_k when pairing with a noncompact divisor: BM curve classes
       alone do not justify that numerical pairing. Record only its
       finite intersection with the compact component. *)
    Table[<|"Type" -> "ExuberantBoundary", "ELine" -> k,
      "Component" -> lines[[k, "OtherComponent"]],
      "Ends" -> {qv[[lines[[k, "OtherComponent"]]]]}, "Compact" -> False,
      "ArithmeticGenus" -> Missing["NonCompact"]|>, {k, boundaryLines}]];
  A0 = Flatten[Position[recs, r_ /; r["Type"] === "Toric", {1}, Heads -> False]];
  Ae = Flatten[Position[recs, r_ /; r["Type"] === "Exuberant", {1}, Heads -> False]];
  Ab = Flatten[Position[recs, r_ /; r["Type"] === "ExuberantBoundary", {1}, Heads -> False]];
  (* labels *)
  dI = data["DivisorIndex"];
  labels = Table[Switch[recs[[a, "Type"]],
     "ExuberantBoundary", "BE" <> ToString[recs[[a, "ELine"]]],
     "Exuberant", If[Lookup[data, "CurveNames", "Numbered"] === "Ends",
       "C" <> ToString[nD0 + recs[[a, "ELine"]]], "C" <> ToString[a]],
     _, If[Lookup[data, "CurveNames", "Numbered"] === "Ends",
       ends = Sort[Lookup[dI, gtpKey /@ recs[[a, "Ends"]]]];
       "C" <> ToString[ends[[1]]] <> "," <> ToString[ends[[2]]], "C" <> ToString[a]]], {a, nC}];
  (* primes on repeated names C_{i,j} *)
  labels = Module[{seen = <||>}, Table[seen[labels[[a]]] = Lookup[seen, labels[[a]], -1] + 1;
     labels[[a]] <> StringRepeat["'", seen[labels[[a]]]], {a, nC}]];
  cur = Table[<|"Index" -> a, "Label" -> labels[[a]], recs[[a]]|>, {a, nC}];
  cIdx = AssociationThread[cur[[A0, "Rep"]] -> A0];
  eIdx = AssociationThread[Range[nE] -> Ae];
  <|"Data" -> data, "ELines" -> lines, "nD0" -> nD0, "nE" -> nE, "nC0" -> nC0,
   "nD" -> nD0 + nE, "nC" -> nC, "Divisors" -> div, "Curves" -> cur,
   "CurveIndexOfRep" -> cIdx, "ExuberantCurveIndex" -> eIdx,
   "BoundaryCurveIndexOfELine" -> AssociationThread[boundaryLines -> Ab],
   "I" -> Range[nD0 + nE], "Ic" -> Select[Range[nD0], div[[#, "Compact"]] &],
   "I0" -> Range[nD0], "Ie" -> Range[nD0 + 1, nD0 + nE],
   "A" -> Range[nC], "Ac" -> Select[Range[nC], cur[[#, "Compact"]] &],
   "A0" -> A0, "Ae" -> Ae, "Ab" -> Ab|>];

(* ================================================================== *)
(*  5.  Step 1 : the central-fiber coefficients                        *)
(* ================================================================== *)

(* alpha^a_ij for toric i,j :   D_i.D_j = sum_a alpha^a_ij C_a         *)
gtpToricAlpha[ix_Association] :=
 Module[{data = ix["Data"], qv, nD0, nC0, ends, out},
  qv = ix["Data"]["VerticesQuotiented"]; nD0 = ix["nD0"]; nC0 = ix["nC0"];
  ends = ix["Curves"][[All, "Ends"]];
  out = Table[
    If[i === j, ConstantArray[0, ix["nC"]],
     Table[If[MemberQ[ix["A0"], a] && Sort[ends[[a]]] === Sort[{qv[[i]], qv[[j]]}], 1, 0], {a, ix["nC"]}]],
    {i, nD0}, {j, nD0}];
  (* X0 = Sum_i D_i is principal. A conductor/self-double curve of one
     component is not an additional component in this sum. *)
  Do[out[[i, i]] = -Total[Delete[out[[i]], i]], {i, nD0}];
  out];

(* beta_ia = D_i.C_a for toric i and toric a ; Undefined when D_i cap C_a
   is non compact *)
gtpToricBeta[ix_Association] :=
 Module[{data = ix["Data"], qv, nD0, nC0, reps, selfInts, vclass, eclass, dIdx,
   triContrib, out, edges, tri},
  data = ix["Data"]; qv = data["VerticesQuotiented"]; nD0 = ix["nD0"]; nC0 = ix["nC0"];
  reps = data["QuotientEdgeRepresentatives"]; selfInts = data["InternalEdgeSelfIntersections"];
  vclass = data["VertexEquivalence"]["ClassOf"]; eclass = data["EdgeEquivalence"]["ClassOf"];
  edges = data["Edges"]; tri = data["Triangles"];
  (* number of (triangle,side) incidences with side in class a and
     opposite vertex mapping to v_i *)
  triContrib[i_, rep_] := Total@Flatten@Table[
     Boole[eclassOfSeg[data, s] === rep &&
       vclass[First[Complement[t, s]]] === qv[[i]]],
     {t, tri}, {s, Subsets[t, {2}]}];
  (* the columns are all the curves C_a of ix; those that are not
     double curves are left Undefined here *)
  out = Table[If[! MemberQ[ix["A0"], a], Undefined,
    Module[{v = qv[[i]], rep = ix["Curves"][[a, "Rep"]], e = ix["Curves"][[a, "Ends"]], val, cmp},
     cmp = MemberQ[data["InternalEdgeRepresentatives"], rep];
     If[! cmp && MemberQ[e, v], Undefined,
      val = triContrib[i, rep];
      If[cmp,
       val += If[e[[1]] === v, selfInts[rep][[1]], 0] + If[e[[2]] === v, selfInts[rep][[2]], 0]];
      val]]],
    {i, nD0}, {a, ix["nC"]}];
  out];

eclassOfSeg[data_, seg_] := Module[{ids},
  ids = Select[Range[Length[data["Edges"]]], edgeKey[data["Edges"][[#, "Ends"]]] === edgeKey[seg] &];
  If[ids === {}, None, data["EdgeEquivalence"]["ClassOf"][First[ids]]]];

(* ================================================================== *)
(*  6.  Steps 2 and 3 : admissibility and the exuberant ansatz         *)
(* ================================================================== *)

Clear[aa, bb];
(* alpha^a_ij and beta_ia as in the draft, with a comma between the
   indices when one of them has two digits: alpha^9_64, beta_6,10      *)
indexRow[l_List] := If[Max[l] >= 10, Row[l, ","], Row[l]];
aa /: MakeBoxes[aa[i_, j_, a_], f_] := With[{r = indexRow[{i, j}]}, MakeBoxes[Subsuperscript["\[Alpha]", r, a], f]];
bb /: MakeBoxes[bb[i_, a_], f_] := With[{r = indexRow[{i, a}]}, MakeBoxes[Subscript["\[Beta]", r], f]];

(* classes of all sides of the triangles swept by the E-line *)
corridorClasses[data_, line_] :=
  DeleteDuplicates@Flatten@Table[eclassOfSeg[data, s], {t, line["Triangles"]}, {s, Subsets[t, {2}]}];

(* Boundary curve intersections on a smooth compact component. Distinct
   boundary curves meet at the retained triangles. Diagonal entries are
   recovered from D_j|D_i = Sum of the curves joining i to j. *)
gtpCompactBoundaryForm[ix_Association, be_, i_Integer] :=
 Module[{data = ix["Data"], v, boundary, neighbors, q, cidx, edges, k2, rho,
   va, a, b, pos, t, genera, euler},
  v = data["VerticesQuotiented"][[i]];
  boundary = Select[ix["A0"], ix["Curves"][[#, "Compact"]] && MemberQ[ix["Curves"][[#, "Ends"]], v] &];
  If[boundary === {} || AnyTrue[boundary, Count[ix["Curves"][[#, "Ends"]], v] =!= 1 &],
   Return[Missing["UnsupportedCompactBoundary"]]];
  neighbors = Table[data["DivisorIndex"][gtpKey[First[DeleteCases[ix["Curves"][[a, "Ends"]], v]]]], {a, boundary}];
  q = ConstantArray[0, {Length[boundary], Length[boundary]}];
  cidx = ix["CurveIndexOfRep"];
  Do[If[data["VertexEquivalence"]["ClassOf"][t[[k]]] === v,
    a = cidx[eclassOfSeg[data, {t[[k]], t[[Mod[k, 3] + 1]]}]];
    b = cidx[eclassOfSeg[data, {t[[k]], t[[Mod[k + 1, 3] + 1]]}]];
    If[a =!= b && MemberQ[boundary, a] && MemberQ[boundary, b],
     pos = First[FirstPosition[boundary, a]]; va = First[FirstPosition[boundary, b]];
     q[[pos, va]]++; q[[va, pos]]++]],
   {t, data["Triangles"]}, {k, 3}];
  Do[q[[k, k]] = be[[neighbors[[k]], boundary[[k]]]] -
     Total@Table[If[l =!= k && neighbors[[l]] === neighbors[[k]], q[[l, k]], 0], {l, Length[boundary]}],
   {k, Length[boundary]}];
  k2 = -Total[be[[i, boundary]]]; rho = 10 - k2;
  genera = ix["Curves"][[boundary, "ArithmeticGenus"]];
  euler = Total[(2 - 2 #["NormalizationGenus"] - #["NodeCount"]) & /@
     ix["Curves"][[boundary]]] - Total[Flatten[UpperTriangularize[q, 1]]];
  <|"Curves" -> boundary, "IntersectionMatrix" -> q, "Rank" -> MatrixRank[q],
    "PicardRank" -> rho, "K2" -> k2, "BoundaryEulerCharacteristic" -> euler,
    "SpansNumericalPicard" -> (IntegerQ[rho] && rho > 0 && MatrixRank[q] === rho)|>];

(* For two E-curves on the same compact surface, boundary pairings
   determine their mutual intersection when the boundary spans its
   numerical Picard group. If Q c_l = b_l, then C_k.C_l = b_k.c_l.
   A pseudoinverse handles redundant boundary representatives exactly.
   Otherwise leave the pairing as an unknown: distinct E-curves are
   not automatically disjoint.
   When some boundary pairings b_k are still unknown, the same identity
   C_k.C_l = b_k.Q^+.b_l is returned as an equation of the compatibility
   system, so that it applies once compatibility has fixed them. *)
gtpNumericalExuberantPairings[ix_Association, al_, beta_] :=
 Module[{be = beta, forms = <||>, records = {}, eqs = {}, lines = ix["ELines"], eC = ix["ExuberantCurveIndex"], s, form,
   boundary, qinv, bk, bl, value, nk, nl, unitK, unitL},
  Do[s = lines[[k, "Start"]];
   If[s === lines[[l, "Start"]] && MemberQ[ix["Ic"], s],
    unitK = UnitVector[ix["nC"], eC[k]];
    unitL = UnitVector[ix["nC"], eC[l]];
    If[al[[s, ix["nD0"] + k]] === unitK && al[[s, ix["nD0"] + l]] === unitL,
     If[! KeyExistsQ[forms, s], forms[s] = gtpCompactBoundaryForm[ix, be, s]];
     form = forms[s];
     If[AssociationQ[form] && TrueQ[form["SpansNumericalPicard"]],
      boundary = form["Curves"]; bk = be[[ix["nD0"] + k, boundary]]; bl = be[[ix["nD0"] + l, boundary]];
      If[AllTrue[Flatten[form["IntersectionMatrix"]], IntegerQ] && ! AllTrue[Join[bk, bl], IntegerQ] &&
         FreeQ[{bk, bl}, Undefined] && MatchQ[be[[ix["nD0"] + k, eC[l]]], _bb],
       value = Expand[bk . PseudoInverse[form["IntersectionMatrix"]] . bl];
       eqs = Join[eqs, {be[[ix["nD0"] + k, eC[l]]] == value, be[[ix["nD0"] + l, eC[k]]] == value}];
       AppendTo[records, <|"ELines" -> {k, l}, "Surface" -> s, "BoundaryForm" -> form,
         "BoundaryPairings" -> {bk, bl}, "Intersection" -> value, "ImposedAsEquation" -> True|>]];
      If[AllTrue[Join[bk, bl, Flatten[form["IntersectionMatrix"]]], IntegerQ],
       qinv = PseudoInverse[form["IntersectionMatrix"]];
       nk = Simplify[bk . qinv . bk]; nl = Simplify[bl . qinv . bl];
       If[Simplify[form["IntersectionMatrix"] . qinv . bk - bk] === ConstantArray[0, Length[bk]] &&
          Simplify[form["IntersectionMatrix"] . qinv . bl - bl] === ConstantArray[0, Length[bl]] &&
          nk === be[[ix["nD0"] + k, eC[k]]] && nl === be[[ix["nD0"] + l, eC[l]]],
        value = Simplify[bk . qinv . bl];
        If[IntegerQ[value],
         be[[ix["nD0"] + k, eC[l]]] = value;
         be[[ix["nD0"] + l, eC[k]]] = value;
         AppendTo[records, <|"ELines" -> {k, l}, "Surface" -> s, "BoundaryForm" -> form,
           "BoundaryPairings" -> {bk, bl}, "Intersection" -> value|>]]]]]]],
   {k, ix["nE"]}, {l, k + 1, ix["nE"]}];
  <|"Beta" -> be, "Derivations" -> records, "Equations" -> eqs|>];

gtpAnsatz[ix_Association] :=
 Module[{data = ix["Data"], nD0 = ix["nD0"], nC0 = ix["nC0"], nE = ix["nE"], nD, nC,
   lines = ix["ELines"], al, be, curves = ix["Curves"], div = ix["Divisors"],
   cIdx = ix["CurveIndexOfRep"], eC = ix["ExuberantCurveIndex"], qv, unk = {}, spanOf, corrOf,
   crossOf, startOf, termOf, pqOf, pairings},
  nD = ix["nD"]; nC = ix["nC"]; qv = data["VerticesQuotiented"];
  (* ---- per E-line bookkeeping ---- *)
  spanOf  = Table[Merge[(Rule[cIdx[#[[1]]], #[[2]]] & /@ lines[[k, "SpanEdges"]]), Total], {k, nE}];
  corrOf  = Table[cIdx /@ corridorClasses[data, lines[[k]]], {k, nE}];
  crossOf = Table[cIdx /@ DeleteDuplicates[lines[[k, "Crossed"]]], {k, nE}];
  startOf = Table[lines[[k, "Start"]], {k, nE}];
  termOf  = Table[lines[[k, "Termination"]], {k, nE}];
  pqOf    = Table[lines[[k, "PQ"]], {k, nE}];

  (* ---------------- alpha ---------------- *)
  al = Table[0, {nD}, {nD}, {nC}];
  (* toric x toric : Step 1 *)
  With[{a0 = gtpToricAlpha[ix]},
   Do[al[[i, j]] = a0[[i, j]], {i, nD0}, {j, nD0}]];
  (* toric x exuberant / exuberant x toric : Step 3, support rules of
     Section 5.3.2 (the curves met by the deformed line l')           *)
  Do[Module[{iE = nD0 + k, vec},
    Do[
     vec = Table[0, {nC}];
     (* (i) the exuberant curve itself *)
     If[j === startOf[[k]], vec[[eC[k]]] = 1];
     If[KeyExistsQ[ix["BoundaryCurveIndexOfELine"], k] && j === lines[[k, "OtherComponent"]],
      vec[[ix["BoundaryCurveIndexOfELine"][k]]] = 1];
     (* (ii) edges crossed transversally by l', including the edges of
        the star of a bypassed vertex on the side of the bypass       *)
     Do[If[MemberQ[curves[[a, "Ends"]], qv[[j]]], vec[[a]] = aa[iE, j, a]], {a, crossOf[[k]]}];
     (* (iii) the E-line stops at a boundary lattice point : E cap D_j is a
        curve on the toric surface D_j, a combination of all the curves
        through v_j -- not only the boundary ones, which miss the curve
        along the last edge when the E-line arrives between two interior
        edges                                                          *)
     If[MatchQ[termOf[[k]], {"Vertex", _}] && j === termOf[[k, 2]],
      Do[If[MemberQ[curves[[a, "Ends"]], qv[[j]]],
        vec[[a]] = aa[iE, j, a]], {a, ix["A0"]}]];
     al[[iE, j]] = vec; al[[j, iE]] = vec,
     {j, nD0}]], {k, nE}];
  (* exuberant x exuberant is left undefined -- see gtpIntersectionValue *)
  Do[al[[nD0 + k, nD0 + l]] = Undefined, {k, nE}, {l, nE}];

  (* ---------------- beta ---------------- *)
  be = Table[Undefined, {nD}, {nC}];
  With[{b0 = gtpToricBeta[ix]},
   Do[be[[i, a]] = b0[[i, a]], {i, nD0}, {a, ix["A0"]}]];
  (* toric divisor . exuberant curve : the exuberant curve rules       *)
  Do[Module[{k = kk, val},
    Do[
     val = 0;
     If[j === startOf[[k]], val -= Total[pqOf[[k]]]];
     Do[If[j === lines[[k, "SpanEdges", r, 3]], val += lines[[k, "SpanEdges", r, 2]]],
      {r, Length[lines[[k, "SpanEdges"]]]}];
     be[[j, eC[k]]] = val, {j, nD0}]], {kk, nE}];
  Do[be[[startOf[[k]], ix["BoundaryCurveIndexOfELine"][k]]] = Total[pqOf[[k]]],
   {k, Keys[ix["BoundaryCurveIndexOfELine"]]}];
  (* exuberant divisor . curve : Steps 2 and 3.  The admissibility
     domain B_e of the draft: a non-compact double curve C_a is
     inadmissible against E_alpha iff the E-line meets the relative
     interior of its edge, i.e. ends on it; meeting an end point of the
     edge does not exclude the pair.  Admissible pairs are fixed by the
     local E-line data, unknown when C_a has an edge in the corridor
     of l', and 0 otherwise (the representatives are disjoint).       *)
  Do[Module[{iE = nD0 + k},
    Do[
     be[[iE, a]] = Which[
       MemberQ[ix["Ab"], a], Undefined,
       (* Mixed pairings require curve-class information or compatibility. *)
       MemberQ[ix["Ae"], a], If[a === eC[k], Times @@ pqOf[[k]] - 1, bb[iE, a]],
       KeyExistsQ[spanOf[[k]], a], spanOf[[k]][a],
       ! curves[[a, "Compact"]] && MatchQ[termOf[[k]], {"Edge", r_} /; cIdx[r] === a], Undefined,
       MemberQ[corrOf[[k]], a], bb[iE, a],
       True, 0], {a, nC}]], {k, nE}];
  pairings = gtpNumericalExuberantPairings[ix, al, be]; be = pairings["Beta"];
  unk = DeleteDuplicates@Cases[{al, be}, _aa | _bb, Infinity];
  <|"IndexSets" -> ix, "Alpha" -> al, "Beta" -> be, "Unknowns" -> unk,
    "ExuberantPairingDerivations" -> pairings["Derivations"],
    "PairingEquations" -> pairings["Equations"]|>];

(* ================================================================== *)
(*  7.  Steps 4 and 5 : compatibility equations, triple intersections  *)
(* ================================================================== *)

(* Always expand the compact curve D_i.D_k first. All its numerical
   pairings are admissible, including two exuberant divisor factors. *)
gtpTripleValue[ans_, i_, j_, k_] :=
 Module[{ix = ans["IndexSets"], al = ans["Alpha"], be = ans["Beta"]},
  Total@Table[al[[i, k, a]] be[[j, a]], {a, ix["Ac"]}]];

gtpCompatibilityEquations[ans_] :=
 Module[{ix = ans["IndexSets"], al = ans["Alpha"], be = ans["Beta"], eqs = {}, ie, lhs, rhs},
  ie = ix["Ie"];
  Do[
   If[MemberQ[ie, j] || MemberQ[ie, k],
    If[MemberQ[ie, j] && MemberQ[ie, k],
     (* symmetry of D_i.E_j.E_k *)
     lhs = Total@Table[al[[i, k, a]] be[[j, a]], {a, ix["Ac"]}];
     rhs = Total@Table[al[[i, j, a]] be[[k, a]], {a, ix["Ac"]}],
     lhs = Total@Table[al[[j, k, a]] be[[i, a]], {a, ix["A"]}];
     rhs = Total@Table[al[[i, k, a]] be[[j, a]], {a, ix["Ac"]}]];
    If[! TrueQ[Simplify[lhs - rhs == 0]], AppendTo[eqs, lhs == rhs]]],
   {i, ix["Ic"]}, {j, ix["I"]}, {k, ix["I"]}];
  DeleteDuplicates[eqs]];

(* purely toric consistency checks (no unknowns involved) *)
gtpToricChecks[ans_] :=
 Module[{ix = ans["IndexSets"], al = ans["Alpha"], be = ans["Beta"], bad = {}, lhs, rhs},
  Do[lhs = Total@Table[al[[j, k, a]] be[[i, a]], {a, ix["A"]}];
   rhs = Total@Table[al[[i, k, a]] be[[j, a]], {a, ix["Ac"]}];
   If[lhs =!= rhs, AppendTo[bad, {i, j, k, lhs, rhs}]],
   {i, ix["Ic"]}, {j, ix["I0"]}, {k, ix["I0"]}];
  bad];

(*  The presentation T_{i;jk} = D_i.(D_j.D_k) = sum_a alpha^a_jk beta_ia
    of Section 5.4 of the draft.  It exists when D_j.D_k has an expansion (j and k
    are not both exuberant) and the expansion is admissible against
    D_i : beta_ia is defined on its whole support.                     *)
gtpPresentation[ans_, i_, j_, k_] :=
 Module[{ix = ans["IndexSets"], al = ans["Alpha"], be = ans["Beta"], ie, terms},
  ie = ix["Ie"];
  If[MemberQ[ie, j] && MemberQ[ie, k], Return[Missing["NoExpansion"]]];
  terms = Table[If[al[[j, k, a]] === 0, 0, {al[[j, k, a]], be[[i, a]]}], {a, ix["A"]}];
  If[MemberQ[terms, {_, Undefined}], Missing["NotAdmissible"],
   Total[Replace[terms, {x_, y_} :> x y, {1}]]]];

(*  All presentations of every triple, compact or not (Step 4): any two
    of T_{i;jk}, T_{j;ik}, T_{k;ij} that exist must agree.  A relation
    free of unknowns among toric divisors is a check of the toric data,
    returned separately; one involving an exuberant divisor stays in the
    system, where it contradicts the ansatz if it fails.               *)
gtpAllPresentationEquations[ans_] :=
 Module[{ix = ans["IndexSets"], n, eqs = {}, bad = {}, vals, d},
  n = Length[ix["I"]];
  Do[
   vals = DeleteCases[{gtpPresentation[ans, i, j, k], gtpPresentation[ans, j, i, k],
      gtpPresentation[ans, k, i, j]}, _Missing];
   If[Length[vals] > 1, Do[d = Expand[vals[[1]] - v];
     Which[d === 0, Null,
      ! FreeQ[d, _aa | _bb] || IntersectingQ[{i, j, k}, ix["Ie"]], AppendTo[eqs, d == 0],
      True, AppendTo[bad, {i, j, k, vals[[1]], v}]], {v, Rest[vals]}]],
   {i, n}, {j, i, n}, {k, j, n}];
  <|"Equations" -> DeleteDuplicates[eqs], "Checks" -> bad|>];

(* ================================================================== *)
(*  8.  Driver                                                         *)
(* ================================================================== *)

(*  "Presentations" -> "All" imposes every relation of Step 4 ;
    "CompactFirst" only those with a compact first divisor.
    The solution branches are kept when they admit an integral solution;
    if the surviving branches disagree on D_i.D_j.D_k the status is
    "MultipleSolutions".  "Rank" is the number of independent equations
    at the solution, "KernelDimension" the number of directions left
    free; "ResidualParameters" lists those entering D_i.D_j.D_k.
    "FreeParameters" are the unknowns that the compatibility equations
    leave free.  When none of them enters D_i.D_j.D_k the status is
    "Solved", the coefficients that depend on them are not determined:
    "SymbolicAlpha" and "SymbolicBeta" keep them, "Alpha" and "Beta"
    use the integral gauge choice "Gauge", "CoefficientsDetermined" is
    False, and the tables of gtpReport show those entries in terms of
    the free unknowns.
    "CurveNames" -> "Numbered" (C_1, C_2, ... as in Example 1 of the
    draft) or "Ends" (C_{i,j}, as in the SU(4) example); "BypassSide"
    -> Automatic, "Left" or "Right" is the side on which l' goes around
    the compact vertices met by an E-line (see traceELine).            *)
Options[gtpIntersections] = {"Verbose" -> False, "Presentations" -> "All",
   "CurveNames" -> "Numbered", "BypassSide" -> Automatic};

gtpIntersections[edges_List, branes_, OptionsPattern[]] :=
 Module[{data, ix, ans, eqs, unk, sol, rules, al, be, tri, free, res, checks, pres, allEq,
   intQ, ints, triOf, tris, rank, fail, gauge, validation, freeUnk, ans0},
  data = <|gtpQuotient[edges, branes], "CurveNames" -> OptionValue["CurveNames"],
    "BypassSide" -> OptionValue["BypassSide"]|>;
  If[data["TriangulatedWedges"] =!= {},
   Return[<|"Status" -> "InvalidGTP", "Data" -> data,
     "Reason" -> "The generalized face of a 7-brane is triangulated: " <>
       StringRiffle[("the triangle(s) " <> ToString[#["Triangles"]] <> " lie in the face of the 7-brane " <>
          ToString[#["Brane"]]) & /@ data["TriangulatedWedges"], "; "] <>
       ".  Remove the edges inside the face."|>]];
  With[{bad = Select[gtpELines[data], ! AssociationQ[#] &]},
   If[bad =!= {}, Return[<|"Status" -> "UnsupportedELine", "Data" -> data,
      "MissingELines" -> bad,
      "Reason" -> "The E-lines have no supported retained-triangle presentation. " <>
       "Fully cut smooth P^2 stars are supported; other empty stars need additional surface data."|>]]];
  ix   = gtpIndexSets[data];
  (* an E-line that stops at a compact vertex -- one that is not
     bypassed by l': a vertex on a cut or the tip of another 7-brane --
     is not covered by the support rules of Step 3 as implemented      *)
  With[{bad = Select[ix["ELines"], MatchQ[#["Termination"], {"Vertex", j_} /; MemberQ[ix["Ic"], j]] &]},
   If[bad =!= {}, Return[<|"Status" -> "UnsupportedELine", "Data" -> data, "IndexSets" -> ix,
      "Reason" -> "The E-line(s) starting at " <> ToString[#["StartVertex"] & /@ bad] <>
        " run into the compact vertex(es) " <> ToString[ix["Divisors"][[#["Termination"][[2]], "Vertex"]] & /@ bad] <>
        ", which l' cannot bypass: a vertex identified by a cut, or the tip of a 7-brane, " <>
        "would need the rule for passing next to another focus-focus singularity, which is not implemented."|>]]];
  ans  = gtpAnsatz[ix];
  checks = gtpToricChecks[ans];
  If[OptionValue["Presentations"] === "CompactFirst",
   eqs = gtpCompatibilityEquations[ans]; pres = {},
   allEq = gtpAllPresentationEquations[ans]; eqs = allEq["Equations"]; pres = allEq["Checks"]];
  eqs = DeleteDuplicates[Join[eqs, ans["PairingEquations"]]];
  unk  = ans["Unknowns"];
  fail[reason_, extra_ : <||>] := <|"Status" -> "Inconsistent", "Reason" -> reason, "Data" -> data,
     "IndexSets" -> ix, "Ansatz" -> ans, "Equations" -> eqs, "ToricChecks" -> checks,
     "PresentationChecks" -> pres, extra|>;
  If[checks =!= {} || pres =!= {}, Return[fail["Central-fiber intersection presentations disagree."]]];
  sol  = Quiet@Solve[eqs, unk];
  If[sol === {}, Return[fail["The compatibility equations of Step 4 have no solution."]]];
  (* integrality: an unverified branch (time out) is kept *)
  intQ[r_] := TimeConstrained[
    Quiet@FindInstance[And @@ Join[eqs, Equal @@@ r], unk, Integers] =!= {}, 60, True];
  ints = Select[sol, intQ];
  If[ints === {}, Return[fail["The compatibility equations have rational solutions but no integral one.",
     <|"Solutions" -> sol|>]]];
  triOf[r_] := With[{a2 = <|ans, "Alpha" -> (ans["Alpha"] /. r), "Beta" -> (ans["Beta"] /. r)|>},
    Association@Table[i ->
      Table[Simplify@gtpTripleValue[a2, i, j, k], {j, ix["I"]}, {k, ix["I"]}], {i, ix["Ic"]}]];
  tris = triOf /@ ints;
  If[Length[DeleteDuplicates[tris]] > 1,
   Return[<|fail["The compatibility equations have " <> ToString[Length[ints]] <>
       " integral solution branches with different triple intersection numbers."],
     "Status" -> "MultipleSolutions", "Solutions" -> ints, "BranchIntersections" -> tris|>]];
  rules = First[ints]; tri = First[tris];
  (* independent equations at the solution; the unknowns left free are
     set to fixed random integers, so that bilinear terms do not vanish *)
  rank = If[eqs === {} || unk === {}, 0,
    BlockRandom[SeedRandom[1];
     MatrixRank[D[(#[[1]] - #[[2]]) & /@ eqs, {unk}] /. rules /.
       Thread[unk -> RandomInteger[{-97, 97}, Length[unk]]]]]];
  free = DeleteDuplicates@Cases[Values[tri], _aa | _bb, Infinity];
  (* the unknowns left free; when they drop out of D_i.D_j.D_k they are
     fixed by a gauge choice : 0 if that is integral, else an integral
     solution                                                          *)
  freeUnk = Complement[unk, Keys[Association[Rule @@@ rules]]];
  gauge = Thread[freeUnk -> 0];
  If[free === {} && gauge =!= {} && ! AllTrue[unk /. rules /. gauge, IntegerQ],
   gauge = Quiet@FindInstance[And @@ Join[eqs, Equal @@@ rules], unk, Integers];
   gauge = If[gauge === {}, {}, Select[First[gauge], MemberQ[gauge[[All, 1]], #[[1]]] &]]];
  If[free =!= {}, gauge = {}];
  al = ans["Alpha"] /. rules /. gauge; be = ans["Beta"] /. rules /. gauge;
  ans0 = ans; ans = <|ans, "Alpha" -> al, "Beta" -> be|>;
  res = <|"Status" -> If[free === {}, "Solved", "Underdetermined"],
    "Data" -> data, "IndexSets" -> ix, "Ansatz" -> ans, "Equations" -> eqs,
    "Presentations" -> OptionValue["Presentations"],
    "Unknowns" -> unk, "Solution" -> rules, "ToricChecks" -> checks,
    "PresentationChecks" -> pres,
    "Rank" -> rank, "KernelDimension" -> Length[unk] - rank,
    "Branches" -> Length[sol], "IntegralBranches" -> Length[ints],
    "Alpha" -> al, "Beta" -> be, "Intersections" -> tri,
    "SymbolicAlpha" -> (ans0["Alpha"] /. rules), "SymbolicBeta" -> (ans0["Beta"] /. rules),
    "CoefficientsDetermined" -> (freeUnk === {}),
    "ExuberantPairingDerivations" -> (ans["ExuberantPairingDerivations"] /. rules /. gauge),
    "ResidualParameters" -> free,
    "FreeParameters" -> freeUnk, "Gauge" -> gauge|>;
  validation = gtpConsistencyChecks[res];
  res = <|res, "ConsistencyChecks" -> validation|>;
  If[AssociationQ[validation] && ! AllTrue[Values[validation], # === {} &],
   res = <|res, "Status" -> "Inconsistent", "Reason" -> "The reconstructed intersections fail geometric consistency checks."|>];
  If[OptionValue["Verbose"], gtpReport[res]];
  res];

(* ------------------------- presentation ---------------------------- *)

gtpDivisorLabels[ix_] := ix["Divisors"][[All, "Label"]];
gtpCurveLabels[ix_]   := ix["Curves"][[All, "Label"]];

(*  "D3" -> D_3, "E15" -> E_15, "C1,7'" -> C_{1,7}' as the front end
    displays them; the labels themselves stay strings                 *)
gtpLabelForm[s_String] := Module[{m = StringCases[s,
     RegularExpression["^([A-Za-z]+)([0-9,]+)('*)$"] -> {"$1", "$2", "$3"}]},
  Which[m === {}, s,
   m[[1, 3]] === "", Subscript[m[[1, 1]], m[[1, 2]]],
   True, Row[{Subscript[m[[1, 1]], m[[1, 2]]], m[[1, 3]]}]]];
gtpLabelForm[x_] := x;

(* sum_a c_a C_a in the order of the curves, e.g.  -C_2 + 2 C_7       *)
gtpCombination[coeffs_List, labels_List] :=
 Module[{t = Select[Transpose[{coeffs, labels}], ! TrueQ[#[[1]] == 0] &], term},
  term[{c_, l_}, first_] := With[{lf = gtpLabelForm[l]}, Which[
     c === 1, If[first, {lf}, {" + ", lf}],
     c === -1, If[first, {"\[Minus]", lf}, {" \[Minus] ", lf}],
     NumberQ[c] && c < 0, If[first, {"\[Minus]", -c, " ", lf}, {" \[Minus] ", -c, " ", lf}],
     NumberQ[c] || MatchQ[c, _aa | _bb], If[first, {c, " ", lf}, {" + ", c, " ", lf}],
     True, If[first, {"(", c, ") ", lf}, {" + (", c, ") ", lf}]]];
  If[t === {}, 0, Row[Join @@ MapIndexed[term[#1, First[#2] === 1] &, t]]]];

(* coefficients left free by the compatibility equations (they drop out
   of D_i.D_j.D_k when the status is Solved) are shown in this colour,
   in terms of the free unknowns                                       *)
gtpUndeterminedColor = RGBColor[0.85, 0.42, 0.0];
gtpUndeterminedQ[res_, x_] := MatchQ[Lookup[res, "FreeParameters", {}], {__}] &&
   ! FreeQ[x, Alternatives @@ res["FreeParameters"]];
gtpSymbolic[res_, key_] := Lookup[res, "Symbolic" <> key, res[key]];

gtpIntersectionTables[res_] :=
 Module[{ix = res["IndexSets"], lab, tri = res["Intersections"]},
  lab = gtpLabelForm /@ gtpDivisorLabels[ix];
  Column@Table[
    Labeled[
     Grid[Join[{Prepend[lab, Style[Row[{lab[[i]], " \[CenterDot]"}], Bold, Blue]]},
       Table[Prepend[Table[If[k < j, "*", tri[i][[j, k]]], {k, ix["nD"]}], Style[lab[[j]], Bold]],
        {j, ix["nD"]}]], Frame -> All, Alignment -> Center,
      Background -> {{LightBlue, None}, {LightBlue, None}}],
     Style[Row[{lab[[i]], " \[CenterDot] ", Subscript["D", "j"], " \[CenterDot] ", Subscript["D", "k"]}], Bold], Top],
    {i, ix["Ic"]}]];

(*  Results that are not solved carry "Status", "Reason" and "Data" (so
    that gtpPlot can still draw them): the tables are replaced by a
    statement of what went wrong.                                      *)
gtpSolvedQ[res_] := AssociationQ[res] && MemberQ[{"Solved", "Underdetermined"}, res["Status"]];
gtpNotSolved[res_] := Style["Not available: gtpIntersections returned the status " <>
    ToString[If[AssociationQ[res], res["Status"], res]] <>
    If[AssociationQ[res] && StringQ[res["Reason"]], ".  " <> res["Reason"], "."], Red];
(gtpIntersectionTables[res_] /; ! gtpSolvedQ[res]) := gtpNotSolved[res];
(gtpDivisorSummary[res_] /; ! gtpSolvedQ[res]) := gtpNotSolved[res];
(gtpCurveSummary[res_] /; ! gtpSolvedQ[res]) := gtpNotSolved[res];
(gtpBetaTable[res_] /; ! gtpSolvedQ[res]) := gtpNotSolved[res];
(gtpAlphaTable[res_] /; ! gtpSolvedQ[res]) := gtpNotSolved[res];
Options[gtpReport] = {"IntersectionTables" -> True};
(gtpReport[res_, OptionsPattern[]] /; ! gtpSolvedQ[res]) := (
   Print[Style["Status: " <> ToString[If[AssociationQ[res], res["Status"], res]], Bold, Red]];
   If[AssociationQ[res] && StringQ[res["Reason"]], Print[res["Reason"]]];
   If[AssociationQ[res] && MatchQ[res["ToricChecks"], {__}],
    Print[Style["toric consistency failures: ", Red], res["ToricChecks"]]];
   If[AssociationQ[res] && MatchQ[res["PresentationChecks"], {__}],
    Print[Style["non-compact presentations that disagree: ", Red], res["PresentationChecks"]]];);

(* one line per E-line, in the terms of Section 6.1 of the draft      *)
gtpELineDescription[ix_, k_] :=
 Module[{l = ix["ELines"][[k]], dl = gtpLabelForm /@ gtpDivisorLabels[ix],
   cl = gtpLabelForm /@ gtpCurveLabels[ix], cIdx = ix["CurveIndexOfRep"],
   dI = ix["Data"]["DivisorIndex"], vcl = ix["Data"]["VertexEquivalence"]["ClassOf"], vec},
  vec[v_] := Row[{"(", Row[v, ","], ")"}];
  If[TrueQ[Lookup[l, "NumericalOnly", False]],
   Return[Row[{"E-line of ", dl[[ix["nD0"] + k]], ": numerical curve class from corner smoothing, ",
      l["ClassDerivation"]}]]];
  Row[{"E-line ", Subscript["\[ScriptL]", ix["nD0"] + k], " of ", dl[[ix["nD0"] + k]], ": starts at ",
    dl[[l["Start"]]], " = ", vec[l["StartVertex"]], " in the direction u = ", vec[l["Direction"]], ", ",
    If[TrueQ[l["AlongEdge"]],
     Row[{"along the edge towards ", dl[[l["SpanEdges"][[1, 3]]]], ", (p,q) = (1,0)"}],
     Row[{"u = p ", Subscript["u", 1], " + q ", Subscript["u", 2], " with ", Subscript["u", 1], " = ",
       vec[l["SpanVectors"][[1]]], " (towards ", dl[[l["SpanEdges"][[1, 3]]]], "), ", Subscript["u", 2], " = ",
       vec[l["SpanVectors"][[2]]], " (towards ", dl[[l["SpanEdges"][[2, 3]]]], "), (p,q) = ", vec[l["PQ"]]}]],
    "; exuberant curve ", cl[[ix["ExuberantCurveIndex"][k]]],
    "; crosses ", If[l["Crossed"] === {}, "no edge", Row[cl[[cIdx /@ DeleteDuplicates[l["Crossed"]]]], ", "]],
    If[Lookup[l, "Bypassed", {}] === {}, "",
     Row[{"; ", Subscript["\[ScriptL]", ix["nD0"] + k], "' passes ",
       Row[MapThread[Row[{"on the ", #2, " of ", #1}] &,
         {dl[[Lookup[dI, gtpKey /@ (vcl /@ l["Bypassed"])]]], l["BypassSides"]}], ", "]}]],
    "; ends ", Switch[l["Termination"],
     {"Edge", _}, Row[{"on the edge of ", cl[[cIdx[l["Termination"][[2]]]]]}],
     {"Vertex", _}, Row[{"at the vertex of ", dl[[l["Termination"][[2]]]]}],
     _, ToString[l["Termination"]]]}]];

gtpReport[res_, OptionsPattern[]] :=
 Module[{ix = res["IndexSets"], dl, cl, dI, row},
  dl = gtpLabelForm /@ gtpDivisorLabels[ix]; cl = gtpLabelForm /@ gtpCurveLabels[ix];
  dI = ix["Data"]["DivisorIndex"]; row[l_] := Row[l, ", "];
  Print[Style["Status: " <> res["Status"], Bold],
   If[res["Status"] === "Solved" && ! TrueQ[Lookup[res, "CoefficientsDetermined", True]],
    Style["   (the triple intersection numbers are determined, some coefficients \[Alpha], \[Beta] are not)",
     gtpUndeterminedColor], ""]];
  Print["divisors  I = {", row[dl], "},   compact ", Subscript["I", "c"], " = {", row[dl[[ix["Ic"]]]], "}"];
  Print["curves    A = {", row[cl], "},   compact ", Subscript["A", "c"], " = {", row[cl[[ix["Ac"]]]], "}"];
  Print["vertices  : ", row[Table[Row[{dl[[i]], " \[RightArrow] ", ix["Data"]["VerticesQuotiented"][[i]]}], {i, ix["I0"]}]]];
  Print["double curves : ", row[Table[Row[{cl[[a]], " = ",
       Row[dl[[Lookup[dI, gtpKey /@ ix["Curves"][[a, "Ends"]]]]], " \[CenterDot] "]}], {a, ix["A0"]}]]];
  Do[Print[gtpELineDescription[ix, k]], {k, ix["nE"]}];
  If[res["ToricChecks"] =!= {}, Print[Style["toric consistency failures: ", Red], res["ToricChecks"]]];
  If[MatchQ[res["PresentationChecks"], {__}],
   Print[Style["non-compact presentations that disagree: ", Red], res["PresentationChecks"]]];
  Print["unknowns: ", Length[res["Unknowns"]], ",  independent equations: ", res["Rank"],
   ",  kernel dimension: ", res["KernelDimension"],
   Which[res["KernelDimension"] === 0, "",
    res["ResidualParameters"] === {}, "  (the free directions drop out of D_i.D_j.D_k)",
    True, "  (the free directions enter D_i.D_j.D_k)"]];
  Print[gtpDivisorSummary[res]];
  Print[gtpCurveSummary[res]];
  Print[Style[Row[{Subscript["D", "i"], " \[CenterDot] ", Subscript["D", "j"], " = \[Sum] ",
      Subsuperscript["\[Alpha]", "ij", "a"], " ", Subscript["C", "a"]}], Bold]]; Print[gtpAlphaTable[res]];
  Print[Style[Row[{Subscript["\[Beta]", "ia"], " = ", Subscript["D", "i"], " \[CenterDot] ", Subscript["C", "a"],
      "   (\[Times] : not in the admissibility domain B)"}], Bold]]; Print[gtpBetaTable[res]];
  If[MatchQ[res["FreeParameters"], {__}],
   Print[Style[Row[{"not fixed by the compatibility equations: ", row[res["FreeParameters"]],
       If[res["ResidualParameters"] === {},
        Row[{".  They drop out of every ", Subscript["D", "i"], " \[CenterDot] ", Subscript["D", "j"], " \[CenterDot] ",
          Subscript["D", "k"], "; the entries of the tables that depend on them are shown in this colour, ",
          "and res[\"Alpha\"], res[\"Beta\"] use the choice ", row[res["Gauge"]], "."}], "."]}], gtpUndeterminedColor]]];
  Print[Style["consistency checks: ", Bold], gtpConsistencyChecks[res]];
  If[TrueQ[OptionValue["IntersectionTables"]], Print[gtpIntersectionTables[res]]];];

(* ------------------------- more tables ----------------------------- *)

gtpDivisorSummary[res_] :=
 Module[{ix = res["IndexSets"], data = res["Data"], qv, lab, nb, chg, tri = res["Intersections"]},
  qv = data["VerticesQuotiented"]; lab = gtpLabelForm /@ ix["Divisors"][[All, "Label"]];
  nb[i_] := Count[ix["A0"], a_ /; MemberQ[ix["Curves"][[a, "Ends"]], qv[[i]]]];
  chg[i_] := Count[ix["ELines"], l_ /; l["Start"] === i];
  Grid[Prepend[
    Table[{lab[[i]], If[i <= ix["nD0"], qv[[i]], "-"],
      If[MemberQ[ix["Ic"], i], "compact", "non-compact"],
      If[i <= ix["nD0"], nb[i], "-"], If[i <= ix["nD0"], chg[i], "-"],
      If[MemberQ[ix["Ic"], i], tri[i][[i, i]], "-"]}, {i, ix["I"]}],
    Style[#, Bold] & /@ {"D", "vertex", "type", "# bdry comps", "charge Q", Superscript["D", 3]}],
   Frame -> All, Alignment -> Center, Background -> {None, {LightBlue, None}}]];

gtpCurveSummary[res_] :=
 Module[{ix = res["IndexSets"], lab, data = res["Data"], dI, type},
  lab = gtpLabelForm /@ ix["Divisors"][[All, "Label"]]; dI = data["DivisorIndex"];
  type = <|"Toric" -> "double curve", "Exuberant" -> "exuberant", "ExuberantBoundary" -> "boundary of E"|>;
  Grid[Prepend[
    Table[{gtpLabelForm[ix["Curves"][[a, "Label"]]], type[ix["Curves"][[a, "Type"]]],
      If[ix["Curves"][[a, "Compact"]], "compact", "non-compact"],
      If[ix["Curves"][[a, "Type"]] === "Toric",
       Row[lab[[Lookup[dI, gtpKey /@ ix["Curves"][[a, "Ends"]]]]], " \[CenterDot] "],
       Row[{lab[[If[ix["Curves"][[a, "Type"]] === "ExuberantBoundary",
           ix["Curves"][[a, "Component"]], ix["ELines"][[ix["Curves"][[a, "ELine"]], "Start"]]]]],
         lab[[ix["nD0"] + ix["Curves"][[a, "ELine"]]]]}, " \[CenterDot] "]],
      If[ix["Curves"][[a, "Type"]] === "Toric", ix["Curves"][[a, "Rep"]], "-"]}, {a, ix["A"]}],
    Style[#, Bold] & /@ {"C", "type", "", Row[{"= ", Subscript["D", "i"], " \[CenterDot] ", Subscript["D", "j"]}], "input edge"}],
   Frame -> All, Alignment -> Center, Background -> {None, {LightBlue, None}}]];

gtpBetaTable[res_] :=
 Module[{ix = res["IndexSets"], be = gtpSymbolic[res, "Beta"]},
  Grid[Join[{Prepend[gtpLabelForm /@ ix["Curves"][[All, "Label"]], Style["\[Beta]", Bold]]},
    Table[Prepend[Table[With[{x = be[[i, a]]}, Which[
         x === Undefined, "\[Times]",
         gtpUndeterminedQ[res, x], Style[x, gtpUndeterminedColor],
         True, x]], {a, ix["nC"]}],
      Style[gtpLabelForm[ix["Divisors"][[i, "Label"]]], Bold]], {i, ix["nD"]}]],
   Frame -> All, Alignment -> Center, Background -> {{LightBlue, None}, {LightBlue, None}}]];

gtpAlphaTable[res_] :=
 Module[{ix = res["IndexSets"], lab, al = gtpSymbolic[res, "Alpha"]},
  lab = ix["Curves"][[All, "Label"]];
  Grid[Join[{Prepend[gtpLabelForm /@ ix["Divisors"][[All, "Label"]], Style["\[Alpha]", Bold]]},
    Table[Prepend[
      Table[Which[al[[i, j]] === Undefined, "\[LongDash]",
        j < i, "*",
        gtpUndeterminedQ[res, al[[i, j]]], Style[gtpCombination[al[[i, j]], lab], gtpUndeterminedColor],
        True, gtpCombination[al[[i, j]], lab]],
       {j, ix["nD"]}], Style[gtpLabelForm[ix["Divisors"][[i, "Label"]]], Bold]], {i, ix["nD"]}]],
   Frame -> All, Alignment -> Center, ItemSize -> Full, Background -> {{LightBlue, None}, {LightBlue, None}}]];

(* --------------------------- picture ------------------------------- *)
(*  gtpPictureData collects everything that is drawn, so that gtpPlot
    and gtpTikZ always show the same thing.                            *)

gtpPictureData[dataIn_Association] :=
 Module[{data, ix, lines, vertices, edges, vclass, eclass, dIdx, cIdx, selfInts, beta,
   faces, faceSegs, whiteDots, blackDots, edgeRecs, vertRecs, eRecs, centroid, outward, allSegs},
  data = If[KeyExistsQ[dataIn, "Data"], dataIn["Data"], dataIn];
  ix = gtpIndexSets[data]; lines = ix["ELines"];
  vertices = data["Vertices"]; edges = data["Edges"][[All, "Ends"]];
  vclass = data["VertexEquivalence"]["ClassOf"]; eclass = data["EdgeEquivalence"]["ClassOf"];
  dIdx = data["DivisorIndex"]; cIdx = ix["CurveIndexOfRep"];
  selfInts = data["InternalEdgeSelfIntersections"]; beta = gtpToricBeta[ix];
  faces     = braneFacePolygon /@ data["Branes"];
  faceSegs  = braneFaceSegment /@ data["Branes"];
  whiteDots = Flatten[braneWhiteDots /@ data["Branes"], 1];
  blackDots = Flatten[braneBlackDots /@ data["Branes"], 1];
  centroid  = Mean[N[vertices]];
  (* where to put the label of a vertex: in the middle of the widest
     angular gap between the edges leaving it *)
  allSegs = Join[edges, faceSegs, Join @@ (#["Path"] & /@ lines)];
  outward[v_] := Module[{nbrs, ang, gaps, k, a},
    nbrs = DeleteCases[Flatten[Select[allSegs, MemberQ[#, v] &] /. v -> Nothing, 1], v];
    nbrs = DeleteDuplicates[nbrs];
    If[nbrs === {}, Return[Normalize[N[v] - centroid + {10^-3, 10^-3}]]];
    ang = Sort[Mod[ArcTan[#[[1]], #[[2]]] & /@ (N[# - v] & /@ nbrs), 2 Pi]];
    gaps = Append[Differences[ang], 2 Pi - Last[ang] + First[ang]];
    k = First[Ordering[gaps, -1]];
    a = ang[[k]] + gaps[[k]]/2;
    {Cos[a], Sin[a]}];
  edgeRecs = Table[
    With[{p1 = edges[[k, 1]], p2 = edges[[k, 2]], rep = eclass[k]},
     <|"P1" -> p1, "P2" -> p2, "Label" -> ix["Curves"][[cIdx[rep], "Label"]],
      "Compact" -> MemberQ[data["InternalEdgeRepresentatives"], rep],
      (* the full beta_{ia} = D_i.C_a: when vertices are identified it
         also counts the triple points where C_a meets D_i once more *)
      "Beta" -> If[KeyExistsQ[selfInts, rep],
        {beta[[dIdx[gtpKey[vclass[p1]]], cIdx[rep]]], beta[[dIdx[gtpKey[vclass[p2]]], cIdx[rep]]]},
        None]|>], {k, Length[edges]}];
  vertRecs = Table[
    <|"P" -> v, "Label" -> ix["Divisors"][[dIdx[gtpKey[vclass[v]]], "Label"]],
     "Compact" -> MemberQ[data["InternalVertices"], vclass[v]], "Out" -> outward[v]|>,
    {v, vertices}];
  eRecs = Table[
    <|"Paths" -> lines[[k, "Path"]], "Label" -> ix["Divisors"][[ix["nD0"] + k, "Label"]],
     "CurveLabel" -> ix["Curves"][[ix["ExuberantCurveIndex"][k], "Label"]],
     "Start" -> lines[[k, "StartVertex"]], "Direction" -> lines[[k, "Direction"]],
     "EndPoint" -> lines[[k, "EndPoint"]], "PQ" -> lines[[k, "PQ"]],
     "Open" -> (First[lines[[k, "Termination"]]] === "Edge")|>,
    {k, Select[Range[Length[lines]], ! TrueQ[Lookup[lines[[#]], "NumericalOnly", False]] &]}];
  <|"Data" -> data, "IndexSets" -> ix, "Faces" -> faces, "FaceSegments" -> faceSegs,
   "WhiteDots" -> whiteDots, "ExtraBlackDots" -> blackDots, "EdgeRecords" -> edgeRecs, "VertexRecords" -> vertRecs,
   "ELineRecords" -> eRecs|>];

(* the label of a double curve sits beside the middle of its edge *)
curveLabelPosition[e_, off_] :=
 Module[{p1 = N[e["P1"]], p2 = N[e["P2"]], d},
  d = Normalize[p2 - p1];
  (p1 + p2)/2 + off {-d[[2]], d[[1]]}];

(* the exuberant curve C_E = D_v . E is labelled next to the beginning
   of the E-line *)
eCurveLabelPosition[l_, off_] :=
 Module[{d = Normalize[N[l["Direction"]]]},
  N[l["Start"]] + 0.42 d - 1.5 off {-d[[2]], d[[1]]}];

(* where to put the label of the exuberant divisor: past the arrow tip
   when the E-line leaves the polygon, else beside its last segment *)
eLineLabelPosition[l_, ov_] :=
 Module[{d = Normalize[N[l["Direction"]]]},
  If[TrueQ[l["Open"]], N[l["EndPoint"]] + (ov + 0.3) d,
   N[l["EndPoint"]] - 0.42 d - 0.26 {-d[[2]], d[[1]]}]];

Options[gtpPlot] = {"ShadeColor" -> Magenta, "ShadeOpacity" -> 0.25,
   "CurveLabels" -> True, "DivisorLabels" -> True, "SelfIntersections" -> True,
   "Shading" -> True, "ELines" -> True, "LabelOffset" -> 0.3, "ELineOvershoot" -> 0.5,
   "CurveLabelOffset" -> 0.17};

gtpPlot[dataIn_Association, opts : OptionsPattern[{gtpPlot, Graphics}]] :=
 Module[{pic, gopt, ov, off, prim},
  pic = gtpPictureData[dataIn];
  gopt = FilterRules[Join[{opts}, Options[gtpPlot]], Options[gtpPlot]];
  ov[k_] := OptionValue[gtpPlot, gopt, k];
  off = ov["LabelOffset"];
  prim = {
    (* generalized faces that have been cut out *)
    If[TrueQ[ov["Shading"]],
     {Directive[ov["ShadeColor"], Opacity[ov["ShadeOpacity"]]], Polygon /@ pic["Faces"]}, {}],
    (* primitive boundary segments closing the generalized faces *)
    {Black, Thick, Line /@ pic["FaceSegments"]},
    (* the triangulation *)
    Table[{If[e["Compact"], Darker[Green], Black], Thick, Line[{e["P1"], e["P2"]}]},
     {e, pic["EdgeRecords"]}],
    (* the E-lines *)
    If[TrueQ[ov["ELines"]],
     Table[{Red, Thick,
       Line /@ l["Paths"],
       If[l["Open"],
        Arrow[{l["EndPoint"],
          l["EndPoint"] + ov["ELineOvershoot"] Normalize[N[l["Direction"]]]}], {}],
       Inset[Style[gtpLabelForm[l["Label"]], 14, Bold, Red], eLineLabelPosition[l, ov["ELineOvershoot"]]]},
      {l, pic["ELineRecords"]}], {}],
    (* lattice points: black dots, and the white dots of the GTP *)
    Table[{If[v["Compact"], Darker[Green], Blue], PointSize[.022], Point[v["P"]]},
     {v, pic["VertexRecords"]}],
    {Blue, PointSize[.022], Point /@ pic["ExtraBlackDots"]},
    {Black, PointSize[.022], Point /@ pic["WhiteDots"]},
    {White, PointSize[.014], Point /@ pic["WhiteDots"]},
    (* divisor labels *)
    If[TrueQ[ov["DivisorLabels"]],
     Table[Inset[Style[gtpLabelForm[v["Label"]], 15, Bold, Red, Background -> White],
       v["P"] + off v["Out"]], {v, pic["VertexRecords"]}], {}],
    (* curve labels, including the exuberant curves C_E = D_v.E *)
    If[TrueQ[ov["CurveLabels"]],
     {Table[Inset[Style[gtpLabelForm[e["Label"]], 11, Bold, If[e["Compact"], Darker[Green], Gray],
         Background -> White], curveLabelPosition[e, ov["CurveLabelOffset"]]],
       {e, pic["EdgeRecords"]}],
      Table[Inset[Style[gtpLabelForm[l["CurveLabel"]], 11, Bold, Red, Background -> White],
        eCurveLabelPosition[l, ov["CurveLabelOffset"]]], {l, pic["ELineRecords"]}]}, {}],
    (* the two normal bundle degrees beta_{ia} = D_i.C_a, next to the
       end of the edge belonging to the divisor D_i they refer to *)
    If[TrueQ[ov["SelfIntersections"]],
     Table[If[e["Beta"] === None, {},
       {Inset[Style[ToString[e["Beta"][[1]]], 10, Background -> White], (3 e["P1"] + e["P2"])/4],
        Inset[Style[ToString[e["Beta"][[2]]], 10, Background -> White], (e["P1"] + 3 e["P2"])/4]}],
      {e, pic["EdgeRecords"]}], {}]};
  Graphics[prim, FilterRules[{opts}, Options[Graphics]], PlotRange -> All, ImageSize -> Large]];


(* --------------------------- TikZ export --------------------------- *)
(*  gtpTikZ[data] returns a String containing a tikzpicture that
    reproduces gtpPlot[data]: the triangulated GTP, the shaded
    generalized faces with their white dots, the divisor labels, the
    E-lines with the exuberant divisor labels, the curve labels and the
    normal bundle degrees beta_{ia}.
    The snippet needs   \usepackage{tikz}  and
    \usetikzlibrary{arrows.meta}  (already loaded in the draft).
    With "DefineStyles"->True (default) all node styles are declared
    locally, so nothing has to be added to the preamble.              *)

(* fixed point, at most three decimals, no trailing zeros, no exponents *)
tkNum[x_] := Module[{r = N[x], sg, ip, fp},
  If[Abs[r - Round[r]] < 10^-9, Return[ToString[Round[r]]]];
  sg = If[r < 0, "-", ""]; r = Abs[r];
  ip = Floor[r]; fp = Round[1000 (r - ip)];
  If[fp >= 1000, ip += 1; fp -= 1000];
  If[fp === 0, Return[sg <> ToString[ip]]];
  sg <> ToString[ip] <> "." <>
   StringReplace[StringPadLeft[ToString[fp], 3, "0"], RegularExpression["0+$"] -> ""]];
tkPt[p_] := "(" <> tkNum[p[[1]]] <> "," <> tkNum[p[[2]]] <> ")";
tkPath[pts_] := StringRiffle[tkPt /@ pts, "--"];
(* D1 -> $D_{1}$ , C12 -> $C_{12}$ , E6 -> $E_{6}$ , C1,7' -> $C_{1,7}'$ *)
tkLabel[s_String] := Module[{m = StringCases[s,
     RegularExpression["^([A-Za-z]+)([0-9,]*)('*)$"] -> {"$1", "$2", "$3"}]},
  Which[m === {}, "$" <> s <> "$",
   m[[1, 2]] === "", "$" <> m[[1, 1]] <> m[[1, 3]] <> "$",
   True, "$" <> m[[1, 1]] <> "_{" <> m[[1, 2]] <> "}" <> m[[1, 3]] <> "$"]];

Options[gtpTikZ] = {"Scale" -> 2, "DefineStyles" -> True, "Shading" -> True,
   "DivisorLabels" -> True, "CurveLabels" -> True, "SelfIntersections" -> True,
   "ELines" -> True, "LabelOffset" -> 0.34, "ELineOvershoot" -> 0.55,
   "CurveLabelOffset" -> 0.17,
   "ShadeColor" -> "magenta", "ShadeOpacity" -> 0.3, "Standalone" -> False,
   "Comments" -> True, "CompactColor" -> "green!55!black", "NonCompactColor" -> "black!45"};

gtpTikZ[dataIn_Association, opts : OptionsPattern[]] :=
 Module[{pic, o, out = {}, add, cmt, off, ov, lbl},
  pic = gtpPictureData[dataIn];
  o[k_] := OptionValue[gtpTikZ, {opts}, k];
  off = o["LabelOffset"]; ov = o["ELineOvershoot"];
  add[l_] := AppendTo[out, l];
  cmt[l_] := If[TrueQ[o["Comments"]], AppendTo[out, "  % " <> l]];

  add["\\begin{tikzpicture}[scale=" <> tkNum[o["Scale"]] <>
    If[TrueQ[o["DefineStyles"]],
     ",\n    bd/.style={circle,fill=black,inner sep=1.5pt}," <>
      "\n    wd/.style={circle,draw=black,fill=white,line width=.6pt,inner sep=1.1pt}," <>
      "\n    dlbl/.style={text=blue,font=\\scriptsize,inner sep=1pt}," <>
      "\n    clbl/.style={font=\\tiny,inner sep=.6pt,fill=white,fill opacity=.7,text opacity=1}," <>
      "\n    nbd/.style={font=\\tiny,inner sep=.4pt,fill=white,fill opacity=.7,text opacity=1}," <>
      "\n    elbl/.style={text=red,font=\\scriptsize,inner sep=1pt}", ""] <> "]"];

  (* shaded generalized faces *)
  If[TrueQ[o["Shading"]] && pic["Faces"] =!= {},
   cmt["generalized faces cut out of the GTP"];
   Do[add["  \\fill[" <> o["ShadeColor"] <> ",opacity=" <> tkNum[o["ShadeOpacity"]] <> "] " <>
      tkPath[f] <> "--cycle;"], {f, pic["Faces"]}]];

  (* primitive boundary segments closing the generalized faces *)
  If[pic["FaceSegments"] =!= {},
   cmt["primitive boundary segments carrying the white dots"];
   Do[add["  \\draw " <> tkPath[sg] <> ";"], {sg, pic["FaceSegments"]}]];

  (* the triangulation *)
  cmt["triangulation of the Symington polygon"];
  Do[add["  \\draw" <> If[e["Compact"], "[" <> o["CompactColor"] <> "]", ""] <> " " <>
     tkPath[{e["P1"], e["P2"]}] <> ";"], {e, pic["EdgeRecords"]}];

  (* the E-lines *)
  If[TrueQ[o["ELines"]] && pic["ELineRecords"] =!= {},
   cmt["E-lines and exuberant divisors"];
   Do[
    Do[add["  \\draw[red,thick] " <> tkPath[p] <> ";"], {p, l["Paths"]}];
    If[l["Open"],
     add["  \\draw[red,thick,-{Stealth}] " <> tkPt[l["EndPoint"]] <> "--" <>
       tkPt[N[l["EndPoint"]] + ov Normalize[N[l["Direction"]]]] <> ";"]];
    add["  \\node[elbl] at " <> tkPt[eLineLabelPosition[l, ov]] <>
      " {" <> tkLabel[l["Label"]] <> "};"],
    {l, pic["ELineRecords"]}]];

  (* lattice points *)
  cmt["lattice points"];
  Do[add["  \\node[bd] at " <> tkPt[v["P"]] <> " {};"], {v, pic["VertexRecords"]}];
  Do[add["  \\node[bd] at " <> tkPt[w] <> " {};"], {w, pic["ExtraBlackDots"]}];
  Do[add["  \\node[wd] at " <> tkPt[w] <> " {};"], {w, pic["WhiteDots"]}];

  (* divisor labels *)
  If[TrueQ[o["DivisorLabels"]],
   cmt["divisors"];
   Do[add["  \\node[dlbl] at " <> tkPt[N[v["P"]] + off v["Out"]] <> " {" <>
      tkLabel[v["Label"]] <> "};"], {v, pic["VertexRecords"]}]];

  (* curve labels *)
  If[TrueQ[o["CurveLabels"]],
   cmt["curves"];
   Do[add["  \\node[clbl,text=" <> If[e["Compact"], o["CompactColor"], o["NonCompactColor"]] <>
      "] at " <> tkPt[curveLabelPosition[e, o["CurveLabelOffset"]]] <>
      " {" <> tkLabel[e["Label"]] <> "};"], {e, pic["EdgeRecords"]}];
   If[TrueQ[o["ELines"]],
    Do[add["  \\node[clbl,text=red] at " <>
       tkPt[eCurveLabelPosition[l, o["CurveLabelOffset"]]] <>
       " {" <> tkLabel[l["CurveLabel"]] <> "};"], {l, pic["ELineRecords"]}]]];

  (* normal bundle degrees *)
  If[TrueQ[o["SelfIntersections"]],
   cmt["normal bundle degrees  beta_{ia} = D_i . C_a  of the compact double curves"];
   Do[If[e["Beta"] =!= None,
     add["  \\node[nbd] at " <> tkPt[(3 N[e["P1"]] + N[e["P2"]])/4] <> " {$" <>
       ToString[e["Beta"][[1]]] <> "$};"];
     add["  \\node[nbd] at " <> tkPt[(N[e["P1"]] + 3 N[e["P2"]])/4] <> " {$" <>
       ToString[e["Beta"][[2]]] <> "$};"]], {e, pic["EdgeRecords"]}]];

  add["\\end{tikzpicture}"];
  lbl = StringRiffle[out, "\n"];
  If[TrueQ[o["Standalone"]],
   lbl = "\\documentclass[border=5pt]{standalone}\n" <>
     "\\usepackage[dvipsnames]{xcolor}\n\\usepackage{tikz}\n" <>
     "\\usetikzlibrary{arrows.meta}\n" <>
     "\\begin{document}\n" <> lbl <> "\n\\end{document}"];
  lbl];

(* ================================================================== *)
(*  9.  Independent consistency checks                                 *)
(* ================================================================== *)
(*  None of these is used by the algorithm; they test the output.
      C1 : adjunction for a compact double curve C_a on the
           normalization of each of the two surfaces D_x containing it:
             (C_a|D_x)^2 + beta_{xa} - Delta_x.C_a = 2 g(C_a) - 2 ,
           with (C_a|D_x)^2 the normal degree of the fan relation, g
           the genus of the normalization of C_a (0, or 1 for a
           recognized smooth elliptic curve), and Delta_x.C_a the number
           of triple points of C_a whose third vertex is again in the
           class of v_x (K of the normalization of D_x is (D_x - Delta_x)|,
           Delta_x its double locus; the nodes of C_a cancel out).  When
           D_x and D_y are normal along C_a it is the relation
           beta_{xa} + beta_{ya} = 2 p_a(C_a) - 2 of the draft; when a
           triangle on C_a has its third vertex in the class of an end,
           beta picks up that transverse branch of D_x and the draft's
           relation does not hold
      C2 : sum over all non exuberant divisors of beta_{ia} = 0 for
           every compact curve C_a
      C3 : D_i^3 = 12 - chi_top(boundary D_i) - Q_i for compact D_i
           (Noether's formula for the log Calabi-Yau surface D_i)
      C4 : D_i.D_j.D_k is totally symmetric
      C5 : sum_k D_i.D_j.D_k = 0 over non exuberant k
      C6 : on a compact D_i whose boundary form is known, the Gram
           matrix of the boundary curves and of the exuberant curves
           C_k = D_i.E_k has rank <= rho = 10 - K^2 and at most one
           positive eigenvalue (Hodge index)                           *)

(gtpConsistencyChecks[res_] /; ! gtpSolvedQ[res]) :=
  Missing["NotSolved", If[AssociationQ[res], res["Status"], res]];

gtpConsistencyChecks[res_] :=
 Module[{ix = res["IndexSets"], data = res["Data"], be = res["Beta"], tri = res["Intersections"],
   qv, out = <||>, bad, lab, nb, chg, eqQ, boundaryForm, expected, gram, transverse},
  qv = data["VerticesQuotiented"]; lab = ix["Divisors"][[All, "Label"]];
  eqQ[x_, y_] := TrueQ[Simplify[x == y]];
  (* C1 *)
  transverse[rep_, x_] := Total@Flatten@Table[
     Boole[eclassOfSeg[data, sd] === rep && data["VertexEquivalence"]["ClassOf"][First[Complement[t, sd]]] === x],
     {t, data["Triangles"]}, {sd, Subsets[t, {2}]}];
  bad = Table[With[{e = ix["Curves"][[a, "Ends"]], rep = ix["Curves"][[a, "Rep"]]},
     If[ix["Curves"][[a, "Compact"]] && e[[1]] =!= e[[2]],
      With[{si = data["InternalEdgeSelfIntersections"][rep]},
       expected = 2 ix["Curves"][[a, "NormalizationGenus"]] - 2;
       Table[With[{x = e[[r]], ix0 = data["DivisorIndex"][gtpKey[e[[r]]]]},
         If[ListQ[si] && eqQ[si[[3 - r]] + be[[ix0, a]] - transverse[rep, x], expected], Nothing,
          {ix["Curves"][[a, "Label"]], lab[[ix0]], If[ListQ[si], si[[3 - r]] + be[[ix0, a]] - transverse[rep, x], si],
           expected}]], {r, 2}]],
      Nothing]], {a, ix["A0"]}];
  out["NormalBundle"] = DeleteCases[Flatten[bad, 1], Nothing];
  (* C2 *)
  out["CurveSumRule"] = DeleteCases[Table[
     If[eqQ[Total[be[[ix["I0"], a]]], 0], Nothing, {ix["Curves"][[a, "Label"]], Total[be[[ix["I0"], a]]]}],
     {a, ix["Ac"]}], Nothing];
  (* C3 *)
  nb[i_] := (boundaryForm = gtpCompactBoundaryForm[ix, be, i];
    If[AssociationQ[boundaryForm], boundaryForm["BoundaryEulerCharacteristic"],
     Count[ix["A0"], a_ /; MemberQ[ix["Curves"][[a, "Ends"]], qv[[i]]]]]);
  chg[i_] := Count[ix["ELines"], l_ /; l["Start"] === i];
  out["Noether"] = DeleteCases[Table[
     If[eqQ[tri[i][[i, i]], 12 - nb[i] - chg[i]], Nothing,
      {lab[[i]], tri[i][[i, i]], 12 - nb[i] - chg[i]}], {i, ix["Ic"]}], Nothing];
  (* C4 *)
  out["Symmetry"] = Join[DeleteCases[Flatten[Table[
      If[eqQ[tri[i][[j, k]], tri[i][[k, j]]], Nothing, {i, j, k}],
      {i, ix["Ic"]}, {j, ix["I"]}, {k, j, ix["nD"]}], 2], Nothing],
    DeleteCases[Flatten[Table[
      If[eqQ[tri[i][[j, k]], tri[j][[i, k]]], Nothing, {i, j, k}],
      {i, ix["Ic"]}, {j, ix["Ic"]}, {k, ix["I"]}], 2], Nothing]];
  (* C5 *)
  out["DivisorSumRule"] = DeleteCases[Flatten[Table[
      If[eqQ[Total[tri[i][[j, ix["I0"]]]], 0], Nothing, {i, j, Total[tri[i][[j, ix["I0"]]]]}],
      {i, ix["Ic"]}, {j, ix["I"]}], 1], Nothing];
  (* C6 *)
  gram[i_] := Module[{form = gtpCompactBoundaryForm[ix, be, i], b, es, g},
    If[! AssociationQ[form], Return[Nothing]];
    b = form["Curves"];
    es = Select[Range[ix["nE"]], ix["ELines"][[#, "Start"]] === i &&
        res["Alpha"][[i, ix["nD0"] + #]] === UnitVector[ix["nC"], ix["ExuberantCurveIndex"][#]] &];
    g = ArrayFlatten[{{form["IntersectionMatrix"], Transpose[be[[ix["nD0"] + es, b]]]},
        {be[[ix["nD0"] + es, b]], be[[ix["nD0"] + es, ix["ExuberantCurveIndex"] /@ es]]}}];
    If[es === {}, g = form["IntersectionMatrix"]];
    If[! MatrixQ[g, IntegerQ], Return[Nothing]];
    If[MatrixRank[g] <= form["PicardRank"] && Count[Eigenvalues[N[g]], x_ /; x > 10^-8] <= 1, Nothing,
     {lab[[i]], "Rank" -> MatrixRank[g], "PicardRank" -> form["PicardRank"],
      "PositiveEigenvalues" -> Count[Eigenvalues[N[g]], x_ /; x > 10^-8]}]];
  out["NumericalLattice"] = Table[gram[i], {i, ix["Ic"]}];
  out];

(* ================================================================== *)
(* 10.  Drawing the triangulation and the 7-brane cuts with the mouse  *)
(* ================================================================== *)
(*  gtpDraw[] asks for the size of the grid, opens it in a dialog and
    returns {edges, branes}, the two arguments of gtpIntersections.
    The editor is a set of pure transitions of a state
      <|"Range", "Mode", "Edges", "Branes", "Sel", "Down", "Draft",
        "Message", "History"|>
    wired to the mouse, so that its logic can be tested without a
    front end.  Edges are kept exactly as drawn: an edge through
    intermediate lattice points is a boundary edge with white dots,
    which are not vertices of the triangulation, and is shown so.
    A 7-brane is entered by drawing the two sides of its wedge from
    the tip; its charge is read off by gtpBraneFromSides.             *)

gtpDraw::cross = "The edges `1` cross, or one of them runs through a vertex of the other; the drawing is not a triangulation.";
gtpDraw::edges = "The initial edges must be a list of pairs of lattice points {{x1,y1},{x2,y2}}.";
gtpDraw::branes = "The initial 7\:2011branes must be a list of 7\:2011branes {{p,q},{{x1,y1},{x2,y2},...}}; a single brane must be wrapped in a list.";
gtpDraw::sides = "The sides `1` of the 7\:2011brane cuts are not edges of the triangulation, so they will not be glued.";
gtpDraw::draft = "An unfinished 7\:2011brane cut was discarded.";
gtpDraw::wedge = "The triangle(s) `1` lie in the generalized face of a 7\:2011brane, which must not be triangulated: remove the edges inside it.";
gtpBraneFromSides::nobrane = "`1`";

(* ---------- the 7-brane folding one side of a wedge onto the other ------
    M = 1 + k N_{p,q} acts on a segment d of the cut as
        M.d = d + k (p d1 + q d2) (q,-p),
    so the difference e - d with the corresponding segment e of the other
    side fixes the invariant direction (q,-p), hence the charge up to
    sign, and then k.  Exchanging the two sides turns k into -k, so
    exactly one of the two orders is a 7-brane: the sides may be given in
    either order, and the one that comes out as the cut is the side that
    cutImage maps onto the other.  k > 1 means k coincident 7-branes.   *)
gtpBraneFromSidesResult[s1_, s2_] := Module[{try, r1, r2},
  Which[
   ! MatchQ[s1, {{_Integer, _Integer} ..}] || ! MatchQ[s2, {{_Integer, _Integer} ..}],
    Return[{"Error", "Each side must be a list of lattice points {x,y}."}],
   Length[s1] < 2 || Length[s2] < 2,
    Return[{"Error", "Each side needs at least two points."}],
   First[s1] =!= First[s2],
    Return[{"Error", "The two sides must start at the same point, the tip of the wedge."}],
   Length[s1] =!= Length[s2],
    Return[{"Error", "The two sides must have the same number of segments, since the monodromy maps each segment of the cut onto a segment of the other side."}],
   MemberQ[Differences[s1], {0, 0}] || MemberQ[Differences[s2], {0, 0}],
    Return[{"Error", "A side repeats a point."}],
   s1 === s2,
    Return[{"Error", "The two sides coincide."}]];
  try[cut_, img_] := Module[{d = Differences[cut], v, i, u, p, q, dot, k, charge},
    v = Differences[img] - d;
    i = SelectFirst[Range[Length[v]], v[[#]] =!= {0, 0} &, None];
    If[i === None, Return[None]];
    u = v[[i]]/(GCD @@ v[[i]]);
    {p, q} = {-u[[2]], u[[1]]};
    dot = {p, q} . d[[i]];
    If[dot == 0, Return[None]];
    k = (v[[i]] . u/(u . u))/dot;
    If[! IntegerQ[k] || k < 1, Return[None]];
    If[p < 0 || (p == 0 && q < 0), {p, q} = -{p, q}];
    charge = If[k == 1, {p, q}, {p, q, k}];
    If[cutImage[{charge, cut}] =!= img, Return[None]];
    {charge, cut}];
  r1 = try[s1, s2]; r2 = try[s2, s1];
  Which[
   r1 =!= None, {"OK", r1},
   r2 =!= None, {"OK", r2},
   True, {"Error", "No 7\:2011brane monodromy maps one of the two sides onto the other."}]];

gtpBraneFromSides[s1_, s2_] := With[{r = gtpBraneFromSidesResult[s1, s2]},
  If[r[[1]] === "OK", r[[2]], Message[gtpBraneFromSides::nobrane, r[[2]]]; $Failed]];

gtpDrawChargeString[{p_, q_}] := "[" <> ToString[p] <> "," <> ToString[q] <> "]";
gtpDrawChargeString[{p_, q_, k_}] := ToString[k] <> "\[Times]" <> gtpDrawChargeString[{p, q}];

(* the two sides of the wedge of a brane : its cut, and the image *)
gtpDrawBraneSides[b_] := {braneCut[b], cutImage[{braneCharge[b], braneCut[b]}]};

(* ---------- state ------------------------------------------------------ *)

gtpDrawNormalize[edges_List] := DeleteDuplicatesBy[Select[edges, #[[1]] =!= #[[2]] &], Sort];

(* the lattice points in the interior of an edge : its white dots *)
gtpDrawWhiteDots[{p_, q_}] := With[{g = GCD @@ Abs[q - p]}, Table[p + k (q - p)/g, {k, 1, g - 1}]];

gtpDrawInit[r_, init_List] := gtpDrawInit[r, init, {}];
gtpDrawInit[{{x0_, x1_}, {y0_, y1_}}, init_List, branes_List] :=
 Module[{es = gtpDrawNormalize[init], pts, r = {{x0, x1}, {y0, y1}}},
  pts = Join[Flatten[es, 1], Flatten[Join @@@ (gtpDrawBraneSides /@ branes), 1]];
  If[pts =!= {},
   r = {{Min[x0, pts[[All, 1]]], Max[x1, pts[[All, 1]]]},
        {Min[y0, pts[[All, 2]]], Max[y1, pts[[All, 2]]]}}];
  <|"Range" -> r, "Mode" -> "Edges", "Edges" -> es, "Branes" -> branes, "Sel" -> None,
    "Down" -> None, "Draft" -> {}, "Message" -> None, "History" -> {}|>];

gtpDrawSnapshot[st_] := Append[st["History"], KeyTake[st, {"Edges", "Branes"}]];

gtpDrawSetMode[st_, m_] :=
  Join[st, <|"Mode" -> m, "Draft" -> {}, "Sel" -> None, "Down" -> None, "Message" -> None|>];

(* the lattice point under the mouse, if there is one *)
gtpDrawSnap[st_, pt_] := If[! VectorQ[pt, NumericQ], None,
  With[{p = Round[pt], r = st["Range"]},
   If[Norm[pt - p] <= 0.3 && r[[1, 1]] <= p[[1]] <= r[[1, 2]] && r[[2, 1]] <= p[[2]] <= r[[2, 2]],
    p, None]]];

gtpDrawSegDist[pt_, {a_, b_}] :=
  With[{t = Clip[(pt - a) . (b - a)/((b - a) . (b - a)), {0, 1}]}, Norm[pt - (a + t (b - a))]];

gtpDrawEdgeAt[st_, pt_] := If[! VectorQ[pt, NumericQ], None,
  SelectFirst[st["Edges"], gtpDrawSegDist[pt, #] < 0.15 &, None]];

(* add the edge [p,q], or remove it if it is already there *)
gtpDrawToggle[st_, p_, q_] := Module[{es = st["Edges"]},
  If[p === q, Return[st]];
  es = If[MemberQ[Sort /@ es, Sort[{p, q}]],
    Select[es, Sort[#] =!= Sort[{p, q}] &],
    Append[es, {p, q}]];
  Join[st, <|"Edges" -> es, "History" -> gtpDrawSnapshot[st]|>]];

gtpDrawDown[st_, pt_] := Join[st, <|"Down" -> gtpDrawSnap[st, pt]|>];

gtpDrawUp[st_, pt_] := If[st["Mode"] === "Cut", gtpDrawCutUp[st, pt], gtpDrawEdgeUp[st, pt]];

(*  edges : a drag between two lattice points toggles that edge; a click
    on a point selects it and a click on a second point toggles the edge
    between them; clicking the selected point again, or empty space,
    deselects; a click on an edge removes it.                          *)
gtpDrawEdgeUp[st0_, pt_] := Module[{p = st0["Down"], q = gtpDrawSnap[st0, pt], st, e},
  st = Join[st0, <|"Down" -> None|>];
  Which[
   p =!= None && q =!= None && p =!= q,
    Join[gtpDrawToggle[st, p, q], <|"Sel" -> None|>],
   p =!= None && q === p,
    Which[st["Sel"] === None, Join[st, <|"Sel" -> p|>],
          st["Sel"] === p, Join[st, <|"Sel" -> None|>],
          True, Join[gtpDrawToggle[st, st["Sel"], p], <|"Sel" -> None|>]],
   p === None && q === None && (e = gtpDrawEdgeAt[st, pt]) =!= None,
    Join[gtpDrawToggle[st, e[[1]], e[[2]]], <|"Sel" -> None|>],
   True, Join[st, <|"Sel" -> None|>]]];

(*  7-brane cuts : clicks on lattice points build the draft {side1} and
    then {side1, side2}.  The first click is the tip; clicking the tip
    again starts the second side, which begins there too.  As soon as
    the second side has as many points as the first, the charge is
    determined and the brane added -- or the reason it cannot be is
    shown, and the draft kept for Undo to correct.                     *)
gtpDrawCutUp[st0_, pt_] := Module[{p = st0["Down"], q = gtpDrawSnap[st0, pt], st},
  st = Join[st0, <|"Down" -> None|>];
  If[p === None || q =!= p, st, gtpDrawCutClick[st, p]]];

gtpDrawCutClick[st_, p_] := Module[{dr = st["Draft"], s1, s2, r},
  Switch[Length[dr],
   0, Join[st, <|"Draft" -> {{p}}, "Message" -> None|>],
   1, s1 = dr[[1]];
    Which[
     p === First[s1], If[Length[s1] >= 2, Join[st, <|"Draft" -> {s1, {p}}|>], st],
     MemberQ[s1, p], st,
     True, Join[st, <|"Draft" -> {Append[s1, p]}, "Message" -> None|>]],
   _, {s1, s2} = dr;
    If[MemberQ[s2, p], Return[st]];
    s2 = Append[s2, p];
    If[Length[s2] < Length[s1], Return[Join[st, <|"Draft" -> {s1, s2}, "Message" -> None|>]]];
    r = gtpBraneFromSidesResult[s1, s2];
    If[r[[1]] === "OK",
     Join[st, <|"Branes" -> Append[st["Branes"], r[[2]]], "Draft" -> {},
       "History" -> gtpDrawSnapshot[st],
       "Message" -> {"OK", "Added a " <> gtpDrawChargeString[braneCharge[r[[2]]]] <> " 7\:2011brane."}|>],
     Join[st, <|"Draft" -> {s1, s2}, "Message" -> {"Error", r[[2]] <> "  Undo removes the last point."}|>]]]];

gtpDrawRemoveBrane[st_, i_Integer] := If[1 <= i <= Length[st["Branes"]],
  Join[st, <|"Branes" -> Delete[st["Branes"], i], "History" -> gtpDrawSnapshot[st], "Message" -> None|>], st];

(* Undo first takes back points of an unfinished cut, then drawing steps *)
gtpDrawDraftBack[{s1_}] := If[Length[s1] <= 1, {}, {Most[s1]}];
gtpDrawDraftBack[{s1_, s2_}] := If[Length[s2] <= 1, {s1}, {s1, Most[s2]}];

gtpDrawUndo[st_] := Which[
  st["Draft"] =!= {}, Join[st, <|"Draft" -> gtpDrawDraftBack[st["Draft"]], "Message" -> None|>],
  st["History"] === {}, st,
  True, Join[st, Last[st["History"]], <|"History" -> Most[st["History"]], "Sel" -> None, "Message" -> None|>]];

gtpDrawClear[st_] := If[st["Edges"] === {} && st["Branes"] === {} && st["Draft"] === {}, st,
  Join[st, <|"Edges" -> {}, "Branes" -> {}, "Draft" -> {}, "History" -> gtpDrawSnapshot[st],
    "Sel" -> None, "Message" -> None|>]];

(* ---------- diagnostics -------------------------------------------------- *)

(* the faces triangulated so far : same criterion as triangleList, but
   through the adjacency so that it stays fast on large grids          *)
gtpDrawTriangles[edges_List] := Module[{nbr, verts, tris},
  If[edges === {}, Return[{}]];
  nbr = GroupBy[Join[edges, Reverse /@ edges], First -> Last];
  verts = DeleteDuplicates[Flatten[edges, 1]];
  tris = DeleteDuplicates[Sort /@ Flatten[Table[
       With[{a = e[[1]], b = e[[2]]},
        {a, b, #} & /@ Intersection[Lookup[nbr, Key[a], {}], Lookup[nbr, Key[b], {}]]],
       {e, edges}], 1]];
  Select[tris, Function[t, det2[t[[2]] - t[[1]], t[[3]] - t[[1]]] =!= 0 &&
     NoneTrue[Complement[verts, t], inClosedTriangle[t, #] &]]]];

(* two edges conflict if they cross at an interior point, or if an end
   of one lies in the interior of the other (a T-junction, which also
   covers collinear overlaps)                                          *)
gtpDrawInsideQ[v_, {a_, b_}] := det2[b - a, v - a] == 0 && 0 < (v - a) . (b - a) < (b - a) . (b - a);

gtpDrawCrossQ[{a_, b_}, {c_, d_}] :=
  (Intersection[{a, b}, {c, d}] === {} &&
    Sign[det2[b - a, c - a]] Sign[det2[b - a, d - a]] < 0 &&
    Sign[det2[d - c, a - c]] Sign[det2[d - c, b - c]] < 0) ||
  AnyTrue[{c, d}, gtpDrawInsideQ[#, {a, b}] &] || AnyTrue[{a, b}, gtpDrawInsideQ[#, {c, d}] &];

gtpDrawCrossings[edges_List] := Select[Subsets[edges, {2}], gtpDrawCrossQ @@ # &];

(* the segments of the sides of the wedges that are not edges : they
   would silently not be glued by gtpQuotient                          *)
gtpDrawMissingSides[edges_List, branes_List] := Module[{keys = Sort /@ edges},
  DeleteDuplicatesBy[Select[Join @@ (Partition[#, 2, 1] & /@ Join @@ (gtpDrawBraneSides /@ branes)),
    ! MemberQ[keys, Sort[#]] &], Sort]];

(* the triangles lying in the face of a 7-brane, which must not be there *)
gtpDrawWedgeTriangles[tris_List, branes_List] := Union @@ (braneWedgeTriangles[tris, #] & /@ branes);

(* ---------- display ------------------------------------------------------ *)

gtpDrawBraneColor = RGBColor[0.85, 0.42, 0.0];

gtpDrawWedge[{s1_, s2_}, label_] := Module[{tip = First[s1], b1, b2, bis},
  b1 = Normalize[s1[[2]] - tip]; b2 = Normalize[s2[[2]] - tip]; bis = b1 + b2;
  If[Norm[bis] < 0.1, bis = {-b1[[2]], b1[[1]]}];
  {{gtpDrawBraneColor, AbsoluteThickness[3], AbsoluteDashing[{9, 3, 2, 3}], Line[{s1, s2}]},
   {gtpDrawBraneColor, AbsolutePointSize[15], Point[tip]}, {White, AbsolutePointSize[6], Point[tip]},
   Text[Style[label, gtpDrawBraneColor, Bold, 12], tip + 0.45 Normalize[bis]]}];

gtpDrawPrimitives[st_] :=
 Module[{r = st["Range"], es = st["Edges"], tris, bad, lattice, miss, dr = st["Draft"]},
  tris = gtpDrawTriangles[es];
  bad = Union @@ gtpDrawCrossings[es];
  miss = gtpDrawMissingSides[es, st["Branes"]];
  lattice = Tuples[{Range @@ r[[1]], Range @@ r[[2]]}];
  {If[tris === {}, {}, {RGBColor[0.80, 0.87, 1.0], EdgeForm[], Polygon[tris]}],
   With[{wt = gtpDrawWedgeTriangles[tris, st["Branes"]]},
    If[wt === {}, {}, {RGBColor[1.0, 0.72, 0.72], EdgeForm[], Polygon[wt]}]],
   {GrayLevel[0.55], AbsolutePointSize[5], Point[lattice]},
   {Black, AbsoluteThickness[3], Line[Complement[es, bad]]},
   If[bad === {}, {}, {Red, AbsoluteThickness[3], Line[bad]}],
   If[miss === {}, {}, {Red, AbsoluteThickness[7], Opacity[0.35], Line[miss]}],
   If[es === {}, {}, {Black, AbsolutePointSize[9], Point[DeleteDuplicates[Flatten[es, 1]]]}],
   With[{w = Join @@ (gtpDrawWhiteDots /@ es)},
    If[w === {}, {}, {{Black, AbsolutePointSize[11], Point[w]}, {White, AbsolutePointSize[7], Point[w]}}]],
   gtpDrawWedge[gtpDrawBraneSides[#], gtpDrawChargeString[braneCharge[#]]] & /@ st["Branes"],
   If[dr === {}, {},
    {{gtpDrawBraneColor, AbsoluteThickness[4], Line[dr]},
     {gtpDrawBraneColor, AbsolutePointSize[11], Point[Join @@ dr]},
     {gtpDrawBraneColor, AbsolutePointSize[16], Point[dr[[1, 1]]]}}],
   If[st["Sel"] === None, {}, {Red, AbsolutePointSize[14], Point[st["Sel"]]}]}];

(* the rubber band from the pending point to the mouse, kept apart so
   that moving the mouse does not recompute the rest of the picture    *)
gtpDrawBand[st_, mouse_] := Which[
  ! VectorQ[mouse, NumericQ], {},
  st["Mode"] =!= "Cut" && st["Down"] =!= None,
   {Red, Dashed, AbsoluteThickness[2], Line[{st["Down"], mouse}]},
  st["Mode"] === "Cut" && st["Draft"] =!= {},
   {gtpDrawBraneColor, Dashed, AbsoluteThickness[2], Line[{Last[Last[st["Draft"]]], mouse}]},
  True, {}];

gtpDrawText[s_String, opts___] := TextCell[s, "Text", FontSize -> 11, opts];

gtpDrawInstructions[st_] := gtpDrawText[If[st["Mode"] === "Cut",
   Switch[Length[st["Draft"]],
    0, "Click the tip of the wedge, then the lattice points along one of its sides.",
    1, "Continue along the first side; when it is complete, click the tip again to start the other side.",
    _, "Now click the points along the other side; the 7\:2011brane is added as soon as both sides have the same length."],
   "Drag from a lattice point to another, or click the two points, to add an edge; do it again to remove it.  Clicking on an edge also removes it."],
  FontColor -> GrayLevel[0.3]];

gtpDrawStatus[st_] := Module[{es = st["Edges"], bad, main, miss, msg = st["Message"]},
  bad = gtpDrawCrossings[es];
  miss = gtpDrawMissingSides[es, st["Branes"]];
  main = ToString[Length[es]] <> " edges, " <> ToString[Length[gtpDrawTriangles[es]]] <> " triangles, " <>
    ToString[Length[st["Branes"]]] <> If[Length[st["Branes"]] == 1, " 7\:2011brane", " 7\:2011branes"] <>
    If[st["Sel"] === None, "", ",   selected " <> ToString[st["Sel"]]];
  Column[{
    Row[{gtpDrawText[main],
      If[bad === {}, "",
       gtpDrawText["      " <> ToString[Length[bad]] <> If[Length[bad] == 1, " pair", " pairs"] <>
         " of conflicting edges, shown in red", FontColor -> Red]]}],
    If[miss === {}, Nothing,
     gtpDrawText["Sides of 7\:2011brane cuts that are not edges, shown in red: " <> ToString[miss], FontColor -> Red]],
    With[{wt = gtpDrawWedgeTriangles[gtpDrawTriangles[es], st["Branes"]]},
     If[wt === {}, Nothing,
      gtpDrawText[ToString[Length[wt]] <> If[Length[wt] == 1, " triangle lies", " triangles lie"] <>
        " in the face of a 7\:2011brane, shown in red; that face must not be triangulated, so remove the edges inside it.",
       FontColor -> Red]]],
    If[ListQ[msg], gtpDrawText[msg[[2]], FontColor -> If[msg[[1]] === "OK", RGBColor[0, 0.5, 0], Red]], Nothing]}]];

gtpDrawPanel[st0_Association] :=
 With[{r = st0["Range"], w = Clip[70 (st0["Range"][[1, 2]] - st0["Range"][[1, 1]] + 1), {320, 700}]},
  DynamicModule[{gtpDrawState = st0, gtpDrawMouse = None},
   Column[{
     Row[{gtpDrawText["Draw   "],
       SetterBar[Dynamic[gtpDrawState["Mode"], (gtpDrawState = gtpDrawSetMode[gtpDrawState, #]) &],
        {"Edges" -> gtpDrawText["edges"], "Cut" -> gtpDrawText["7\:2011brane cuts"]}]}],
     Pane[Dynamic[gtpDrawInstructions[gtpDrawState]], ImageSize -> w],
     EventHandler[
      Deploy@Graphics[{Dynamic[gtpDrawPrimitives[gtpDrawState]], Dynamic[gtpDrawBand[gtpDrawState, gtpDrawMouse]]},
       PlotRange -> {r[[1]] + {-1/2, 1/2}, r[[2]] + {-1/2, 1/2}},
       Frame -> True, FrameTicks -> {{Range @@ r[[2]], None}, {Range @@ r[[1]], None}},
       GridLines -> {Range @@ r[[1]], Range @@ r[[2]]}, GridLinesStyle -> GrayLevel[0.9],
       AspectRatio -> Automatic, ImageSize -> w],
      {"MouseDown" :> (gtpDrawState = gtpDrawDown[gtpDrawState, MousePosition["Graphics"]]),
       "MouseDragged" :> (gtpDrawMouse = MousePosition["Graphics"]),
       "MouseMoved" :> (gtpDrawMouse = MousePosition["Graphics"]),
       "MouseUp" :> (gtpDrawState = gtpDrawUp[gtpDrawState, MousePosition["Graphics"]]; gtpDrawMouse = None)}],
     Pane[Dynamic[gtpDrawStatus[gtpDrawState]], ImageSize -> w],
     Dynamic[If[gtpDrawState["Branes"] === {}, "",
       Row[Prepend[Table[With[{i = i},
          Button[gtpDrawText[gtpDrawChargeString[braneCharge[gtpDrawState["Branes"][[i]]]] <> " at " <>
             ToString[First[braneCut[gtpDrawState["Branes"][[i]]]]] <> "  \[Times]"],
           gtpDrawState = gtpDrawRemoveBrane[gtpDrawState, i], Appearance -> "Palette"]],
          {i, Length[gtpDrawState["Branes"]]}], gtpDrawText["Remove a 7\:2011brane:  "]], " "]]],
     Row[{Button["Undo", gtpDrawState = gtpDrawUndo[gtpDrawState]],
       Button["Clear", gtpDrawState = gtpDrawClear[gtpDrawState]],
       Spacer[40], CancelButton[DialogReturn[$Canceled]],
       DefaultButton["Done", DialogReturn[{gtpDrawState["Edges"], gtpDrawState["Branes"],
          gtpDrawState["Draft"] =!= {}}]]}, " "]},
    Spacings -> 1],
   SaveDefinitions -> True]];

(* ---------- entry points ------------------------------------------------- *)

gtpDraw[] := Module[{n},
  n = DialogInput[DynamicModule[{k = 5},
     Column[{"Size of the grid (number of lattice points on each side):",
       InputField[Dynamic[k], Number, FieldSize -> 6],
       Row[{CancelButton[DialogReturn[$Canceled]], DefaultButton[DialogReturn[k]]}, " "]}]],
    WindowTitle -> "GTP grid size"];
  If[IntegerQ[n] && n >= 2, gtpDraw[n], $Canceled]];

gtpDraw[n_Integer?(# >= 2 &), rest___] := gtpDraw[{{0, n - 1}, {0, n - 1}}, rest];

(* Cancel returns what one started from, so that
   {gtp, branes} = gtpDraw[n, gtp, branes] leaves both untouched; from an
   empty grid it returns $Canceled                                     *)
gtpDraw[r : {{_Integer, _Integer}, {_Integer, _Integer}}, init_List : {}, branes_List : {}] :=
 Module[{res, es, bs, left, bad, miss},
  If[! MatchQ[init, {{{_Integer, _Integer}, {_Integer, _Integer}} ...}],
   Message[gtpDraw::edges]; Return[$Failed]];
  If[! AllTrue[branes, braneQ], Message[gtpDraw::branes]; Return[$Failed]];
  res = DialogInput[gtpDrawPanel[gtpDrawInit[Sort /@ r, init, branes]],
    WindowTitle -> "Draw the triangulated GTP"];
  If[res === $Canceled, Return[If[init === {} && branes === {}, $Canceled, {init, branes}]]];
  {es, bs, left} = res;
  bad = gtpDrawCrossings[es]; If[bad =!= {}, Message[gtpDraw::cross, bad]];
  miss = gtpDrawMissingSides[es, bs]; If[miss =!= {}, Message[gtpDraw::sides, miss]];
  With[{wt = gtpDrawWedgeTriangles[gtpDrawTriangles[es], bs]}, If[wt =!= {}, Message[gtpDraw::wedge, wt]]];
  If[left, Message[gtpDraw::draft]];
  {es, bs}];


(* ================================================================== *)
(* 11.  The whole computation in one command                          *)
(* ================================================================== *)
(*  gtpRun draws (or takes) the GTP, runs gtpIntersections and prints
    every step: the input, the picture, the report, the intersection
    matrices and the TikZ code.                                        *)

Options[gtpRun] = {"Report" -> True, "PlotOptions" -> {}, "TikZOptions" -> {},
   "TikZFile" -> None, "Presentations" -> "All", "CurveNames" -> "Numbered", "BypassSide" -> Automatic};

gtpRunEdgesQ[x_] := MatchQ[x, {{{_Integer, _Integer}, {_Integer, _Integer}} ..}];

gtpRun[opts : OptionsPattern[]] := gtpRunDrawn[gtpDraw[], opts];
gtpRun[edges_?gtpRunEdgesQ, branes_List, opts : OptionsPattern[]] := gtpRunCompute[edges, branes, opts];
gtpRun[r : (_Integer | {{_Integer, _Integer}, {_Integer, _Integer}}), init_List : {}, branes_List : {},
   opts : OptionsPattern[]] := gtpRunDrawn[gtpDraw[r, init, branes], opts];

(* gtpDraw returns $Canceled when the editor is closed with Cancel *)
gtpRunDrawn[{edges_List, branes_List}, opts___] := gtpRunCompute[edges, branes, opts];
gtpRunDrawn[other_, ___] := other;

gtpRunHeading[s_String] := Print[Style[s, Bold, 15, RGBColor[0.12, 0.29, 0.55]]];

gtpRunCompute[edges_List, branes_List, opts : OptionsPattern[gtpRun]] :=
 Module[{res, pic, code, standalone, file, dir},
  res = gtpIntersections[edges, branes, "Presentations" -> OptionValue[gtpRun, {opts}, "Presentations"],
    "CurveNames" -> OptionValue[gtpRun, {opts}, "CurveNames"], "BypassSide" -> OptionValue[gtpRun, {opts}, "BypassSide"]];
  gtpRunHeading["GTP"];
  Print[Row[{Style["input:  ", Gray],
     "gtpRun[" <> ToString[edges, InputForm] <> ", " <> ToString[branes, InputForm] <> "]"}]];
  pic = Quiet@Check[gtpPlot[res, Sequence @@ OptionValue[gtpRun, {opts}, "PlotOptions"]], $Failed];
  Print[If[pic === $Failed, Style["The GTP could not be drawn.", Red], pic]];
  If[TrueQ[OptionValue[gtpRun, {opts}, "Report"]],
   gtpRunHeading["Computation"]; gtpReport[res, "IntersectionTables" -> False],
   If[! gtpSolvedQ[res], Print[gtpNotSolved[res]]]];
  If[gtpSolvedQ[res],
   gtpRunHeading["Triple intersection numbers"];
   Print[If[res["IndexSets"]["Ic"] === {}, "There is no compact divisor.", gtpIntersectionTables[res]]]];
  code = Quiet@Check[gtpTikZ[res, Sequence @@ OptionValue[gtpRun, {opts}, "TikZOptions"]], $Failed];
  If[StringQ[code],
   standalone = Quiet@Check[gtpTikZ[res, Sequence @@ OptionValue[gtpRun, {opts}, "TikZOptions"],
       "Standalone" -> True], code];
   gtpRunHeading["TikZ"];
   (* the strings are inlined: the buttons outlive the local variables *)
   With[{c = code, sa = standalone},
    Print[Row[{Button["Copy the TikZ code", CopyToClipboard[c]],
       Button["Copy as a standalone document", CopyToClipboard[sa]]}, "  "]]];
   file = OptionValue[gtpRun, {opts}, "TikZFile"];
   If[StringQ[file],
    dir = Quiet@Check[NotebookDirectory[], Directory[]];
    If[! StringQ[dir], dir = Directory[]];
    file = If[StringStartsQ[file, "/" | "~" | RegularExpression["[A-Za-z]:"]],
      ExpandFileName[file], FileNameJoin[{dir, file}]];
    Export[file, code, "Text"];
    Print[Style["written to " <> file, Gray]]];
   CellPrint[Cell[code, "Program"]]];
  Append[res, "Input" -> {edges, branes}]];

End[];
EndPackage[];
