#######################################################################
##  Hecke - symmcomb.gd : Combinatorial functions on partitions.     ##
##                                                                   ##
##     This file contains most of the combinatorial functions used   ##
##     by Hecke. Most are standard operations on Young diagrams or   ##
##     partitions.                                                   ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-package SPECHT 2.4 by Andrew Mathas)  ##
##                                                                   ##
#######################################################################

## Hecke 1.0: June 2010:
##   - initial

MakeDispatcherFunc("Lexicographic", [[IsList]],[2],[2]);
MakeDispatcherFunc("LengthLexicographic", [[IsList]],[2],[2]);
MakeDispatcherFunc("ReverseDominance", [[IsList]],[2],[2]);
MakeDispatcherFunc("Dominates", [[IsList]],[2],[2]);

MakeDispatcherFunc("ConjugatePartition", [[]],[1],[1]);

MakeDispatcherFunc("LittlewoodRichardsonRule", [[IsList]],[2],[2]);
MakeDispatcherFunc("LittlewoodRichardsonCoefficient", [[IsList,IsList]],[3],[3]);
MakeDispatcherFunc("InverseLittlewoodRichardsonRule", [[]],[1],[1]);

MakeDispatcherFunc("SpechtDimension", [[],[IsHeckeSpecht]],[1,0],[1,1]);

MakeDispatcherFunc("BetaNumbers", [[]],[1],[1]);
MakeDispatcherFunc("PartitionBetaSet", [[]],[1],[1]);

MakeDispatcherFunc("EAbacusRunners", [[IsInt]],[2],[2]);
MakeDispatcherFunc("ECore",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("IsECore",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("EWeight",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("EQuotient",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("EAbacus",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("CombineEQuotientECore",[[IsInt,IsList],[IsAlgebraObj,IsList]],[3,3],[3,3]);
MakeDispatcherFunc("IsERegular",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("ERegularPartitions",[[IsInt,IsInt],[IsAlgebraObj,IsInt]],[0,0],[2,2]);

MakeDispatcherFunc("EResidueDiagram",[[IsInt],[IsAlgebraObj],[IsHeckeSpecht]],[2,2,0],[2,2,1]);

MakeDispatcherFunc("ETopLadder",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("EHookDiagram",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);

MakeDispatcherFunc("HookLengthDiagram",[[]],[1],[1]);

MakeDispatcherFunc("NormalNodes",[[IsInt],[IsAlgebraObj],[IsInt,IsInt],[IsAlgebraObj,IsInt]],[2,2,2,2],[2,2,3,3]);
MakeDispatcherFunc("RemoveNormalNodes",[[IsInt,IsInt],[IsAlgebraObj,IsInt]],[2,2],[3,3]);
MakeDispatcherFunc("GoodNodes",[[IsInt],[IsAlgebraObj],[IsInt,IsInt],[IsAlgebraObj,IsInt]],[2,2,2,2],[2,2,3,3]);
MakeDispatcherFunc("GoodNodeSequence",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("GoodNodeSequences",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("PartitionGoodNodeSequence",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("GoodNodeLatticePath",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("GoodNodeLatticePaths",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("LatticePathGoodNodeSequence",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);

MakeDispatcherFunc("MullineuxSymbol",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);
MakeDispatcherFunc("PartitionMullineuxSymbol",[[IsInt],[IsAlgebraObj]],[2,2],[2,2]);

MakeDispatcherFunc("RemoveRimHook",[[IsInt,IsInt],[IsList,IsInt,IsInt]],[1,4],[3,4]);
MakeDispatcherFunc("AddRimHook",[[IsInt,IsInt]],[1],[3]);

