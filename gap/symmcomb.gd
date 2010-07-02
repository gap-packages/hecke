
DeclareOperation("Lexicographic", [IsList,IsList]);
DeclareOperation("LengthLexicographic", [IsList,IsList]);
DeclareOperation("ReverseDominance", [IsList,IsList]);
DeclareOperation("Dominance", [IsList,IsList]);

DeclareOperation("ConjugatePartition", [IsList]);

DeclareOperation("LittlewoodRichardsonRule", [IsList,IsList]);
DeclareOperation("LittlewoodRichardsonCoefficient", [IsList,IsList,IsList]);
DeclareOperation("InverseLittlewoodRichardsonRule", [IsList]);

DeclareOperation("SpechtDimension", [IsHeckeSpecht]);
DeclareOperation("SpechtDimension", [IsList]);

DeclareOperation("BetaNumbers", [IsList]);
## ALREADY AVAILABLE IN GAP4
## DeclareOperation("BetaSet", [IsList]);
DeclareOperation("PartitionBetaSet", [IsList]);

DeclareOperation("EAbacusRunners", [IsInt,IsList]);
DeclareOperation("ECore",[IsInt,IsList]);
DeclareOperation("ECore",[IsAlgebraObj,IsList]);
DeclareOperation("IsECore",[IsInt,IsList]);
DeclareOperation("IsECore",[IsAlgebraObj,IsList]);
DeclareOperation("EWeight",[IsInt,IsList]);
DeclareOperation("EWeight",[IsAlgebraObj,IsList]);
DeclareOperation("EQuotient",[IsInt,IsList]);
DeclareOperation("EQuotient",[IsAlgebraObj,IsList]);
DeclareOperation("EAbacus",[IsInt,IsList]);
DeclareOperation("EAbacus",[IsAlgebraObj,IsList]);
DeclareOperation("CombineEQuotientECore",[IsInt,IsList,IsList]);
DeclareOperation("CombineEQuotientECore",[IsAlgebraObj,IsList,IsList]);
DeclareOperation("IsERegular",[IsInt,IsList]);
DeclareOperation("IsERegular",[IsAlgebraObj,IsList]);
DeclareOperation("ERegularPartitions",[IsInt,IsInt]);
DeclareOperation("ERegularPartitions",[IsAlgebraObj,IsInt]);

DeclareOperation("EResidueDiagram",[IsInt,IsList]);
DeclareOperation("EResidueDiagram",[IsAlgebraObj,IsList]);
DeclareOperation("EResidueDiagram",[IsHeckeSpecht]);

DeclareOperation("ETopLadder",[IsInt,IsList]);
DeclareOperation("ETopLadder",[IsAlgebraObj,IsList]);
DeclareOperation("EHookDiagram",[IsInt,IsList]);
DeclareOperation("EHookDiagram",[IsAlgebraObj,IsList]);

DeclareOperation("HookLengthDiagram",[IsList]);

DeclareOperation("NormalNodes",[IsInt,IsList]);
DeclareOperation("NormalNodes",[IsAlgebraObj,IsList]);
DeclareOperation("NormalNodes",[IsInt,IsList,IsInt]);
DeclareOperation("NormalNodes",[IsAlgebraObj,IsList,IsInt]);
DeclareOperation("RemoveNormalNodes",[IsInt,IsList,IsInt]);
DeclareOperation("RemoveNormalNodes",[IsAlgebraObj,IsList,IsInt]);
DeclareOperation("GoodNodes",[IsInt,IsList]);
DeclareOperation("GoodNodes",[IsAlgebraObj,IsList]);
DeclareOperation("GoodNodes",[IsInt,IsList,IsInt]);
DeclareOperation("GoodNodes",[IsAlgebraObj,IsList,IsInt]);
DeclareOperation("GoodNodeSequence",[IsInt,IsList]);
DeclareOperation("GoodNodeSequence",[IsAlgebraObj,IsList]);
DeclareOperation("GoodNodeSequences",[IsInt,IsList]);
DeclareOperation("GoodNodeSequences",[IsAlgebraObj,IsList]);
DeclareOperation("PartitionGoodNodeSequence",[IsInt,IsList]);
DeclareOperation("PartitionGoodNodeSequence",[IsAlgebraObj,IsList]);
DeclareOperation("GoodNodeLatticePath",[IsInt,IsList]);
DeclareOperation("GoodNodeLatticePath",[IsAlgebraObj,IsList]);
DeclareOperation("GoodNodeLatticePaths",[IsInt,IsList]);
DeclareOperation("GoodNodeLatticePaths",[IsAlgebraObj,IsList]);
DeclareOperation("LatticePathGoodNodeSequence",[IsInt,IsList]);
DeclareOperation("LatticePathGoodNodeSequence",[IsAlgebraObj,IsList]);

DeclareOperation("MullineuxSymbol",[IsInt,IsList]);
DeclareOperation("MullineuxSymbol",[IsAlgebraObj,IsList]);
DeclareOperation("PartitionMullineuxSymbol",[IsInt,IsList]);
DeclareOperation("PartitionMullineuxSymbol",[IsAlgebraObj,IsList]);

DeclareOperation("RemoveRimHook",[IsList,IsInt,IsInt,IsList]);
DeclareOperation("RemoveRimHook",[IsList,IsInt,IsInt]);
DeclareOperation("AddRimHook",[IsList,IsInt,IsInt]);

