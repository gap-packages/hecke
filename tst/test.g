## Specht documentation test
H:=Specht(4);
RInducedModule(MakePIM(H,12,2));
#############################
H:=Specht(3,3);
d:=DecompositionMatrix(H,5);
for n in [6..14] do d:=InducedDecompositionMatrix(d); SaveDecompositionMatrix(d);  od;
#############################
H:=Specht(5);
MakeSimple(H,3,2,1);
MakeSpecht(last);
RInducedModule(MakePIM(H,3,2,1));
MakeSpecht(last);
MakeSimple(H,3,1) * MakeSimple(H,3);
RRestrictedModule(last);
MakeSpecht(last);
MakePIM(last);
#############################
MakePIM(H,4,3,2);
MakeSimple(MakePIM(H,4,3,2));
MakeSpecht(MakeSimple(MakeSpecht(H,1,1,1,1,1)));
#############################
H:=Specht(3,3);
d:=InducedDecompositionMatrix(DecompositionMatrix(H,14));
SaveDecompositionMatrix(d);
MakePIM(d,4,3,3,2,2,1);
MakeSpecht(d,7,3,3,2);
MakeSimple(d,14,1);
MakeSpecht(d, MakeSimple(d,10,5) );
#############################
H:=Specht(5,5);; SimpleDimension(H,6);
#############################
val:=function(x) local v;
x:=Sum([0..x-1],v->4^v);  # x-${>}$[x]\_q
v:=0; while x mod 5=0 do x:=x/5; v:=v+1; od;
return v;
end;;
H:=Specht(2,5,val,"e2q4");
#############################
H:=Specht(4);
MakeFockPIM(H,6,2);
RRestrictedModule(last);
MakePIM(last);
Specialized(last);
MakeFockSpecht(H,5,3,2);
RInducedModule(last,0);
#############################
DecompositionMatrix(Specht(3),6,LengthLexicographic);
#############################
CrystalDecompositionMatrix(Specht(3), 6);
Specialized(last);
#############################
H:=Specht(6);;
DecompositionNumber(H,[6,4,2],[6,6]);
#############################
H:=Specht(2,2);
RInducedModule(MakeSpecht(H,7,4,3,1));
RInducedModule(MakePIM(H,5,3,1));
RInducedModule(MakeSimple(H,11,2,1));
#############################
H:=Specht(4);
RInducedModule(MakeSpecht(H,5,2,1));
RInducedModule(MakeSpecht(H,5,2,1),0);
RInducedModule(MakeSpecht(H,5,2,1),1);
RInducedModule(MakeSpecht(H,5,2,1),2);
RInducedModule(MakeSpecht(H,5,2,1),3);
EResidueDiagram(H,5,2,1);
#############################
H:=Specht(3);
RInducedModule(MakeFockPIM(H,4,2),1,2);
MakePIM(last);
#############################
H:=Specht(4);
SInducedModule(MakePIM(H,5,2,1),3);
SInducedModule(MakePIM(H,5,2,1),3,1);
RInducedModule(MakePIM(H,5,2,1),1,1,1);
#############################
H:=Specht(6);
RRestrictedModule(MakePIM(H,5,3,2,1),4);
RRestrictedModule(MakeSimple(H,5,3,2),1);
#############################
H:=Specht(6);
SRestrictedModule(MakeSpecht(H,4,3,2),3);
SRestrictedModule(MakePIM(H,5,4,1),2,4);
#############################
d:=DecompositionMatrix(Specht(3,3),14);
InducedDecompositionMatrix(d);
#############################
H:=Specht(2,2);
d:=InducedDecompositionMatrix(DecompositionMatrix(H,9));
x:=RInducedModule(MakePIM(H,9),1);
IsNewIndecomposable(d,x);
x:=x-MakePIM(d,6,3,1);
IsNewIndecomposable(d,x,6,3,1);
AddIndecomposable(d,x);
SaveDecompositionMatrix(d);
#############################
H:=Specht(4);
d:=CrystalDecompositionMatrix(H,5);
InvertDecompositionMatrix(d);
#############################
H:=Specht(2);; Hp:=Specht(2,2);;
d:=DecompositionMatrix(H,13);; dp:=DecompositionMatrix(Hp,13);;
a:=AdjustmentMatrix(dp,d);
MatrixDecompositionMatrix(dp)=MatrixDecompositionMatrix(d)*MatrixDecompositionMatrix(a);
#############################
H:=Specht(5,5);;
d:=DecompositionMatrix(H,9);;
for r in [10..20] do
d:=InducedDecompositionMatrix(d);
SaveDecompositionMatrix(d);
od;
#############################
H:=Specht(2,2);
d:=DecompositionMatrix(H,15);
d:=CalculateDecompositionMatrix(H,15);
MissingIndecomposables(d);
#############################
MatrixDecompositionMatrix(DecompositionMatrix(Specht(3),5));
#############################
H:=Specht(3);
m:=[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 1, 0, 1, 0 ], [ 0, 0, 0, 1 ], [ 0, 0, 1, 0 ] ];
DecompositionMatrixMatrix(H,m,4);
#############################
H:=Specht(6);
SimpleDimension(H,11,3);
d:=DecompositionMatrix(H,5);
SimpleDimension(d,3,2);
SimpleDimension(d);
#############################
SpechtDimension(6,3,2,1);
#############################
H:=Specht(2);
Schaper(H,9,5,3,2,1);
Schaper(H,9,6,5,2);
#############################
H:=Specht(3);
IsSimpleModule(H,45,31,24);
#############################
MullineuxMap(Specht(2),12,5,2);
MullineuxMap(Specht(4),12,5,2);
MullineuxMap(Specht(6),12,5,2);
MullineuxMap(Specht(8),12,5,2);
MullineuxMap(Specht(10),12,5,2);
#############################
MullineuxSymbol(5,[8,6,5,5]);
#############################
PartitionMullineuxSymbol(5, MullineuxSymbol(5,[8,6,5,5]) );
#############################
GoodNodes(5,[5,4,3,2]);
GoodNodes(5,[5,4,3,2],0);
GoodNodes(5,[5,4,3,2],4);
#############################
NormalNodes(5,[6,5,4,4,3,2,1,1,1]);
NormalNodes(5,[6,5,4,4,3,2,1,1,1],0);
#############################
H:=Specht(4);
GoodNodeSequence(H,4,3,1);
GoodNodeSequence(H,4,3,2);
GoodNodeSequence(H,4,4,2);
GoodNodeSequence(H,5,4,2);
GoodNodeSequences(H,5,2,1);
#############################
H:=Specht(4);
PartitionGoodNodeSequence(H,0, 3, 1, 0, 2, 2, 1, 3, 3, 2);
#############################
GoodNodeLatticePath(3,3,2,1);
GoodNodeLatticePaths(3,3,2,1);
GoodNodeSequence(4,6,3,2);
LatticePathGoodNodeSequence(4,last);
#############################
H:=Specht(0);
LittlewoodRichardsonRule([3,2,1],[4,2]);
MakeSpecht(H,3,2,1)*MakeSpecht(H,4,2);
LittlewoodRichardsonCoefficient([3,2,1],[4,2],[5,4,2,1]);
#############################
InverseLittlewoodRichardsonRule([3,2,1]);
#############################
H:=Specht(2);
EResidueDiagram(MakeSpecht(MakePIM(H,7,5)));
#############################
RemoveRimHook([6,5,4],1,2);
RemoveRimHook([6,5,4],2,3);
HookLengthDiagram(6,5,4);
#############################
AddRimHook([6,4,3],1,3);
AddRimHook([6,4,3],2,3);
AddRimHook([6,4,3],3,3);
AddRimHook([6,4,3],4,3);
AddRimHook([6,4,3],5,3);
#############################
H:=Specht(6);
ECore(H,16,8,6,5,3,1);
#############################
H:=Specht(8);
EQuotient(H,22,18,16,12,12,1,1);
#############################
H:=Specht(11);
mu:=[100,98,57,43,12,1];
Q:=EQuotient(H,mu);
C:=ECore(H,mu);
CombineEQuotientECore(H,Q,C);
#############################
EWeight(6,[16,8,6,5,3,1]);
#############################
H:=Specht(3);
ERegularPartitions(H,6);
#############################
ConjugatePartition(6,4,3,2);
#############################
BetaSet([5,4,2,2]);
PartitionBetaSet([ 2, 3, 6, 8 ]);
#############################
H:=Specht(4);
ETopLadder(H,1,1,1,1,1,1,1,1,1,1);
ETopLadder(6,1,1,1,1,1,1,1,1,1,1);
#############################
Dominates([5,4],[4,4,1]);
p:=Partitions(6);;Sort(p,LengthLexicographic); p;
p:=Partitions(6);;Sort(p,Lexicographic); p;
p:=Partitions(6);;Sort(p,ReverseDominance); p;
#############################
H:=Specht(2);
x:=MakeFockPIM(H,6,2);
Specialized(x);
Specialized(x,2);
#############################
H:=Specht(8);
x:=MakeSpecht(RInducedModule(MakePIM(H,8,5,3)) );
ERegulars(x);
MakePIM(x);
#############################
H:=Specht(2);
SplitECores(RInducedModule(MakeSpecht(H,5,3,1)));
RInducedModule(MakeSpecht(H,5,3,1),0);
RInducedModule(MakeSpecht(H,5,3,1),1);
#############################
H:=Specht(3);
x:=MakeSpecht(MakePIM(H,7,3));
Coefficient(x,5,2,2,1);
#############################
H:=Specht(4);
InnerProduct(MakeSpecht(H,2,2,2,1), MakePIM(H,4,3));
DecompositionNumber(H,[2,2,2,1],[4,3]);
#############################
H:=Specht(2);;
x:=MakeSpecht(MakePIM(H,6));
Display(x);
Print(x,"\n");
#############################
SemiStandardTableaux([4,3],1,1,1,2,2);
StandardTableaux(4,2);
ConjugateTableau(Tableau([ [ 1, 3, 5, 6 ], [ 2, 4 ] ]));
ShapeTableau(Tableau( [ [ 1, 1, 2, 3 ], [ 4, 5 ] ] ));
List(SemiStandardTableaux([5,4,2],[4,3,0,1,3]),TypeTableau);

