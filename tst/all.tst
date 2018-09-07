## Specht documentation test
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> RInducedModule(MakePIM(H,12,2));
<direct sum of 5 P-modules>

## Specht documentation test
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> RInducedModule(MakePIM(H,12,2));
<direct sum of 5 P-modules>

#############################
gap> H:=Specht(3,3);
<Hecke algebra with e = 3>
gap> d:=DecompositionMatrix(H,5);
<7x5 decomposition matrix>
gap> for n in [6..14] do d:=InducedDecompositionMatrix(d); od;
# Inducing..
# Inducing..
# Inducing...
# Inducing...
# Inducing...
# Inducing....

#############################
gap> H:=Specht(5);
<Hecke algebra with e = 5>
gap> MakeSimple(H,3,2,1);
<direct sum of 1 D-modules>
gap> MakeSpecht(last);
<direct sum of 3 S-modules>
gap> RInducedModule(MakePIM(H,3,2,1));
<direct sum of 4 P-modules>
gap> MakeSpecht(last);
<direct sum of 6 S-modules>
gap> MakeSimple(H,3,1) * MakeSimple(H,3);
<direct sum of 7 D-modules>
gap> RRestrictedModule(last);
<direct sum of 6 D-modules>
gap> MakeSpecht(last);
<direct sum of 6 S-modules>
gap> MakePIM(last);
<direct sum of 5 P-modules>

#############################
gap> MakePIM(H,4,3,2);
<direct sum of 1 P-modules>
gap> MakeSimple(MakePIM(H,4,3,2));
<direct sum of 3 D-modules>
gap> MakeSpecht(MakeSimple(MakeSpecht(H,1,1,1,1,1)));
<direct sum of 4 S-modules>

#############################
gap> H:=Specht(3,3);
<Hecke algebra with e = 3>
gap> d:=InducedDecompositionMatrix(DecompositionMatrix(H,14));
# Inducing....
The following projectives are missing from <d>:
    [ 15 ]  [ 8, 7 ]
<176x70 decomposition matrix>
gap> MakePIM(d,4,3,3,2,2,1);
<direct sum of 4 S-modules>
gap> MakeSpecht(d,7,3,3,2);
<direct sum of 6 D-modules>
gap> MakeSimple(d,14,1);
fail
gap> MakeSpecht(d, MakeSimple(d,10,5) );
<direct sum of 2 S-modules>

#############################
gap> H:=Specht(5,5);; SimpleDimension(H,6);
6       : 1
5,1     : 5
4,2     : 8
4,1^2   : 10
3^2     : 5
3,2,1   : 8
3,1^3   : 10
2^3     : 5
2^2,1^2 : 1
2,1^4   : 5
true

#############################
gap> val:=function(x) local v;
> x:=Sum([0..x-1],v->4^v);  # x-${>}$[x]\_q
> v:=0; while x mod 5=0 do x:=x/5; v:=v+1; od;
> return v;
> end;;
gap> H:=Specht(2,5,val,"e2q4");
<Hecke algebra with e = 2>

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> MakeFockPIM(H,6,2);
<direct sum of 2 Sq-modules>
gap> RRestrictedModule(last);
<direct sum of 3 Sq-modules>
gap> MakePIM(last);
<direct sum of 2 Pq-modules>
gap> Specialized(last);
<direct sum of 2 P-modules>
gap> MakeFockSpecht(H,5,3,2);
<direct sum of 1 Sq-modules>
gap> RInducedModule(last,0);
<direct sum of 1 Sq-modules>

#############################
gap> DecompositionMatrix(Specht(3),6,LengthLexicographic);
<11x7 decomposition matrix>

#############################
gap> CrystalDecompositionMatrix(Specht(3), 6);
<11x7 decomposition matrix>
gap> Specialized(last);
<11x7 decomposition matrix>

#############################
gap> H:=Specht(6);;
gap> DecompositionNumber(H,[6,4,2],[6,6]);
0

#############################
gap> H:=Specht(2,2);
<Hecke algebra with e = 2>
gap> RInducedModule(MakeSpecht(H,7,4,3,1));
<direct sum of 5 S-modules>
gap> RInducedModule(MakePIM(H,5,3,1));
<direct sum of 3 P-modules>
gap> RInducedModule(MakeSimple(H,11,2,1));
# D(<x>), unable to rewrite <x> as a sum of simples
<direct sum of 4 S-modules>

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> RInducedModule(MakeSpecht(H,5,2,1));
<direct sum of 4 S-modules>
gap> RInducedModule(MakeSpecht(H,5,2,1),0);
<direct sum of 1 S-modules>
gap> RInducedModule(MakeSpecht(H,5,2,1),1);
<direct sum of 3 S-modules>
gap> RInducedModule(MakeSpecht(H,5,2,1),2);
<direct sum of 1 S-modules>
gap> RInducedModule(MakeSpecht(H,5,2,1),3);
<direct sum of 1 S-modules>
gap> EResidueDiagram(H,5,2,1);
   0   1   2   3   0
   3   0
   2
true

#############################
gap> H:=Specht(3);
<Hecke algebra with e = 3>
gap> RInducedModule(MakeFockPIM(H,4,2),1,2);
<direct sum of 3 Sq-modules>
gap> MakePIM(last);
<direct sum of 1 Pq-modules>

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> SInducedModule(MakePIM(H,5,2,1),3);
<direct sum of 8 P-modules>
gap> SInducedModule(MakePIM(H,5,2,1),3,1);
<direct sum of 1 P-modules>
gap> RInducedModule(MakePIM(H,5,2,1),1,1,1);
<direct sum of 1 P-modules>

#############################
gap> H:=Specht(6);
<Hecke algebra with e = 6>
gap> RRestrictedModule(MakePIM(H,5,3,2,1),4);
<direct sum of 1 P-modules>
gap> RRestrictedModule(MakeSimple(H,5,3,2),1);
<direct sum of 1 D-modules>

#############################
gap> H:=Specht(6);
<Hecke algebra with e = 6>
gap> SRestrictedModule(MakeSpecht(H,4,3,2),3);
<direct sum of 5 S-modules>
gap> SRestrictedModule(MakePIM(H,5,4,1),2,4);
<direct sum of 1 P-modules>

#############################
gap> d:=DecompositionMatrix(Specht(3,3),14);
<135x57 decomposition matrix>
gap> InducedDecompositionMatrix(d);
# Inducing....
The following projectives are missing from <d>:
    [ 15 ]  [ 8, 7 ]
<176x70 decomposition matrix>

#############################
gap> H:=Specht(2,2);
<Hecke algebra with e = 2>
gap> d:=InducedDecompositionMatrix(DecompositionMatrix(H,9));
# Inducing.
<42x10 decomposition matrix>
gap> x:=RInducedModule(MakePIM(H,9),1);
<direct sum of 2 P-modules>
gap> IsNewIndecomposable(d,x);
# This module is a sum of known indecomposables.
false
gap> x:=x-MakePIM(d,6,3,1);
<direct sum of 32 S-modules>
gap> IsNewIndecomposable(d,x,6,3,1);
# This module is a sum of known indecomposables.
false
gap> AddIndecomposable(d,x);
# AddIndecomposable: overwriting old value of P(10) in <d>

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> d:=CrystalDecompositionMatrix(H,5);
<7x6 decomposition matrix>
gap> InvertDecompositionMatrix(d);
<6x6 decomposition matrix>

#############################
gap> H:=Specht(2);; Hp:=Specht(2,2);;
gap> d:=DecompositionMatrix(H,13);; dp:=DecompositionMatrix(Hp,13);;
gap> a:=AdjustmentMatrix(dp,d);
<18x18 decomposition matrix>
gap> MatrixDecompositionMatrix(dp)=MatrixDecompositionMatrix(d)*MatrixDecompositionMatrix(a);
true

#############################
gap> H:=Specht(5,5);;
gap> d:=DecompositionMatrix(H,9);;
gap> for r in [10..20] do
> d:=InducedDecompositionMatrix(d);
> od;
# Inducing...
# Inducing....
# Inducing....
# Inducing.....
# Inducing......
# Inducing.......
# Inducing........
# Inducing..........
# Inducing............
# Inducing..............
# Inducing.................

#############################
gap> H:=Specht(2,2);
<Hecke algebra with e = 2>
gap> d:=DecompositionMatrix(H,15);
# This decomposition matrix is not known; use CalculateDecompositionMatrix()
# or InducedDecompositionMatrix() to calculate with this matrix.
fail
gap> d:=CalculateDecompositionMatrix(H,15);
# Projective indecomposable P(6,4,3,2) not known.
# Projective indecomposable P(6,5,3,1) not known.
# Projective indecomposable P(6,5,4) not known.
# Projective indecomposable P(7,4,3,1) not known.
# Projective indecomposable P(7,5,2,1) not known.
# Projective indecomposable P(7,5,3) not known.
# Projective indecomposable P(7,6,2) not known.
# Projective indecomposable P(8,4,2,1) not known.
# Projective indecomposable P(8,4,3) not known.
# Projective indecomposable P(8,5,2) not known.
# Projective indecomposable P(8,6,1) not known.
# Projective indecomposable P(8,7) not known.
# Projective indecomposable P(9,3,2,1) not known.
# Projective indecomposable P(9,4,2) not known.
# Projective indecomposable P(9,5,1) not known.
# Projective indecomposable P(9,6) not known.
# Projective indecomposable P(10,3,2) not known.
# Projective indecomposable P(10,4,1) not known.
# Projective indecomposable P(10,5) not known.
# Projective indecomposable P(11,3,1) not known.
# Projective indecomposable P(11,4) not known.
# Projective indecomposable P(12,2,1) not known.
# Projective indecomposable P(12,3) not known.
# Projective indecomposable P(13,2) not known.
# Projective indecomposable P(14,1) not known.
# Projective indecomposable P(15) not known.
<176x27 decomposition matrix>
gap> MissingIndecomposables(d);
The following projectives are missing from <d>:
    [ 15 ]  [ 14, 1 ]  [ 13, 2 ]  [ 12, 3 ]  [ 12, 2, 1 ]  [ 11, 4 ]  
[ 11, 3, 1 ]  [ 10, 5 ]  [ 10, 4, 1 ]  [ 10, 3, 2 ]  [ 9, 6 ]  [ 9, 5, 1 ]  
[ 9, 4, 2 ]  [ 9, 3, 2, 1 ]  [ 8, 7 ]  [ 8, 6, 1 ]  [ 8, 5, 2 ]  [ 8, 4, 3 ]  
[ 8, 4, 2, 1 ]  [ 7, 6, 2 ]  [ 7, 5, 3 ]  [ 7, 5, 2, 1 ]  [ 7, 4, 3, 1 ]  
[ 6, 5, 4 ]  [ 6, 5, 3, 1 ]  [ 6, 4, 3, 2 ]

#############################
gap> MatrixDecompositionMatrix(DecompositionMatrix(Specht(3),5));
[ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 1, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
  [ 1, 0, 0, 0, 1 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 0 ] ]

#############################
gap> H:=Specht(3);
<Hecke algebra with e = 3>
gap> m:=[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 1, 0, 1, 0 ], [ 0, 0, 0, 1 ], [ 0, 0, 1, 0 ] ];
[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 1, 0, 1, 0 ], [ 0, 0, 0, 1 ], 
  [ 0, 0, 1, 0 ] ]
gap> DecompositionMatrixMatrix(H,m,4);
<5x4 decomposition matrix>

#############################
gap> H:=Specht(6);
<Hecke algebra with e = 6>
gap> SimpleDimension(H,11,3);
272
gap> d:=DecompositionMatrix(H,5);
<7x7 decomposition matrix>
gap> SimpleDimension(d,3,2);
5
gap> SimpleDimension(d);
5     : 1
4,1   : 4
3,2   : 5
3,1^2 : 6
2^2,1 : 5
2,1^3 : 4
1^5   : 1
true

#############################
gap> SpechtDimension(6,3,2,1);
5632

#############################
gap> H:=Specht(2);
<Hecke algebra with e = 2>
gap> Schaper(H,9,5,3,2,1);
<direct sum of 10 S-modules>
gap> Schaper(H,9,6,5,2);
<direct sum of 1 S-modules>

#############################
gap> H:=Specht(3);
<Hecke algebra with e = 3>
gap> IsSimpleModule(H,45,31,24);
false

#############################
gap> MullineuxMap(Specht(2),12,5,2);
[ 12, 5, 2 ]
gap> MullineuxMap(Specht(4),12,5,2);
[ 4, 4, 4, 2, 2, 1, 1, 1 ]
gap> MullineuxMap(Specht(6),12,5,2);
[ 4, 3, 2, 2, 2, 2, 2, 1, 1 ]
gap> MullineuxMap(Specht(8),12,5,2);
[ 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1 ]
gap> MullineuxMap(Specht(10),12,5,2);
[ 3, 3, 3, 3, 2, 1, 1, 1, 1, 1 ]

#############################
gap> MullineuxSymbol(5,[8,6,5,5]);
[ [ 10, 6, 5, 3 ], [ 4, 4, 3, 2 ] ]

#############################
gap> PartitionMullineuxSymbol(5, MullineuxSymbol(5,[8,6,5,5]) );
[ 8, 6, 5, 5 ]

#############################
gap> GoodNodes(5,[5,4,3,2]);
[ fail, fail, 2, fail, 1 ]
gap> GoodNodes(5,[5,4,3,2],0);
fail
gap> GoodNodes(5,[5,4,3,2],4);
1

#############################
gap> NormalNodes(5,[6,5,4,4,3,2,1,1,1]);
[ [ 1, 4 ], [  ], [  ], [ 2, 5 ], [  ] ]
gap> NormalNodes(5,[6,5,4,4,3,2,1,1,1],0);
[ 1, 4 ]

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> GoodNodeSequence(H,4,3,1);
[ 0, 3, 1, 0, 2, 2, 1, 3 ]
gap> GoodNodeSequence(H,4,3,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3 ]
gap> GoodNodeSequence(H,4,4,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3, 2 ]
gap> GoodNodeSequence(H,5,4,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3, 2, 0 ]
gap> GoodNodeSequences(H,5,2,1);
[ [ 0, 1, 2, 3, 3, 2, 0, 0 ], [ 0, 3, 1, 2, 2, 3, 0, 0 ], 
  [ 0, 1, 3, 2, 2, 3, 0, 0 ], [ 0, 1, 2, 3, 3, 0, 2, 0 ], 
  [ 0, 1, 2, 3, 0, 3, 2, 0 ], [ 0, 1, 2, 3, 3, 0, 0, 2 ], 
  [ 0, 1, 2, 3, 0, 3, 0, 2 ] ]

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> PartitionGoodNodeSequence(H,0, 3, 1, 0, 2, 2, 1, 3, 3, 2);
[ 4, 4, 2 ]

#############################
gap> GoodNodeLatticePath(3,3,2,1);
[ [ 1 ], [ 1, 1 ], [ 2, 1 ], [ 2, 1, 1 ], [ 2, 2, 1 ], [ 3, 2, 1 ] ]
gap> GoodNodeLatticePaths(3,3,2,1);
[ [ [ 1 ], [ 1, 1 ], [ 2, 1 ], [ 2, 1, 1 ], [ 2, 2, 1 ], [ 3, 2, 1 ] ], 
  [ [ 1 ], [ 1, 1 ], [ 2, 1 ], [ 2, 2 ], [ 2, 2, 1 ], [ 3, 2, 1 ] ] ]
gap> GoodNodeSequence(4,6,3,2);
[ 0, 3, 1, 0, 2, 2, 3, 3, 0, 1, 1 ]
gap> LatticePathGoodNodeSequence(4,last);
[ [ 1 ], [ 1, 1 ], [ 2, 1 ], [ 2, 2 ], [ 3, 2 ], [ 3, 2, 1 ], [ 4, 2, 1 ], 
  [ 4, 2, 2 ], [ 5, 2, 2 ], [ 6, 2, 2 ], [ 6, 3, 2 ] ]

#############################
gap> H:=Specht(0);
<Hecke algebra with e = 0>
gap> LittlewoodRichardsonRule([3,2,1],[4,2]);
[ [ 4, 3, 2, 2, 1 ], [ 4, 3, 3, 1, 1 ], [ 4, 3, 3, 2 ], [ 4, 4, 2, 1, 1 ], 
  [ 4, 4, 2, 2 ], [ 4, 4, 3, 1 ], [ 5, 2, 2, 2, 1 ], [ 5, 3, 2, 1, 1 ], 
  [ 5, 3, 2, 2 ], [ 5, 4, 2, 1 ], [ 5, 3, 2, 1, 1 ], [ 5, 3, 3, 1 ], 
  [ 5, 4, 1, 1, 1 ], [ 5, 4, 2, 1 ], [ 5, 5, 1, 1 ], [ 5, 3, 2, 2 ], 
  [ 5, 3, 3, 1 ], [ 5, 4, 2, 1 ], [ 5, 4, 3 ], [ 5, 5, 2 ], [ 6, 2, 2, 1, 1 ],
  [ 6, 3, 1, 1, 1 ], [ 6, 3, 2, 1 ], [ 6, 4, 1, 1 ], [ 6, 2, 2, 2 ], 
  [ 6, 3, 2, 1 ], [ 6, 4, 2 ], [ 6, 3, 2, 1 ], [ 6, 3, 3 ], [ 6, 4, 1, 1 ], 
  [ 6, 4, 2 ], [ 6, 5, 1 ], [ 7, 2, 2, 1 ], [ 7, 3, 1, 1 ], [ 7, 3, 2 ], 
  [ 7, 4, 1 ] ]
gap> MakeSpecht(H,3,2,1)*MakeSpecht(H,4,2);
<direct sum of 27 S-modules>
gap> LittlewoodRichardsonCoefficient([3,2,1],[4,2],[5,4,2,1]);
3

#############################
gap> InverseLittlewoodRichardsonRule([3,2,1]);
[ [ [  ], [ 3, 2, 1 ] ], [ [ 1 ], [ 3, 2 ] ], [ [ 1 ], [ 2, 2, 1 ] ], 
  [ [ 1 ], [ 3, 1, 1 ] ], [ [ 1, 1 ], [ 2, 2 ] ], [ [ 1, 1 ], [ 3, 1 ] ], 
  [ [ 1, 1 ], [ 2, 1, 1 ] ], [ [ 1, 1, 1 ], [ 2, 1 ] ], [ [ 2 ], [ 2, 2 ] ], 
  [ [ 2 ], [ 3, 1 ] ], [ [ 2 ], [ 2, 1, 1 ] ], [ [ 2, 1 ], [ 3 ] ], 
  [ [ 2, 1 ], [ 2, 1 ] ], [ [ 2, 1 ], [ 2, 1 ] ], [ [ 2, 1 ], [ 1, 1, 1 ] ], 
  [ [ 2, 1, 1 ], [ 2 ] ], [ [ 2, 1, 1 ], [ 1, 1 ] ], [ [ 2, 2 ], [ 2 ] ], 
  [ [ 2, 2 ], [ 1, 1 ] ], [ [ 2, 2, 1 ], [ 1 ] ], [ [ 3 ], [ 2, 1 ] ], 
  [ [ 3, 1 ], [ 2 ] ], [ [ 3, 1 ], [ 1, 1 ] ], [ [ 3, 1, 1 ], [ 1 ] ], 
  [ [ 3, 2 ], [ 1 ] ], [ [ 3, 2, 1 ], [  ] ] ]

#############################
gap> H:=Specht(2);
<Hecke algebra with e = 2>
gap> EResidueDiagram(MakeSpecht(MakePIM(H,7,5)));
[ 7, 5 ]
   0   1   0   1   0   1   0
   1   0   1   0   1
[ 6, 5, 1 ]
   0   1   0   1   0   1
   1   0   1   0   1
   0
[ 5, 4, 2, 1 ]
   0   1   0   1   0
   1   0   1   0
   0   1
   1
# There are 3 2-regular partitions.
true

#############################
gap> RemoveRimHook([6,5,4],1,2);
[ 4, 3, 1 ]
gap> RemoveRimHook([6,5,4],2,3);
[ 6, 3, 2 ]
gap> HookLengthDiagram(6,5,4);
   8   7   6   5   3   1
   6   5   4   3   1
   4   3   2   1
true

#############################
gap> AddRimHook([6,4,3],1,3);
[ [ 9, 4, 3 ], 0 ]
gap> AddRimHook([6,4,3],2,3);
fail
gap> AddRimHook([6,4,3],3,3);
[ [ 6, 5, 5 ], 1 ]
gap> AddRimHook([6,4,3],4,3);
[ [ 6, 4, 3, 3 ], 0 ]
gap> AddRimHook([6,4,3],5,3);
fail

#############################
gap> H:=Specht(6);
<Hecke algebra with e = 6>
gap> ECore(H,16,8,6,5,3,1);
[ 4, 3, 1, 1 ]

#############################
gap> H:=Specht(8);
<Hecke algebra with e = 8>
gap> EQuotient(H,22,18,16,12,12,1,1);
[ [ 1, 1 ], [  ], [  ], [  ], [  ], [ 2, 2 ], [  ], [ 1 ] ]

#############################
gap> H:=Specht(11);
<Hecke algebra with e = 11>
gap> mu:=[100,98,57,43,12,1];
[ 100, 98, 57, 43, 12, 1 ]
gap> Q:=EQuotient(H,mu);
[ [ 9 ], [  ], [  ], [  ], [  ], [  ], [ 3 ], [ 1 ], [ 9 ], [  ], [ 5 ] ]
gap> C:=ECore(H,mu);
[ 7, 2, 2, 1, 1, 1 ]
gap> CombineEQuotientECore(H,Q,C);
[ 100, 98, 57, 43, 12, 1 ]

#############################
gap> EWeight(6,[16,8,6,5,3,1]);
5

#############################
gap> H:=Specht(3);
<Hecke algebra with e = 3>
gap> ERegularPartitions(H,6);
[ [ 2, 2, 1, 1 ], [ 3, 2, 1 ], [ 3, 3 ], [ 4, 1, 1 ], [ 4, 2 ], [ 5, 1 ], 
  [ 6 ] ]

#############################
gap> ConjugatePartition(6,4,3,2);
[ 4, 4, 3, 2, 1, 1 ]

#############################
gap> BetaSet([5,4,2,2]);
[ 2, 3, 6, 8 ]
gap> PartitionBetaSet([ 2, 3, 6, 8 ]);
[ 5, 4, 2, 2 ]

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> ETopLadder(H,1,1,1,1,1,1,1,1,1,1);
[ 4, 3, 3 ]
gap> ETopLadder(6,1,1,1,1,1,1,1,1,1,1);
[ 2, 2, 2, 2, 2 ]

#############################
gap> Dominates([5,4],[4,4,1]);
true
gap> p:=Partitions(6);;Sort(p,LengthLexicographic); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 3, 3 ], [ 4, 1, 1 ], [ 3, 2, 1 ], [ 2, 2, 2 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1 ] ]
gap> p:=Partitions(6);;Sort(p,Lexicographic); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 4, 1, 1 ], [ 3, 3 ], [ 3, 2, 1 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 2 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], 
  [ 1, 1, 1, 1, 1, 1 ] ]
gap> p:=Partitions(6);;Sort(p,ReverseDominance); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 3, 3 ], [ 4, 1, 1 ], [ 3, 2, 1 ], [ 2, 2, 2 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1 ] ]

#############################
gap> H:=Specht(2);
<Hecke algebra with e = 2>
gap> x:=MakeFockPIM(H,6,2);
<direct sum of 13 Sq-modules>
gap> Specialized(x);
<direct sum of 13 S-modules>
gap> Specialized(x,2);
<direct sum of 13 S-modules>

#############################
gap> H:=Specht(8);
<Hecke algebra with e = 8>
gap> x:=MakeSpecht(RInducedModule(MakePIM(H,8,5,3)) );
<direct sum of 7 S-modules>
gap> ERegulars(x);
[ 9, 5, 3 ]  [ 8, 6, 3 ]  [ 8, 5, 4 ]  [ 8, 5, 3, 1 ]  
[ 6, 5, 3, 3 ]  [ 5, 5, 4, 3 ]  [ 5, 5, 3, 3, 1 ]  
gap> MakePIM(x);
<direct sum of 4 P-modules>

#############################
gap> H:=Specht(2);
<Hecke algebra with e = 2>
gap> SplitECores(RInducedModule(MakeSpecht(H,5,3,1)));
[ <direct sum of 3 S-modules>, <direct sum of 1 S-modules> ]
gap> RInducedModule(MakeSpecht(H,5,3,1),0);
<direct sum of 1 S-modules>
gap> RInducedModule(MakeSpecht(H,5,3,1),1);
<direct sum of 3 S-modules>

#############################
gap> H:=Specht(3);
<Hecke algebra with e = 3>
gap> x:=MakeSpecht(MakePIM(H,7,3));
<direct sum of 8 S-modules>
gap> Coefficient(x,5,2,2,1);
1

#############################
gap> H:=Specht(4);
<Hecke algebra with e = 4>
gap> InnerProduct(MakeSpecht(H,2,2,2,1), MakePIM(H,4,3));
1
gap> DecompositionNumber(H,[2,2,2,1],[4,3]);
1

#############################
gap> H:=Specht(2);;
gap> x:=MakeSpecht(MakePIM(H,6));
<direct sum of 6 S-modules>
gap> Display(x);
S(6) + S(5,1) + S(4,1^2) + S(3,1^3) + S(2,1^4) + S(1^6)
gap> Print(x,"\n");
S(6) + S(5,1) + S(4,1,1) + S(3,1,1,1) + S(2,1,1,1,1) + S(1,1,1,1,1,1)

#############################
gap> SemiStandardTableaux([4,3],1,1,1,2,2);
[ <tableau of shape [ 4, 3 ]>, <tableau of shape [ 4, 3 ]>, 
  <tableau of shape [ 4, 3 ]>, <tableau of shape [ 4, 3 ]>, 
  <tableau of shape [ 4, 3 ]>, <tableau of shape [ 4, 3 ]> ]
gap> StandardTableaux(4,2);
[ <tableau of shape [ 4, 2 ]>, <tableau of shape [ 4, 2 ]>, 
  <tableau of shape [ 4, 2 ]>, <tableau of shape [ 4, 2 ]>, 
  <tableau of shape [ 4, 2 ]>, <tableau of shape [ 4, 2 ]>, 
  <tableau of shape [ 4, 2 ]>, <tableau of shape [ 4, 2 ]>, 
  <tableau of shape [ 4, 2 ]> ]
gap> ConjugateTableau(Tableau([ [ 1, 3, 5, 6 ], [ 2, 4 ] ]));
<tableau of shape [ 2, 2, 1, 1 ]>
gap> ShapeTableau(Tableau( [ [ 1, 1, 2, 3 ], [ 4, 5 ] ] ));
[ 4, 2 ]
gap> List(SemiStandardTableaux([5,4,2],[4,3,0,1,3]),TypeTableau);
[ [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], 
  [ 4, 3, 0, 1, 3 ] ]
