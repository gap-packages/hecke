# hecke, chapter 3
#
# DO NOT EDIT THIS FILE - EDIT EXAMPLES IN THE SOURCE INSTEAD!
#
# This file has been generated by AutoDoc. It contains examples extracted from
# the package documentation. Each example is preceded by a comment which gives
# the name of a GAPDoc XML file and a line range from which the example were
# taken. Note that the XML file in turn may have been generated by AutoDoc
# from some other input.
#
gap> START_TEST("hecke02.tst");

# doc/functionality.xml:142-149
gap> H:=Specht(5);
<Hecke algebra with e = 5>
gap> Display(last);
Specht(e=5, S(), P(), D())
gap> IsZeroCharacteristic(last);
true

# doc/functionality.xml:204-207
gap> H:=Specht(5);; MakePIM(H,4,3,2);; Display(last);
P(4,3,2)

# doc/functionality.xml:226-231
gap> Display( MakeSimple( MakePIM(H,4,3,2) ) );
D(5,3,1) + 2D(4,3,2) + D(2^4,1)
gap> Display( MakeSpecht( MakeSimple( MakeSpecht(H,1,1,1,1,1) ) ) );
 - S(5) + S(4,1) - S(3,1^2) + S(2,1^3)

# doc/functionality.xml:257-269
gap> H:=Specht(3,3);;   # e = 3, p = 3 = characteristic of 'R'
gap>  d:=InducedDecompositionMatrix(DecompositionMatrix(H,14));;
# Inducing....
The following projectives are missing from <d>:
    [ 15 ]  [ 8, 7 ]
gap> Display(MakePIM(d,4,3,3,2,2,1));
S(4,3^2,2^2,1) + S(4,3^2,2,1^3) + S(4,3,2^3,1^2) + S(3^3,2^2,1^2)
gap> Display(MakeSpecht(d,7, 3, 3, 2));
D(11,2,1^2) + D(10,3,1^2) + D(8,5,1^2) + D(8,3^2,1) + D(7,6,1^2) + D(7,3^2,2)
gap> Display(MakeSimple(d,14,1));
fail

# doc/functionality.xml:290-293
gap> Display(MakeSpecht(d, MakeSimple(d,10,5) ));
 - S(13,2) + S(10,5)

# doc/functionality.xml:310-323
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

# doc/functionality.xml:351-359
gap> val:=function(x) local v;
>       x:=Sum([0..x-1],v->4^v);  # x->[x]_q
>       v:=0; while x mod 5=0 do x:=x/5; v:=v+1; od;
>       return v;
>     end;;
gap> H:=Specht(2,5,val,"e2q4");; Display(H);
Specht(e=2, p=5, S(), P(), D(), HeckeRing="e2q4")

# doc/functionality.xml:428-442
gap> H:=Specht(4);; MakeFockPIM(H,6,2);; Display(last);
Sq(6,2) + vSq(5,3)
gap> RRestrictedModule(last); Display(last);
<direct sum of 3 Sq-modules>
Sq(6,1) + (v+v^-1)Sq(5,2) + vSq(4,3)
gap> MakePIM(last);; Display(last);
Pq(6,1) + (v+v^-1)Pq(5,2)
gap> Specialized(last);; Display(last);
P(6,1) + 2P(5,2)
gap> MakeFockSpecht(H,5,3,2);; Display(last);
Sq(5,3,2)
gap> RInducedModule(last,0);; Display(last);
v^-1Sq(5,3^2)

# doc/functionality.xml:503-519
gap> S:=Schur(2);
<Schur algebra with e = 2>
gap> InducedDecompositionMatrix(DecompositionMatrix(S,3));
The following projectives are missing from <d>:
    [ 2, 2 ]
<5x5 decomposition matrix>
gap> Display(last);
4    | 1         
3,1  | 1 1       
2^2  | . 1 .     
2,1^2| 1 1 . 1   
1^4  | 1 . . 1 1

# DecompositionMatrix(S,4) returns the full decomposition matrix. The point of
# this example is to emphasize the current limitations of Schur.

# doc/functionality.xml:583-598
gap> DecompositionMatrix(Specht(3),6,LengthLexicographic);
<11x7 decomposition matrix>
gap> Display(last);
6      | 1             
5,1    | 1 1           
4,2    | . . 1         
3^2    | . 1 . 1       
4,1^2  | . 1 . . 1     
3,2,1  | 1 1 . 1 1 1   
2^3    | 1 . . . . 1   
3,1^3  | . . . . 1 1   
2^2,1^2| . . . . . . 1
2,1^4  | . . . 1 . 1 . 
1^6    | . . . 1 . . . 

# doc/functionality.xml:639-668
gap> CrystalDecompositionMatrix(Specht(3), 6);
<11x7 decomposition matrix>
gap> Display(last);
6      |   1                         
5,1    |   v   1                     
4,2    |   .   .   1                 
4,1^2  |   .   v   .   1             
3^2    |   .   v   .   .   1         
3,2,1  |   v v^2   .   v   v   1     
3,1^3  |   .   .   . v^2   .   v     
2^3    | v^2   .   .   .   .   v     
2^2,1^2|   .   .   .   .   .   .   1
2,1^4  |   .   .   .   .   v v^2   . 
1^6    |   .   .   .   . v^2   .   . 
gap> Specialized(last);   # set v equal to 1.
<11x7 decomposition matrix>
gap> Display(last);
6      | 1             
5,1    | 1 1           
4,2    | . . 1         
4,1^2  | . 1 . 1       
3^2    | . 1 . . 1     
3,2,1  | 1 1 . 1 1 1   
3,1^3  | . . . 1 . 1   
2^3    | 1 . . . . 1   
2^2,1^2| . . . . . . 1
2,1^4  | . . . . 1 1 . 
1^6    | . . . . 1 . . 

# doc/functionality.xml:688-691
gap> H:=Specht(6);; DecompositionNumber(H,[6,4,2],[6,6]);
0

# doc/functionality.xml:705-710
gap> H:=Specht(4);; Print(MakeSpecht(MakePIM(H,6,4)),"\n");
S(6,4) + S(6,3,1) + S(5,3,1,1) + S(3,3,2,1,1) + S(2,2,2,2,2)
gap> Print(MakeSpecht(MakePIM(H,[6,4])),"\n");
S(6,4) + S(6,3,1) + S(5,3,1,1) + S(3,3,2,1,1) + S(2,2,2,2,2)

# doc/functionality.xml:715-724
gap> ECore(3, [6,4,2]);
[ 6, 4, 2 ]
gap> ECore(3, 6,4,2);
[ 6, 4, 2 ]
gap> GoodNodes(3, 6,4,2);
[ fail, fail, 3 ]
gap> GoodNodes(3, [6,4,2]);
[ fail, fail, 3 ]

# doc/functionality.xml:772-781
gap> H:=Specht(2,2);;
gap> Display(RInducedModule(MakeSpecht(H,7,4,3,1)));
S(8,4,3,1) + S(7,5,3,1) + S(7,4^2,1) + S(7,4,3,2) + S(7,4,3,1^2)
gap> Display(RInducedModule(MakePIM(H,5,3,1)));
P(6,3,1) + 2P(5,4,1) + P(5,3,2)
gap> Display(RInducedModule(MakeSimple(H,11,2,1)));
# D(<x>), unable to rewrite <x> as a sum of simples
S(12,2,1) + S(11,3,1) + S(11,2^2) + S(11,2,1^2)

# doc/functionality.xml:807-818
gap> H:=Specht(4);; Display(RInducedModule(MakeSpecht(H,5,2,1)));
S(6,2,1) + S(5,3,1) + S(5,2^2) + S(5,2,1^2)
gap> Display(RInducedModule(MakeSpecht(H,5,2,1),0));
0S()
gap> Display(RInducedModule(MakeSpecht(H,5,2,1),1));
S(6,2,1) + S(5,3,1) + S(5,2,1^2)
gap> Display(RInducedModule(MakeSpecht(H,5,2,1),2));
0S()
gap> Display(RInducedModule(MakeSpecht(H,5,2,1),3));
S(5,2^2)

# doc/functionality.xml:824-830
gap> EResidueDiagram(H,5,2,1);
   0   1   2   3   0
   3   0
   2
true

# doc/functionality.xml:844-849
gap> H:=Specht(3);; x:=RInducedModule(MakeFockPIM(H,4,2),1,2);;
gap> Display(x); Display(MakePIM(x));
Sq(6,2) + vSq(4^2) + v^2Sq(4,2^2)
Pq(6,2)

# doc/functionality.xml:867-876
gap> SizeScreen([80,20]);;
gap> H:=Specht(4);; Display(SInducedModule(MakePIM(H,5,2,1),3));
P(8,2,1) + 3P(7,3,1) + 2P(7,2^2) + 6P(6,3,2) + 6P(6,3,1^2) + 3P(6,2,1^3) + 2P(\
5,3^2) + P(5,2^2,1^2)
gap> Display(SInducedModule(MakePIM(H,5,2,1),3,1));
P(6,3,1^2)
gap> Display(RInducedModule(MakePIM(H,5,2,1),1,1,1));
6P(6,3,1^2)

# doc/functionality.xml:918-923
gap> H:=Specht(6);; Display(RRestrictedModule(MakePIM(H,5,3,2,1),4));
2P(4,3,2,1)
gap> Display(RRestrictedModule(MakeSimple(H,5,3,2),1));
D(5,2^2)

# doc/functionality.xml:951-956
gap> H:=Specht(6);; Display(SRestrictedModule(MakeSpecht(H,4,3,2),3));
3S(4,2) + 2S(4,1^2) + 3S(3^2) + 6S(3,2,1) + 2S(2^3)
gap> Display(SRestrictedModule(MakePIM(H,5,4,1),2,4));
P(4^2)

# doc/functionality.xml:992-996
gap> H:=Specht(4);; for n in [8..20] do
>      SaveDecompositionMatrix(DecompositionMatrix(H,n));
>    od;

# doc/functionality.xml:1019-1027
gap> d:=DecompositionMatrix(Specht(3,3),14);
<135x57 decomposition matrix>
gap> InducedDecompositionMatrix(d);
# Inducing....
The following projectives are missing from <d>:
    [ 15 ]  [ 8, 7 ]
<176x70 decomposition matrix>

# doc/functionality.xml:1161-1170
gap> H:=Specht(4);; d:=CrystalDecompositionMatrix(H,5);;
gap> Display(InvertDecompositionMatrix(d));
5    |   1                     
4,1  |   .   1                 
3,2  |  -v   .   1             
3,1^2|   .   .   .   1         
2^2,1| v^2   .  -v   .   1     
2,1^3|   .   .   .   .   .   1

# doc/functionality.xml:1189-1216
gap> H:=Specht(2);; Hp:=Specht(2,2);;
gap> d:=DecompositionMatrix(H,13);; dp:=DecompositionMatrix(Hp,13);;
gap> a:=AdjustmentMatrix(dp,d);
<18x18 decomposition matrix>
gap> Display(a);
13     | 1                                   
12,1   | . 1                                 
11,2   | 1 . 1                               
10,3   | . . . 1                             
10,2,1 | . . . . 1                           
9,4    | 1 . 1 . . 1                         
9,3,1  | 2 . . . . . 1                       
8,5    | . 1 . . . . . 1                     
8,4,1  | 1 . . . . . . . 1                   
8,3,2  | . 2 . . . . . 1 . 1                 
7,6    | 1 . . . . 1 . . . . 1               
7,5,1  | . . . . . . 1 . . . . 1             
7,4,2  | 1 . 1 . . 1 . . . . 1 . 1           
7,3,2,1| . . . . . . . . . . . . . 1         
6,5,2  | . 1 . . . . . 1 . 1 . . . . 1       
6,4,3  | 2 . . . 1 . . . . . . . . . . 1     
6,4,2,1| . 2 . 1 . . . . . . . . . . . . 1   
5,4,3,1| 4 . 2 . . . . . . . . . . . . . . 1
gap> MatrixDecompositionMatrix(dp)=
>           MatrixDecompositionMatrix(d)*MatrixDecompositionMatrix(a);
true

# doc/functionality.xml:1247-1265
gap> H:=Specht(5,5);;
gap> d:=DecompositionMatrix(H,9);;
gap> for r in [10..20] do
>      d:=InducedDecompositionMatrix(d);
>      SaveDecompositionMatrix(d);
>    od;
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

# doc/functionality.xml:1301-1340
gap> H:=Specht(2,2);; d:=DecompositionMatrix(H,15);
# This decomposition matrix is not known; use CalculateDecompositionMatrix()
# or InducedDecompositionMatrix() to calculate with this matrix.
fail
gap> d:=CalculateDecompositionMatrix(H,15);;
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
gap> SizeScreen([80,20]);; MissingIndecomposables(d);
The following projectives are missing from <d>:
    [ 15 ]  [ 14, 1 ]  [ 13, 2 ]  [ 12, 3 ]  [ 12, 2, 1 ]  [ 11, 4 ]  
[ 11, 3, 1 ]  [ 10, 5 ]  [ 10, 4, 1 ]  [ 10, 3, 2 ]  [ 9, 6 ]  [ 9, 5, 1 ]  
[ 9, 4, 2 ]  [ 9, 3, 2, 1 ]  [ 8, 7 ]  [ 8, 6, 1 ]  [ 8, 5, 2 ]  [ 8, 4, 3 ]  
[ 8, 4, 2, 1 ]  [ 7, 6, 2 ]  [ 7, 5, 3 ]  [ 7, 5, 2, 1 ]  [ 7, 4, 3, 1 ]  
[ 6, 5, 4 ]  [ 6, 5, 3, 1 ]  [ 6, 4, 3, 2 ]

# doc/functionality.xml:1359-1364
gap> SizeScreen([80,20]);;
gap> MatrixDecompositionMatrix(DecompositionMatrix(Specht(3),5));
[ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 1, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
  [ 1, 0, 0, 0, 1 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 0 ] ]

# doc/functionality.xml:1383-1393
gap> H:=Specht(3);;
gap> m:=[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 1, 0, 1, 0 ],
>         [ 0, 0, 0, 1 ], [ 0, 0, 1, 0 ] ];;
gap> Display(DecompositionMatrixMatrix(H,m,4));
4    | 1       
3,1  | . 1     
2^2  | 1 . 1   
2,1^2| . . . 1
1^4  | . . 1 . 

# doc/functionality.xml:1470-1485
gap> H:=Specht(6);;
gap> SimpleDimension(H,11,3);
272
gap> d:=DecompositionMatrix(H,5);; SimpleDimension(d,3,2);
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

# doc/functionality.xml:1500-1503
gap> SpechtDimension(6,3,2,1);
5632

# doc/functionality.xml:1545-1552
gap> SizeScreen([80,20]);; H:=Specht(2);;
gap> Display(Schaper(H,9,5,3,2,1));
S(17,2,1) - S(15,2,1^3) + S(13,2^3,1) - S(11,3^2,2,1) + S(10,4,3,2,1) - S(9,8,\
3) - S(9,8,1^3) + S(9,6,3,2) + S(9,6,3,1^2) + S(9,6,2^2,1)
gap> Display(Schaper(H,9,6,5,2));
0S()

# doc/functionality.xml:1575-1579
gap> H:=Specht(3);;
gap> IsSimpleModule(H,45,31,24);
false

# doc/functionality.xml:1613-1624
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

# doc/functionality.xml:1660-1663
gap> MullineuxSymbol(5,[8,6,5,5]);
[ [ 10, 6, 5, 3 ], [ 4, 4, 3, 2 ] ]

# doc/functionality.xml:1675-1678
gap> PartitionMullineuxSymbol(5, MullineuxSymbol(5,[8,6,5,5]) );
[ 8, 6, 5, 5 ]

# doc/functionality.xml:1705-1712
gap> GoodNodes(5,[5,4,3,2]);
[ fail, fail, 2, fail, 1 ]
gap> GoodNodes(5,[5,4,3,2],0);
fail
gap> GoodNodes(5,[5,4,3,2],4);
1

# doc/functionality.xml:1729-1734
gap> NormalNodes(5,[6,5,4,4,3,2,1,1,1]);
[ [ 1, 4 ], [  ], [  ], [ 2, 5 ], [  ] ]
gap> NormalNodes(5,[6,5,4,4,3,2,1,1,1],0);
[ 1, 4 ]

# doc/functionality.xml:1753-1762
gap> H:=Specht(4);; GoodNodeSequence(H,4,3,1);
[ 0, 3, 1, 0, 2, 2, 1, 3 ]
gap> GoodNodeSequence(H,4,3,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3 ]
gap> GoodNodeSequence(H,4,4,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3, 2 ]
gap> GoodNodeSequence(H,5,4,2);
[ 0, 3, 1, 0, 2, 2, 1, 3, 3, 2, 0 ]

# doc/functionality.xml:1767-1773
gap> H:=Specht(4);; GoodNodeSequences(H,5,2,1);
[ [ 0, 1, 2, 3, 3, 2, 0, 0 ], [ 0, 3, 1, 2, 2, 3, 0, 0 ], 
  [ 0, 1, 3, 2, 2, 3, 0, 0 ], [ 0, 1, 2, 3, 3, 0, 2, 0 ], 
  [ 0, 1, 2, 3, 0, 3, 2, 0 ], [ 0, 1, 2, 3, 3, 0, 0, 2 ], 
  [ 0, 1, 2, 3, 0, 3, 0, 2 ] ]

# doc/functionality.xml:1787-1791
gap> H:=Specht(4);;
gap> PartitionGoodNodeSequence(H,0, 3, 1, 0, 2, 2, 1, 3, 3, 2);
[ 4, 4, 2 ]

# doc/functionality.xml:1810-1821
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

# doc/functionality.xml:1849-1870
gap> SizeScreen([80,20]);;
gap> H:=Specht(0);; # the generic Hecke algebra with R=C[q]
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
gap> Display(MakeSpecht(H,3,2,1)*MakeSpecht(H,4,2));
S(7,4,1) + S(7,3,2) + S(7,3,1^2) + S(7,2^2,1) + S(6,5,1) + 2S(6,4,2) + 2S(6,4,\
1^2) + S(6,3^2) + 3S(6,3,2,1) + S(6,3,1^3) + S(6,2^3) + S(6,2^2,1^2) + S(5^2,2\
) + S(5^2,1^2) + S(5,4,3) + 3S(5,4,2,1) + S(5,4,1^3) + 2S(5,3^2,1) + 2S(5,3,2^\
2) + 2S(5,3,2,1^2) + S(5,2^3,1) + S(4^2,3,1) + S(4^2,2^2) + S(4^2,2,1^2) + S(4\
,3^2,2) + S(4,3^2,1^2) + S(4,3,2^2,1)
gap> LittlewoodRichardsonCoefficient([3,2,1],[4,2],[5,4,2,1]);
3

# doc/functionality.xml:1890-1901
gap> SizeScreen([80,20]);; InverseLittlewoodRichardsonRule(3,2,1);
[ [ [  ], [ 3, 2, 1 ] ], [ [ 1 ], [ 3, 2 ] ], [ [ 1 ], [ 2, 2, 1 ] ], 
  [ [ 1 ], [ 3, 1, 1 ] ], [ [ 1, 1 ], [ 2, 2 ] ], [ [ 1, 1 ], [ 3, 1 ] ], 
  [ [ 1, 1 ], [ 2, 1, 1 ] ], [ [ 1, 1, 1 ], [ 2, 1 ] ], [ [ 2 ], [ 2, 2 ] ], 
  [ [ 2 ], [ 3, 1 ] ], [ [ 2 ], [ 2, 1, 1 ] ], [ [ 2, 1 ], [ 3 ] ], 
  [ [ 2, 1 ], [ 2, 1 ] ], [ [ 2, 1 ], [ 2, 1 ] ], [ [ 2, 1 ], [ 1, 1, 1 ] ], 
  [ [ 2, 1, 1 ], [ 2 ] ], [ [ 2, 1, 1 ], [ 1, 1 ] ], [ [ 2, 2 ], [ 2 ] ], 
  [ [ 2, 2 ], [ 1, 1 ] ], [ [ 2, 2, 1 ], [ 1 ] ], [ [ 3 ], [ 2, 1 ] ], 
  [ [ 3, 1 ], [ 2 ] ], [ [ 3, 1 ], [ 1, 1 ] ], [ [ 3, 1, 1 ], [ 1 ] ], 
  [ [ 3, 2 ], [ 1 ] ], [ [ 3, 2, 1 ], [  ] ] ]

# doc/functionality.xml:1925-1941
gap> H:=Specht(2);; EResidueDiagram(MakeSpecht(MakePIM(H,7,5)));
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

# doc/functionality.xml:1954-1961
gap> HookLengthDiagram(11,6,3,2);
  14  13  11   9   8   7   5   4   3   2   1
   8   7   5   3   2   1
   4   3   1
   2   1
true

# doc/functionality.xml:1973-1983
gap> RemoveRimHook([6,5,4],1,2);
[ 4, 3, 1 ]
gap> RemoveRimHook([6,5,4],2,3);
[ 6, 3, 2 ]
gap> HookLengthDiagram(6,5,4);
   8   7   6   5   3   1
   6   5   4   3   1
   4   3   2   1
true

# doc/functionality.xml:1999-2010
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

# doc/functionality.xml:2032-2035
gap> H:=Specht(6);; ECore(H,16,8,6,5,3,1);
[ 4, 3, 1, 1 ]

# doc/functionality.xml:2062-2065
gap> H:=Specht(8);; EQuotient(H,22,18,16,12,12,1,1);
[ [ 1, 1 ], [  ], [  ], [  ], [  ], [ 2, 2 ], [  ], [ 1 ] ]

# doc/functionality.xml:2081-2089
gap> H:=Specht(11);; mu:=[100,98,57,43,12,1];;
gap> Q:=EQuotient(H,mu);
[ [ 9 ], [  ], [  ], [  ], [  ], [  ], [ 3 ], [ 1 ], [ 9 ], [  ], [ 5 ] ]
gap> C:=ECore(H,mu);
[ 7, 2, 2, 1, 1, 1 ]
gap> CombineEQuotientECore(H,Q,C);
[ 100, 98, 57, 43, 12, 1 ]

# doc/functionality.xml:2104-2107
gap> EWeight(6,[16,8,6,5,3,1]);
5

# doc/functionality.xml:2123-2127
gap> H:=Specht(3);; ERegularPartitions(H,6);
[ [ 2, 2, 1, 1 ], [ 3, 2, 1 ], [ 3, 3 ], [ 4, 1, 1 ], [ 4, 2 ], [ 5, 1 ], 
  [ 6 ] ]

# doc/functionality.xml:2149-2152
gap> ConjugatePartition(6,4,3,2);
[ 4, 4, 3, 2, 1, 1 ]

# doc/functionality.xml:2165-2168
gap> PartitionBetaSet([ 2, 3, 6, 8 ]);
[ 5, 4, 2, 2 ]

# doc/functionality.xml:2187-2193
gap> H:=Specht(4);;
gap> ETopLadder(H,1,1,1,1,1,1,1,1,1,1);
[ 4, 3, 3 ]
gap> ETopLadder(6,1,1,1,1,1,1,1,1,1,1);
[ 2, 2, 2, 2, 2 ]

# doc/functionality.xml:2209-2212
gap> Dominates([5,4],[4,4,1]);
true

# doc/functionality.xml:2225-2229
gap> p:=Partitions(6);;Sort(p,LengthLexicographic); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 3, 3 ], [ 4, 1, 1 ], [ 3, 2, 1 ], [ 2, 2, 2 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1 ] ]

# doc/functionality.xml:2241-2246
gap> p:=Partitions(6);;Sort(p,Lexicographic); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 4, 1, 1 ], [ 3, 3 ], [ 3, 2, 1 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 2 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], 
  [ 1, 1, 1, 1, 1, 1 ] ]

# doc/functionality.xml:2260-2264
gap> p:=Partitions(6);;Sort(p,ReverseDominance); p;
[ [ 6 ], [ 5, 1 ], [ 4, 2 ], [ 3, 3 ], [ 4, 1, 1 ], [ 3, 2, 1 ], [ 2, 2, 2 ], 
  [ 3, 1, 1, 1 ], [ 2, 2, 1, 1 ], [ 2, 1, 1, 1, 1 ], [ 1, 1, 1, 1, 1, 1 ] ]

# doc/functionality.xml:2294-2306
gap> SizeScreen([80,20]);; H:=Specht(2);; x:=MakeFockPIM(H,6,2);; Display(x);
Sq(6,2) + vSq(6,1^2) + vSq(5,3) + v^2Sq(5,1^3) + vSq(4,3,1) + v^2Sq(4,2^2) + (\
v^3+v)Sq(4,2,1^2) + v^2Sq(4,1^4) + v^2Sq(3^2,1^2) + v^3Sq(3,2^2,1) + v^3Sq(3,1\
^5) + v^3Sq(2^3,1^2) + v^4Sq(2^2,1^4)
gap> Display(Specialized(x));
S(6,2) + S(6,1^2) + S(5,3) + S(5,1^3) + S(4,3,1) + S(4,2^2) + 2S(4,2,1^2) + S(\
4,1^4) + S(3^2,1^2) + S(3,2^2,1) + S(3,1^5) + S(2^3,1^2) + S(2^2,1^4)
gap> Display(Specialized(x,2));
S(6,2) + 2S(6,1^2) + 2S(5,3) + 4S(5,1^3) + 2S(4,3,1) + 4S(4,2^2) + 10S(4,2,1^2\
) + 4S(4,1^4) + 4S(3^2,1^2) + 8S(3,2^2,1) + 8S(3,1^5) + 8S(2^3,1^2) + 16S(2^2,\
1^4)

# doc/functionality.xml:2325-2335
gap> H:=Specht(8);;
gap> x:=MakeSpecht(RInducedModule(MakePIM(H,8,5,3)));; Display(x);
S(9,5,3) + S(8,6,3) + S(8,5,4) + S(8,5,3,1) + S(6,5,3^2) + S(5^2,4,3) + S(5^2,\
3^2,1)
gap> ERegulars(x);
[ 9, 5, 3 ]  [ 8, 6, 3 ]  [ 8, 5, 4 ]  [ 8, 5, 3, 1 ]  
[ 6, 5, 3, 3 ]  [ 5, 5, 4, 3 ]  [ 5, 5, 3, 3, 1 ]  
gap> Display(MakePIM(x));
P(9,5,3) + P(8,6,3) + P(8,5,4) + P(8,5,3,1)

# doc/functionality.xml:2367-2375
gap> H:=Specht(2);;
gap> Display(SplitECores(RInducedModule(MakeSpecht(H,5,3,1))));
[ S(6,3,1) + S(5,3,2) + S(5,3,1,1), S(5,4,1) ]
gap> Display(RInducedModule(MakeSpecht(H,5,3,1),0));
S(5,4,1)
gap> Display(RInducedModule(MakeSpecht(H,5,3,1),1));
S(6,3,1) + S(5,3,2) + S(5,3,1^2)

# doc/functionality.xml:2388-2395
gap> SizeScreen([80,20]);;
gap> H:=Specht(3);; x:=MakeSpecht(MakePIM(H,7,3));; Display(x);
S(7,3) + S(7,2,1) + S(6,2,1^2) + S(5^2) + S(5,2^2,1) + S(4^2,1^2) + S(4,3^2) +\
 S(4,3,2,1)
gap> Coefficient(x,5,2,2,1);
1

# doc/functionality.xml:2410-2415
gap> H:=Specht(2);; InnerProduct(MakeSpecht(H,2,2,2,1), MakePIM(H,4,3));
1
gap> DecompositionNumber(H,[2,2,2,1],[4,3]);
1

# doc/functionality.xml:2444-2452
gap> SizeScreen([80,20]);; Display(SemiStandardTableaux([4,3],[1,1,1,2,2]));
[ Tableau( [ [ 1, 2, 3, 4 ], [ 4, 5, 5 ] ] ), 
  Tableau( [ [ 1, 2, 3, 5 ], [ 4, 4, 5 ] ] ), 
  Tableau( [ [ 1, 2, 4, 4 ], [ 3, 5, 5 ] ] ), 
  Tableau( [ [ 1, 2, 4, 5 ], [ 3, 4, 5 ] ] ), 
  Tableau( [ [ 1, 3, 4, 4 ], [ 2, 5, 5 ] ] ), 
  Tableau( [ [ 1, 3, 4, 5 ], [ 2, 4, 5 ] ] ) ]

# doc/functionality.xml:2463-2474
gap> SizeScreen([80,20]);; Display(StandardTableaux(4,2));
[ Tableau( [ [ 1, 2, 3, 4 ], [ 5, 6 ] ] ), 
  Tableau( [ [ 1, 2, 3, 5 ], [ 4, 6 ] ] ), 
  Tableau( [ [ 1, 2, 3, 6 ], [ 4, 5 ] ] ), 
  Tableau( [ [ 1, 2, 4, 5 ], [ 3, 6 ] ] ), 
  Tableau( [ [ 1, 2, 4, 6 ], [ 3, 5 ] ] ), 
  Tableau( [ [ 1, 2, 5, 6 ], [ 3, 4 ] ] ), 
  Tableau( [ [ 1, 3, 4, 5 ], [ 2, 6 ] ] ), 
  Tableau( [ [ 1, 3, 4, 6 ], [ 2, 5 ] ] ), 
  Tableau( [ [ 1, 3, 5, 6 ], [ 2, 4 ] ] ) ]

# doc/functionality.xml:2485-2492
gap> Display(ConjugateTableau(Tableau([ [ 1, 3, 5, 6 ], [ 2, 4 ] ])));
Standard Tableau:
1	2	
3	4	
5	
6	

# doc/functionality.xml:2501-2504
gap> ShapeTableau( Tableau([ [ 1, 1, 2, 3 ], [ 4, 5 ] ]) );
[ 4, 2 ]

# doc/functionality.xml:2516-2521
gap> SizeScreen([80,20]);;
gap> List(SemiStandardTableaux([5,4,2],[4,3,0,1,3]),TypeTableau);
[ [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], [ 4, 3, 0, 1, 3 ], 
  [ 4, 3, 0, 1, 3 ] ]

#
gap> STOP_TEST("hecke02.tst", 1);
