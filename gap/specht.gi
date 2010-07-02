#######################################################################
##  SPECHT - specht.g : the kernel of SPECHT                         ##
##                                                                   ##
##     A GAP package for calculating the decomposition numbers of    ##
##     Hecke algebras of type A (over fields of characteristic       ##
##     zero). The functions provided are primarily combinatorial in  ##
##     nature. Many of the combinatorial tools from the (modular)    ##
##     representation theory of the symmetric groups appear in the   ##
##     package.                                                      ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-version by Andrew Mathas)             ##
##                                                                   ##
#######################################################################

## 3.0: June 2010:
##   - Translated to GAP4

## Change log
## 2.4:
##  - fixed more bugs in H.valuation; returned incorrect answers before
##    when e=0 or e=p (symmetric group case).
##  - fixed bug in Dq(), via Sasha Kleshchev.

## 2.3:
##  - fixed bug in H.valuation  reported by Johannes Lipp
##  - fixed bug in Sq() reported by Johannes Lipp.
##  - corrected FindDecompositionMatrix() so that it updates the matrices
##    CrystalMatrices[] and DecompositionMatrices[] after calculating a
##    crystalized decomposition matrix.

## 2.2: June 1996: various changes requested by the referee.
##  - mainly changing function names.
##  - DecompositionMatrix changed so that it no longer attempts to
##    calculate decomposition matrices in the finite field case; added
##    CalculateDecompositionMatrix() to do this.
##  - replaced Matrix() function with MatrixDecompositionMatrix() and
##    added the function DecompositionMatrixMatix().

## 2.1: April 1996:
##  - Added a filename argument to SaveDecompositionMatrix and made it
##    save the non-decomposition matrices under (more) sensible names; the 
##    later is done using the existence of a record component d.matname.
##    Need to do something here about reading such matrices in (this can be
##    done using DecompositionMatrix)..
##  - Changed ReadDecompositionmatrix so that it automatically reads the
##    format of version 1.0 files (and deleted ReadOldDecompositionMatrix).
##    Also fixed a bug here which was causing confusion between Hecke algebra
##    and Schur algebra matrices.
##  - Renamed FindDecompositionMatrix as KnownDecompositionMatrix because it
##    doesn't try to calculate a crytalized decomposition matrix; the new
##    FindDecompositionMatrix will return the decomposition matrix if at all
##    possible. In particular, this fixes a 'bug' in SimpleDimension.
##  - Rewrote AdjustmentMatrix so that it actually works.
##  - Changed P()->S() module conversions so that it no longer requires
##    the projective module to have positive coefficients (and fixed bug in
##    matrix ops).

## 2.0: March 1996: 
##   - LLT algorithm implemented for calculating crystal basis of
##     the Fock space and hence by specialization the decomposition
##     matrices for all Hecke algebras over fields of characteristic
##     zero; most of the work is done by the internal function Pq(). 
##     This required a new set of 'module' types ("Pq", "Sq", and "Dq"),
##     with correspondining operation sets H.operations.X. In particular,
##     these include q-induction and q-restriction, effectively allowing
##     computations with $U_q(\hat sl_e)$ on the Fock space.
##   - crystallized decomposition matrices added and decomposition
##     matrix type overhauled. d.d[c] is now a list of records with
##     partitions replaced by references to d.rows etc. Changed format
##     of decomposition matrix library files so as to handle polynomial
##     entries in the matrices..
##   - printing of Specht records changed allowing more compact and
##     readable notation. eg. S(1,1,1,1)->S(1^4). (see SpechtPrintFn).
##   - reversed order of parts and coeffs records in H.S(), d.d etc;
##     these lists are now sets which improves use of Position().
##   - reorganised the Specht function and the record H() which it returns; 
##     in particular added H.info.
##   - extended SInducedModule() and SRestrict to allow multiple inducing and
##     restricting (all residues).

## 1.0: December 1995: initial release.

######################################################################

## Here is a description of the structure of the main records used 
## in SPECHT

## 1. Specht()
## Specht() is the main function in specht.g, it returns a record 'H'
## which represents the family of Hecke algebras, or Schur algebras, 
## corresponding to some fixed field <R> and parameter <q>. 'H' has 
## the components:
##   IsSpecht      : is either 'true' or 'false' depending on whether
##                   'H' is a Hecke algebra or Schur algebra resp.
##   S(), P(), D() : these three functions return records which 
##                   represent Specht modules, PIMs, and simple
##                   'H' modules repectively. These functions exist
##                   only when 'H' is a Hecke algebra record.
##   W(), P(), F() : these are the corresponding functions for Weyl
##                   modules, PIMs, and simple modules of Schur algbras.
##   info          : this is a record with components
##                     version: SPECHT version number,
##                     Library: path to SPECHT library files
##                     SpechtDirectory: addition directory searched by
##                            SPECHT (defaults to current directory)
##                     Indeterminate: the indedeterminate used by SPECHT
##                            in the LLT algorithm when H.p=0.
##   operations    : apart from the obvious things like the Print()
##                   function for 'H' this record also contains the
##                   operation records for the modules S(), P() etc.
##                   as well as functions for manipulating these records
##                   and accessing decomposition matrix files. The most
##                   most important of these are:
##                     S, P, D, Pq, Sq, Dq : operations records for modules
##                     New : creation function for modules. Internally
##                       modules are created via
##                         H.operations.new(module,coeffs,parts)
##                       where module is one of "S", "P", "D", "Sq", "Pq",
##                       or "Dq" ("S" and "D" are used even for Schur
##                       algebras), coeffs is an integer or a *set* of
##                       integers and parts is a partition or a *set* of
##                       partitions. In any programs the use of New is 
##                       better than H.S(mu), for example, because the 
##                       function names are different for Hecke and Schur 
##                       algebras. Note that coeffs and parts must be 
##                       ordered reverse lexicographically (ie. they are 
##                       *sets*).
##                     Collect: like New() except that  coeffs and parts
##                       need not be sets (and may contain repeats).
##                     NewDecompositionMatrix : creates a decomposition
##                       matrix.
##                     ReadDecompositionMatrix : reads, and returns, a 
##                       decomposition matrix file.
##                     KnownDecompositionMatrix : returns a decomposition
##                       matrix of a given size; will either extract this
##                       matrix from Specht's internal lists or call
##                       ReadDecompositionMatrix(), or try to calculate 
##                       the decomposition matrix (without using the
##                       crystalized decomposition matrices).
##                     FindDecompositionMatrix : like KnownDM except that
##                       it will calculate the crystalized DM if needed.
##   Ordering      : a function for ordering partitions; controls how
##                   decomposition matrices for H are printed.
##   e             : order of <q> in <R>
##   p             : characteristic of <R>
##   valuation     : the valuation map of [JM2]; used primarily by
##                   the q-Schaper theorem.D
##   HeckeRing     : bookkeeping string used primarily in searching for
##                   library files.
##   Pq(), Sq()    : Functions for computing elements of the Fock space
##                   when H.p=0 (used in LLT algorithm). Note that there is
##                   no Dq; also unlike their counter parts S(), P(), and
##                   D() they accept only partitions as arguments.
##
## 2. The module functions S(), P() and D() (and Schur equivalents)
## These functions return record 'x' which represents some 'H'--module.
## 'x' is a record with the following components:
##   H      : a pointer back to the corresponding algebra
##   module : one of "S", "P", "D", "Sq", "Pq", or "Dq", (not "W", or "F").
##   coeffs : a *set* of coefficients
##   parts  : the corresponding *set* of partitions
##   operations :
##       + - * / : for algebric manipulations
##       Print : calls PrintModule
##       Coefficient : returns the coefficient of a given partition
##       PositiveCoefficients : true if all coefficients are non-negative
##       IntegralCoefficients : true if all coefficients are integral
##       InnerProduct : computes the 'Kronecker' inner product
##       Induce, Restrict, SInduce, SRestrict : induction and restriction
##                 functions taylored to 'x'. These functions convert 'x'
##                 to a linear combination of Specht modules, induce and
##                 then convert back to the type of 'x' (if possible).
##                 Quantized versions are applied as appropriate.
##       S, P, D : functions for rewriting 'x' into the specified type
##                 (so, for example, S('x') rewrites 'x' as a linear
##                 combination of Specht modules).
## 'x'.operations is a pointer to 'H'.operations.('x'.module).
## 
## 3. DecompositionMatrices
## Decomposition matrices 'd' in Specht are represented as records with the
## following components:
##   
##   d : a list, indexed by d.cols, where each entry is a record
##       corresponding to a column of 'd'; this record has components
##       two * sets* coeffs and parts, where parts is the index of the 
##       corresponding partition in d.rows.
##   rows : the *set* of the partitions which make up the rows or 'd'.
##   cols : the *set* of the partitions which make up the rows or 'd'.
##   inverse : a list of records containing the inverse of 'd'. These
##          records are computed only as needed.
##   dimensions : a list of the dimenions of the simle modules; again
##          comuted only as needed.
##   IsDecompositionMatrix : false if 'd' is a crystallized decomposition
##          matrix, and true otherwise.
##   H :a pointer back to the corresponding algebra
##   operations :
##     = : equality.
##     Print, TeX, Matrix : printing, TeX, and a GAP Matrix.
##     AddIndecomposable : for adding a PIM to 'd'.
##     Store : for updating Specht's internal record of 'd'.
##     S, P, D: for accessing the entries of 'd' and using 'd' to
##              convert between the various types of 'H'--modules.
##              These are actually records, each containing three
##              functions S(), P(), and D(); so X.Y() tells 'd' how
##              to write an X-module as a linear comination of Y-modules.
##     Invert : calculates D(mu) using 'd'.
##     IsNewIndecomposable : the heart of the 'IsNewIndecomposable'
##              function.
##     Induce : for inducing decomposition matrices (non--crystallized).
##   P : a short-hand for d.H.P('d',<mu>).

InstallMethod(
  Specht,
  "generates a Hecke-Algebra object",
  [IsInt],
  function(e)
    local H;
    
    H := rec(
      e:=e,
      ## bits and pieces about H
      info:=rec(version:=PackageInfo("specht")[1].Version),

      ## ordering used when printing decomposition matrices
      Ordering:=Lexicographic,
    );

    H.S:=function(arg) return NewModule(H,"S",arg); end;
    H.P:=function(arg) return NewModule(H,"P",arg); end;
    H.D:=function(arg) return NewModule(H,"D",arg); end;

    Objectify(HeckeType,H);

    return H;
  end
);

InstallMethod(OrderOfQ,"reading access to H.e",[IsAlgebraObj],
  function(H) return H!.e; end
);

InstallMethod(OrderOfQ,"reading access to x.H.e",[IsAlgebraObjModule],
  function(x) return x!.H!.e; end
);

InstallMethod(SetOrdering,"writing access to H.Ordering",
  [IsAlgebraObj,IsFunction],
  function(H,ord) H!.Ordering := ord; end
);

InstallMethod(SpechtCoefficients,"reading access to S.coeffs",[IsHeckeSpecht],
  function(S) return S!.coeffs; end
);

InstallMethod(SpechtPartitions,"reading access to S.parts",[IsHeckeSpecht],
  function(S) return S!.parts; end
);

InstallMethod(ListERegulars,[IsAlgebraObjModule],
  function(x) local e,parts,coeffs,p;
    e:=x!.H!.e;
    parts:=x!.parts;
    coeffs:=x!.coeffs;
    if e=0 then return parts;
    elif x=0*x then return [];
    else return List(Filtered([Length(parts),Length(parts)-1..1],
           p->IsERegular(e,parts[p])),p->[coeffs[p], parts[p]]);
    fi;
  end
); # ListERegulars
