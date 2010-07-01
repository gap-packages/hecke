#######################################################################
##  SPECHT - symmcomb.g : Combinatorial functions on partitions.     ##
##                                                                   ##
##     This file contains most of the combinatorial functions used   ##
##     by Specht. Most are standard operations on Young diagrams     ##
##     or partitions.                                                ##
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

## 2.4: October 1997:
##  - added funtions MullineuxSymbol, PartitionMullineuxSymbol,
##    NormalNodes (plus undocumented friends), BetaSet, PartitionBetaSet.

## 2.2: June 1996:
##  - mainly change of function names to make it more compatible with
##    GAP conventions.

## 2.1: April 1996:
##  - added functions for finding paths in the good node partition
##    lattice.

## 2.0: March 1996 : symmcomb.g file created, breaking by specht.g
##  - added functions for finding Kleshchev's "good nodes" and implemented
##    his (and James') algorithm for the Mullineux map.

## 1.0: December 1995: initial release.

######################################################################

#F Lexicographic ordering on partitions
## Lexicographic(mu,nu); -mu and nu are lists
InstallMethod(
  Lexicographic,
  "for two partitions",
  [IsList,IsList],      
  function(lambda,mu) return lambda=mu or lambda>mu; end
);

#F LengthLexicographic(mu,nu); -mu and nu are lists
## By default this is used by DecompositionMatrix().
InstallMethod(
  LengthLexicographic,
  "for two partitions",
  [IsList,IsList],      
  function(mu,nu)
    if Length(mu)=Length(nu) then
      return mu=nu or mu>nu;
    else return Length(mu)<Length(nu);
    fi;
  end
);

#F Yet another total order. *** undocumented
InstallMethod(
  ReverseDominance,
  "for two partitions",
  [IsList,IsList],      
  function(nu,mu) local i, Mu, Nu;
    if Length(nu)=Length(mu) then
      i:=Length(mu);
      Mu:=0; Nu:=0;
      while i > 0 do
        Mu:=Mu + mu[i]; Nu:=Nu + nu[i];
        if Nu < Mu then return true;
        elif Nu > Mu then return false;
        fi;
        i:=i - 1;
      od;
    else return Length(nu)<Length(mu);
    fi;
  end
);

## dominance ordering: returns true is mu dominates, or equals, nu
#F Dominance(mu,nu);  -mu and nu are lists
InstallMethod(
  Dominance,
  "for two partitions",
  [IsList,IsList],      
  function(mu, nu) local i, m, n;
    if nu=mu then return true;
    elif Sum(nu)=0 then return true;
    elif Sum(nu)>Sum(mu) then return false;
    fi;
    m:=0; n:=0; # partial sums
    i:=1;
    while i <=Length(nu) and i<=Length(mu) do
      m:=m + mu[i]; n:=n + nu[i];
      if n > m then return false; fi;
      i:=i + 1;
    od;
    return true;
  end
);

## The coonjugate partition to arg.
#F ConjugatePartition(mu);  -mu is a sequence or a list
InstallMethod(
  ConjugatePartition,  
  "for partition",
  [IsList],      
  function(arg) local part, d, l, dl, x;
    part:=Flat(arg);
    d:=[];
    l:=Length(part);
    dl:=0;
    while l > 0 do
      if part[l] > dl then
        Append(d, List([dl+1..part[l]], x->l));
        dl:=part[l];
      fi;
      l:=l - 1;
    od;
    return d;
  end
);

#F The Littlewood-Richardson Rule.
## the algorithm has (at least), one obvious improvement in that it should
## collect like terms using something like H.operations.Collect after wrapping 
## on each row of beta.
InstallMethod(
  LittlewoodRichardsonRule,
  "for two partitions",
  [IsList,IsList],
  function(alpha, beta)
    local lrr, newlrr, x, i, j, row, Place, max, newbies;

    # place k nodes in row r>1 and above; max is the maximum number
    # of new nodes which may be added on this row and below (so
    # this is dependent upon p.new=the number of nodes added to a
    # given row from the previous row of beta).
    Place:=function(p, k, r, max) local newp, np, m, i, M;

      if r > Length(p.lam) then  # top of the partition
        Add(p.lam, k); Add(p.new, k); return [ p ]; else
        if r > 1 and p.lam[r]=p.lam[r-1] then
          max:=max + p.new[r];
          p.new[r]:=0;
          return Place(p, k, r+1, max);
        else
          if r=1 then            # m number of nodes that can be new
            m:=Minimum(k, max);  # to row r
          else m:=Minimum(p.lam[r-1]-p.lam[r], k, max);
          fi;
          if m >=0 and k-m <=p.lam[r] then  # something may fit
            newp:=[];
            for i in [0..m] do  # i nodes on row r
              if k-i <=p.lam[r] then      # remaining nodes will fit on top
                M:=max - i + p.new[r];    # uncovered nodes in previous rows
                np:=StructuralCopy(p);
                if k-i=0 and m > 0 then   # all gone
                  np.lam[r]:=np.lam[r] + i;
                  np.new[r]:=i;
                  Add(newp, np);
                else                      # more nodes can still be placed
                  np.new[r]:=i;
                  for np in Place(np, k-i, r+1, M) do
                    np.lam[r]:=np.lam[r] + i;
                    Add(newp, np);
                  od;
                fi;
              fi;
            od;
            return newp;
          fi;
          return [];
        fi;
      fi;
    end;  # end of Place; LRR internal

    if alpha=[] or alpha=[0] then return [ beta ];
    elif beta=[] or beta=[0] then return [ alpha ];
    elif Length(beta)*Sum(beta) > Length(alpha)*Sum(alpha) then
      return LittlewoodRichardsonRule(beta, alpha); 
    else
      lrr:=Place(rec(lam:=StructuralCopy(alpha),# partition
                   new:=List(alpha, i->0)),  # new nodes added from this row
                   beta[1], 1, beta[1]);
      for i in [2..Length(beta)] do
        newlrr:=[];
        for x in lrr do
          row:=1;
          while x.new[row]=0 do row:=row + 1; od;
          max:=x.new[row];
          x.new[row]:=0;
          Append(newlrr, Place(x, beta[i], row+1, max));
        od;
        lrr:=newlrr;
      od;
      return List(lrr, x->x.lam);
    fi;
  end  # LittlewoodRichardsonRule
);

#F Not used anywhere, but someone might want it. It wouldn't be too hard 
## to write something more efficient, but...
InstallMethod(
  LittlewoodRichardsonCoefficient,
  "for three partitions",
  [IsList,IsList,IsList],  
  function(lambda,mu,nu) 
    local x;
  
    if Sum(nu)<>Sum(mu)+Sum(lambda) then return 0;
    else return Length(Filtered(LittlewoodRichardsonRule(lambda,mu),x->x=nu));
    fi;
  end
);

## returns a set of the beta numbers for the partition mu
InstallMethod(
  BetaNumbers,
  "for partitions",
  [IsList],  
  function(mu)
    return mu + [Length(mu)-1, Length(mu)-2..0];
  end
);

## ALREADY AVAILABLE IN GAP4
## returns a set of the beta numbers for the partition mu
## InstallMethod(
##   BetaSet,
##   [IsList], 
##   function(mu)
##     if mu=[] then return [0];
##     else return Reversed(mu) + [0..Length(mu)-1];
##     fi;
##   end
## );

## given a beta set return the corresponding partition
InstallMethod(
  PartitionBetaSet,
  "for beta sets",
  [IsList], 
  function(beta) local i;
    if beta[Length(beta)]=Length(beta)-1 then return []; fi;
    beta:=beta-[0..Length(beta)-1];
    if beta[1]=0 then 
      beta:=beta{[First([1..Length(beta)],i->beta[i]>0)..Length(beta)]};
    fi;
    return Reversed(beta);
  end
);

