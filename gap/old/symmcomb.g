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
##     Andrew Mathas                                                 ##
##                                                                   ##
#######################################################################

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
Lexicographic:=function(lambda,mu) return lambda=mu or lambda>mu; end;

#F LengthLexicographic(mu,nu); -mu and nu are lists
## By default this is used by DecompositionMatrix().
LengthLexicographic:=function(mu,nu)
  if Length(mu)=Length(nu) then
    return mu=nu or mu>nu;
  else return Length(mu)<Length(nu);
  fi;
end;

#F Yet another total order. *** undocumented
ReverseDominance:=function(nu,mu) local i, Mu, Nu;
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
end; #ReverseDominance

## dominance ordering: returns true is mu dominates, or equals, nu
#F Dominates(mu,nu);  -mu and nu are lists
Dominates:=function(mu, nu) local i, m, n;
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
end;  # Dominates

## The coonjugate partition to arg.
#F ConjugatePartition(mu);  -mu is a sequence or a list
ConjugatePartition:=function(arg) local part, d, l, dl, x;
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
end; # ConjugatePartition

#F The Littlewood-Richardson Rule.
## the algorithm has (at least), one obvious improvement in that it should
## collect like terms using something like H.operations.Collect after wrapping 
## on each row of beta.
LittlewoodRichardsonRule:=function(alpha, beta)
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
              np:=Copy(p);
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
    lrr:=Place(rec(lam:=Copy(alpha),       # partition
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
end;  # LittlewoodRichardsonRule

#F Not used anywhere, but someone might want it. It wouldn't be too hard 
## to write something more efficient, but...
LittlewoodRichardsonCoefficient:=function(lambda,mu,nu) local x;
  if Sum(nu)<>Sum(mu)+Sum(lambda) then return 0;
  else return Length(Filtered(LittlewoodRichardsonRule(lambda,mu),x->x=nu));
  fi;
end;

#F the inverse Littlewood-Richardson Rule
InverseLittlewoodRichardsonRule:=function(arg)
  local initialise, fill, alpha, n, l, invlr, p, r, npp, newp, row, max, x;

  initialise:=function(p, r) local M, np, newp, i;
    if r=1 then newp:=[ ]; M:=alpha[1];
    else newp:=[ Copy(p) ]; M:=Minimum(alpha[r], p[r-1]);
    fi;
    for i in [1..M] do
      np:=Copy(p);
      np[r]:=i;
      if r < Length(alpha) then Append(newp, initialise(np, r+1));
      else Add(newp, np);
      fi;
    od;
    return newp;
  end;

  fill:=function(p, row, r, max) local m, M, np, newp, i, x;
    newp:=[];
    m:=Minimum(Minimum(p.total[r-1],alpha[r])-p.total[r], max);
    if row > 1 then m:=Minimum(m, p.mu[row-1]-p.mu[row]); fi;
    max:=max + p.new[r];
    for i in [0..m] do
      np:=Copy(p);
      np.new[r]:=i;
      np.mu[row]:=np.mu[row] + i;
      if r=Length(alpha) then np.total[r]:=np.total[r] + i; Add(newp, np);
      else
        for x in fill(np, row, r+1, max-i) do
           x.total[r]:=x.total[r] + i; Add(newp, x);
        od;
      fi;
    od;
    return newp;
  end;

  alpha:=Flat(arg);
  n:=Sum(alpha);
  invlr:=[ [ [], alpha ] ];
  for l in initialise([], 1) do
    npp:=[rec(total:=Copy(l), new:=List(alpha, r -> 0), mu:=[])];
    for r in [Length(npp[1].total)+1..Length(alpha)] do
      npp[1].total[r]:=0;
    od;
    row:=1;
    while npp<>[] do
      newp:=[];
      for p in npp do
        if row > 1 then r:=row - 1;
        else r:=1;
        fi;
        max:=0;
        while r < Length(p.total) and p.total[r]=alpha[r] do
          max:=max + p.new[r]; p.new[r]:=0; r:=r + 1;
        od;
        p.mu[row]:=alpha[r] - p.total[r];
        if row=1 or p.mu[row] <=max then
          if row=1 then max:=p.total[1];
          else max:=max + p.new[r] - p.mu[row];
          fi;
          p.new[r]:=p.mu[row];
          if r < Length(alpha) then
            for x in fill(p, row, r+1, max) do
              x.total[r]:=x.total[r] + p.mu[row];
              if Sum(x.total)=n then Add(invlr, [l, x.mu]);
              else Add(newp, x);
              fi;
            od;
          else Add(invlr, [l, p.mu]);
          fi;
        fi;
      od;
      row:=row + 1;
      npp:=newp;
    od;
  od;
  invlr[Length(invlr)]:=[alpha,[]];   ## rough hack...
  return invlr;
end;  # Inverse Littlewood-Richardson Rule

#F dimension of a Specht module
SpechtDimension:=function(arg) local Dim,y;

  Dim:=function(mu) local mud, i,j,d;
    mud:=ConjugatePartition(mu);
    d:=Factorial(Sum(mu));
    for i in [1..Length(mu)] do
      for j in [1..mu[i]] do
        d:=d/(mu[i] + mud[j] - i - j + 1);
      od;
    od;
    return d;
  end;

  if Length(arg)=1 and IsSpecht(arg[1]) then
    return Sum([1..Length(arg[1].coeffs)],
               y->arg[1].coeffs[y]*Dim(arg[1].parts[y]));
  else return Dim(Flat(arg));
  fi;
end; #SpechtDimension

## returns a set of the beta numbers for the partition mu
BetaNumbers:=function(mu)
  return mu + [Length(mu)-1, Length(mu)-2..0];
end;

## returns a set of the beta numbers for the partition mu
BetaSet:=function(mu)
  if mu=[] then return [0];
  else return Reversed(mu) + [0..Length(mu)-1];
  fi;
end;

## given a beta set return the corresponding partition
PartitionBetaSet:=function(beta) local i;
  if beta[Length(beta)]=Length(beta)-1 then return []; fi;
  beta:=beta-[0..Length(beta)-1];
  if beta[1]=0 then 
    beta:=beta{[First([1..Length(beta)],i->beta[i]>0)..Length(beta)]};
  fi;
  return Reversed(beta);
end;

## **** undocumented
## The runners for a partition on an abacus; a multiple of e-runners
## is returned  
#F EAbacusRunners(mu);  -mu is a list
EAbacusRunners:=function(e,mu) local i, j, k, aba, beta;
  aba:=List([1..e], i->[]);
  if mu=[] or mu=[0] then return aba; fi;

  ## first we find a set of beta numbers for mu; we want an e-multiple
  ## of (strictly) decreasing beta numbers for mu.
  beta:=BetaNumbers(mu);
  
  if Length(beta) mod e <> 0 then ## now add beta numbers back to get 
    i:=-Length(beta) mod e;       ## an e-multiple of beta numbers
    beta:=beta+i;
    Append(beta,[i-1,i-2..0]);
  fi;

  for i in beta do
    Add(aba[ (i mod e)+1 ], Int(i/e) );
  od;
  return aba;
end; # EAbacusRunners

#F ECore(e,mu), ECore(H,mu); -mu is a sequence or a list
##   Find the core of a partition (all partitions are 0-cores).
ECore:=function(arg) local core, beta, i, j, e;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, ECore(<H>,<mu>), ECore(<e>,<mu>)");
  fi;

  if e=0 then return Flat(arg{[2..Length(arg)]}); fi;  
  beta:=List(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})), i->Length(i));
  beta:=beta - Minimum(beta);  # remove all fully occupied rows
  if Maximum(beta)=0 then return [];
  else
    ## at present beta contains the number of beads on each runner of
    ## the abacus. next we get the beta numbers for all of the beads
    core:=[];
    for i in [1..e] do
      Append(core, List([1..beta[i]], j->e*(j-1)+i-1));
    od;
    Sort(core);
    if core[1]=0 then ## remove the beads which don't affect the beta numbers
      if core[Length(core)]=Length(core)-1 then return []; fi; #empty
      i:=First([1..Length(core)], i->core[i]<>i-1);
      core:=core{[i..Length(core)]}-i+1;
    fi;
    
    ## finally, we unravel the beta numbers of our core
    core:=List([1..Length(core)],i->core[i]-i+1);
    return core{[Length(core),Length(core)-1..1]};
  fi;
end;  # ECore

#F True is mu is an e-more. slightly better than the test mu=ECore(e,mu)
IsECore:=function(arg) local e, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, IsECore(<H>,<mu>), EWeight(<e>,<mu>)");
  fi;

  return ForAll(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})), 
                          r->r=[] or Length(r)=r[1]+1);
end;

## returns the e-weight of a partition
#F EWeight(e,mu);  -mu is a sequence or a list
## again, a slight improvement on (Sum(mu)-Sum(ECore(e,mu))/e
EWeight:=function(arg) local e, r, wt;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EWeight(<H>,<mu>), EWeight(<e>,<mu>)");
  fi;

  if e=0 then return 0;
  else return Sum(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})),
                  r->Sum(r)-Length(r)*(Length(r)-1)/2);
  fi;
end;

#F EQuotient(e,mu);  -mu is a sequence or a list
## e-quotient of a partition. algorithm based on the "star diagram"
EQuotient:=function(arg) local e, q, mu, d, i, j, qj, x;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EQuotient(<H>,<mu>), EQuotient(<e>,<mu>)");
  fi;

  if e=0 then return []; fi;
  mu:=Flat(arg{[2..Length(arg)]});
  d:=ConjugatePartition(mu);
  q:=List([1..e], j->[]);
  for i in [1..Length(mu)] do
    x:=0;
    qj:=(mu[i]-i) mod e;
    for j in [1..mu[i]] do
      if (j-d[j]-1) mod e=qj then x:=x + 1; fi;
    od;
    if x<>0 then Add(q[qj+1], x); fi;
  od;
  return q;
end;  # EQuotient

## Prints the e-abacus for the partition arg (the number of beads is
## divisible by e, and it is the smallest abacus for arg with this
## property). Pretty to look at, but useful?
#P EAbacus(mu);  -mu is a sequence or a list
EAbacus:=function(arg) local e, alpha, i, j, m;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EAbacus(<H>,<mu>), EAbacus(<e>,<mu>)");
  fi;

  alpha:=Flat(arg{[2..Length(arg)]});
  if e=0 then Print(List([1..Sum(alpha)],i->'.'),"\n");
  elif alpha=[] or alpha=[0] then
    for j in [1..e] do Print("  ."); od;
    Print("\n\n");
  else
    alpha:=EAbacusRunners(e,alpha);
    m:=Maximum(Flat(alpha)) + 1;
    for i in [0..m] do
      for j in [1..e] do
        if  i in alpha[j] then Print("  0");
        else Print("  .");
        fi;
      od;
      Print("\n");
    od;
    Print("\n");
  fi;
end;  # EAbacus

## combine a quotient and core to give a partition using abacuses.
#F CombineEQuotientECore(e,quot,core);  
##   <quot> is a list of e-partitions and <core> is a partition.
CombineEQuotientECore:=function(e, q, c) local aba, m, beta, i, j;
  if not IsInt(e) and IsBound(e.e) then e:=e.e; fi;
  if not IsInt(e) or e<>Length(q) then
    Error("usage, CombineEQuotientECore(<e>,<q>,<c>) or ",
                 "CombineEQuotientECore(<e>,<q>,<c>)");
  fi;

  aba:=EAbacusRunners(e,c); # abacus with an e-multiple of runners to which
                            # we need to add m beads to fit the quotient
  m:=Maximum(List([1..e], i->Length(q[i])-Length(aba[i])));
  m:=Maximum(m, 0);
  beta:=[];
  for i in [1..e] do
    if q[i]<>[] then
      q[i]:=q[i] + Length(aba[i]) + m;
      for j in [1..Length(q[i])] do Add(beta, (q[i][j]-j)*e + i - 1); od;
    fi;
    for j in [1..Length(aba[i])+m-Length(q[i])] do
      Add(beta, (j-1)*e + i - 1);
    od;
  od;
  Sort(beta);
  if beta[1]=0 then  ## remove irrelevant beta numbers; see ECore()
    if beta[Length(beta)]=Length(beta)-1 then return []; fi;
    m:=First([1..Length(beta)],i->beta[i]<>i-1);
    beta:=beta{[m..Length(beta)]}-m+1;
  fi;
  beta:=List([1..Length(beta)], i->beta[i]-i+1);
  return beta{[Length(beta),Length(beta)-1..1]};
end;  # CombineEQuotientECore

## true is arg is a ERegular partition
#F IsERegular(e,mu), IsERegular(H,mu)  -mu is a sequence or a list
IsERegular:=function(arg) local mu, e, i;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then  e:=arg[1].e;
  else Error("usage, IsERegular(<H>,<mu>) or IsERegular(<e>,<mu>)");
  fi;

  mu:=Flat(arg{[2..Length(arg)]});
  if e=0 then return false;
  else ## assume that mu is ordered
    e:=e-1;
    return ForAll([1..Length(mu)-e], i->mu[i]<>mu[i+e]);
  fi;
end;

#F list of the ERegular partitions of n
## ??? add support for e-regular partitions of length k ???
ERegularPartitions:=function(e,n) 
  if IsRec(e) and IsBound(e.e) then e:=e.e; 
  elif not IsInt(e) then
    Error("usage, ERegularPartitions(<H>,<mu>) or ",
          "ERegularPartitions(<e>,<mu>)");
  fi;

  if n<2 then return [ [n] ];
  elif e=0 then return Partitions(n);
  fi;

  e:=e-1;
  return Filtered(Partitions(n),
                  p->ForAll([1..Length(p)-e], i->p[i]<>p[i+e]) );
end;

#P usage: EResidueDiagram(e,mu) or EResidueDiagram(x), the second form
## returns the residue daigrams of the e-regular partitions in x
EResidueDiagram:=function(arg) local e, rs, r, x, PrintEResidueDiagram;

  PrintEResidueDiagram:=function(e,mu) local i, j;
    if mu=[] then Print("\n");
    else
      for i in [1..Length(mu)] do
        for j in [1..mu[i]] do 
          Print(String((j-i) mod e,4)); 
        od;
        Print("\n");
      od;
    fi;
  end;

  if IsSpecht(arg[1]) or (Length(arg)=2 and IsSpecht(arg[2])) then
    if IsSpecht(arg[1]) then x:=arg[1]; else x:=arg[2]; fi;
    rs:=ListERegulars(x);
    if rs=[] or IsInt(rs[1]) then PrintEResidueDiagram(x.H.e, rs);
    else
      for r in rs do
        if r[1]<>1 then Print(r[1],"*"); fi;
        Print(r[2],"\n");
        PrintEResidueDiagram(x.H.e, r[2]);
      od;
      if Length(rs) > 1 then
        Print("# There are ", Length(rs), " ", x.H.e, 
                "-regular partitions.\n");
      fi;
    fi;
  else 
    if IsInt(arg[1]) then e:=arg[1];
    elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
    else Error("usage, EResidueDiagram(<H>,<mu>), EResidueDiagram(<e>,<mu>)");
    fi;

    PrintEResidueDiagram(e,Flat(arg{[2..Length(arg)]}));
  fi;
end; # EResidueDiagram

#F Returns the partion obtained from mu by pushing nodes to the top
## of their e-ladders (see [JK]; there the notation is mu^R).
ETopLadder:=function(arg) local e, mu, ladder, r, c, C, k;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, ETopLadder(<H>,<mu>), ETopLadder(<e>,<mu>)");
  fi;

  mu:=Flat(arg{[2..Length(arg)]});
  ladder:=List(mu,r->List([1..r],c->0));

  for r in [2..Length(mu)] do  
    for c in [1..mu[r]] do
      k:=r-(e-1)*Int(r/(e-1));
      if k<1 then k:=k+e-1; fi;
      while k<r do
        C:=c+(r-k)/(e-1);
        if IsBound(ladder[k][C]) then k:=k+e-1;
        else 
          ladder[k][C]:=0; 
          Unbind(ladder[r][c]);
          k:=r;
        fi;
      od;
    od;
  od;
  return List(Filtered(ladder,r->Length(r)>0), r->Length(r));
end;  ## ETopLadder
  
#P hook lengths in a diagram mod e 
## *** undocumented: useful when lookng at the q-Schaper theorem
EHookDiagram:=function(arg) local e, mu, mud, i, j;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, EHookDiagram(<H>,<mu>), EHookDiagram(<e>,<mu>)");
  fi;

  mu:=Flat(arg{[2..Length(arg)]});
  mud:=ConjugatePartition(mu);
  for i in [1..Length(mu)] do
    for j in [1..mu[i]] do
        Print("  ", (mu[i]+mud[j]-i-j+1) mod e);
    od;
    Print("\n");
  od;
end;

#P hook length diagram
HookLengthDiagram:=function(arg) local mu, mud, i, j;
  mu:=Flat(arg);
  mud:=ConjugatePartition(mu);
  for i in [1..Length(mu)] do
    for j in [1..mu[i]] do
      Print(String(mu[i]+mud[j]-i-j+1, 4));
    od;
    Print("\n");
  od;
end; #HookLengthDiagram

#F Returns the numbers of the rows which end in one of Kleshchev's 
## "normal nodes" (see [LLT] or one of Kleshchev's papers for a description).
##   usage: NormalNodes(H|e, mu [,i]);
NormalNodes:=function(arg) local e, mu, I, normalnodes, res, i, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)=3 and IsList(arg[2]) and IsInt(arg[3]) then
    mu:=arg[2]; I:=arg[3];
  else mu:=Flat(arg{[2..Length(arg)]});
  fi;
  if not IsBound(e) or (IsBound(I) and (I<0 or I>=e) ) then
    Error("usage: NormalNodes(<e|H>, mu [, I])\n");
  fi;

  normalnodes:=List([1..e],i->[]);    ## will hold the normal nodes
  res:=List([1..e], i->0);          ## tally of #removable-#addable r-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    r:=r+1; 
    if i=Length(mu) or mu[i]>mu[i+1] then  ## removable r-node
      if res[r]=0 then Add(normalnodes[r],i); 
      else res[r]:=res[r]+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if i=1 or mu[i]<mu[i-1] then           ## addable r-node
      res[r]:=res[r]-1;
    fi;
  od;
  if IsBound(I) then return normalnodes[I+1];
  else return normalnodes;
  fi;
end;

## usage: RemoveNormalNodes(H|e, mu, i)
## returnsthe partition obtained from <mu> by removing all the normal
## nodes of residue <i>.
RemoveNormalNodes:=function(e, mu, I) local res,i,r;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not ( IsInt(e) and IsList(mu) and IsInt(I) and I<e and I>=0) then 
    Error("usage, RemoveNormalNodes(e|H, mu, I)"); 
  fi;

  mu:=Copy(mu);               ## we are going to change this so...
  res:=0;                     ## tally of #removable-#addable I-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    if r=I and (i=Length(mu) or mu[i]>mu[i+1]) then  ## removable I-node
      if res=0 then mu[i]:=mu[i]-1;                  ## normal I-node
      else res:=res+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if r=I and (i=1 or mu[i]<mu[i-1]) then           ## addable I-node
      res:=res-1;
    fi;
  od;
  return mu;
end;

#F Returns the numbers of the rows which end in one of Kleshchev's 
## "good nodes" (see [LLT] or one of Kleshchev's papers for a description).
## Basically, reading from the top down, count +1 for a *removable* node
## of residue r and -1 for an *addable* node of residue r. The last
## removable r-node with all of these tallies strictly positive is
## the (unique) good node of residue r - should it exist.
##   usage: GoodNodes(H|e, mu [,I]);
## If <I> is supplied the number of the row containing the unique good node
## of residue I is return.
GoodNodes:=function(arg) local e, mu, I, goodnodes, res, i, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)=3 and IsList(arg[2]) and IsInt(arg[3]) then
    mu:=arg[2]; I:=arg[3];
  else mu:=Flat(arg{[2..Length(arg)]});
  fi;
  if not IsBound(e) or (IsBound(I) and (I<0 or I>=e) ) then
    Error("usage: GoodNodes(<e|H>, mu [, I])\n");
  fi;

  goodnodes:=List([1..e],i->false); ## will hold the good nodes
  res:=List([1..e], i->0);          ## tally of #removable-#addable r-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    r:=r+1; 
    if i=Length(mu) or mu[i]>mu[i+1] then  ## removable r-node
      if res[r]=0 then goodnodes[r]:=i; 
      else res[r]:=res[r]+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if i=1 or mu[i]<mu[i-1] then           ## addable r-node
      res[r]:=res[r]-1;
    fi;
  od;
  if IsBound(I) then return goodnodes[I+1];
  else return goodnodes;
  fi;
end;
      
#F Given an e-regular partition mu this function returns the corresponding
## good node sequence (= path is Kleshchev's e-good partition lattice).
##   usage: GoodNodeSequence(e|H, mu);
GoodNodeSequence:=function(arg) local e, mu, goodnodeseq,  row, res, r;
  if IsInt(arg[1])  then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>,<mu>) or ",
              "GoodNodeSequence(<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]}); 
  if not IsERegular(e,mu) then
    Error("GoodNodeSequence(<e>,<mu>): <mu> must be <e>-regular\n");
  fi;
  goodnodeseq:=[];
  while mu<>[] do
    row:=1;
    while row<Length(mu) and mu[row]=mu[row+1] do
      row:=row+1;      ## there is a good node with the same residue as
    od;                ## the first removable node
    r:=(mu[row]-row) mod e;
    res:=0;
    repeat
      if r=(mu[row]-row) mod e and (row=Length(mu) or mu[row]>mu[row+1])
      then
         if res=0 then 
           if mu[row]=1 then Unbind(mu[row]);
           else mu[row]:=mu[row]-1;
           fi;
           Add(goodnodeseq, r);
         else res:=res+1; 
         fi;
       elif r=(mu[row]+1-row) mod e and mu[row]<mu[row-1] then ## addable
         res:=res-1;
       fi;
       row:=row+1;
    until row>Length(mu);
  od;
  return goodnodeseq{[Length(goodnodeseq),Length(goodnodeseq)-1..1]};
end;
 
#F Returns the list of all good node sequences for the partition <mu>
GoodNodeSequences:=function(arg) local e, mu, r, gnss, nu, s, res;
  if IsInt(arg[1]) then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>|<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  if not IsERegular(e,mu) then
    Error("GoodNodeSequence(<e>,<mu>): <mu> must be <e>-regular\n");
  fi;

  if mu=[1] then gnss:=[ [0] ]; 
  else
    gnss:=[];
    for r in GoodNodes(e,mu) do
      if r<>false then
        nu:=Copy(mu);
        nu[r]:=nu[r]-1;
        if nu[r]=0 then Unbind(nu[r]); fi;
        res:=(mu[r]-r) mod e;
        for s in GoodNodeSequences(e,nu) do
          Add(s,res); 
          Add(gnss, s);
        od;
      fi;
    od;
  fi;
  return gnss;
end;

#F Given a good node sequence this function returns the corresponding
## partition, or false if the sequence is not a good node sequence.
##   usage: GoodNodeSequence(H|e, mu)
PartitionGoodNodeSequence:=function(arg) local e, gns, mu, r, i, res, row;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, PartitionGoodNodeSequence(<H>|<e>, <gns>)");
  fi;
  gns:=Flat(arg{[2..Length(arg)]});

  mu:=[];
  for r in gns do
    row:=0;
    res:=0;
    for i in [1..Length(mu)] do
      if r=(mu[i]-i) mod e and (i=Length(mu) or mu[i]>mu[i+1]) and res<0 
      then res:=res+1;
      elif r=(mu[i]+1-i) mod e and (i=1 or mu[i]<mu[i-1]) then
        if res=0 then row:=i; fi;
        res:=res-1;
      fi;
    od;
    if res=0 and r=(-Length(mu))mod e then mu[Length(mu)+1]:=1;
    elif row>0 then mu[row]:=mu[row]+1;
    else return false;  ## bad sequence
    fi;
  od;
  return mu;
end;
      
#F GoodNodeLatticePath: returns a path in the good partition lattice
## from the empty partition to <mu>.
GoodNodeLatticePath:=function(arg) local e, gns;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, GoodNodeLatticePath(<H>|<e>, <mu>)");
  fi;
  gns:=GoodNodeSequence(e,Flat(arg{[2..Length(arg)]}));
  return List([1..Length(gns)],i->PartitionGoodNodeSequence(e,gns{[1..i]}));
end;


#F GoodNodeLatticePath: returns the list of all paths in the good partition 
## lattice from the empty partition to <mu>.
GoodNodeLatticePaths:=function(arg) local e, gns,g;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, GoodNodeLatticePaths(<H>|<e>, <mu>)");
  fi;
  gns:=GoodNodeSequences(e,Flat(arg{[2..Length(arg)]}));
  return List(gns, g->List([1..Length(g)],
              i->PartitionGoodNodeSequence(e,g{[1..i]})));
end;

#F LatticePathGoodNodeSequence()
## Returns the path in the e-good partition lattice corresponding
## to the good node sequence <gns>.
LatticePathGoodNodeSequence:=function(arg) local e, gns, i;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, GoodNodeSequenceLattice(<H>|<e>, <gns>)");
  fi;
  gns:=Flat(arg{[2..Length(arg)]});
  gns:=List([1..Length(gns)],i->PartitionGoodNodeSequence(e,gns{[1..i]}));
  if false in gns then return gns{[1..Position(gns,false)]};
  else return gns;
  fi;
end;

#F returns the Mullineux symbol of the <e>-regular partition <mu>
##   usage, MullineuxSymbol(<H>|<e>, <mu>)
## the algorithms is basically to shuffle the first column hooks lengths;
## this is a reformulation of Mullineux's approach.
## e.g. if e=3 and mu=[4,3,2] then we do the following:
##    betanums =  [6, 4, 2, 0] 
##             -> [4, 3, 1, 0] :6->4, 4->3 ( we want 2 but can only 
##                                           remove 1 more node as e=3 )
##             -> [3, 2, 1, 0].
## To get the Mullineux symbols we record the number of beads removed at
## each stage and also the number of signiciant numbers in the previous 
## beta number (i.e. the numebr of rows); here we get
##                  5, 3, 1
##                  3, 2, 1
MullineuxSymbol:=function(arg) 
  local e, mu, betaset, newbetaset, tally, difference,i,ms;

  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, MullineuxSymbol(<H>|<e>, <mu>)");
  fi;

  mu:=arg{[2..Length(arg)]};
  if mu=[] or mu=[0] then return [ [0],[0] ];
  elif IsList(mu[1]) then mu:=mu[1]; 
  fi;
  betaset:=BetaSet(mu);
  ms:=[ [],[] ];
  while betaset<>[] do
    newbetaset:=Copy(betaset);
    RemoveSet(newbetaset, newbetaset[Length(newbetaset)]);
    AddSet(newbetaset,0);
    difference:=betaset-newbetaset;
    tally:=0;
    Add(ms[1], 0);
    Add(ms[2], Length(betaset));
    for i in [Length(betaset),Length(betaset)-1..1] do
      tally:=tally+difference[i];
      if tally>=e then
        newbetaset[i]:=newbetaset[i]+tally-e;
        ms[1][Length(ms[1])]:=ms[1][Length(ms[1])]+e;
        tally:=0;
      fi;
    od;
    ms[1][Length(ms[1])]:=ms[1][Length(ms[1])]+tally;
    betaset:=newbetaset;
    if not IsSet(betaset) then return false; fi; ## can happen?
    if betaset[1]=0 then
      if betaset[Length(betaset)]=Length(betaset)-1 then 
        betaset:=[];
      else
        i:=First([1..Length(betaset)], i->betaset[i]<>i-1);
        betaset:=betaset{[i..Length(betaset)]}-i+1;
      fi;
    fi;
  od;
  return ms;
end;

#F given a Mullineux Symbol <ms> and an integer <e>, return the corresponding 
## <e>-regular partition.
PartitionMullineuxSymbol:=function(arg) local e, ms, betaset, i,tally,betaN;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)<>2 or not IsBound(e) then
    Error("usage, MullineuxSymbol(<H>|<e>, <mu>)");
  fi;
  ms:=Copy(arg[2]);

  betaset:=[0..ms[2][1]-1];
  ms[2]:=ms[2][1]-ms[2]+1;  # significant numbers in betaset
  i:=Length(ms[1]);
  while i>0 do
    tally:=0;
    betaN:=ms[2][i];
    repeat
      if tally=0 then
        tally:=ms[1][i] mod e;
        if tally=0 then tally:=e; fi;
        ms[1][i]:=ms[1][i]-tally;
      fi;
      if betaN=Length(betaset) then 
        betaset[betaN]:=betaset[betaN]+tally;
        tally:=0;
      else
        if betaset[betaN+1]-betaset[betaN]>tally then
          betaset[betaN]:=betaset[betaN]+tally;
          tally:=0;
        else
          tally:=tally-betaset[betaN+1]+betaset[betaN];
          betaset[betaN]:=betaset[betaN+1];
        fi;
      fi;
      betaN:=betaN+1;
    until (tally=0 and ms[1][i]=0) or betaN>Length(betaset);
    if tally>0 or ms[1][i]>0 then return false; fi;
    i:=i-1;
  od; ## while
  return PartitionBetaSet(betaset);
end;

#F removes the rim hook from mu which corresponding to the 
## (row,cols)-th hook.
RemoveRimHook:=function(arg) local mu, row, col, mud, r, c, x, nx;
  mu:=Copy(arg[1]);
  row:=arg[2];
  col:=arg[3];
  if Length(arg)=3 then mud:=ConjugatePartition(mu);
  else mud:=arg[4];
  fi;
  r:=mud[col];
  x:=col;
  while r >=row do
    nx:=mu[r];
    if x=1 then Unbind(mu[r]);
    else mu[r]:=x - 1;
    fi;
    x:=nx;
    r:=r - 1;
  od;
  return mu;
end;

#F Returns the partition obtained from mu by adding a rim hook with 
## foot in row <row>, of length of length <h>. The empty partition []
## is returned if the resulting diagram is not a partition.
AddRimHook:=function(nu, row, h) local r;
  nu:=Copy(nu);
  r:=row;
  if r=Length(nu) + 1 then nu[r]:=0;
  elif r > Length(nu) then h:=0;
  fi;
  while r > 1 and h > 0 do
    h:=h-nu[r-1]+nu[r]-1;
    if h > 0 then 
      nu[r]:=nu[r-1]+1; 
      r:=r-1;
    elif h < 0 then 
      nu[r]:=h+nu[r-1]+1;
    fi;
  od;
  if h > 0 then nu[1]:=nu[1] + h; r:=1;
  elif h=0 then return false;
  fi;
  return [nu, row - r];
end;
