############################################################
##                                                        ##
##   SPECHT This file contains functions for generating   ##
##     semistandard tableaux                              ##
##                                                        ##
##   Andrew Mathas                                        ##
##                                                        ##
############################################################
 
## July 1997
##   o fixed a bug (reported by Schmuel Zelikson) in 
##     SemiStandardTableau(); added Type and Shape
##     functions for tableaux and allowed the type of
##     a tableau to be a composition which has parts 
##     which are zero.

## March 1996

#F the dual tableaux of t = [ [t_1, t_2, ...], [t_k, ...], ... ]
ConjugateTableau:=function(t)
  local  d, r, c;

  d:=List([1..Length(t[1])], r->[]);
  for r in [1..Length(t)] do
    for c in [1..Length(t[r])] do
      d[c][r]:=t[r][c];
    od;
  od;
  return d;
end;   ## ConjugateTableau

#F true if t is a semi-standard tableaux
IsSemiStandardTableau:=function(t) local r, sort;
  sort:=function(x) local X; 
    X:=ShallowCopy(x); Sort(X); return X; 
  end;

  return ForAll(t, r->r=sort(r)) and ForAll(ConjugateTableau(t), r->IsSet(r));
end;   ## IsSemiStandardTableau

#F Returns the list of semi-standard nu-tableaux of content mu.
##   Usage: SemiStandardTableaux(nu, mu) or SemiStandardTableaux(nu).
## If mu is omitted the list of all semistandard nu-tableaux is return;
## otherwise only those semistandard nu-tableaux of type mu is returned.
##   Nodes are placed recursively via FillTableau in such a way so as
## so avoid repeats.
SemiStandardTableaux:=function(arg) local FillTableau, ss, i, nu, mu;

  ## FillTableau adds upto <n>nodes of weight <i> into the tableau <t> on 
  ## or below its <row>-th row (similar to the Littlewood-Richardson rule).
  FillTableau:=function(t, i, n, row)
    local row, max, nn, nodes, nt, r;

    if row>Length(mu) then return;
    elif n=0 then          # nodes from i-th row of mu have all been placed
      if i=Length(mu) then # t is completely full
        Add(ss, t);
        return;
      else
        while i<Length(mu) and n=0 do
          i:=i + 1;          # start next row of mu
          n:=mu[i];          # mu[i] nodes to go into FillTableau
        od;
        row:=1;            # starting from the first row
      fi;
    fi;
    for r in [row..Length(t)] do
      if Length(t[r]) < nu[r] then
        if r = 1 then max:=nu[1]-Length(t[1]);
        elif Length(t[r-1]) > Length(t[r]) and t[r-1][Length(t[r]+1)]<i 
             and Length(t[r-1])>=n then
          max:=Position(t[r-1], i);
          if max=false then max:=Length(t[r-1])-Length(t[r]); #no i in t[r-1]
          else max:=max-1-Length(t[r]);
          fi;
        else max:=0;
        fi;
        max:=Minimum(max, n, nu[r]-Length(t[r]));
        nodes:=[];
        for nn in [1..max] do
          Add(nodes, i);
          nt:=Copy(t);
          Append(nt[r], nodes);
          FillTableau(nt, i, n - nn, r + 1);
        od;
      fi;
    od; 
    r:=Length(t);
    if r < Length(nu)  and n <= nu[r+1] and n <= Length(t[r]) 
    and t[r][n] < i  then
      Add(t, List([1..n], nn->i));
      if i < Length(mu) then FillTableau(t, i+1, mu[i+1], 1);
      else Add(ss, t);
      fi;
    fi;
  end;

  if Length(arg)=2 and IsList(arg[1]) and IsList(arg[2]) then
    nu:=arg[1]; mu:=arg[2];
    if Sum(nu) <> Sum(mu) then
      Error("<nu> and <mu> must be partitions of the same integer.\n");
    fi;

    ## no semi-standard nu-tableau with content mu
    if Dominates(mu, nu) then    
      if mu<>nu then return [];
      else return [ List([1..Length(mu)], i->List([1..mu[i]], ss->i)) ];
      fi;
    fi;

    ss:=[];              ## will hold the semi-standard tableaux
    FillTableau([ List([1..mu[1]],i->1) ], 2, mu[2], 1);

    return Set(ss); 
  else
    nu:=Flat(arg);
    ss:=[];
    for mu in Partitions(Sum(nu)) do
      if Dominates(mu, nu) then
        if mu = nu then
          Append(ss,[ List([1..Length(mu)], i->List([1..mu[i]], ss->i)) ]);
          return Set(ss);
        fi;
      else FillTableau([List([1..mu[1]], i->1)], 2, mu[2], 1);
      fi;
    od;
  fi;
end;   ## SemiStandardTableaux

#F Standard tableau of shape ([nu,] mu)
StandardTableaux:=function(arg) local lam, i;
  lam:=Flat(arg);
  return SemiStandardTableaux(lam, List([1..Sum(lam)], i->1) );
end;

#F return the type of the tableau <tab>; note the slight complication
## because the composition which is the type of <tab> may contain zeros.
TypeTableau:=function(tab)
  tab:=Flat(tab);
  Append(tab, [1..Maximum(tab)]);
  tab:=Collected(tab);
  return tab{[1..Length(tab)]}[2]-1;
end;

#F return the shape of the tableau <tab>
ShapeTableau:=tab->List(tab, Length);
