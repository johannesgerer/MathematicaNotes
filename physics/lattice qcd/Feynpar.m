(* FeynmanParameter and trace were written by Todd West, Theory Group,
  Department of Physics, University of Texas, Austin, Texas 78712.
  e-mail: toddwest@utaphy.ph.utexas.edu
  Please direct all questions, comments, or bug reports to this address.
  This program was last updated on 12/16/93 and is compatible with
  Mathematica versions 1.2, 2.X.  Any publication using results
  from these programs should contain a reference to the paper 
  "FeynmanParameter and trace-programs for expressing Feynman amplitudes
  as integrals over Feynman parameters" in Computer Physics Communications
  77, 286-298 (1993).  These programs may be distributed freely.  More
  details of their use may be found in the file tutorial.txt and in the
  previously mentioned reference.*)

trace::usage = "trace[expr] or tr[expr] computes the trace of a product of gamma matrices
 where expr is a sequence of terms delimited by commas.  
Each term in expr may contain at most one gamma matrix which is represented by
 an index or by a 'five' ( for gamma 5 ).  By default, trace works in four dimensions.  
If the user wishes to work in other than four dimensions, they should type:
 dimension = d [Enter] , where d is the new number of dimensions.";

tr::usage = "trace[expr] or tr[expr] computes the trace of a product of gamma matrices
 where expr is a sequence of terms delimited by commas.  
Each term in expr may contain at most one gamma matrix which is represented by
 an index or by a 'five' ( for gamma 5 ).  By default, trace works in four dimensions.  
If the user wishes to work in other than four dimensions, they should type:
 dimension = d [Enter] , where d is the new number of dimensions.";

FeynmanParameter::usage = "FeynmanParameter[expr,k(,l(,j))] or 
fp[expr,k(,l(,j))] computes the integral of expr over k 
(divided by (2 Pi)^d), then computes the integral of the result over l 
(divided by (2 Pi)^d) if the third argument is provided and, finally, 
computes the integral of this result over j (divided by (2 Pi)^d), 
if the fourth argument is given.
  By default, FeynmanParameter works in four dimensions.  
If the user wishes to work in other than four dimensions, they should type:
dimension = d [Enter] , where d is the new number of dimensions.";

fp::usage = "FeynmanParameter[expr,k(,l(,j))] or 
fp[expr,k(,l(,j))] computes the integral of expr over k 
(divided by (2 Pi)^d), then computes the integral of the result over l 
(divided by (2 Pi)^d) if the third argument is provided and, finally, 
computes the integral of this result over j (divided by (2 Pi)^d), 
if the fourth argument is given.
  By default, FeynmanParameter works in four dimensions.  
If the user wishes to work in other than four dimensions, they should type:
dimension = d [Enter] , where d is the new number of dimensions.";

  
five::usage= "five represents the matrix gamma5.";

dimension::usage="dimension = the number of dimensions in which trace[] and
 FeynmanParameter[] will work.
   By default, dimension=4.";

g::usage = "g is the metric tensor.";

ghat::usage = "ghat is the tensor that appears in the anticommutator of gamma5 
and a gamma matrix in other than four dimensions.";

e::usage = "e is the Levi-Civita tensor.";

x::usage = "x[1], x[2], etc. are Feynman parameters.";

y::usage = "y[1], y[2], etc. are Feynman parameters.";

z::usage = "z[1], z[2], etc. are Feynman parameters.";

$Pre = ReplaceAll[#,{tr -> trace, Tr -> trace,
      fp -> FeynmanParameter, FP -> FeynmanParameter}]&

Unprotect[g,ghat,e];

(* Give the computer properties of the metric tensor, g. *)

Attributes[g] = {Orderless}

 g/: g[x_, y_]*g[x_, z_] := g[y, z] /; !NumberQ[x]
 
 g/: g[x_, y_]*(z_)[w___,x_,v___] := z[w,y,v] /; (z =!= m) &&
                              (z =!= Complex) && (!NumberQ[x])
 
 g[x_, x_] := dimension /; !NumberQ[x]

 g/: g[x_, y_]^2 := dimension /; !NumberQ[x]

(* Tell the computer about ghat, a 2 index tensor that appears in 
    traces outside of 4 dimensions. *)

Attributes[ghat] = {Orderless}

 ghat/: ghat[x_,y_]*ghat[x_,z_] := ghat[y,z] /; !NumberQ[x]

 ghat[x_, x_] := dimension - 4 /; !NumberQ[x]

 ghat/: ghat[x_, y_]^2 := dimension - 4 /; !NumberQ[x]

 ghat/: ghat[x_,y_]*e[xxx___,x_,yyy___] := 0 /; !NumberQ[x]

(* Here are some of the properties of the Levi-Civita tensor. *)

 e[w_, x_, y_, z_] := 
   Signature[{w, x, y, z}]*Apply[e, Sort[{w, x, y, z}]] /; 
    !OrderedQ[{w,x,y,z}]
 
 e/: (x_)[a_]*(x_)[b_]*e[y___, a_, z___, b_, w___] := 0
 
 e[x___, y_, z___, y_, w___] := 0

 e/: e[xx__]*e[xx__] := -24

 e/: e[w_,x_,y_,z_]*e[wp_,xp_,yp_,zp_] := (
 -g[w, zp]*g[wp, z]*g[x, yp]*g[xp, y] + 
  g[w, yp]*g[wp, z]*g[x, zp]*g[xp, y] + 
  g[w, zp]*g[wp, y]*g[x, yp]*g[xp, z] - 
  g[w, yp]*g[wp, y]*g[x, zp]*g[xp, z] + 
  g[w, zp]*g[wp, z]*g[x, xp]*g[y, yp] - 
  g[w, xp]*g[wp, z]*g[x, zp]*g[y, yp] - 
  g[w, zp]*g[wp, x]*g[xp, z]*g[y, yp] + 
  g[w, wp]*g[x, zp]*g[xp, z]*g[y, yp] - 
  g[w, yp]*g[wp, z]*g[x, xp]*g[y, zp] + 
  g[w, xp]*g[wp, z]*g[x, yp]*g[y, zp] + 
  g[w, yp]*g[wp, x]*g[xp, z]*g[y, zp] - 
  g[w, wp]*g[x, yp]*g[xp, z]*g[y, zp] - 
  g[w, zp]*g[wp, y]*g[x, xp]*g[yp, z] + 
  g[w, xp]*g[wp, y]*g[x, zp]*g[yp, z] + 
  g[w, zp]*g[wp, x]*g[xp, y]*g[yp, z] - 
  g[w, wp]*g[x, zp]*g[xp, y]*g[yp, z] - 
  g[w, xp]*g[wp, x]*g[y, zp]*g[yp, z] + 
  g[w, wp]*g[x, xp]*g[y, zp]*g[yp, z] + 
  g[w, yp]*g[wp, y]*g[x, xp]*g[z, zp] - 
  g[w, xp]*g[wp, y]*g[x, yp]*g[z, zp] - 
  g[w, yp]*g[wp, x]*g[xp, y]*g[z, zp] + 
  g[w, wp]*g[x, yp]*g[xp, y]*g[z, zp] + 
  g[w, xp]*g[wp, x]*g[y, yp]*g[z, zp] - 
  g[w, wp]*g[x, xp]*g[y, yp]*g[z, zp])
 
(* trace computes traces of products of gamma matrices. *)

 trace[argument1__] :=
   Block[{finalresult},

 (* Use different algorithms depending on the number of dimensions. *)

      finalresult = If[dimension === 4, privateouterloop[{argument1}],
         privateouterloopdim[{argument1}]];

 (* The following steps simplify the result. *)

finalresult = ExpandAll[finalresult];
finalresult = finalresult /. x_[ae___,be_Symbol,ce___] y_[de___,be_Symbol
 ,ee___] :> x[ae,alpha,ce] y[de,alpha,ee] /; x =!= m && y =!= m;
finalresult = finalresult /. x_[ae___,be_Symbol,ce___] y_[de___,be_Symbol
 ,ee___] :> x[ae,beta,ce] y[de,beta,ee] /; be =!= alpha && x =!= m && y =!= m;
finalresult = finalresult /. x_[ae___,be_Symbol,ce___] y_[de___,be_Symbol
 ,ee___] :> x[ae,gamma,ce] y[de,gamma,ee] /; (be =!= alpha) && 
  (be =!= beta) && x =!= m && y =!= m ;    
finalresult = finalresult /. x_[ae___,be_Symbol,ce___] y_[de___,be_Symbol
 ,ee___] :> x[ae,delta,ce] y[de,delta,ee] /; (be =!= alpha) && 
  (be =!= beta) && (be =!= gamma) && x =!= m && y =!= m;
finalresult = finalresult /. x_[ae_]^2 :> x^2 /; x =!= m 
	]

(* Implement linearity of traces *)

privateouterloop[{z___,c_ x_Symbol,y___}]:= c privateouterloop[{z,x,y}]

privateouterloop[{z___,c_ (1+five),y___}] :=c privateouterloop[{z,1+five,y}]

privateouterloop[{z___,c_ (1-five),y___}] :=c privateouterloop[{z,1-five,y}]

privateouterloop[{z___,x_+y_,w___}] := 
     privateouterloop[{z,x,w}]+privateouterloop[{z,y,w}] /; 
          (x+y =!= 1+five)&&(x+y =!= 1-five)

privateouterloop[{x___,c_,y___}]:=
  c privateouterloop[{x,y}] /;(Head[c] =!= Plus)&&(Head[c] =!= Symbol)

(* The trace of an odd number of gamma matrices is 0. *)

privateouterloop[gammalist_]:= 0 /; OddQ[Length[gammalist] -
     Count[gammalist,five,Infinity]]
  
(* Use a different procedure if the list of gamma matrices contains 
   a gammafive. *)

privateouterloop[gammalist_] := If[Count[gammalist,five,Infinity] == 0,
      privatetrsub[gammalist], privatetrsub5[gammalist]]

(* Here are some simple traces.  The goal of privatetrsub is to get to
   one of these limiting cases. *)

privatetrsub/: privatetrsub[{}] = 4

privatetrsub[{w_,x_}] := 4 g[w,x] 

privatetrsub[{w_,x_,y_,z_}] := 4 g[w,x] g[y,z] - 4 g[w,y] g[x,z] + 
                                4 g[w,z] g[x,y]

(* Apply a recursive formula to reduce the number of gamma matrices
   in the trace. *)

privatetrsub[gammalist_] := Block[{gammanumber}, 
  Sum[(-1)^gammanumber*g[gammalist[[1]], gammalist[[gammanumber]]]*
     privatetrsub[Drop[Rest[gammalist], {gammanumber - 1, gammanumber - 1}]], 
     {gammanumber, 2, Length[gammalist]}]]  

(* Use {gammafive,gamma}=0 to simplify traces containing gammafive or 
1 +/- gammafive. *)

privatetrsub5[{xxx___,1+five,yyy___,1+five,zzz___}] := If[
  OddQ[Length[{yyy}]],0,2 privatetrsub5[{xxx,1+five,yyy,zzz}]] /; 
    FreeQ[{yyy},five]

privatetrsub5[{xxx___,1-five,yyy___,1-five,zzz___}] :=If[
  OddQ[Length[{yyy}]],0,2 privatetrsub5[{xxx,1-five,yyy,zzz}]] /; 
    FreeQ[{yyy},five]

privatetrsub5[{xxx___,1+five,yyy___,1-five,zzz___}] :=If[
  OddQ[Length[{yyy}]],2 privatetrsub5[{xxx,1+five,yyy,zzz}],0] /; 
    FreeQ[{yyy},five]

privatetrsub5[{xxx___,1-five,yyy___,1+five,zzz___}] :=If[
  OddQ[Length[{yyy}]],2 privatetrsub5[{xxx,1-five,yyy,zzz}],0] /; 
    FreeQ[{yyy},five]

privatetrsub5[{xxx___,five,yyy___,five,zzz___}] :=If[
  OddQ[Length[{yyy}]],-1,1] privatetrsub5[{xxx,yyy,zzz}] /;
    FreeQ[{yyy},five]

privatetrsub5[{xxx___,five,yyy___,1+five,zzz___}] := If[
  OddQ[Length[{yyy}]],-privatetrsub5[{xxx,1-five,yyy,zzz}] ,
    privatetrsub5[{xxx,1+five,yyy,zzz}]] /; FreeQ[{yyy},five]

privatetrsub5[{xxx___,five,yyy___,1-five,zzz___}] :=If[
  OddQ[Length[{yyy}]],privatetrsub5[{xxx,1+five,yyy,zzz}],
    -privatetrsub5[{xxx,1-five,yyy,zzz}] ] /; FreeQ[{yyy},five]

privatetrsub5[{xxx___,1+five,yyy___,five,zzz___}] := If[
  OddQ[Length[{yyy}]],-privatetrsub5[{xxx,1+five,yyy,zzz}] ,
    privatetrsub5[{xxx,1+five,yyy,zzz}]] /; FreeQ[{yyy},five]

privatetrsub5[{xxx___,1-five,yyy___,five,zzz___}] :=If[
  OddQ[Length[{yyy}]],privatetrsub5[{xxx,1-five,yyy,zzz}],
    -privatetrsub5[{xxx,1-five,yyy,zzz}] ] /; FreeQ[{yyy},five]

privatetrsub5[{xxx___,five,yyy___}] := If[
  OddQ[Length[{xxx}]],-1,1]*privatefinish5[{five,xxx,yyy}]

privatetrsub5[{xxx___,1+five,yyy___}] := privatetrsub[
  {xxx,yyy}] + If[OddQ[Length[{xxx}]],-1,1]*privatefinish5[{five,xxx,yyy}]

privatetrsub5[{xxx___,1-five,yyy___}] := privatetrsub[
{xxx,yyy}] - If[OddQ[Length[{xxx}]],-1,1]*privatefinish5[{five,xxx,yyy}]

(* If the trace no longer contains any gammafives, switch to the procedures
   for traces that do not involve gammafives. *)

privatetrsub5[{xxx___}] := 
   privatetrsub[{xxx}] /; FreeQ[{xxx},five]

(* Use some simple limiting cases for traces of gammafive times 4 
    or fewer gamma matrices. *)

privatefinish5[{five,xx_,yy_,zz_,ww_}] := -4I e[xx,yy,zz,ww]

privatefinish5[{five,xxx___}] := 0 /; Length[{xxx}] < 4

(* If there are more than 4 gamma matrices multiplying the gammafive,
    apply a recursive formula which reduces the number of gamma matrices
     by 2. *)

privatefinish5[{five,gamma2_,gamma3_,gamma4_,gamman__}] :=
  Block[{sigma1}, g[gamma2,gamma3]*
       privatefinish5[{five,gamma4,gamman}] + g[gamma3,gamma4]* 
       privatefinish5[{five,gamma2,gamman}] -
       g[gamma2,gamma4] privatefinish5[{five,gamma3,gamman}] - 
       I e[gamma2,gamma3,gamma4,sigma1=Unique["zzzzzzz"]]* 
       privatetrsub[{sigma1,gamman}]]

(* The next several lines perform traces in other than four dimensions. *)

privateouterloopdim[{z___,c_ x_Symbol,y___}]:= 
  c privateouterloopdim[{z,x,y}]

privateouterloopdim[{z___,c_ (1+five),y___}] := 
  c privateouterloopdim[{z,1+five,y}]

privateouterloopdim[{z___,c_ (1-five),y___}] := 
  c privateouterloopdim[{z,1-five,y}]

privateouterloopdim[{z___,x_+y_,w___}] := privateouterloopdim[{z,x,w}]+
  privateouterloopdim[{z,y,w}] 

privateouterloopdim[{x___,c_,y___}]:=
  c privateouterloopdim[{x,y}] /; Head[c] =!= Symbol

privateouterloopdim[xmass_]:= 0 /; (OddQ[Length[xmass] - 
                                     Count[xmass,five,Infinity]])

(* If there are no gammafives, we can use the dimension 4 procedures. *)

 privateouterloopdim[gammalist_]:= If[
  Count[gammalist,five,Infinity] == 0,
    privatetrsub[gammalist],privatetrsubdim5[gammalist]]

(* Apply some limiting cases. *)

 privatetrsubdim5[{five,xx_,yy_,zz_,ww_}] := -4I e[xx,yy,zz,ww]

 privatetrsubdim5[gammafivelist_] := 0 /; 
       (Length[gammafivelist] < 5 && Count[gammafivelist,five] == 1)

(* Use {gammafive,gamma[a]} = 2 ghat[b,a] five gamma[b] and 
   gammafive^2 = 1 to reduce the number of gammafive's. *)

 privatetrsubdim5[{xxx___,five,yyy___,five,zzz___}]:=
    Block[{temp1={yyy},index,unique},
	    Product[-g[ind={yyy}[[index]],
		  (unique[index]=temp1[[index]]=
		  Unique[ToString[ind]])] + 2 ghat[ind,unique[index]],
		  {index,1,Length[{yyy}]}] privatetrsubdim5[Flatten[
                                             {xxx,temp1,zzz}]]]

(* Bring any remaining gammafive to the beginning of the trace. *)

 privatetrsubdim5[{xxx__,five,yyy___}] := privatetrsubdim5[{five,yyy,xxx}]

(* Use a recursive relation to reduce the number of gamma matrices 
    in the trace until some limiting case is reached. *)

 privatetrsubdim5[{five,fivefree__}] :=
   Block[{ct,iter,jter},
     ct=Length[{fivefree}];
     Sum[(-1)^(iter+jter+1) 2 g[{fivefree}[[jter]],{fivefree}[[iter]]]*
	privatetrsubdim5[Flatten[{five,Drop[Drop[{fivefree},{iter,iter}],
	{jter,jter}]}]]/(ct-4), {iter,2,ct}, {jter,1,iter-1}]]

(* If privatetrsubdim5 makes it to this point, then there are no gammafive's 
   remaining in the trace.  Therefore we can use the dimension 4 procedures at
   this point. *)

 privatetrsubdim5[nofives_] := privatetrsub[nofives]

(* This is the end of the section which computes traces of products 
   of gamma matrices.  Next we'll convert integrals over momenta 
   space to integrals over Feynman parameters. *)


(* FeynmanParameter is the procedure which converts integrals 
   over momentum space into integrals over Feynman parameters.  This 
   procedure calls privatefeynmanparamter once for each momentum 
   integral to be so converted. *)


 FeynmanParameter[input1_,input2_,input3_:Null,input4_:Null]:=
     Block[
     {integrator, intermediateresult1, intermediateresult2},
       integrator = If[FreeQ[input1,x],
         x, If[FreeQ[input1,y], y, z]];
       intermediateresult1=privatefeyparamter[input1,input2];
       integrator = If[integrator === x, y, z];
       intermediateresult2=intermediateresult1;
       If[input3 =!= Null,intermediateresult2=
         privatefeyparamter[intermediateresult1,input3]];
       integrator=z;
       intermediateresult1=intermediateresult2;
       If[input4 =!= Null,intermediateresult1=
         privatefeyparamter[intermediateresult2,input4]];
       intermediateresult1 /.  vector1_[_]^2 :> vector1^2 /;
	(vector1 =!= m && vector1 =!= x && vector1 =!= y && vector1 =!= z)
            
    ]


(* The following procedure is called once for each momentum 
     integral that we are converting. *)

privatefeyparamter[ex_, kay_] := 
 Block[{denom = 1, nokay, count, coefofksquared, enn, finaldenominator,
        power, coefficientlist, ell, ellnumerator, elldenominator, msquared,
        msquaredminusellsquared, numerator, eyezero, alphaenn}, 
   denominator = Denominator[ex]; 

(* Pick out terms in the denominator which depend on the momentum
over which we're integrating (kay). *)

   nokay = Map[If[FreeQ[#1,kay],#1,(denom *= #1 ; 1)]&,denominator];
    enn=Length[denom];

(* Take care of terms in the denominator which are of the form
(kay^2+2p kay+M^2)^power , where power is not equal to 1. *)

    Block[{iter},
      Do[If[Head[denom[[iter]]] ===
        Power,(power[iter]=denom[[iter,2]];denom[[iter]]=
	denom[[iter,1]];),power[iter]=1
	 ],
      {iter,enn}
      ]];

(* Expand each term in the denominator. *)

    denom = Map[Expand[#1] & , denom, 2];                

(* If a term in the denominator is of the form (a kay^2+2p kay+M^2), 
with a not equal to 1, divide this term by a. *)

  Block[{iter},
    Do[If[(coefofksquared=Coefficient[denom[[iter]],kay,2])=!=1,(nokay *= 
       coefofksquared^power[iter];denom[[iter]]=Apart[denom[[iter]]/
       coefofksquared];)] , 
    {iter,enn}]];

(* Introduce Feynman parameters (x[i], y[i], z[i]). *)

    finaldenominator= denom[[1]]+Sum[(denom[[iter]]-denom[[iter-1]])*
      integrator[iter-1],{iter,2,enn}]; 
    coefficientlist=CoefficientList[finaldenominator,kay];

(* The new denominator can be written in the form
(kay^2+2 ell kay+msquared). *)

 ell=Factor[coefficientlist[[2]]/2];
 ellnumerator=Expand[Numerator[ell]];
 elldenominator=Denominator[ell];
 msquared=coefficientlist[[1]];
 msquaredminusellsquared = msquared - Expand[ell ell];
 numerator=Expand[Numerator[ex]];

(* eyezero is approximately equal to the result we get if the 
numerator of the expression which we are converting does
not contain any factors of kay. *)

 eyezero= (I(-1)^(dimension/2)/(2^dimension Pi^(dimension/2)))*
   Gamma[(alphaenn=Block[{iterator,firstletter},
     Sum[power[iterator],{iterator,enn}]])-dimension/2]/
     Block[{iterator},Product[(power[iterator]-1)!,{iterator,enn}]]*
     integrator[enn-1]^(power[enn]-1)*
     Block[{iterator},Product[(integrator[iterator-1]-
             integrator[iterator])^(power[iterator]-1),
	   {iterator,2,enn-1}]] (1-integrator[1])^
     (power[1]-1) msquaredminusellsquared^(dimension/2-alphaenn);

(* The rest of this procedure handles any factors of kay which 
may appear in the numerator (up to terms of the form kay^5). *)

numerator = numerator //. kay^2 :> kay[Unique["dumb$"]]^2;
numerator = numerator //. kay mom_Symbol xxxxxx___ :> 
              kay[dumb=Unique["dumb$"]] mom[dumb] xxxxxx;
numerator=numerator /. kay[xx1_] kay[yy2_] kay[zz3_] kay[ww4_] kay[vv5_]:> 
 -(privateel[xx1] privateel[yy2] privateel[zz3] privateel[ww4]privateel[vv5]+
  (g[xx1,yy2]privateel[zz3]privateel[ww4]privateel[vv5]+
   g[xx1,zz3]privateel[yy2]privateel[ww4]privateel[vv5]+
   g[yy2,zz3]privateel[xx1]privateel[ww4]privateel[vv5]+
   g[xx1,ww4]privateel[zz3]privateel[yy2]privateel[vv5]+
   g[yy2,ww4]privateel[xx1]privateel[zz3]privateel[vv5]+
   g[zz3,ww4]privateel[xx1]privateel[yy2]privateel[vv5]+
   g[xx1,vv5]privateel[zz3]privateel[ww4]privateel[yy2]+
   g[yy2,vv5]privateel[xx1]privateel[ww4]privateel[zz3]+
   g[zz3,vv5]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[ww4,vv5]privateel[xx1]privateel[yy2]privateel[zz3])*
     msquaredminusellsquared/(2 alphaenn-2-dimension)+
  (g[xx1,yy2]g[zz3,ww4]privateel[vv5]+g[xx1,zz3]g[yy2,ww4]privateel[vv5]+
   g[xx1,ww4]g[yy2,zz3]privateel[vv5]+g[xx1,yy2]g[zz3,vv5]privateel[ww4]+
   g[xx1,yy2]g[ww4,vv5]privateel[zz3]+g[xx1,zz3]g[yy2,vv5]privateel[ww4]+
   g[xx1,zz3]g[ww4,vv5]privateel[yy2]+g[yy2,zz3]g[xx1,vv5]privateel[ww4]+
   g[yy2,zz3]g[ww4,vv5]privateel[xx1]+g[xx1,ww4]g[yy2,vv5]privateel[zz3]+
   g[xx1,ww4]g[zz3,vv5]privateel[yy2]+g[yy2,ww4]g[xx1,vv5]privateel[zz3]+
   g[yy2,ww4]g[zz3,vv5]privateel[xx1]+g[zz3,ww4]g[xx1,vv5]privateel[yy2]+
   g[zz3,ww4]g[yy2,vv5]privateel[xx1])msquaredminusellsquared^2/
    ((2 alphaenn-2-dimension)(2 alphaenn-4-dimension)));

numerator=numerator /. kay[xx1_] kay[yy2_] kay[zz3_] kay[ww4_]^2 :> 
  -(privateel[xx1] privateel[yy2] privateel[zz3] privateel[ww4]*
    privateel[ww4]+(g[xx1,yy2]privateel[zz3]privateel[ww4]privateel[ww4]+
   g[xx1,zz3]privateel[yy2]privateel[ww4]privateel[ww4]+
   g[yy2,zz3]privateel[xx1]privateel[ww4]privateel[ww4]+
   g[xx1,ww4]privateel[zz3]privateel[yy2]privateel[ww4]+
   g[yy2,ww4]privateel[xx1]privateel[zz3]privateel[ww4]+
   g[zz3,ww4]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[xx1,ww4]privateel[zz3]privateel[ww4]privateel[yy2]+
   g[yy2,ww4]privateel[xx1]privateel[ww4]privateel[zz3]+
   g[zz3,ww4]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[ww4,ww4]privateel[xx1]privateel[yy2]privateel[zz3])*
     msquaredminusellsquared/(2 alphaenn-2-dimension)+
  (g[xx1,yy2]g[zz3,ww4]privateel[ww4]+g[xx1,zz3]g[yy2,ww4]privateel[ww4]+
   g[xx1,ww4]g[yy2,zz3]privateel[ww4]+g[xx1,yy2]g[zz3,ww4]privateel[ww4]+
   g[xx1,yy2]g[ww4,ww4]privateel[zz3]+g[xx1,zz3]g[yy2,ww4]privateel[ww4]+
   g[xx1,zz3]g[ww4,ww4]privateel[yy2]+g[yy2,zz3]g[xx1,ww4]privateel[ww4]+
   g[yy2,zz3]g[ww4,ww4]privateel[xx1]+g[xx1,ww4]g[yy2,ww4]privateel[zz3]+
   g[xx1,ww4]g[zz3,ww4]privateel[yy2]+g[yy2,ww4]g[xx1,ww4]privateel[zz3]+
   g[yy2,ww4]g[zz3,ww4]privateel[xx1]+g[zz3,ww4]g[xx1,ww4]privateel[yy2]+
   g[zz3,ww4]g[yy2,ww4]privateel[xx1])msquaredminusellsquared^2/
      ((2 alphaenn-2-dimension)(2 alphaenn-4-dimension)));

numerator=numerator /. kay[xx1_] kay[yy2_]^2 kay[ww4_]^2 :> 
-(privateel[xx1] privateel[yy2] privateel[yy2] privateel[ww4]privateel[ww4]+
  (g[xx1,yy2]privateel[yy2]privateel[ww4]privateel[ww4]+
   g[xx1,yy2]privateel[yy2]privateel[ww4]privateel[ww4]+
   g[yy2,yy2]privateel[xx1]privateel[ww4]privateel[ww4]+
   g[xx1,ww4]privateel[yy2]privateel[yy2]privateel[ww4]+
   g[yy2,ww4]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[yy2,ww4]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[xx1,ww4]privateel[yy2]privateel[ww4]privateel[yy2]+
   g[yy2,ww4]privateel[xx1]privateel[ww4]privateel[yy2]+
   g[yy2,ww4]privateel[xx1]privateel[yy2]privateel[ww4]+
   g[ww4,ww4]privateel[xx1]privateel[yy2]privateel[yy2])*
      msquaredminusellsquared/(2 alphaenn-2-dimension)+
  (g[xx1,yy2]g[yy2,ww4]privateel[ww4]+g[xx1,yy2]g[yy2,ww4]privateel[ww4]+
   g[xx1,ww4]g[yy2,yy2]privateel[ww4]+g[xx1,yy2]g[yy2,ww4]privateel[ww4]+
   g[xx1,yy2]g[ww4,ww4]privateel[yy2]+g[xx1,yy2]g[yy2,ww4]privateel[ww4]+
   g[xx1,yy2]g[ww4,ww4]privateel[yy2]+g[yy2,yy2]g[xx1,ww4]privateel[ww4]+
   g[yy2,yy2]g[ww4,ww4]privateel[xx1]+g[xx1,ww4]g[yy2,ww4]privateel[yy2]+
   g[xx1,ww4]g[yy2,ww4]privateel[yy2]+g[yy2,ww4]g[xx1,ww4]privateel[yy2]+
   g[yy2,ww4]g[yy2,ww4]privateel[xx1]+g[yy2,ww4]g[xx1,ww4]privateel[yy2]+
   g[yy2,ww4]g[yy2,ww4]privateel[xx1])*msquaredminusellsquared^2/
       ((2 alphaenn-2-dimension)(2 alphaenn-4-dimension)));

numerator=numerator /. kay[xx1_] kay[yy2_] kay[zz3_] kay[ww4_]:>
(privateel[xx1] privateel[yy2] privateel[zz3] privateel[ww4]+
  (g[xx1,yy2]privateel[zz3]privateel[ww4]+
   g[xx1,zz3]privateel[yy2]privateel[ww4]+
   g[yy2,zz3]privateel[xx1]privateel[ww4]+
   g[xx1,ww4]privateel[zz3]privateel[yy2]+
   g[yy2,ww4]privateel[xx1]privateel[zz3]+
   g[zz3,ww4]privateel[xx1]privateel[yy2])*
      msquaredminusellsquared/(2 alphaenn-2-dimension)+
  (g[xx1,yy2]g[zz3,ww4]+g[xx1,zz3]g[yy2,ww4]+g[xx1,ww4]g[yy2,zz3])*
      msquaredminusellsquared^2/((2 alphaenn-2-dimension)*
      (2 alphaenn-4-dimension)));
   
numerator=numerator /. kay[xx1_] kay[yy2_] kay[zz3_]^2 :>
(privateel[xx1] privateel[yy2] privateel[zz3] privateel[zz3]+
  (g[xx1,yy2]privateel[zz3]privateel[zz3]+
   g[xx1,zz3]privateel[yy2]privateel[zz3]+
   g[yy2,zz3]privateel[xx1]privateel[zz3]+
   g[xx1,zz3]privateel[zz3]privateel[yy2]+
   g[yy2,zz3]privateel[xx1]privateel[zz3]+
   g[zz3,zz3]privateel[xx1]privateel[yy2])*
      msquaredminusellsquared/(2 alphaenn-2-dimension)+
  (g[xx1,yy2]g[zz3,zz3]+g[xx1,zz3]g[yy2,zz3]+g[xx1,zz3]g[yy2,zz3])*
     msquaredminusellsquared^2/((2 alphaenn-2-dimension)*
     (2 alphaenn-4-dimension)));

numerator=numerator /. kay[xx1_]^2 kay[zz3_]^2 :> 
(privateel[xx1] privateel[xx1] privateel[zz3]privateel[zz3]+
  (g[xx1,xx1]privateel[zz3]privateel[zz3]+
   g[xx1,zz3]privateel[xx1]privateel[zz3]+
   g[xx1,zz3]privateel[xx1]privateel[zz3]+
   g[xx1,zz3]privateel[zz3]privateel[xx1]+
   g[xx1,zz3]privateel[xx1]privateel[zz3]+
   g[zz3,zz3]privateel[xx1]privateel[xx1])msquaredminusellsquared/
      (2 alphaenn-2-dimension)+
  (g[xx1,xx1]g[zz3,zz3]+g[xx1,zz3]g[xx1,zz3]+g[xx1,zz3]g[xx1,zz3])*
      msquaredminusellsquared^2/((2 alphaenn-2-dimension)*
      (2 alphaenn-4-dimension)));

numerator=numerator /. kay[xx1_] kay[yy2_] kay[zz3_]:> 
-(privateel[xx1] privateel[yy2] privateel[zz3]+
  (1/2)(g[xx1,yy2]privateel[zz3]+g[xx1,zz3]privateel[yy2]+
  g[yy2,zz3]privateel[xx1])msquaredminusellsquared/(alphaenn-1-dimension/2));

numerator = numerator /. kay[xx1_]^2 kay[zz3_] :> 
-(privateel[xx1]^2privateel[zz3]+(1/2)(g[xx1,xx1]privateel[zz3]+
 privateel[zz3]+privateel[zz3])msquaredminusellsquared/
    (alphaenn-1-dimension/2));

numerator=numerator /. kay[xx1_] kay[yy2_] :> (privateel[xx1] privateel[yy2]+
   g[xx1,yy2] msquaredminusellsquared/(2(alphaenn-1-dimension/2))) ;

numerator=numerator /. kay[xx1_]^2 :> (privateel[xx1]^2+g[xx1,xx1]*
   msquaredminusellsquared/(2(alphaenn-1-dimension/2))) ;

numerator=numerator /. kay[xx1_] -> -privateel[xx1];
   
numerator*eyezero/nokay]

(* privateel[x] appends the Lorentz index x to momenta in ell *)

privateel[xxx2_] := Block[{temporary},temporary=ellnumerator /.
       xxx3___ yyy3_Symbol :> xxx3 yyy3[xxx2];
       Map[If[Head[#] === Symbol,#[xxx2],#]&,temporary,{1}]/elldenominator]


Protect[g,ghat,e];

(* By default, the above procedures work in 4 dimensions. *)

dimension = 4;


