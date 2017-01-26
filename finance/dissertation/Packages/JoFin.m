(* ::Package:: *)

(*This generates the Feynman-Kac equation in n_ Dimensions 
for the martingale V_ dependend on the vector process S_, 
which solves the following system of stoch. diff. eqns:
dS[i] = a[i] dt + b[i,j]dW[j]
(Where the sum is to be performed over all
independend Brownian Motions dW)*)
FK[V_,S_,a_,b_,n_,r_]:=DFK[V,S,a,b,n]-r V@@Prepend[Table[S[i],{i,n}],t];
DFK[V_,S_,a_,b_,n_]:=D[#,t]+
Sum[D[#,S[i] ]a[i],{i,n}]+
1/2 Sum[D[#,S[i] ,S[j]]b[i,k]b[j,k],{k,n},{j,n},{i,n}]&[
V@@Prepend[Table[S[i],{i,n}],t]
];



FK[V_,MM_,r_]:=DFK[V,MM]-r V@@Prepend[Table[S[i],{i,Length[MM[[1]]]}],t];
DFK[V_,MM_]:=D[#,t]+
Sum[D[#,S[i] ]MM[[1,i]],{i,Length[MM[[1]]]}]+
1/2 Sum[D[#,S[i] ,S[j]]MM[[2,i,k]]MM[[2,j,m]]MM[[3,k,m]],
{m,Length[MM[[1]]]},{k,Length[MM[[1]]]},
{j,Length[MM[[1]]]},{i,Length[MM[[1]]]}]&[
V@@Prepend[Table[S[i],{i,Length[MM[[1]]]}],t]
]

(* This generates the SDE coefficients 
a[], b[][] and the correlation matrix \[Rho][][] for a 
m_ dimensional market model within a larger n_ dimensional model.
This can be passed to GFK*)
MM[n_,m_]:= {
Normal[SparseArray[{i_}/;i<=m->S[i] r ,{n}]],(*a*)
Normal[SparseArray[{i_,i_}/;i<=m->S[i] \[Sigma][i],{n,n}]],(*b*)
Normal[SparseArray[{{i_,i_}->1,{i_,j_}/;i<j&&j<=m->\[Rho][i,j],(*\[Rho]*)
{i_,j_}/;i>j&&i<=m->\[Rho][j,i]},{n,n}]]
}

(*This function generates the coefficients assuming 
a mutlidimensional correlated market model and requires only the 
temporal and brownian motion coefficients for the additional 
processes*)
MMc[a_,c_]:=MM[Length[Transpose[c]],Length[Transpose[c]]-Length[c]]+{
Join[Table[0,{Length[Transpose[c]]-Length[a]}],a],
Join[Table[0,{Length[Transpose[c]]-Length[a]},{Length[Transpose[c]]}],c],
0}


(*Transforms the SDE into discounted processes*)
MMdisc[MM0_]:=Module[
{
n=Length[MM0[[1]]],
MM
},
MM=Simplify[{Exp[-r t],Exp[-r t],1}*
(MM0/.Table[S[i]->S[i] Exp[r t],{i,n}])-
{r Table[S[i],{i,n}],0,0}];
Print[Text["The resulting system of SDEs:"],MatrixForm/@MM];
MM
]

(* generates SDEs for the Logs of the processes
(and removes constant time coefficients)*)
MMlog[MM0_]:=Module[
{
n=Length[MM0[[1]]],
MM=MM0/.Table[S[i]->Exp[S[i]],{i,Length[MM0[[1]]]}],
f=Table[Exp[-S[i]],{i,Length[MM0[[1]]]}],
a,
i
},
MM=Simplify[{f,f,1}*MM-{
f^2*Table[Sum[MM[[2,i,j]]MM[[2,i,k]]MM[[3,j,k]],{j,n},{k,n}],{i,n}]/2,0,0}];
a=MM[[1]];
For[i=1,i<=n,i++,If[Sum[Abs[D[a[[i]],S[j]]],{j,n}]==0,MM[[1,i]]=0]];
Print[Text["The following time coefficients have been set zero:"],MatrixForm[a-MM[[1]]]];
Print[Text["The resulting system of SDEs:"],MatrixForm/@MM];
MM
]



covariancMatrix[n_,\[Rho]_]:=Normal[SparseArray[{{i_,i_}->1,{i_,j_}/;i<j->\[Rho][i,j],(*\[Rho]*)
{i_,j_}/;i>j->\[Rho][j,i]},{n,n}]];
covolatilityMatrix[n_,\[Sigma]_,\[Rho]_]:=Normal[SparseArray[{{i_,i_}->\[Sigma][i]^2,{i_,j_}/;i<j->\[Rho][i,j]\[Sigma][i]\[Sigma][j],(*\[Rho]*)
{i_,j_}/;i>j->\[Rho][j,i]\[Sigma][i]\[Sigma][j]},{n,n}]]


ncdf[x_]:=1/2 Erfc[-(x/Sqrt[2])];
npdf[x_]:=E^(-(x^2/2))/Sqrt[2 \[Pi]];

(* Black-Scholes values and greeks with parameters:
S - Stock Price
K - Strike Price
T - Time to expiry
r - interest rate
\[Sigma] - volatility
d - dividend rate
*)
BlackScholesDPlus[S_,K_,T_,r_,\[Sigma]_,d_]:=(Log[S/K]+(r-d+\[Sigma]^2/2)T)/\[Sigma]/Sqrt[T];
BlackScholesDMinus[S_,K_,T_,r_,\[Sigma]_,d_]:=(Log[S/K]+(r-d-\[Sigma]^2/2)T)/\[Sigma]/Sqrt[T];
BlackScholesCall[S_,K_,T_,r_,\[Sigma]_,d_]:=
S Exp[-d T]ncdf[BlackScholesDPlus[S,K,T,r,\[Sigma],d]]-K Exp[-r T] ncdf[BlackScholesDMinus[S,K,T,r,\[Sigma],d]];
BlackScholesPut[S_,K_,T_,r_,\[Sigma]_,d_]:=
K Exp[-r T] ncdf[-BlackScholesDMinus[S,K,T,r,\[Sigma],d]]-S Exp[-d T] ncdf[-BlackScholesDPlus[S,K,T,r,\[Sigma],d]];
BlackScholesCallDelta[S_,K_,T_,r_,\[Sigma]_,d_]:=
Exp[-d T]ncdf[BlackScholesDPlus[S,K,T,r,\[Sigma],d]];
BlackScholesPutDelta[S_,K_,T_,r_,\[Sigma]_,d_]:=
Exp[-d T](ncdf[BlackScholesDPlus[S,K,T,r,\[Sigma],d]]-1);
BlackScholesCallGamma[S_,K_,T_,r_,\[Sigma]_,d_]:=
Exp[-d T]npdf[BlackScholesDPlus[S,K,T,r,\[Sigma],d]]/\[Sigma]/S/Sqrt[T];
BlackScholesPutGamma[S_,K_,T_,r_,\[Sigma]_,d_]:=
Exp[-d T]npdf[BlackScholesDPlus[S,K,T,r,\[Sigma],d]]/\[Sigma]/S/Sqrt[T];

(* Useage:

Quelle: siehe Kovalov, P., Linetsky, V., & Marcozzi, M. (2007). Pricing Multi-Asset American Options: A Finite Element Method-of-Lines with Smooth Penalty. Journal of Scientific Computing, 33(3), 209-237. doi:10.1007/s10915-007-9150-z
*)
GeometricAverageParams[S_,K_,T_,r_,\[Sigma]_,d_,\[Rho]_]:=Module[{\[Sigma]I2=\[Sigma].\[Rho].\[Sigma]/Length[\[Sigma]]^2},{S,K,T,r,Sqrt[\[Sigma]I2],
1/Length[\[Sigma]] Total[d+1/2 \[Sigma]^2]-1/2 \[Sigma]I2}];

IsometricGeometricAverageParams[S_,K_,T_,r_,\[Sigma]_,d_,\[Rho]_,n_]:=Assuming[\[Sigma]>0,
Simplify[
	With[{c=ConstantArray[1,n]},
		GeometricAverageParams[S,K,T,r,\[Sigma] c,d,covariancMatrix[n,\[Rho]]/.\[Rho][_,_]:>\[Rho]]]
	]
];

(*Mathematica Style*)
MM[C_,{S_,K_,T_,r_,\[Sigma]_,d_},O_]:=FinancialDerivative[C,
{"StrikePrice"-> K, "Expiration"->T},
{"InterestRate"-> r, "Volatility" ->\[Sigma], "CurrentPrice"-> S, "Dividend"->d},
##]&@@Flatten[{O}];

IsometricGeometricAverageParameters[\[Sigma]_,d_,\[Rho]_,n_]:={Sqrt[ (1- \[Rho])/n+\[Rho] ]\[Sigma],
d+(1 -\[Rho])(1- 1/n)\[Sigma]^2/2};

(*Test:
Simplify[IsometricGeometricAverageParams[S,K,T,r,\[Sigma],d,\[Rho],#][[-2;;]]-IsometricGeometricAverageParameters[\[Sigma],d,\[Rho],#]]&/@Range[8]
*)


CentralDifferences[a_]:=Differences[a][[2;;]]+Differences[a][[;;-2]]



