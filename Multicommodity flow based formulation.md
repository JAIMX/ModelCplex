Model
===================





Parameters
-------------
**卡车相关参数**

 - $\alpha:$ fixed cost per day per truck
 - $\beta$ :   transportation cost per package per unit distance
 - $C_{k}$:   capacity of truck $k $
 - $\gamma_{k}$ :  fixed cost for the use of truck $k $
 - $L:$ max number of legs allowed to be traveled by a truck
 -  $D:$ max distance allowed to be traveled by a truck
 -  $Speed:$ average speed of trucks, if necessary it can be truck specific
 -  $DrivingTimePerDay:$ driving time per day allowed for trucks


**节点相关参数**

 - $q^{p}$: quantity of pickup and delivery demand $p $
 - $l_{i,j}$: distance of arc$(i,j) $
 - $\overline{M}:$ a sufficiently large value

***Auxiliary graph***  $\qquad G^{'}(V^{'},A^{'})$

 - $V_{0}:$ depot of all vehicles
 - $V_{st}:$ for each $u \in V \setminus V_{0}$, associate T+1 vertices: $u_{0},u_{1},...,u_{T}$
 - $A_{T}=\left \{ (u_{t},u_{t+1})|u\in  V \setminus V_{0}, t\in\left \{ 0,1,...,T-1 \right \}\right \}$
 - $\tilde{A}=\left \{ (u_{t},w_{t}+t(u,w))|(u,w)\in  A \setminus \delta (V_{0}), t\in\left \{ 0,1,...,T-t(u,w) \right \}\right \}$
 - $O,D: $  origin and destination of vehicles(depot)
 - $A^{O}=\left \{ (o_{k},u_{0})|u\in  V \setminus V_{0}, k\in K\right \}$
 - $A^{D}=\left \{ (u_{T},d_{k})|u\in  V \setminus V_{0}, k\in K\right \}$
 -  $V^{'}=V_{st}\cup \left \{ O,D \right \}$
 - $A^{'}=A_{t}\cup \tilde{A}\cup A^{O} \cup A^{D}$
 - cost: $A_{T}=0  \qquad\tilde{A}=l(u,w)$     
 

<br/>
<br/>
<br/>
<br/>

 
Decision variables
-------------

- $X_{i,j}^{k}$:=1 if arc $(i,j)$belongs to the route of vehicle $k$, otherwise 0
- $y_{i,j}^{p}$: a split of demand $q^{p}$ shipped on arc $(i,j)\in \tilde{A}\cup A_{T}$


Sets
-------------
- $V$:  set of nodes
- $A$:  set of arcs
- $K$:  set of tracks
- $P$:  set of demand O-D pairs

Indices
-------------
- $i,j$:  index of nodes
- $(i,j)$:  index of arcs
- $p$:  index of O-D pairs
- $k$:  index of tracks

Const
-------------


$$
b_{u_{t}}^{p}=
 \begin{cases}
   q^{p}  &\mbox{$u=o^{p},t=0$}\\
   -q^{p}  &\mbox{$u=d^{p},t=T$}\\
   0&\mbox{otherwise}
   \end{cases}
$$

<br/>
<br/>
<br/>
<br/>

***Minimize***
$\sum_{k\in K} \sum_{(i,j)\in \tilde{A}\cup A_{T}}\frac{\alpha l_{ij}X_{ij}^{k}}{Speed*DrivingTimePerDay}+\sum_{k\in K} \sum_{(o_{k},i)\in A^{'}}\gamma_{k}X_{o_{k},i}^{k}+\sum_{(i,j)\in \tilde{A}\cup A_{T}}\sum_{p\in P} \beta l_{ij}y_{ij}^{p}$

***Subject to:***


   $ \sum_{(j,i)\in A^{'}}X_{ji}^{k}=\sum_{(i,j)\in A^{'}}X_{ij}^{k}  \qquad\forall i\in V_{st},k\in K   \qquad$(1)
   
   $\sum_{(o_{k},i)\in A^{'}}X_{o_{k},i}^{k}\leqslant 1  \qquad \forall k\in K  \qquad$(2)

   $\sum_{(o_{k^{'}},i)\in A^{'},k^{'}\neq k}X_{o_{k^{'}},i}^{k}=0  \qquad \forall k\in K  \qquad$(3)


   $\sum_{(i,d_{k^{'}})\in A^{'},d_{k^{'}}\neq d_{k}}X_{i,d_{k^{'}}}^{k}=0  \qquad \forall k\in K  \qquad$(4)  
   
   $\sum_{(i,j)\in\tilde{A}}X_{ij}^{k}\leqslant L  \qquad\forall k\in K  \qquad$(5)
   
   $\sum_{(i,j)\in \tilde{A}}l_{ij}X_{ij}^{k}\leqslant D  \qquad\forall k\in K\qquad$(6)
   
  
   
   $\sum_{(i,j)\in \delta ^{+}(i)\setminus A^{D}} y_{ij}^{p}-\sum_{(j,i)\in \delta ^{-}(i)\setminus A^{O}} y_{ji}^{p}=b_{i}^{p} \qquad \forall p\in P, i\in V_{st}\qquad$(7)
   
$\sum_{p\in P}y_{ij}^{p} \leqslant \sum_{k\in K} C_{k}X_{ij}^{k} \qquad \forall(i,j)\in \tilde{A}\qquad$(8)
   
   $X_{ij}^{k}\in \left \{ 0,1 \right \}\qquad \forall k\in K,(i,j)\in A^{'}\qquad$(9)
   
   $y_{ij}^{p}\geqslant 0\qquad \forall p\in P,(i,j) \in \tilde{A}\cup A_{T}\qquad$(10)