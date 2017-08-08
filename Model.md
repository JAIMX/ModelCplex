Model
===================

Parameters
-------------
**卡车相关参数**

 - $\alpha:$ fixed cost per day per truck
 - $\beta$ :   transportation cost per package per unit distance
 - $C_{k}$:   capacity of truck $k $
 - $L:$ max number of legs allowed to be traveled by a truck
 -  $D:$ max distance allowed to be traveled by a truck
 -  $Speed:$ average speed of trucks, if necessary it can be truck specific
 -  $DrivingTimePerDay:$ driving time per day allowed for trucks


**节点相关参数**

 - $a_{i}$:   arrival time lower bound of node $i $
 - $b_{i}$:  departure time upper bound of node $i $
 - $HL_{i}$ :  open time  of node $i $
 - $HU_{i}$:  close time  of node $i $
 - $P_{i}$:  processing time at node $i $
 - $d_{s,t}$: demand quantity from node $s $ to node $t $
 - $l_{i,j}$: distance of arc$(i,j) $
 - $\overline{M}:$ a sufficiently large value


Decision variables
-------------

- $X_{i,j,k}$:=1 if arc $(i,j)$belongs to the route of truck $k$, otherwise 0
- $y_{i,j}^{k,s,t}$: a split of demand $d_{s,t}$ shipped on arc $(i,j)$ by truck  $k$
- $x_{i,j}^{k,s,t}$:=1 if arc$(i,j)$ is part of the route to ship demand from $s$ to $t$ by truck $k$
- $t_{k,i}^{a}$: arrival time at node $i$ by truck $k$
- $t_{k,i}^{d}$: departure time at node $i$ by truck $k$

Sets
-------------
- $N$:  set of nodes
- $A$:  set of arcs
- $K$:  set of tracks
- $DemandODs$:  set of demand O-D pairs

Indices
-------------
- $i,j$:  index of nodes
- $(i,j)$:  index of arcs
- $(s,t)$:  index of O-D pairs
- $k$:  index of tracks


----------


***Minimize***
$\sum_{j\in N:(j,i)\in A}\sum_{k\in K} \frac{\alpha l_{ij}X_{ijk}}{Speed*DrivingTimePerDay}+ \sum_{k\in K}\sum_{(s,t)\in demandODs}\sum_{(i,j)\in A}\beta l_{ij}y_{ij}^{kst}$

***Subject to:***


   $ \sum_{(j,i)\in A}X_{jik}=\sum_{(j,i)\in A}X_{ijk}  \qquad\forall i\in N,k\in K   \qquad$(1)
   
   $\sum_{(i,j)\in A}X_{ijk}\leqslant 1  \qquad\forall i\in N:i=StartNode_{k},k\in K  \qquad$(2)
   
   $\sum_{(i,j)\in A}X_{ijk}\leqslant L  \qquad\forall k\in K  \qquad$(3)
   
   $\sum_{(i,j)\in A}l_{ij}X_{ijk}\leqslant D  \qquad\forall k\in K\qquad$(4)
   
   $t_{ki}^{d}\leqslant b_{i}+\overline{M}(1-\sum_{(j,i)\in A}X_{jik}) \qquad\forall k\in K,i\in N\qquad$(5) 
   (t % 24)
   $t_{ki}^{a}\geqslant  a_{i}-\overline{M}(1-\sum_{(j,i)\in A}X_{jik}) \qquad\forall k\in K,i\in N\qquad$(6)
   
   $t_{ki}^{d}\leqslant HU_{i}+\overline{M}(1-\sum_{(j,i)\in A}X_{jik})\qquad\forall k\in K,i\in N \qquad$(7)
   
   $t_{ki}^{a}\geqslant HL_{i}-\overline{M}(1-\sum_{(j,i)\in A}X_{jik})\qquad\forall k\in K,i\in N \qquad$(8)
   
   $t_{kj}^{a}\leqslant t_{ki}^{d}+\frac{l_{ij}}{Speed}+\overline{M}(1-X_{ijk})\qquad\forall k\in K,(i,j)\in A  \qquad$(9)
   
   $t_{kj}^{a}\geqslant t_{ki}^{d}+\frac{l_{ij}}{Speed}-\overline{M}(1-X_{ijk})\qquad\forall k\in K,(i,j)\in A  \qquad$(10)
   
   $t_{ki}^{d}\geqslant t_{ki}^{a}+P_{i} \qquad \forall k\in K, i\in N:i\neq StartNode_{k} \qquad$(11)
   
   $t_{k1,i}^{a}+P_{i}\leqslant t_{k2,i}^{d}+\overline{M}(2-x_{i1,i}^{k1,st}-x_{i,i2}^{k2,st}) \qquad \forall k1,k2\in K: k1\neq k2;(i1,i),(i,i2)\in A;(s,t)\in DemandODs \qquad$(12)
   
   $\sum_{k\in K}\sum_{(s,j)\in A}y_{sj}^{kst}=d_{st} \qquad \forall (s,t)\in DemandODs \qquad$(13)
   
   $\sum_{k\in K}\sum_{(j,t)\in A}y_{jt}^{kst}=d_{st} \qquad \forall (s,t)\in DemandODs \qquad$(14)
   
   $\sum_{k\in K}\sum_{(i,j)\in A}y_{ij}^{kst}-\sum_{k\in K}\sum_{(j,i)\in A}y_{ji}^{kst}=0 \qquad \forall k\in K;(s,t)\in DemandODs;i\in N:i\neq s,i\neq t\qquad$(15)
   
   $y_{ij}^{kst}\leqslant \overline{M}x_{ij}^{kst} \qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(16)
   
   $x_{ij}^{kst}\leqslant X_{ijk} \qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(17)
   
   $\sum_{(s,j)\in A}x_{sj}^{kst}\leqslant 1 \qquad \forall k\in K,(s,t)\in DemandODs\qquad$(18)
   
   $\sum_{(s,t)\in DemandODs}y_{ij}^{kst}\leqslant C_{k} \qquad \forall k\in K,(i,j)\in A\qquad$(19)
   
$\sum_{(i,j)\in A}X_{ijk}\leqslant M\sum_{(s,j)\in A}X_{sjk}\qquad \forall k\in K:s=StartNode_{k}\qquad$(20)

$\sum_{(s,j)\in A}x_{sj}^{kst}+ \sum_{(i,s)\in A}x_{is}^{kst}\leqslant 1\qquad \forall k\in K:s=StartNode_{k},\forall(s,t)\in DemandODs\qquad$(21)
   
   $X_{ijk}\in \left \{ 0,1 \right \}\qquad \forall k\in K,(i,j)\in A\qquad$(22)
   
   $x_{ijk}^{kst}\in \left \{ 0,1 \right \}\qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(23)









