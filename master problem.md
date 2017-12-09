Model
===================





Parameters
-------------
**卡车相关参数**

 - $\alpha:$ fixed cost per day per truck
 - $\beta$ :   transportation cost per package per unit distance
 - $s$ :  a sequence of index of different capacities
 - $C_{s}$:   capacity of $s th$ type of truck
 - $L:$ max number of legs allowed to be traveled by a truck
 -  $D:$ max distance allowed to be traveled by a truck
 -  $Speed:$ average speed of trucks, if necessary it can be truck specific
 -  $DrivingTimePerDay:$ driving time per day allowed for trucks
 - $b_{ij}^{kos}=1$ if  $\theta_{k}^{os}$ contains arc(i,j)
 - $L_{os}:$  the number of trucks available starting from origin o with capacity of  $C_{s}$


**节点相关参数**

 - $q^{p}$: quantity of pickup and delivery demand $p $
 - $l_{i,j}$: distance of arc$(i,j) $

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

- $\theta_{k}^{os}$: the number of trucks with capacity $C_{s}$ start from depot $o$ choses the $kth$ route which $r_{k}\in \Omega _{o}$
- $y_{i,j}^{p}$: a split of demand $q^{p}$ shipped on arc $(i,j)\in \tilde{A}\cup A_{T}$


Sets
-------------
- $V$:  set of nodes
- $A$:  set of arcs
- $P$:  set of demand O-D pairs
- $O$:  set of origins
- $S$:  set of index of different capacities of trucks

Indices
-------------
- $i,j$:  index of nodes
- $(i,j)$:  index of arcs
- $p$:  index of O-D pairs

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
$\sum_{o\in O}\sum_{s\in S}\sum_{r_{k}\in \Omega _{o}} \sum_{(i,j)\in \tilde{A}\cup A_{T}} \frac{\alpha l_{ij}b_{ij}^{kos}\theta_{k}^{os}}{Speed*DrivingTimePerDay}+ \sum_{o\in O}\sum_{s\in S}\sum_{r_{k}\in \Omega _{o}} \gamma _{s}\theta _{k}^{os}+\sum_{(i,j)\in \tilde{A}\cup A_{T}}\sum_{p\in P} \beta l_{ij}y_{ij}^{p}$

***Subject to:***
  
   
   $\sum_{(i,j)\in \delta ^{+}(i)\setminus A^{D}} y_{ij}^{p}-\sum_{(j,i)\in \delta ^{-}(i)\setminus A^{O}} y_{ji}^{p}=b_{i}^{p} \qquad \forall p\in P, i\in V_{st}\qquad (1)$
   
$\sum_{p\in P}y_{ij}^{p} \leqslant \sum_{o\in O}\sum_{s\in S}\sum_{r_{k}\in \Omega _{o}} C_{s}b_{ij}^{kos}\theta_{k}^{os} \qquad \forall(i,j)\in \tilde{A}\qquad(2)$
   
   $\sum_{r_{k}\in \Omega _{o}}\theta_{k}^{os} \leqslant L_{os}\qquad \forall o\in O, \forall s\in S\qquad(3)$
   
   $y_{ij}^{p}\geqslant 0\qquad \forall p\in P,(i,j) \in \tilde{A}\cup A_{T}\qquad(4)$
   $\theta_{k}^{os}\in Z\qquad \forall o\in O, s\in S \qquad(5)$
   