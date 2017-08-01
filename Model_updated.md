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
   
   $t_{k1,i}^{a}+P_{i}\leqslant t_{k2,i}^{d}+\overline{M}(2-x_{i1,i}^{k1,st}-x_{i,i2}^{k2,st}) \qquad \forall k1,k2\in K: k1\neq k2;(i1,i),(i,i2)\in A;(s,t)\in DemandODs \qquad$(12-1)
      $t_{k2,i}^{a}+P_{i}\leqslant t_{k1,i}^{d}+\overline{M}(2-x_{i1,i}^{k1,st}-x_{i,i2}^{k2,st}) \qquad \forall k1,k2\in K: k1\neq k2;(i1,i),(i,i2)\in A, i\neq StartNode_{k1} ,i\neq StartNode_{k2} ;(s,t)\in DemandODs \qquad$(12-2)
   
   $\sum_{k\in K}\sum_{(s,j)\in A}y_{sj}^{kst}=d_{st} \qquad \forall (s,t)\in DemandODs \qquad$(13)
   
   $\sum_{k\in K}\sum_{(j,t)\in A}y_{jt}^{kst}=d_{st} \qquad \forall (s,t)\in DemandODs \qquad$(14)
   
   $\sum_{k\in K}\sum_{(i,j)\in A}y_{ij}^{kst}-\sum_{k\in K}\sum_{(j,i)\in A}y_{ji}^{kst}=0 \qquad \forall k\in K;(s,t)\in DemandODs;i\in N:i\neq s,i\neq t\qquad$(15)
   
   $y_{ij}^{kst}\leqslant \overline{M}x_{ij}^{kst} \qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(16)
   
   $x_{ij}^{kst}\leqslant X_{ijk} \qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(17)
   
   $\sum_{(s,j)\in A}x_{sj}^{kst}\leqslant 1 \qquad \forall k\in K,(s,t)\in DemandODs\qquad$(18)
   
   $\sum_{(s,t)\in DemandODs}y_{ij}^{kst}\leqslant C_{k} \qquad \forall k\in K,(i,j)\in A\qquad$(19)
   
   $X_{ijk}\in \left \{ 0,1 \right \}\qquad \forall k\in K,(i,j)\in A\qquad$(20)
   
   $x_{ijk}^{kst}\in \left \{ 0,1 \right \}\qquad \forall k\in K,(s,t)\in DemandODs,(i,j)\in A\qquad$(21)




模型的目标是最小化运输成本，目前所用的目标函数中包含两项，第一项是卡车的固定费用，第二项是与卡车具体运输路径里程相关的可变费用。模型的约束可以分为四类，第一类约束是(1)-(4)，用于为每一辆卡车构造完整的路径，值得注意的是，这些约束并不能保证剔除sub-tour。在此模型中，我们通过卡车访问节点的时间的定义以及添加卡车在各点的到达与出发时刻的关系防止了sub-tour的出现。第二类是约束(5)-(12)，用于体现每个节点的与时间相关的限制以及时间计算的关系。第三类是约束(13)-(19)，用于将运输需求分配到具体的卡车上并建立货物路径和卡车路径的关联。最后一类则是关于决策变量的0/1约束。

具体地说，约束(1)表示对每一辆卡车，其在每一点的进出弧的数目相等，结合约束(2)即卡车在出发点最多只能选择一条链接出发，所以卡车在每一点的出入弧数均为1。约束(3)和(4)是针对卡车所能行驶的最大弧数和最长距离所做的限制。约束(5)和(6)分别体现了卡车在各点的发车时刻和到达时刻的限制。与此类似，约束(7)和(8)表示各节点的操作时间窗的限制。约束(9)和(10)用于计算卡车在跑完一段弧后的时刻，而约束(11)则用于表达卡车在各点的到达时刻、发车时刻和处理所需时间的关系。约束(12)是用于保证当货物在节点处换车时的前后车的时间关系。结合约束(13)-(19)我们可以理解成在卡车的总路径（一辆或者多辆卡车）上寻找相应的能满足运输需求的最短子路径。更具体一点的解释是，约束(16)表明货物只能在子路径上运输，约束(17)表示子路径的任一组成弧必须是相应卡车的路径的一部分，约束(18)则是用于表明任一卡车最多只能运行在一条子路径上（通过限制出发子路径数量不超过1）。约束(19)体现的则是卡车的载荷约束。剩余的约束是关于决策变量的0/1约束。





