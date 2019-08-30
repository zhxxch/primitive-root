原根 - primitive root
======================

**********************

## 寻找GF(p)的原根

如果g是原根，群的阶数为m，其他原根g<sup>s</sup>满足gcd(s, m-1)=1.

p是素数，有限域GF(p)的原根有ϕ(ϕ(p))个。

ϕ(p)=p-1.

如果元素w使w<sup>ϕ(p)/p<sub>i</sub></sup> = 1(mod p)对ϕ(p)的某个因子p<sub>i</sub>成立则w不是原根。

**********************

## 寻找GF(p<sup>n</sup>)的本原多项式

本原多项式的根为原根。

GF(p<sup>n</sup>)的本原多项式f(x)满足：

1. x<sup>p<sup>n</sup>-1</sup>=1(mod f(x))
2. x<sup>m</sup>≠1(mod f(x))

其中多项式的系数模p，1≤m≤p<sup>n</sup>-1

如果多项式w(x)有实根则w不是本原多项式。

如果w(0)不是GF(p)的原根则w不是本原多项式。

如果w使w<sup>(p<sup>n</sup>-1)/p<sub>i</sub></sup> = 1(mod p)对(p<sup>n</sup>-1)的某个因子p<sub>i</sub>成立则w不是本原多项式。

**********************

### GF(Mp<sub>31</sub><sup>2</sup>)

Mp<sub>31</sub>=2<sup>31</sup>-1是梅森素数。GF(Mp<sub>31</sub><sup>2</sup>)原根的阶为(2<sup>31</sup>-1)<sup>2</sup>-1=2<sup>32</sup>(2<sup>30</sup>-1)

本原多项式：

x<sup>2</sup>+3x+7（系数最小）

x<sup>2</sup>-14x+53（原根7±2i）

x<sup>2</sup>-16x+73（原根8±3i）

对本原多项式w(x)=x<sup>2</sup>+a<sub>1</sub>x+a<sub>0</sub>，x<sup>2</sup>-a<sub>1</sub>x+a<sub>0</sub>也是本原多项式。

### GF((3*2^30+1)^2)

### GF((15*2^27+1)^2)