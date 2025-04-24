# 相互作用ヒストグラムの生成

vdWP 理論の式のなかで、もっとも計算コストが高いのは、次の分配関数に現れる積分である。

$$
f_0=-kT\ln\int_{\bf \Omega}\int_{\bf r}\exp\left[-\beta w({\bf r},{\bf \Omega})\right]{\rm d}{\bf r}{\rm d}{\bf \Omega}\tag{0.6}
$$

ケージ内を細かい格子に分割し、各点でさらに分子配向を数百通り選んでそれぞれの配置でのゲスト分子と周囲の水分子の相互作用を計算する必要がある。計算すべき配置の総数は数百万に及び、計算には 1 分子種あたり 30 分(THF の場合)程度かかる。温度が変わるたびにこの計算を行うのは非常にコストが高いので、あらかじめこれらの数百万点での配置でのポテンシャルエネルギーを計算し、ヒストグラムを作って保存しておく。具体的には、分配関数

$$
\begin{aligned}
&&\int_{\bf \Omega}\int_{v_k}\exp\left[-\beta w({\bf r},{\bf \Omega})\right]{\rm d}{\bf r}{\rm d}{\bf \Omega}\\
&=&\int\int\int_{x,y,z\in v_k}\int_{\phi=0}^{2\pi}\int_{\psi=0}^{2\pi}\int_{\theta=0}^{\pi}\sin\theta\exp\left[-\beta w_k\right]{\rm d}x{\rm d}y{\rm d}z{\rm d}\phi{\rm d}\psi{\rm d}\theta\\
&=&\sum_{x,y,z\in v_k}\sum\sum\sum\exp\left[-\beta w_k\right]\sin\theta{\rm \Delta}x{\rm \Delta}y{\rm \Delta}z{\rm \Delta}\phi{\rm \Delta}\psi{\rm \Delta}\theta
\end{aligned}
\tag{4.1}
$$

の総和を計算する代わりに、$w_k$の値を横軸とし、縦軸には$\sin\theta{\rm \Delta}x{\rm \Delta}y{\rm \Delta}z{\rm \Delta}\phi{\rm \Delta}\psi{\rm \Delta}\theta$を積算したヒストグラムを作成する。任意の温度での分配関数はこのヒストグラムに$\exp[-\beta w_k]$をかけながら積分すれば得られる。また、平均ポテンシャル$\bar{w}_k$は、このヒストグラムに$w_k\exp[-\beta w_k]$をかけながら積分し、分配関数で割れば得られる。
