---
title: "田中秀樹教授講義メモ2025-11"
emoji: "🔥"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["math", "science"]
published: true
---

#包接水和物 #vdWP

# 2025-11-04 正準集団と大正準集団

カノニカル分配関数(あるいは状態和)は次のように書ける。

$$Z_N=\frac{1}{N!}\int\int\exp\left(-\frac{H}{kT}\right)d\mathbf r^N d\mathbf p^N\tag{\text{0.1, Classical}}$$

$$Z_N=\sum_jw_j\exp\left(-\frac{E_j}{kT}\right)\tag{0.1, \text{Quantum}}$$

これから、Helmholtz エネルギー$A$が次のように書ける。

$$A=-kT\ln Z_N        \tag{0.2}$$

一方、グランドカノニカル分配関数は次のように書ける。

$$\Xi=\sum_{n=0}^\infty \exp\left(\frac{\mu n}{kT}\right)Z_n        \tag{0.3}$$

これを化学ポテンシャル$\mu$で微分すると平均分子数$\left<n \right>$が得られる。

$$kT\frac{\partial\ln\Xi}{\partial\mu}=\frac{1}{\Xi}\sum_{n=0}^\infty n\exp\left(\frac{\mu n}{kT}\right)Z_n = \frac{\sum_{n=0}^\infty n\exp\left(\frac{\mu n}{kT}\right)Z_n}{\sum_{n=0}^\infty \exp\left(\frac{\mu n}{kT}\right)Z_n}=\left<n\right>        \tag{0.4}$$

# 1. 包接水和物

包接水和物は、水分子が形成するカゴ状の結晶構造(ホスト格子)の中にゲスト分子がとりこまれたもの。ホスト格子は丈夫な完全な結晶で、分子の個数が変動しないのでカノニカル集団として扱いたい。一方、ゲスト分子の入る割合は化学ポテンシャルによって左右されるので、グランドカノニカル(大正準)集団として扱いたい。

前者と後者が互いに独立であると仮定する。つまり、ホスト格子の状態に関わらずゲスト分子はケージに入るかどうかを自由に決定でき、またゲスト配置によってホスト格子の熱力学状態が影響を受けないと仮定すると、これらの状態は独立に選べ、分配関数はホストの分配関数と、ゲストの分配関数の積で表せる。

$$Q=\Xi\cdot Z_H        \tag{1.1}$$

$\Xi(\mu_g,V,T)$、$Z_H(N_w,V,T)$はそれぞれゲスト、ホストの分配関数で、$V$と$T$は共通である。$\mu$はゲストの化学ポテンシャル、$N_w$は水の分子数である。

## 1.1 ホスト格子

ホスト格子の Helmholtz エネルギー$A_H$は、一般的な固体と同じように近似できる。

$$A_H=U_q+F_H+F_A-TS_R        \tag{1.1.1}$$

ただし、$U_q$は絶対零度での凝集エネルギー、$F_H$と$F_A$は調和振動と非調和振動の自由エネルギー、$S_R$は残余エントロピー。これらの計算方法は後述する(かも)

## 1.2 ゲスト分子

ゲスト分子 1 つが、結晶の中にいる場合の分配関数は

$$\int\int\exp\left(-\frac{h_g}{kT}\right)d\mathbf r d\mathbf p        \tag{1.2.1}$$

運動量での積分は温度のみに依存するので先に計算して除外できる。

$$\int\exp\left(-\frac{h_g}{kT}\right)d\mathbf r        \tag{1.2.2}$$

実際には、結晶内部は多数のケージに区分されていて、それぞれのケージの中は等価なので、結晶全体での積分を、ケージ内だけの積分に限定しそれをケージ数$N_c$倍する。

$$N_c\int_\text{cage}\exp\left(-\frac{h_g}{kT}\right)d\mathbf r        \tag{1.2.3}$$

ゲストが 2 個ある場合。

$$\frac{N_c(N_c-1)}{2!}\int_1\int_2\exp\left(-\frac{h_1+h_2}{kT}\right)d\mathbf r_1d\mathbf r_2=\frac{N_c(N_c-1)}{2!}\left(\int_\text{cage}\exp\left(-\frac{h_g}{kT}\right)d\mathbf r\right)^2        \tag{1.2.4}$$

ゲストが$n$個ある場合。

$$Z_g=\left(\begin{array}{c}N_c\\n\end{array}\right)\left(\int_\text{cage}\exp\left(-\frac{h_g}{kT}\right)d\mathbf r\right)^n        \tag{1.2.5}$$

右のカッコ内の 1 分子分配関数を$z$と書き、$f=-kT\ln z=-kT\ln \int_\text{cage}\exp\left(-\frac{h_g}{kT}\right)d\mathbf r$をケージ占有自由エネルギーと呼ぶことにすると、

$$Z_g=\left(\begin{array}{c}N_c\\n\end{array}\right)\left[\exp\left(-\frac{f}{kT}\right)\right]^n        \tag{1.2.6}$$

これを使い、大分配関数$\Xi$を書き下す。(二項定理を用いた)$\mu_g$はゲスト分子の化学ポテンシャル。

> $\mu_g$に関しては次節に解説がある。

$$
\begin{array}{rcl}\Xi&=&\sum_{n=0}^{N_c}\exp\left(\frac{\mu_g n}{kT}\right)\cdot Z_g(T,V,n)\\
&=&\sum_{n=0}^{N_c}\left(\begin{array}{c}N_c\\n\end{array}\right)\exp\left(-\frac{n}{kT}(\mu_g-f)\right)\\
&=&\left[1+\exp\left(\frac{\mu_g-f}{kT}\right)\right]^{N_c}
\end{array}        \tag{1.2.7}
$$

実際の典型的な包接水和物の場合、ケージには 2 種類(12 面体と 14 面体)あり、それぞれにケージ内包自由エネルギーが異なる。それぞれのケージの個数を$N_{12}, N_{14}$とすると、大分配関数は以下のようになる。

$$\Xi=\left[1+\exp\left(\frac{\mu_g-f_{12}}{kT}\right)\right]^{N_{12}}\left[1+\exp\left(\frac{\mu_g-f_{14}}{kT}\right)\right]^{N_{14}}        \tag{1.2.8}$$

ケージの個数は結晶構造によってあらかじめ決まっており、I 型ハイドレート(メタンハイドレートなど)の場合、水分子数$N_w$との比は次のようになっている。

$$\alpha_{12}=N_{12}/N_w=2/46        \tag{1.2.9}$$
$$\alpha_{14}=N_{14}/N_w=6/46        \tag{1.2.10}$$

これを用いると、

$$\Xi=\left[1+\exp\left(\frac{\mu_g-f_{12}}{kT}\right)\right]^{\alpha_{12}N_w}\left[1+\exp\left(\frac{\mu_g-f_{14}}{kT}\right)\right]^{\alpha_{14}N_w}        \tag{1.2.11}$$

一般的にケージの種類を$j$で表すと、

$$\Xi=\prod_j\left[1+\exp\left(\frac{\mu_g-f_j}{kT}\right)\right]^{\alpha_jN_w}        \tag{1.2.12}$$

これを冒頭の式$Q=\Xi\cdot Z_H$に代入し、偏微分することでさまざまな物理量を得ることができる。

例えば、$\mu_g$で微分すると、ゲスト分子の個数$\left<N_g\right>$が得られる。$Z_H$はゲスト分子に関係がないので、

$$
\begin{array}{rcl}
\left<N_g\right>&=&kT\left(\frac{\partial\ln Q}{\partial \mu_g}\right)_{T,V,N_w}\\
&=&kT\left(\frac{\sum_j\alpha_jN_w\partial\ln\left(1+\exp\left(\frac{\mu_g-f_j}{kT}\right)\right)}{\partial\mu_g}\right)\\
&=&N_w\sum_j\alpha_j\frac{\exp\left(\frac{\mu_g-f_j}{kT}\right)}{1+\exp\left(\frac{\mu_g-f_j}{kT}\right)}
        \tag{1.2.13}
\end{array}
$$

ケージ種$j$の占有率を$x_j$とすると、全ゲスト分子の個数は$N_g=\sum_jN_jx_j=N_w\sum_j\alpha_jx_j$と書ける。これと上の式を見比べると、占有率は次のように表せる。

$$x_j=\frac{\exp\left(\frac{\mu_g-f_j}{kT}\right)}{1+\exp\left(\frac{\mu_g-f_j}{kT}\right)}        \tag{1.2.14}$$

また、$Q$を$N_w$で微分すると、水の化学ポテンシャルが求められる。

$$
\begin{array}{rcl}
\mu_w&=&\mu_c\\
&=&-kT\left(\frac{\partial\ln Q}{\partial N_w}\right)_{T,V,\mu_g}\\
&=&-kT\left(\frac{\partial\ln Z_H}{\partial N_w}\right)-kT\left(\frac{\partial\ln \Xi}{\partial N_w}\right)\\
&=&\mu_h^0-kT\sum_j\alpha_j\ln\left(1+\exp\left(\frac{\mu_g-f_j}{kT}\right)\right)\\
&=&\mu_h^0+kT\sum_j\alpha_j\ln\left(1-x_j\right)\\
        \tag{1.2.15}
\end{array}
$$

この式から、ケージ占有率$x_j$が大きいほど、水の化学ポテンシャル$\mu_w$が低くなることがわかる。この式が、van der Waals--Platteeuw 理論の要諦である。

三相平衡では、水(氷)、ハイドレート、ゲスト流体(ガス)の三者が平衡となる。この時、水とハイドレートの中の水の化学ポテンシャルは等しく、またハイドレートとゲスト流体の中のゲストの化学ポテンシャルは等しい。

(2025-11-11)

ハイドレート内外にあるゲスト分子の化学ポテンシャル$\mu_g$は互いに等しい。ハイドレート外のゲスト分子が気体である場合の化学ポテンシャルを一般的な形で書き下す。

$$
Z(N,T,V)=\frac{1}{N!h^{3N}}\int\int\exp\left(-\frac{\phi}{kT}\right)\exp\left(-\frac{\sum\mathbf{p}^2}{2mkT}\right)d\mathbf{r}^N\mathbf{p}^N
\tag{1.2.16}
$$

## 1.2.1 理想気体

理想気体なら相互作用$\phi$はない。

$$Z=\frac{1}{N!h^{3N}}\int_Vd\mathbf{r}^N\prod_N\int\exp\left(-\frac{\mathbf{p}^2}{2mkT}\right)\mathbf{p}^N     \tag{1.2.1.1}$$

体積全域での空積分はただの体積、運動量の部分はガウス積分。

$$Z=\frac{1}{N!h^{3N}}\left(2\pi mkT\right)^\frac{3N}{2}     \tag{1.2.1.2}$$

Helmholtz エネルギーは

$$A=-kT\ln Z=-kTN\left(\ln V-3\ln h+\frac{3}{2}\ln(2\pi mkT)-\ln N+1\right)     \tag{1.2.1.3}$$

ここで de Broglie 波長$\lambda=\frac{h}{\sqrt{2\pi mkT}}$を導入すると、

$$A=kTN\left(\ln \frac{N}{V}+3\ln\lambda-1\right)=NkT(\ln\rho\lambda^3-1)     \tag{1.2.1.4}$$

$$G=A+pV=NkT\ln\rho\lambda^3     \tag{1.2.1.5}$$

気体の化学ポテンシャルは、

$$\mu_g=\frac{\partial G}{\partial N}=kT\ln\rho\lambda^3     \tag{1.2.1.6}$$

## 1.2.2 実在気体

一般に、

$$\left(\frac{\partial A}{\partial V}\right)_{T,N}=-p     \tag{1.2.2.1}$$

Helmholtz エネルギーから理想気体の寄与$A_0$を引く。

$$\left(\frac{\partial (A-A_0)}{\partial V}\right)_{T,N}=-(p-p_0)     \tag{1.2.2.2}$$

両辺を体積で積分すると、

$$A-A_0=-\int_\infty^V (p-p_0)dV     \tag{1.2.2.3}$$

移項して

$$A=A_0+\int_V^\infty (p-p_0)dV     \tag{1.2.2.4}$$

つまり、理想気体からの圧力のずれを、体積無限大(ずれなし)から積分すれば、実在気体の Helmholtz エネルギーが求められる。

## 1.3 三相平衡

三相平衡(水、ハイドレート、ゲスト)の場合、ギブズの相律より、自由に選べるパラメータは 1 つだけであり、例えば温度$T$を指定すると、体積$V$、圧力$p$、組成$x$が定まる。

> 実際の計算手順は、
>
> 1. 温度$T$を指定する。
> 1. ガスの圧力$p$あるいはその温度でのゲストの化学ポテンシャル$\mu_g$(1.2.1.6)を指定する。
> 1. ハイドレートの体積$V_c$は圧力でほとんど変化しないので、あらかじめ与える。
> 1. その温度でのケージ内包自由エネルギー$f_j$を計算する。(方法は後述)
> 1. 空のホスト格子の自由エネルギーを計算し、$\mu_h^0$(1.1.1)を求める。
> 1. van der Waals-Platteuw の式(1.2.15)により、$\mu_c$を得る。
> 1. 水の化学ポテンシャル$\mu_w$を別途計算し、$\mu_c$と等しくなるように、温度あるいは圧力を調節する。

三相平衡の場合、ゲストのモル分率は結果として得られるので、大正準集団的に扱う。(水の分子数は固定されているので、こちらは正準集団的)

## 1.4 二相共存

2 相共存には次の 2 種類が考えられる。

1. 水とハイドレートの共存
2. ハイドレートとゲスト流体の共存

以下では 2 の場合を議論する。

二相共存の場合は、ギブズの相律により、2 つのパラメータが指定できる。ここでは温度$T$とゲストのモル分率$y$を選ぶ。

---

> この部分の話は、計算には関係ない気がする

ハイドレートの Helmholtz エネルギーは次の式で表現される。

$$A_c=A_c^0+kTN_w\sum_j\alpha_j\left[\beta x_jf_j+x_j\ln x_j+(1-x_j)\ln(1-x_j)\right]     \tag{1.4.1}$$

ただし、$A_c^0$は空のホスト格子の Helmholtz エネルギー、$\alpha_j=N_j/N_w$は$j$ケージの数と水分子数の比(1.2.9)、$f_j$は$j$ケージへのゲストの内包自由エネルギー、$x_j$は$j$ケージの占有率(1.2.14)。

これを変形すると、次の式に至る。

$$A_c=A_c^0+N_wkT\sum_j\alpha_j\left[\beta x_j\mu_g - \ln\left(1+\exp\left\{\beta(\mu_g-f_j)\right\}\right)\right]     \tag{1.4.2}$$

$j$ケージ内のゲスト分子数は$j$ケージの個数$N_j$と水分子数$N_w$と占有率$x_j$などから計算できる。

$$n_j=x_j\alpha_j N_w     \tag{1.4.3}$$

これを上の式に入れると、

$$
A_c=A_c^0
+\sum_jn_j\mu_g
-N_wkT\sum_j\alpha_j\ln\left(1+\exp\left\{\beta(\mu_g-f_j)\right\}\right)     \tag{1.4.4}
$$

水分子数$N_w$で微分すると、水分子の化学ポテンシャル$\mu_w$が得られる。(1.2.15)

$$\mu_w=\frac{\partial A_c}{\partial N_w}=\mu_c^0-kT\sum_j\alpha_j\ln\left(1+\exp\left\{\beta(\mu_g-f_j)\right\}\right)     \tag{1.4.5}$$

---

(1.4.3)を用いて、ハイドレートを構成するゲストのモル分率は次の式で書ける。

$$y=\frac{\sum_jn_j}{N_w+\sum_jn_j}     \tag{1.4.6}$$

この式の中で、$j$ケージの占有率$x_j$は$\mu_g$に依存している。だから、これを逆向きに計算すると、与えられたモル分率$y$になるような、ゲストの化学ポテンシャル$\mu_g$が計算できるはずだ。

計算を簡単にするために、$C=\exp\beta\mu_g, F_j=\exp(-\beta f_j)$とおく。また、一般に$a=\frac{b}{b+c}$なら$b=\frac{ac}{1-a}$である。

> イメージしやすいように、メタンハイドレートを想定しよう 。メタンハイドレートには 12 面体と 14 面体の 2 種類のケージがある。これらの個数をそれぞれ$N_{12}, N_{14}$と書く。$\alpha_{12}=N_{12}/N_w=1/23, \alpha_{14}=3/23$である。

$$n_{12}=x_{12}\alpha_{12} N_w$$

$$x_{12}=\frac{\exp\beta(\mu_g-f_{12})}{1+\exp\beta(\mu_g-f_{12})}=\frac{CF_{12}}{1+CF_{12}}$$

これらを、(1.4.6)に代入すると、

$$\sum_jn_j=\frac{yN_w}{1-y}=x_{12}\alpha_{12}N_w+x_{14}\alpha_{14}N_w$$

$$\frac{y}{1-y}=\frac{\alpha_{12}CF_{12}}{1+CF_{12}}+\frac{\alpha_{14}CF_{14}}{1+CF_{14}}$$

$\frac{y}{1-y}\equiv Y$と書くと、
$$Y(1+CF_{12})(1+CF_{14})=\alpha_{12}CF_{12}(1+CF_{14})+\alpha_{14}CF_{14}(1+CF_{12})$$
整理すると、
$$Y+YCF_{12}+YCF_{14}+YC^2F_{12}F_{14}=\alpha_{12}CF_{12}+\alpha_{14}CF_{14}+(\alpha_{12}+\alpha_{14})C^2F_{12}F_{14}$$
$$C^2F_{12}F_{14}(Y-\alpha_{12}-\alpha_{14})+C((Y-\alpha_{12})F_{12}+(Y-\alpha_{14})F_{14})+Y=0$$

$C$に関するこの二次方程式を解くことで、$y$から$\mu_g$が得られる。
