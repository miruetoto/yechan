---
title: "연구 ▷ HST ▷ Toward an Occupation SLLN for the Heavy-Snow Walker"
author: 신록예찬
date: 05/15/2026
draft: false
output-file: 260515_dc486d.html
---

> **Abstract.** We investigate whether the empirical occupation frequencies of the heavy-snow walker converge almost surely: $\frac{1}{t}\sum_{s=1}^t\mathbf{1}\{X_s=v_i\}\to\rho_i$ a.s. A martingale decomposition (Section 2) reduces the problem to the convergence of a Cesàro average of conditional visit probabilities; this step is unconditional and rigorous. The remaining step---showing that the Cesàro average converges---is open in general. We identify the precise obstruction (Section 4), prove a conditional result under an explicit *downstream-stabilization* hypothesis (Section 5), and verify the hypothesis for specific graph families (Section 6).

## 1. Setup

Algorithm 1' (random-step variant) of [Choi--Oh 2026]: connected weighted graph $\mathcal{G}=(V,\mathbf{E},\mathbf{W})$ with $|V|=n$, i.i.d. increments $b'_t\sim\mathrm{Unif}(0,b)$, block timer $Z_t\in\{0,\ldots,T_{\max}\}$, fall distribution $\boldsymbol{\mu}_0$ with $\mu_{\min}:=\min_i\mu_0(v_i)>0$.

**Order of operations at step $t$.**

`1.` Determine the walker destination $X_t$:

- If $Z_{t-1}=T_{\max}$ (block-flag): draw $X_t\sim\boldsymbol{\mu}_0$, independently of $\mathcal{F}_{t-1}$.
- Otherwise (flow step): $X_t$ is drawn by the flow rule applied to heights $\{h(v,t-1)\}$ (see below).

`2.` Draw $b'_t\sim\mathrm{Unif}(0,b)$, independently of $(\mathcal{F}_{t-1},X_t)$.

`3.` Update $h(X_t,t):=h(X_t,t-1)+b'_t$.

**Flow rule (precise).** At a non-block step with current node $X_{t-1}=u$: the set of downstream neighbors is $\mathcal{D}_t:=\{w\in\mathcal{N}(u):h(w,t-1)\leq h(u,t-1)\}$, i.e. those neighbors whose height (after all step-$(t-1)$ updates) is at or below that of $u$. Note that $h(u,t-1)$ already includes the increment $b'_{t-1}$ if $u=X_{t-1}$ (by the update rule in step 3 of the previous iteration), so the comparison is between *post-snowfall* heights. If $\mathcal{D}_t\neq\emptyset$, the walker moves to $X_t\sim\mathrm{Cat}\bigl(\{w_{uw}/\sum_{w'\in\mathcal{D}_t}w_{uw'}:w\in\mathcal{D}_t\}\bigr)$---a *weight-proportional random choice* among downstream neighbors. If $\mathcal{D}_t=\emptyset$, the walker falls back to $X_t\sim\boldsymbol{\mu}_0$.

> **Remark** (Flow is stochastic, not deterministic). Even with a permanently fixed height ranking and all pairwise gaps exceeding $(T_{\max}+1)b$, the flow step from node $u$ picks among the downstream neighbors of $u$ by weight-proportional random choice, not by "pick the single lowest neighbor." On graphs with branching (degree $>1$ at some nodes), the per-round visit sequence is *not* a deterministic function of the starting node and the ranking.

**Filtration.** $\mathcal{F}_t:=\sigma\bigl(S(0),X_1,U_1,b'_1,\ldots,X_t,U_t,b'_t\bigr)$, where $S(t):=(\boldsymbol{\delta}_t,X_t,Z_t)$ and $U_s$ denotes the randomness used to determine $X_s$ at step $s$: at a block-flag step, $U_s$ is the $\boldsymbol{\mu}_0$-coin; at a flow step, $U_s$ is the weight-proportional coin among downstream neighbors. Note:

- $X_t$ is $\sigma(\mathcal{F}_{t-1},U_t)$-measurable;
- $b'_t$ is independent of $\sigma(\mathcal{F}_{t-1},U_t)$.

**Conditional visit probability.** Define

$$p_i(s)\;:=\;\mathbb{P}(X_s=v_i\mid\mathcal{F}_{s-1}).$$

This is $\mathcal{F}_{s-1}$-measurable, with $\sum_i p_i(s)=1$ and $p_i(s)\geq 0$. At a block-flag time $t_r$: $p_i(t_r+1)=\mu_0(v_i)$. At a flow step: $p_i(s)$ is the weight-proportional probability determined by heights at time $s-1$ and the previous position $X_{s-1}$.

## 2. Martingale decomposition

**Lemma** (Martingale noise removal). For each $v_i\in V$, define

$$M_t^{(i)}\;:=\;\sum_{s=1}^{t}\bigl[\mathbf{1}\{X_s=v_i\}-p_i(s)\bigr].$$

::: {.callout-tip collapse="true" title="보충설명"}
$M_t^{(i)}$는 "실제 방문 횟수"와 "기대 방문 횟수"의 누적 차이. 매 step에서 $v_i$에 실제로 갔으면 $+1$, 안 갔으면 $0$을 기록하고, 거기서 "갈 확률" $p_i(s)$를 빼준 것을 $t$번 더한 것이다. 일종의 "예측 오차의 합산".
:::

Then:

**(a)** $\{M_t^{(i)}\}_{t\geq 0}$ is a martingale w.r.t. $\{\mathcal{F}_t\}$ with bounded increments $|\Delta M_t^{(i)}|\leq 1$.

::: {.callout-tip collapse="true" title="보충설명: (a)가 말하는 것"}
$M_t^{(i)}$가 마팅게일이라는 건: 과거 정보를 다 알아도, 다음 예측 오차의 기댓값은 0이라는 뜻. 즉 예측 오차가 체계적으로 한 쪽으로 쏠리지 않는다. 그리고 한 step에서 오차는 최대 1 (갔거나 안 갔거나).
:::

**(b)** $M_t^{(i)}/t\to 0$ a.s.

::: {.callout-tip collapse="true" title="보충설명: (b)가 말하는 것"}
예측 오차의 합 $M_t^{(i)}$가 $t$보다 훨씬 느리게 자란다. $t$로 나누면 0으로 간다. 직관: 동전 던지기의 누적 편차가 $\sqrt{t}$ 정도로 자라는 것과 비슷. $t$로 나누면 사라진다.
:::

**(c)** Therefore

$$\frac{1}{t}\sum_{s=1}^{t}\mathbf{1}\{X_s=v_i\}\;=\;\frac{1}{t}\sum_{s=1}^{t}p_i(s)\;+\;o(1)\qquad\text{a.s.}$$

::: {.callout-important collapse="true" title="핵심 결론"}
**실제 방문 비율** $=$ **조건부 확률의 평균** $+$ **사라지는 noise**.

즉 $\rho_i$가 존재하는지는 오직 $\frac{1}{t}\sum p_i(s)$가 수렴하는지에 달렸다. 실제 방문의 랜덤성(noise)은 $1/t$ 스케일에서 자동으로 사라진다.
:::

*Proof.*

**(a)** The increment is $\Delta M_s^{(i)}=\mathbf{1}\{X_s=v_i\}-p_i(s)$. Since $p_i(s)$ is $\mathcal{F}_{s-1}$-measurable, and $\mathbb{E}[\mathbf{1}\{X_s=v_i\}\mid\mathcal{F}_{s-1}]=p_i(s)$ by definition, we have $\mathbb{E}[\Delta M_s^{(i)}\mid\mathcal{F}_{s-1}]=0$.

Since $X_s$ and $U_s$ are included in $\mathcal{F}_s$, the increment $\Delta M_s^{(i)}$ is $\mathcal{F}_s$-measurable. The bound $|\Delta M_s^{(i)}|\leq 1$ is immediate.

**(b)** By the strong law for square-integrable martingales (Hall--Heyde, Thm 2.18): if $\sum_{s=1}^\infty\mathbb{E}[(\Delta M_s^{(i)})^2\mid\mathcal{F}_{s-1}]/s^2<\infty$ a.s., then $M_t^{(i)}/t\to 0$ a.s. Since $\mathbb{E}[(\Delta M_s^{(i)})^2\mid\mathcal{F}_{s-1}]=p_i(s)(1-p_i(s))\leq 1$, the sum is $\leq\sum_{s\geq 1}1/s^2<\infty$.

**(c)** Immediate from (b). $\square$

> **Remark** (Scope). This lemma is *unconditional*: it holds for every realization, on every connected graph, without any assumption on $\rho_i$, ergodicity, or regime. The problem of proving $\rho_i$ exists is now equivalent to showing that the Cesàro average $\frac{1}{t}\sum_{s=1}^t p_i(s)$ converges.

## 3. Height-difference process

**Lemma** (Height-difference decomposition). For any pair $(i,j)$, define $D_{ij}(t):=h(v_i,t)-h(v_j,t)$. Then

$$D_{ij}(t)\;=\;D_{ij}(0)+\bar{b}\sum_{s=1}^t\bigl[\mathbf{1}\{X_s=v_i\}-\mathbf{1}\{X_s=v_j\}\bigr]+\widetilde{M}_{ij}(t),$$

where $\bar{b}:=b/2$ and $\widetilde{M}_{ij}(t):=\sum_{s=1}^t(b'_s-\bar{b})[\mathbf{1}\{X_s=v_i\}-\mathbf{1}\{X_s=v_j\}]$ satisfies $|\widetilde{M}_{ij}(t)|=O(\sqrt{t\log\log t})$ a.s.

*Proof.* Write $b'_s=\bar{b}+(b'_s-\bar{b})$. Define the intermediate filtration $\mathcal{H}_s:=\sigma(\mathcal{F}_{s-1},U_s)\supseteq\sigma(\mathcal{F}_{s-1},X_s)$. Since $b'_s$ is independent of $\mathcal{H}_s$ and has mean $\bar{b}$, the increment $(b'_s-\bar{b})[\mathbf{1}\{X_s=v_i\}-\mathbf{1}\{X_s=v_j\}]$ has conditional mean zero given $\mathcal{H}_s$; since $\mathcal{H}_s\supseteq\mathcal{F}_{s-1}$, the tower property gives $\mathbb{E}[\,\cdot\mid\mathcal{F}_{s-1}]=0$ as well. Since $b'_s,X_s,U_s\in\mathcal{F}_s$, the increment is $\mathcal{F}_s$-measurable. Hence $\widetilde{M}_{ij}$ is a martingale w.r.t. $\{\mathcal{F}_t\}$.

The conditional variance of each increment is $\leq\mathrm{Var}(b')\cdot 1=b^2/12$, so $\langle\widetilde{M}_{ij}\rangle_t\leq(b^2/12)\,t$. The LIL for martingales gives $|\widetilde{M}_{ij}(t)|=O(\sqrt{t\log\log t})$ a.s. $\square$

**Lemma** (Height--occupation link). Define $\hat\rho_i(t):=\frac{1}{t}\sum_{s=1}^t\mathbf{1}\{X_s=v_i\}$. Then

$$D_{ij}(t)\;=\;\bar{b}\,t\,[\hat\rho_i(t)-\hat\rho_j(t)]+O(\sqrt{t\log\log t})\qquad\text{a.s.}$$

*Proof.* Immediate from the decomposition and $D_{ij}(0)=O(1)$. $\square$

## 4. The obstruction

The martingale decomposition reduces the problem to: does $\frac{1}{t}\sum p_i(s)$ converge? A sufficient condition (by Cesàro's lemma) is $p_i(s)\to p_i^\ast$, but this is *not necessary*---the Cesàro average can converge even if $p_i(s)$ oscillates.

The core difficulty is the *feedback loop*:

$$\text{heights}\;\xrightarrow{\text{flow rule}}\;\text{walker visits}\;\xrightarrow{\text{accumulation}}\;\text{heights}.$$

The conditional visit probability $p_i(s)$ depends on heights at time $s-1$, which depend on all previous visits.

**Why subsequential limits are insufficient.** The running occupation $\hat\rho(t)=(\hat\rho_1(t),\ldots,\hat\rho_n(t))$ lies in the compact simplex $\Delta^{n-1}$, so subsequential limits exist. One might hope to argue:

`(i)` if a subsequential limit has $\hat\rho_i\neq\hat\rho_j$, the height difference $D_{ij}$ diverges, forcing ranking stabilization.

This fails because:

- **Oscillation.** $\hat\rho_i(t)-\hat\rho_j(t)$ can have subsequences converging to $+\gamma$ and others converging to $-\gamma$. Then $D_{ij}(t)$ alternates between $+\Theta(t)$ and $-\Theta(t)$, and the ranking never freezes.
- **No a priori bound on inter-subsequence gaps.** Between a subsequence where $D_{ij}(t_k)\approx+\bar{b}\gamma t_k$ and the next where $D_{ij}(t'_k)\approx-\bar{b}\gamma t'_k$, only $O(t_k)$ steps are needed (bounded increments $\leq b$). The linear growth rate is exactly matched by the linear time scale.

**What would close the argument.** The occupation SLLN follows if we can show that the set of subsequential limits of $\hat\rho(t)$ is a *singleton*. Possible strategies:

**(A)** *Self-reinforcement*: once $D_{ij}$ is large and positive, the dynamics reinforce $\hat\rho_i>\hat\rho_j$. Graph-structure-dependent.

**(B)** *Harris recurrence* of the augmented chain. Works in the balanced regime; in the drift regime, $M(t)\to\infty$ prevents standard recurrence.

**(C)** *Stochastic approximation* (ODE method): view $\hat\rho(t)$ as tracking a mean-field ODE and show the ODE has a unique attractor.

## 5. Conditional result

We formulate a hypothesis strong enough to close the SLLN. A mere ranking stabilization (strict height ordering eventually fixed) is *not* sufficient: even with a fixed ranking, the downstream set $\mathcal{D}_s$ at a flow step from $u$ depends on the actual height *gap* $h(u,s-1)-h(w,s-1)$, not just the sign. If adjacent gaps remain bounded (say $O(b)$), within-round snowfall accumulation can change the downstream membership from step to step. We therefore require that the gaps diverge.

**Definition** (Downstream stabilization, **DS**). Say that a realization of Algorithm 1' satisfies *downstream stabilization* (DS) if there exist a finite time $T_0$ and, for every node $u\in V$, a fixed set $\mathcal{D}_\infty(u)\subseteq\mathcal{N}(u)$ (possibly empty) such that for all $t\geq T_0$ and every flow step from $u$ at time $t$:

$$\mathcal{D}_t\;=\;\mathcal{D}_\infty(u).$$

A sufficient condition (used in all our applications) is: for every pair $(u,w)$ with $w\in\mathcal{N}(u)$, eventually either $h(u,t)-h(w,t)>(T_{\max}+1)b$ permanently (so $w$ is always downstream of $u$), or $h(w,t)-h(u,t)>(T_{\max}+1)b$ permanently (so $w$ is never downstream of $u$). This gap condition implies DS because the within-round height fluctuations are at most $(T_{\max}+1)b$, so the downstream membership of $w$ w.r.t. $u$ cannot change once the gap exceeds this threshold.

> **Remark** (DS is strictly stronger than ranking stabilization). A fixed strict height ranking does *not* imply DS. With a fixed ranking $v_1>v_2>\cdots>v_n$, adjacent nodes $v_k,v_{k+1}$ might have $0<h(v_k,t)-h(v_{k+1},t)<(T_{\max}+1)b$ for all $t$---the ranking never flips, but the gap stays bounded, allowing the downstream set to include or exclude $v_{k+1}$ depending on the within-round fluctuations.
>
> DS requires all *adjacent* height gaps to diverge past the $(T_{\max}+1)b$ threshold, which happens in particular when every pair of neighbors has distinct asymptotic occupation rates.

**Theorem** (Occupation SLLN under DS). Under Algorithm 1' with connected $\mathcal{G}$ and $\mu_{\min}>0$, on the event $\{\text{DS holds}\}$, the limit

$$\rho_i\;:=\;\lim_{t\to\infty}\frac{1}{t}\sum_{s=1}^t\mathbf{1}\{X_s=v_i\}\quad\text{exists a.s.}$$

Moreover, $\rho_i$ is determined by the stabilized downstream structure $\{\mathcal{D}_\infty(u)\}_{u\in V}$, $\mathcal{G}$, $\boldsymbol{\mu}_0$, and $T_{\max}$.

*Proof.*

**Step 1: stabilized flow transition.** After time $T_0$, at every flow step from node $u$, the downstream set is $\mathcal{D}_\infty(u)$. The flow transition is therefore the fixed kernel

$$P_\infty(u,w)\;=\;\frac{w_{uw}}{\sum_{w'\in\mathcal{D}_\infty(u)}w_{uw'}}\qquad\text{for }w\in\mathcal{D}_\infty(u),$$

if $\mathcal{D}_\infty(u)\neq\emptyset$, or a $\boldsymbol{\mu}_0$-fallback if $\mathcal{D}_\infty(u)=\emptyset$. This is a fixed stochastic kernel, independent of $t$.

**Step 2: round structure and regeneration.** A *round* begins at a block-flag time $t_r$ (i.e. $Z_{t_r}=T_{\max}$) and ends at the next block-flag time $t_{r+1}$. Within a round, the block timer increments at each flow step ($Z_t\leftarrow Z_{t-1}+1$) and the round terminates when $Z=T_{\max}$ is reached again ($m_r:=t_{r+1}-t_r\leq T_{\max}+1$). The first step of each round is a fresh $\boldsymbol{\mu}_0$-draw: $X_{t_r+1}\sim\boldsymbol{\mu}_0$, independently of $\mathcal{F}_{t_r}$. This is a *regeneration*---the walker's position at the start of round $r+1$ carries no memory of the previous round's internal trajectory, so the round structure restarts identically.

**Step 3: i.i.d. rounds.** For rounds $r\geq R_1$ (starting after $T_0$):

- Step 1 of each round: $\tilde X_1^{(r)}\sim\boldsymbol{\mu}_0$, i.i.d. across rounds (by regeneration in Step 2).
- Steps 2 through $m_r$: the walker follows the Markov chain with transition $P_\infty$ (stochastic---the walker randomly picks among the fixed downstream set---but with a time-invariant kernel).

Across rounds, the tuples $(\boldsymbol{N}^{(r)},m_r)_{r\geq R_1}$ are i.i.d., where $\boldsymbol{N}^{(r)}=(N_1^{(r)},\ldots,N_n^{(r)})$ are visit counts and $m_r$ is the round length. Independence across rounds holds because: (i) each round starts from a fresh $\boldsymbol{\mu}_0$-draw (regeneration); (ii) the flow kernel $P_\infty$ is fixed (by DS); and (iii) the $b'_t$ increments (which affect heights but *not* the downstream sets after stabilization) are i.i.d. and independent of walker positions.

**Step 4: renewal-reward SLLN.** By the ratio SLLN for i.i.d. numerator/denominator (renewal-reward theorem):

$$\frac{\sum_{r=R_1}^{R}N_i^{(r)}}{\sum_{r=R_1}^{R}m_r}\;\xrightarrow{R\to\infty}\;\frac{\mathbb{E}[N_i^{(r)}]}{\mathbb{E}[m_r]}\;=:\;\rho_i^\infty\qquad\text{a.s.}$$

The first $R_1$ rounds contribute $O(1)$ to both numerator and denominator. Converting from round-index $R$ to time-index $t$: since $m_r\leq T_{\max}+1$, the number of rounds $R(t)$ up to time $t$ satisfies $R(t)\to\infty$, and the ratio $\hat\rho_i(t)=\bigl(\sum_{r=0}^{R(t)}N_i^{(r)}\bigr)/t$ inherits the same limit $\rho_i^\infty$ a.s. (the pre-stabilization rounds contribute $O(1)/t\to 0$). Hence $\rho_i=\rho_i^\infty$ exists. $\square$

> **Remark** ($\mathbf{y}$-dependence). The stabilized kernel $P_\infty$ depends on $\{\mathcal{D}_\infty(u)\}$, which in turn depends on the asymptotic ranking. The asymptotic ranking is determined by the $\rho_i$ values themselves (height--occupation link), creating a fixed-point structure: $\rho_i=\rho_i^\infty(\sigma^\ast)$ where $\sigma^\ast=\sigma^\ast(\rho)$. The initial signal $\mathbf{y}$ can in principle select among different fixed points (different rankings), so $\rho_i$ may depend on $\mathbf{y}$ through $\sigma^\ast$.

## 6. Verified examples

### Complete graph $K_n$ with uniform $\boldsymbol{\mu}_0$

**Proposition.** On $K_n$ with $\boldsymbol{\mu}_0=(1/n,\ldots,1/n)$, $\rho_i=1/n$ for all $i$.

*Proof.* By symmetry of $K_n$ and uniform $\boldsymbol{\mu}_0$, every permutation of node labels preserves the law of the process. Hence $\rho_i$ (if it exists) must equal $1/n$ for all $i$. Existence of $\rho_i$ follows from the positive Harris recurrence established in the companion file (balanced-regime machinery).

Note: DS does *not* hold in this case (the balanced regime has no diverging height gaps), so the Occupation SLLN Theorem is not the mechanism. The SLLN comes instead from the ergodic theorem for Harris chains. $\square$

### Star graph $S_n$ with degree-proportional $\boldsymbol{\mu}_0$

**Proposition.** On the star $S_n$ (hub $v_0$, leaves $v_1,\ldots,v_{n-1}$) with $\boldsymbol{\mu}_0\propto\deg$: $\mu_0(v_0)=1/2$, $\mu_0(v_k)=1/(2n-2)$ for $k\geq 1$. Then DS holds a.s. for the hub-vs-leaf pairs (the hub eventually stays permanently above all leaves), and $\rho_0>1/n$ for $n\geq 3$.

*Proof sketch.* On the star, each leaf $v_k$ has a single neighbor: the hub $v_0$. Consider the regime where the hub is already the highest node ($h(v_0,t)>h(v_k,t)+(T_{\max}+1)b$ for all leaves $k$).

*Flow from the hub $v_0$.* All leaves are below $v_0$, so the downstream set is $\mathcal{D}=\{v_1,\ldots,v_{n-1}\}$. The flow sends the walker to a uniformly random leaf (unweighted star).

*Flow from a leaf $v_k$.* The only neighbor is the hub $v_0$, which has $h(v_0)>h(v_k)+(T_{\max}+1)b$. After adding $b'<b$ to $v_k$, the hub is still higher: $h(v_0)>h(v_k)+b'$. So the hub is *not* downstream of $v_k$ (the hub is above, and flow goes *downhill*). Hence $\mathcal{D}=\emptyset$, and the walker falls back to $\boldsymbol{\mu}_0$.

*Per-round structure.* Each round starts with a $\boldsymbol{\mu}_0$-draw. With probability $1/2$, the draw lands on the hub; the walker then flows to a leaf, and from the leaf falls back to $\boldsymbol{\mu}_0$ (block-flag or empty downstream). With probability $1/(2n-2)$, the draw lands on leaf $v_k$; the walker immediately falls back to $\boldsymbol{\mu}_0$.

In all cases, the hub receives visits at rate $\geq\mu_0(v_0)=1/2$ (the block-flag draws alone give the hub $1/2$ of all block-flag visits). Since each leaf gets only $1/(2n-2)$ of block-flag draws plus occasional flow visits, we get $\rho_0>1/n$ for $n\geq 3$.

The height gap $D_{0k}(t)\sim\bar{b}(\rho_0-\rho_k)t\to+\infty$ confirms that the hub-above-leaves regime is self-reinforcing, verifying DS for the hub-leaf pairs.

Leaf-leaf pairs have $\rho_k=\rho_\ell=(1-\rho_0)/(n-1)$ by symmetry among leaves, so their height gaps do not diverge. However, the flow rule between leaves is trivial (leaves are not neighbors on the star), so DS for leaf-leaf pairs is vacuously satisfied (no edge, no downstream set to stabilize). $\square$

> **Remark** (Partial DS). On graphs where some adjacent pairs have equal asymptotic rates, the full DS may not hold for those pairs. In the star example, this is not an issue because equal-rate pairs (leaf-leaf) are not adjacent.
>
> On general graphs, equal-rate adjacent pairs can exist, and their downstream sets may continue to fluctuate. Whether the resulting oscillation in $p_i(s)$ still allows the Cesàro average to converge is an open question; it depends on the graph structure and cannot be resolved by exchangeability arguments alone (equal occupation rates do not imply exchangeable roles in the flow dynamics).

## 7. Application to Theorem B

**Corollary** (Drift-regime convergence under DS). Under Algorithm 1' with connected $\mathcal{G}$, $\mu_{\min}>0$, and on the event $\{\text{DS holds}\}$:

$$\frac{SD^2_{ij}(t)}{t^3}\;\xrightarrow{t\to\infty}\;\frac{\bar{b}^{\,2}(\rho_i-\rho_j)^2}{3}\qquad\text{a.s.},$$

where $\rho_i$ is given by the Occupation SLLN Theorem.

*Proof.* The theorem gives $\rho_i$ exists a.s. on $\{\text{DS}\}$. The proof of Theorem B (height decomposition + Toeplitz summation) applies verbatim. $\square$

> **Remark** (Status of Theorem B). Theorem B of the companion file assumed $\rho_i$ exists. The corollary replaces that assumption with the structural hypothesis DS, which is checkable on specific graphs (it is a property of the sample path, not a conclusion about occupation statistics).
>
> The logical relationship between DS and "$\rho_i$ exists" is: $\text{DS}\Rightarrow\rho_i\text{ exists}$, but the converse is not true in general ($\rho_i$ can exist in the balanced regime without any gap divergence).

**Conjecture.** For any connected $\mathcal{G}$ and any $\boldsymbol{\mu}_0$ with $\mu_{\min}>0$, $\rho_i$ exists a.s. for every $i\in V$.

## 8. Summary of what is and is not proved

| Statement | Status |
|-----------|--------|
| Martingale noise removal | **proved** |
| $\rho_i$ exists $\Leftrightarrow$ $\frac{1}{t}\sum p_i(s)$ converges (each $i$) | **proved** |
| Height--occupation link | **proved** |
| $\rho_i$ exists under DS | **proved** |
| $\rho_i=1/n$ on $K_n$ + uniform $\mu_0$ | **proved** (via Harris) |
| DS for star + deg-prop $\mu_0$ | **sketch** |
| $\rho_i$ exists for general graphs | **open** |
| Thm B without "$\rho_i$ exists" | **proved under DS** |

## References

- Choi, G. and Oh, H.-S. (2026). Heavy-Snow Transform for Analysis of Data on Graphs. *Manuscript*.
- Hall, P. and Heyde, C. C. (1980). *Martingale Limit Theory and Its Application*. Academic Press.
- Meyn, S. and Tweedie, R. L. (2009). *Markov Chains and Stochastic Stability* (2nd ed.). Cambridge University Press.
