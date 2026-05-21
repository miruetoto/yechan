---
title: 연구 ▷ HST ▷ 마팅게일 SLLN — Durrett 4.5.3 책공부 + HST Theorem C 의 출발점
author: 리이
date: 05/22/2026
draft: false
output-file: 260522_f2fd96.html
fontsize: 0.85em
---
<!-- 소유권자: 리이 | 사용자: 리이 -->

```{=html}
<style>
.math.display { text-align: left; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
</style>
```

> 도서관 `BOOK_Durrett_Probability_TheoryAndExamples_5e.pdf` (490p) §4.5 Square Integrable Martingales 의 Theorem 4.5.3 (martingale SLLN 의 표준 form) 을 한 줄씩 따라가고, 그 결과가 에이의 ABC증명 (`paper/260514_guebin/에이/abc오희석교수님공유용_v2.tex`) §5 Theorem C 의 Step 1 (Height 분해) 에서 어떻게 쓰이는지 짚는다.

# §0 한 줄 직관

::: {.callout-important collapse="false" title="한 줄 직관"}
Martingale 의 누적량이 — 증분이 유계이고 누적 분산이 발산할 때 — 시간에 비해 무시할 만큼 작다. 즉 *"공정한 게임의 누적 손익은 시간의 ratio 척도에서 평균 0"*.
:::

# §1 큰 그림

::: {.callout-important collapse="true" title="해석: 왜 마팅게일 SLLN 이 필요한가"}

i.i.d. case 의 Kolmogorov SLLN $\overline X_n \to \mu$ 는 합의 항이 독립이라는 가정을 본질로 쓴다. 그런데 HST 의 height 갱신처럼 한 step 의 증분이 walker 위치 $X_s$ (과거 history $\mathcal F_{s-1}$ 에 의존) 와 increment $b'_s$ (과거와 독립) 의 곱으로 표현되는 경우 — 곱 자체는 independent 아님 — 일반 SLLN 이 직접 안 먹는다.

해법: 증분의 조건부 평균이 0 인 *마팅게일 difference* 로 재포장하면 마팅게일 SLLN 한 줄로 ratio $\to 0$ 가 나온다. Durrett 4.5.3 이 정확히 그 일을 한다.

본 글의 핵심 한 줄: **bounded increment + predictable quadratic variation 발산 → $M_n / n \to 0$ a.s.**

:::

# §2 Statement — Durrett Theorem 4.5.3

$\{X_n\}_{n \geq 0}$ 이 filtration $\{\mathcal F_n\}$ 에 대해 square integrable martingale (즉 $E[X_n^2] < \infty$, $E[X_n \mid \mathcal F_{n-1}] = X_{n-1}$, $X_0 \in \mathcal F_0$) 이라 하자. 그 **predictable quadratic variation** 을

$$A_n := \sum_{m=1}^{n} E\!\bigl[(X_m - X_{m-1})^2 \,\big|\, \mathcal F_{m-1}\bigr]$$

으로 정의하면 $A_n \in \mathcal F_{n-1}$ 이며 $n$ 에 대해 증가한다.

::: {.callout-note collapse="false" title="Theorem 4.5.3 (Durrett 5e)"}

Let $f : [0, \infty) \to [1, \infty)$ be an increasing function with

$$\int_0^\infty f(t)^{-2}\, dt < \infty.$$

Then

$$\frac{X_n}{f(A_n)} \xrightarrow{n \to \infty} 0 \quad \text{a.s. on } \{A_\infty = \infty\}.$$

:::

`-` 가정 $A_\infty = \infty$ 는 "누적 분산이 발산한다" 는 가설. HST 에서는 거의 자동 성립.
`-` $f$ 의 전형: $f(t) = (t \log^{1+\epsilon} t)^{1/2} \vee 1$ 또는 $f(t) = t^{1/2 + \epsilon}$. 어느 쪽이든 $f$ 가 $\sqrt{t}$ 보다 살짝 큰 정도.
`-` 결론은 ratio $X_n / f(A_n) \to 0$. $f(A_n) \sim \sqrt{A_n \cdot \log^{1+\epsilon} A_n}$ 정도이므로 $X_n$ 이 $\sqrt{A_n}$ 보다 ε 만큼 천천히 자란다는 강한 진술.

# §3 증명 — 한 줄씩

증명 핵심은 *martingale transform* (Thm 4.2.8) → *square integrable martingale 의 finite limit 보조정리* (Thm 4.5.2) → *Kronecker's lemma* (Thm 2.5.9) 세 가지의 합성. 한 줄씩 풀어쓰면 다음과 같다.

::: {.callout-note collapse="true" title="Proof of Theorem 4.5.3"}

**Step 1 — Define a bounded predictable process.**

Set $H_m := f(A_m)^{-1}$. Since $A_m \in \mathcal F_{m-1}$, the process $\{H_m\}$ is predictable. Since $f \geq 1$, $H_m \leq 1$ a.s.

**Step 2 — Form the martingale transform.**

By Theorem 4.2.8 (bounded predictable transform preserves martingale property),

$$Y_n := (H \cdot X)_n := \sum_{m=1}^{n} \frac{X_m - X_{m-1}}{f(A_m)}$$

is again a martingale with respect to $\{\mathcal F_n\}$.

**Step 3 — Compute the predictable quadratic variation of $Y_n$.**

$$
\begin{aligned}
B_{n+1} - B_n
&:= E\!\bigl[(Y_{n+1} - Y_n)^2 \,\big|\, \mathcal F_n\bigr] \\
&= E\!\left[\frac{(X_{n+1} - X_n)^2}{f(A_{n+1})^2} \,\bigg|\, \mathcal F_n\right]
   && (\because f(A_{n+1}) \in \mathcal F_n) \\
&= \frac{1}{f(A_{n+1})^2} \cdot E\!\bigl[(X_{n+1} - X_n)^2 \,\big|\, \mathcal F_n\bigr]
   && (\because \text{factor out predictable term}) \\
&= \frac{A_{n+1} - A_n}{f(A_{n+1})^2}
   && (\because \text{definition of }A).
\end{aligned}
$$

**Step 4 — Bound $B_\infty$ using the integral assumption.**

Since $f$ is increasing, for $t \in [A_n, A_{n+1})$ we have $f(t) \leq f(A_{n+1})$ and hence $f(t)^{-2} \geq f(A_{n+1})^{-2}$, so

$$
\begin{aligned}
B_\infty
&= \sum_{n=0}^{\infty} \frac{A_{n+1} - A_n}{f(A_{n+1})^2} \\
&\leq \sum_{n=0}^{\infty} \int_{[A_n,\, A_{n+1})} f(t)^{-2}\, dt
   && (\because f \text{ increasing} \Rightarrow \text{lower-Riemann bound}) \\
&\leq \int_0^\infty f(t)^{-2}\, dt < \infty
   && (\because \text{hypothesis}).
\end{aligned}
$$

So $B_\infty < \infty$ a.s.

**Step 5 — Apply Theorem 4.5.2.**

Theorem 4.5.2 says that a square integrable martingale converges a.s. to a finite limit on $\{B_\infty < \infty\}$. By Step 4 this event is all of $\Omega$ (modulo null sets), so

$$Y_n \xrightarrow{n \to \infty} Y_\infty \quad \text{(finite)} \quad \text{a.s.}$$

**Step 6 — Conclude by Kronecker's lemma.**

Kronecker's lemma (Theorem 2.5.9): if $a_n \uparrow \infty$ and $\sum_n x_n / a_n$ converges, then $a_n^{-1} \sum_{m=1}^n x_m \to 0$.

Apply with $a_n := f(A_n)$, $x_m := X_m - X_{m-1}$. The hypothesis is

$$\sum_{m=1}^{n} \frac{x_m}{a_m} = \sum_{m=1}^{n} \frac{X_m - X_{m-1}}{f(A_m)} = Y_n \xrightarrow{n \to \infty} Y_\infty \;(\text{finite, by Step 5}),$$

and the assumption $A_\infty = \infty$ together with $f$ increasing gives $a_n = f(A_n) \uparrow \infty$. Hence Kronecker applies:

$$\frac{X_n}{f(A_n)} \;=\; \frac{1}{f(A_n)} \sum_{m=1}^{n} (X_m - X_{m-1}) \xrightarrow{n \to \infty} 0 \quad \text{a.s.} \qquad \square$$

:::

::: {.callout-tip collapse="true" title="보충: 세 도구의 역할 분담"}

- **Martingale transform (Step 1–2)**: 원래 $X_n$ 의 증분을 $f(A_m)$ 으로 가중해서 변동이 작아진 $Y_n$ 으로 변환. 증분이 작아지면 분산도 작아진다.
- **Theorem 4.5.2 (Step 3–5)**: predictable quadratic variation 이 유한하면 square integrable martingale 은 거의 확실히 finite limit 로 수렴한다. 즉 $Y_n$ 의 *진동*이 사라진다.
- **Kronecker (Step 6)**: 가중합 $\sum x_m/a_m$ 의 수렴을 비가중 합 $a_n^{-1} \sum x_m$ 의 0 수렴으로 옮긴다. 가중치 $a_n = f(A_n)$ 가 분모로 가는 게 핵심.

세 도구가 따로 놀면 의미가 없고, 같이 묶일 때만 ratio $X_n / f(A_n) \to 0$ 라는 강한 결론이 나온다.

:::

# §4 Example: 고전적 Kolmogorov SLLN 의 일반화

::: {.callout-note collapse="false" title="Example 4.5.4 (Durrett)"}

Let $\xi_1, \xi_2, \dots$ be independent with $E\xi_i = 0$ and $E\xi_i^2 = \sigma_i^2 < \infty$. Set $X_n := \xi_1 + \cdots + \xi_n$. Then $\{X_n\}$ is a square integrable martingale with respect to its natural filtration, and

$$A_n = \sum_{i=1}^{n} E[\xi_i^2 \mid \mathcal F_{i-1}] = \sum_{i=1}^{n} \sigma_i^2 \quad (\because \xi_i \perp \mathcal F_{i-1}).$$

Choose $f(t) = (t \log^{1+\epsilon} t)^{1/2} \vee 1$. Then $\int_0^\infty f(t)^{-2}\, dt < \infty$ (the integrand behaves like $1/(t \log^{1+\epsilon} t)$ near $\infty$, which is integrable). If $\sum_i \sigma_i^2 = \infty$, Theorem 4.5.3 gives

$$\frac{X_n}{\sqrt{A_n\, \log^{1+\epsilon} A_n}} \xrightarrow{n \to \infty} 0 \quad \text{a.s.}$$

In the i.i.d.\ case $\sigma_i^2 = \sigma^2$, we get $A_n = n\sigma^2$, hence

$$\frac{X_n}{\sqrt{n\, \log^{1+\epsilon} n}} \to 0 \;\Rightarrow\; \frac{X_n}{n} \to 0 \quad \text{a.s.,}$$

recovering Kolmogorov's classical SLLN with a logarithmic rate refinement.

:::

::: {.callout-important collapse="true" title="해석: independence 가 사라져도 마팅게일 구조만 살아있으면 된다"}

위 example 은 $\xi_i$ 의 *독립성* 을 본질로 쓰지 않았다. 우리가 쓴 것은:

- $E[\xi_i \mid \mathcal F_{i-1}] = 0$ (마팅게일 difference)
- $E[\xi_i^2 \mid \mathcal F_{i-1}]$ 가 잘 정의됨 (square integrable)

이 두 가지뿐이다. 따라서 *correlated 가 있어도 conditional 평균이 0 이면 SLLN 이 살아난다*. 이게 마팅게일 SLLN 이 i.i.d. SLLN 보다 광범위한 이유.

HST 에서 increment $b'_s$ 는 walker $X_s$ 와 *독립* 하지만, $b'_s \cdot \mathbf{1}\{X_s = v_i\}$ 자체는 walker 가 $v_i$ 에 자주 들리느냐에 따라 *highly correlated* 가 된다. 그래도 *conditional 평균이 0* 이라는 마팅게일 구조는 깨지지 않는다 — §5 에서 그 적용을 본다.

:::

# §5 HST 본업 적용 — 에이 ABC증명 Theorem C Step 1

에이의 `abc오희석교수님공유용_v2.tex` §5 Theorem C (drift regime) 의 증명 Step 1 (Height 분해) 에서 마팅게일 SLLN 이 정확히 한 번 호출된다. 그 자리를 들여다본다.

## 5.1 HST height 의 마팅게일 분해

HST 의 height 정의를 그대로 가져오자.

$$h(v_i, t) := \sum_{s=1}^{t} b'_s \cdot \mathbf{1}\{X_s = v_i\}, \qquad h(v_i, 0) = 0.$$

여기서 $b'_s \sim \mathrm{Unif}(0, b)$ 는 i.i.d. 이고 walker $\{X_s\}$ 와 독립. $\bar b := E[b'_s] = b/2$ 로 둔다. 증분에서 평균을 빼내면

$$
\begin{aligned}
h(v_i, t)
&= \sum_{s=1}^{t} b'_s \cdot \mathbf{1}\{X_s = v_i\} \\
&= \underbrace{\bar b \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\}}_{\text{(drift 항)}}
 + \underbrace{\sum_{s=1}^{t} (b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\}}_{=:\, M_t^{(i)}\;\text{(martingale 항)}}.
\end{aligned}
$$

`-` Drift 항: 점유율 $\rho_i(t) := t^{-1} \sum_{s \leq t} \mathbf{1}\{X_s = v_i\}$ 의 정의에서 $\sum_{s} \mathbf{1}\{X_s = v_i\} = t \cdot \rho_i(t) \to t \cdot \rho_i$ (장기 적립률 가정).
`-` Martingale 항 $M_t^{(i)}$: 다음 5.2 에서 마팅게일임을 확인.

## 5.2 $M_t^{(i)}$ 가 마팅게일임 — 한 줄

::: {.callout-note collapse="true" title="Proof: $M_t^{(i)}$ is a martingale"}

Let $\mathcal F_s := \sigma(X_1, b'_1, \dots, X_s, b'_s)$ be the natural filtration. The increment

$$\Delta M_s^{(i)} := M_s^{(i)} - M_{s-1}^{(i)} = (b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\}.$$

By the HST algorithm order, $X_s$ is determined from $\mathcal F_{s-1}$ *before* $b'_s$ is drawn; in particular $b'_s \perp (\mathcal F_{s-1}, X_s)$. Therefore

$$
\begin{aligned}
E\!\bigl[\Delta M_s^{(i)} \,\big|\, \mathcal F_{s-1}\bigr]
&= E\!\bigl[(b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\} \,\big|\, \mathcal F_{s-1}\bigr] \\
&= E\!\left[\mathbf{1}\{X_s = v_i\} \cdot \underbrace{E[b'_s - \bar b \,|\, \mathcal F_{s-1},\, X_s]}_{=\, 0} \,\bigg|\, \mathcal F_{s-1}\right]
   && (\because \text{tower property}) \\
&= 0
   && (\because b'_s \perp (\mathcal F_{s-1}, X_s)).
\end{aligned}
$$

So $\{M_t^{(i)}\}_{t \geq 0}$ is a martingale w.r.t.\ $\{\mathcal F_t\}$. The increment is also bounded:

$$|\Delta M_s^{(i)}| \;\leq\; |b'_s - \bar b| \cdot 1 \;\leq\; \tfrac{b}{2} \quad \text{a.s.} \qquad \square$$

:::

## 5.3 마팅게일 SLLN 적용 — $M_t^{(i)} / t \to 0$

이제 본문의 핵심. Durrett 4.5.3 의 가설을 점검한다.

`-` **Square integrable**: $|\Delta M_s^{(i)}| \leq b/2 \Rightarrow E[(M_t^{(i)})^2] \leq t \cdot (b/2)^2 < \infty$. ✓
`-` **Predictable quadratic variation**:

$$
\begin{aligned}
A_t
&= \sum_{s=1}^{t} E\!\bigl[(\Delta M_s^{(i)})^2 \,\big|\, \mathcal F_{s-1}\bigr] \\
&= \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\}^2 \cdot E\!\bigl[(b'_s - \bar b)^2 \,\big|\, \mathcal F_{s-1},\, X_s\bigr]
   && (\because \text{indicator squared} = \text{indicator}) \\
&= \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\} \cdot \operatorname{Var}(b'_s)
   && (\because b'_s \perp (\mathcal F_{s-1}, X_s),\; \mathrm{Var}(b'_s) = b^2/12) \\
&= \frac{b^2}{12} \cdot t \cdot \rho_i(t).
\end{aligned}
$$

`-` **$A_\infty = \infty$**: 장기 적립률 $\rho_i > 0$ 가정 하 $A_t \sim (b^2/12) \cdot \rho_i \cdot t \to \infty$. ✓

`-` $f(t) := t^{1/2 + \epsilon}$ 로 두면 $\int_0^\infty f(t)^{-2} dt = \int_0^\infty t^{-1-2\epsilon} dt < \infty$ (near $\infty$ 에서 수렴, near $0$ 에서는 $f \geq 1$ 의 clipping 으로 처리). ✓

Theorem 4.5.3 적용:

$$\frac{M_t^{(i)}}{f(A_t)} \to 0 \quad \text{a.s.}$$

$f(A_t) \sim (A_t)^{1/2+\epsilon} \sim t^{1/2+\epsilon}$ 이므로

$$\frac{M_t^{(i)}}{t^{1/2+\epsilon}} \to 0 \;\Rightarrow\; \frac{M_t^{(i)}}{t} \to 0 \quad \text{a.s.}$$

## 5.4 결과를 Theorem C 의 출발점으로 옮기기

5.1 의 분해에 5.3 의 결과를 대입하면

$$
\begin{aligned}
\frac{h(v_i, t)}{t}
&= \bar b \cdot \rho_i(t) + \frac{M_t^{(i)}}{t} \\
&\xrightarrow{t \to \infty} \bar b \cdot \rho_i + 0
   && (\because \rho_i(t) \to \rho_i,\; M_t^{(i)}/t \to 0) \\
&= \bar b\, \rho_i \quad \text{a.s.}
\end{aligned}
$$

즉

$$h(v_i, t) = \bar b\, \rho_i\, t + o(t) \quad \text{a.s.}$$

이게 Theorem C 증명의 Step 1 결과. Step 2 부터는 이 점근식의 차

$$\Delta h_{ij}(s) := h(v_i, s) - h(v_j, s) = \bar b\, (\rho_i - \rho_j)\, s + o(s)$$

를 제곱·합산해서 $\mathrm{SD}^2_{ij}(t)/t^3 \to \bar b^2(\rho_i - \rho_j)^2/3$ 를 끌어낸다 — Cesàro / Toeplitz 적용. 마팅게일 SLLN 이 들어가는 자리는 Step 1 단 한 번.

# §6 마무리

::: {.callout-important collapse="true" title="한 줄 요약 — 다섯 거리 survey 와의 연결"}

본 글의 메시지를 한 줄로 압축하면 — *"독립성을 양보하고 마팅게일 구조만 보존하면 SLLN 은 끝까지 살아남는다"* — 이다. HST 의 height 갱신은 walker 와 increment 가 따로 i.i.d. 인 *곱의 합* 이라 직접 i.i.d. SLLN 은 안 먹지만, 마팅게일 분해 한 번이면 마팅게일 SLLN 이 곧장 들어맞는다.

본인의 다섯 거리 survey (`260516_연구_HST_그래프도메인에서의거리`) 에서 다룬 다른 거리들 — diffusion, $R_{\mathrm{eff}}$, biharmonic, hitting time — 은 *정적* (snapshot) 그래프 거리라 SLLN 자체가 등장하지 않는다. 반면 HST 의 $\mathrm{SD}^2$ 만이 동적 (dynamic) 거리이고, 그 점근 분석의 *모든* 출발점에 마팅게일 SLLN 이 한 번 호출된다. 이게 HST 가 다섯 거리 family 와 본질적으로 다른 이유의 일부.

:::

# 참고

- R. Durrett, *Probability: Theory and Examples* (5th ed., 2019), Cambridge UP. §2.5.9 Kronecker's lemma, §4.5 Square Integrable Martingales (Thm 4.5.1/4.5.2/**4.5.3**, Example 4.5.4).
- 에이, *Three-Regime Convergence of the Snow Distance — Unified Proof under Model A*, `paper/260514_guebin/에이/abc오희석교수님공유용_v2.tex` (2026), §5 Theorem C Step 1.
- 유진, *Three-Regime Convergence 자세히 따라가기*, blog `260515_8b8636` (2026), §5 Theorem C 의 본인 한국어 해설.
