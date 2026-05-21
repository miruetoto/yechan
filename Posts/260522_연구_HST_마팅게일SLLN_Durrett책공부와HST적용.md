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

도서관 `BOOK_Durrett_Probability_TheoryAndExamples_5e.pdf` (490p) §4.5 Square Integrable Martingales 의 Theorem 4.5.3 (마팅게일 SLLN 의 표준 form) 을 한 번 따라가고, 그 결과가 에이의 ABC증명 (`paper/260514_guebin/에이/abc오희석교수님공유용_v2.tex`) §5 Theorem C 의 Step 1 (Height 분해) 에서 어떻게 쓰이는지를 짚는다. 한 줄 직관은 이것이다 — 마팅게일의 누적량은 증분이 유계이고 누적 조건부 분산이 발산할 때 시간에 비해 무시할 만큼 작다. 즉 "공정한 게임의 누적 손익은 시간 ratio 척도에서 평균 0". HST 의 height 갱신처럼 walker $X_s$ (과거 history 에 의존) 와 increment $b'_s$ (과거와 독립) 의 곱으로 표현되는 합은 항이 independent 가 아니어서 Kolmogorov i.i.d. SLLN 이 직접 안 먹지만, 마팅게일 분해 한 번이면 마팅게일 SLLN 이 곧장 들어맞는다.

# 1. Durrett Theorem 4.5.3

$\{X_n\}_{n \geq 0}$ 이 filtration $\{\mathcal F_n\}$ 에 대해 square integrable martingale (즉 $E[X_n^2] < \infty$, $E[X_n \mid \mathcal F_{n-1}] = X_{n-1}$, $X_0 \in \mathcal F_0$) 이라 하자. 그 predictable quadratic variation 을

$$A_n := \sum_{m=1}^{n} E\!\bigl[(X_m - X_{m-1})^2 \,\big|\, \mathcal F_{m-1}\bigr]$$

으로 정의하면 $A_n \in \mathcal F_{n-1}$ 이며 $n$ 에 대해 증가한다. Durrett 의 Theorem 4.5.3 은 — $f : [0, \infty) \to [1, \infty)$ 가 증가함수이고 $\int_0^\infty f(t)^{-2}\, dt < \infty$ 를 만족하면 —

$$\frac{X_n}{f(A_n)} \xrightarrow{n \to \infty} 0 \quad \text{a.s. on } \{A_\infty = \infty\}$$

라고 말한다. 가정 $A_\infty = \infty$ 는 "누적 분산이 발산한다" 는 가설이고, $f$ 의 전형은 $f(t) = (t \log^{1+\epsilon} t)^{1/2} \vee 1$ 또는 $f(t) = t^{1/2 + \epsilon}$ — 어느 쪽이든 $\sqrt{t}$ 보다 살짝 큰 정도다. 결론은 $X_n$ 이 $\sqrt{A_n}$ 보다 ε 만큼 천천히 자란다는 강한 진술. 증명은 세 도구의 합성으로 끝난다 — martingale transform (Thm 4.2.8), square integrable martingale 의 finite limit 보조정리 (Thm 4.5.2), Kronecker's lemma (Thm 2.5.9).

증명 자체를 한 줄씩 풀어쓴다. $H_m := f(A_m)^{-1}$ 로 두면 $A_m \in \mathcal F_{m-1}$ 이므로 $\{H_m\}$ 은 predictable, $f \geq 1$ 이라 $H_m \leq 1$ 로 bounded. Theorem 4.2.8 (bounded predictable transform 이 마팅게일을 보존) 적용으로

$$Y_n := (H \cdot X)_n := \sum_{m=1}^{n} \frac{X_m - X_{m-1}}{f(A_m)}$$

이 다시 마팅게일이 된다. $Y_n$ 의 predictable quadratic variation $B_n$ 을 한 step 계산하면

$$
\begin{aligned}
B_{n+1} - B_n
&= E\!\bigl[(Y_{n+1} - Y_n)^2 \,\big|\, \mathcal F_n\bigr] \\
&= E\!\left[\frac{(X_{n+1} - X_n)^2}{f(A_{n+1})^2} \,\bigg|\, \mathcal F_n\right]
   && (\because f(A_{n+1}) \in \mathcal F_n) \\
&= \frac{1}{f(A_{n+1})^2}\, E\!\bigl[(X_{n+1} - X_n)^2 \,\big|\, \mathcal F_n\bigr]
   && (\because \text{predictable factor out}) \\
&= \frac{A_{n+1} - A_n}{f(A_{n+1})^2}
   && (\because A \text{ 정의}).
\end{aligned}
$$

$f$ 가 증가하므로 $t \in [A_n, A_{n+1})$ 에서 $f(t) \leq f(A_{n+1}) \Rightarrow f(t)^{-2} \geq f(A_{n+1})^{-2}$, 따라서

$$B_\infty = \sum_{n=0}^{\infty} \frac{A_{n+1} - A_n}{f(A_{n+1})^2} \;\leq\; \sum_{n=0}^{\infty} \int_{[A_n,\, A_{n+1})} f(t)^{-2}\, dt \;\leq\; \int_0^\infty f(t)^{-2}\, dt \;<\; \infty$$

이 가설로 $B_\infty < \infty$ a.s. 가 보장된다. Theorem 4.5.2 (predictable quadratic variation 이 유한하면 square integrable martingale 은 finite limit 로 a.s. 수렴) 에 의해 $Y_n \to Y_\infty$ 가 유한한 극한으로 수렴. 끝으로 Kronecker's lemma — $a_n \uparrow \infty$ 이고 $\sum_n x_n/a_n$ 이 수렴하면 $a_n^{-1} \sum_{m=1}^n x_m \to 0$ — 에 $a_n := f(A_n)$, $x_m := X_m - X_{m-1}$ 을 대입하면 가설은 $\sum_m x_m/a_m = Y_n \to Y_\infty$ (유한) 으로 이미 확보되었고 $A_\infty = \infty$ 가설이 $a_n \uparrow \infty$ 를 준다. 따라서

$$\frac{X_n}{f(A_n)} \;=\; \frac{1}{f(A_n)} \sum_{m=1}^{n} (X_m - X_{m-1}) \xrightarrow{n \to \infty} 0 \quad \text{a.s.}$$

이걸로 증명 끝. 세 도구의 역할 분담을 한 줄로 보면 — transform 이 증분을 $f$ 로 가중해서 변동을 줄이고, Theorem 4.5.2 가 변동을 다 빼낸 $Y_n$ 의 finite limit 을 보장하고, Kronecker 가 가중합의 수렴을 비가중 합의 ratio 0 수렴으로 옮긴다.

이 정리의 가장 유명한 specialisation 은 Durrett 의 Example 4.5.4 인데, $\xi_1, \xi_2, \dots$ 가 독립이고 $E\xi_i = 0$, $E\xi_i^2 = \sigma_i^2 < \infty$ 일 때 $X_n := \xi_1 + \cdots + \xi_n$ 이 square integrable martingale 이고 $A_n = \sum_i \sigma_i^2$ ($\xi_i \perp \mathcal F_{i-1}$ 이라 조건부 분산이 무조건부 분산과 일치). $f(t) = (t \log^{1+\epsilon} t)^{1/2} \vee 1$ 로 두면 $\int_0^\infty f^{-2} dt$ 가 $\infty$ 근방에서 $\int dt/(t \log^{1+\epsilon} t)$ 처럼 수렴해서 가설을 채운다. $\sum_i \sigma_i^2 = \infty$ 이면 $X_n / \sqrt{A_n \log^{1+\epsilon} A_n} \to 0$ a.s.; i.i.d. 케이스 $\sigma_i^2 = \sigma^2$ 에서는 $A_n = n\sigma^2$ 이라 $X_n/n \to 0$ a.s. 가 따라와 고전 Kolmogorov SLLN 을 회복한다 (게다가 로그 보정까지). 이 example 이 보여주는 것은 — independence 가 사라져도 conditional 평균이 0 이라는 마팅게일 구조만 살아있으면 SLLN 이 끝까지 살아남는다는 사실이다.

# 2. HST Theorem C 의 출발점 적용

에이의 `abc오희석교수님공유용_v2.tex` §5 Theorem C (drift regime) 의 증명은 두 step 으로 갈리는데, Step 1 이 height 의 마팅게일 분해 + 마팅게일 SLLN 으로 점근식 $h(v_i, t) = \bar b \rho_i t + o(t)$ 를 끌어내고, Step 2 가 그걸 제곱·합산해서 $\mathrm{SD}^2_{ij}(t)/t^3 \to \bar b^2(\rho_i - \rho_j)^2/3$ 를 받아낸다. 마팅게일 SLLN 이 들어가는 자리는 Step 1 단 한 번. 그 한 번을 들여다본다.

HST height 의 정의 $h(v_i, t) := \sum_{s=1}^{t} b'_s \cdot \mathbf{1}\{X_s = v_i\}$ 에서 ($b'_s \sim \mathrm{Unif}(0, b)$ i.i.d., walker $\{X_s\}$ 와 독립, $\bar b := E[b'_s] = b/2$), 증분에서 평균을 빼내면

$$
\begin{aligned}
h(v_i, t)
&= \sum_{s=1}^{t} b'_s \cdot \mathbf{1}\{X_s = v_i\} \\
&= \underbrace{\bar b \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\}}_{\text{drift 항}}
 + \underbrace{\sum_{s=1}^{t} (b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\}}_{=:\, M_t^{(i)} \;\text{(martingale 항)}}
\end{aligned}
$$

로 갈라진다. Drift 항은 점유율 $\rho_i(t) := t^{-1} \sum_{s \leq t} \mathbf{1}\{X_s = v_i\}$ 의 정의에서 $\bar b \cdot t \cdot \rho_i(t)$, 장기 적립률 $\rho_i(t) \to \rho_i$ 가정 하 $\bar b \rho_i t + o(t)$. 남은 일은 $M_t^{(i)} / t \to 0$ a.s. 를 마팅게일 SLLN 으로 받아내는 것이다.

$M_t^{(i)}$ 가 마팅게일임은 한 줄. natural filtration $\mathcal F_s := \sigma(X_1, b'_1, \dots, X_s, b'_s)$ 에서 증분 $\Delta M_s^{(i)} = (b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\}$ 의 조건부 평균은

$$
\begin{aligned}
E\!\bigl[\Delta M_s^{(i)} \,\big|\, \mathcal F_{s-1}\bigr]
&= E\!\bigl[(b'_s - \bar b) \cdot \mathbf{1}\{X_s = v_i\} \,\big|\, \mathcal F_{s-1}\bigr] \\
&= E\!\left[\mathbf{1}\{X_s = v_i\} \cdot \underbrace{E[b'_s - \bar b \,|\, \mathcal F_{s-1}, X_s]}_{=\, 0} \,\bigg|\, \mathcal F_{s-1}\right]
   && (\because \text{tower property}) \\
&= 0
   && (\because b'_s \perp (\mathcal F_{s-1}, X_s)).
\end{aligned}
$$

HST 의 algorithm 순서가 핵심인데, 한 step 에서 $X_s$ 가 $\mathcal F_{s-1}$ 로부터 *먼저* 결정되고 $b'_s$ 는 그 *다음에* 추첨되므로 $b'_s \perp (\mathcal F_{s-1}, X_s)$ 가 성립. 그 다음 tower property 로 안쪽 조건부 평균이 0 이 되어 외부 평균도 0. 증분 크기는 $|\Delta M_s^{(i)}| \leq |b'_s - \bar b| \leq b/2$ 로 유계.

이제 Durrett 4.5.3 의 가설 점검. square integrable 은 $E[(M_t^{(i)})^2] \leq t(b/2)^2 < \infty$ 로 자동. predictable quadratic variation 은

$$
\begin{aligned}
A_t
&= \sum_{s=1}^{t} E\!\bigl[(\Delta M_s^{(i)})^2 \,\big|\, \mathcal F_{s-1}\bigr] \\
&= \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\}^2 \cdot E\!\bigl[(b'_s - \bar b)^2 \,\big|\, \mathcal F_{s-1}, X_s\bigr]
   && (\because \text{indicator}^2 = \text{indicator}) \\
&= \sum_{s=1}^{t} \mathbf{1}\{X_s = v_i\} \cdot \operatorname{Var}(b'_s)
   && (\because b'_s \perp (\mathcal F_{s-1}, X_s),\; \operatorname{Var}(b'_s) = b^2/12) \\
&= \frac{b^2}{12} \cdot t \cdot \rho_i(t).
\end{aligned}
$$

$\rho_i > 0$ 가정 하 $A_t \sim (b^2 \rho_i / 12) \cdot t \to \infty$ 이므로 $A_\infty = \infty$ ✓. $f(t) := t^{1/2 + \epsilon}$ 로 두면 $\int_0^\infty f^{-2} dt$ 는 $\infty$ 근방에서 $\int t^{-1-2\epsilon} dt$ 로 수렴 (near $0$ 에서는 $f \geq 1$ 의 clipping 으로 처리) ✓. Theorem 4.5.3 적용 결과는

$$\frac{M_t^{(i)}}{f(A_t)} \to 0 \;\text{a.s.} \quad\Rightarrow\quad \frac{M_t^{(i)}}{t^{1/2 + \epsilon}} \to 0 \;\Rightarrow\; \frac{M_t^{(i)}}{t} \to 0 \;\text{a.s.}$$

이걸 원래 분해에 대입하면

$$\frac{h(v_i, t)}{t} = \bar b \cdot \rho_i(t) + \frac{M_t^{(i)}}{t} \xrightarrow{t \to \infty} \bar b \rho_i \;\text{a.s.} \quad\Leftrightarrow\quad h(v_i, t) = \bar b \rho_i t + o(t) \;\text{a.s.}$$

가 되고, 이게 Theorem C 증명 Step 1 의 결론이다. Step 2 부터는 $\Delta h_{ij}(s) = h(v_i, s) - h(v_j, s) = \bar b(\rho_i - \rho_j) s + o(s)$ 의 제곱을 시간에 대해 합산해서 $\mathrm{SD}^2_{ij}(t)/t^3 \to \bar b^2(\rho_i - \rho_j)^2/3$ 를 Cesàro / Toeplitz 형 평가로 받아낸다 — 그 단계엔 마팅게일이 다시 등장하지 않는다. 즉 Theorem C 전체에서 마팅게일 SLLN 은 정확히 한 번, height 의 평균제거된 잔차 $M_t^{(i)}$ 가 $t$ 보다 천천히 자란다는 사실을 받는 데에만 호출된다.

본인이 작성한 다섯 거리 survey (`260516_연구_HST_그래프도메인에서의거리`) 에서 다룬 다른 거리들 — diffusion, $R_{\mathrm{eff}}$, biharmonic, hitting time — 은 모두 정적 (snapshot) 그래프 거리라 시간 점근의 SLLN 자체가 들어갈 자리가 없다. HST 의 $\mathrm{SD}^2$ 만이 누적 동적 (cumulative dynamic) 거리이고, 그 점근 분석의 모든 출발점에 마팅게일 SLLN 이 한 번 호출된다. 이게 HST 가 다섯 거리 family 와 본질적으로 다른 위치를 차지하는 이유 중 하나.

# 참고

- R. Durrett, *Probability: Theory and Examples* (5th ed., 2019), Cambridge UP. §2.5.9 Kronecker's lemma, §4.5 Square Integrable Martingales (Thm 4.5.1 / 4.5.2 / **4.5.3**, Example 4.5.4).
- 에이, *Three-Regime Convergence of the Snow Distance — Unified Proof under Model A*, `paper/260514_guebin/에이/abc오희석교수님공유용_v2.tex` (2026), §5 Theorem C Step 1.
- 유진, *Three-Regime Convergence 자세히 따라가기*, blog `260515_8b8636` (2026), §5 Theorem C 의 한국어 해설.
