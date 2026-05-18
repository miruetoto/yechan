---
title: 책공부 ▷ Brémaud Ch.4 ▷ Long-Run Behavior — 팟캐스트
author: 연
date: 05/18/2026
draft: false
output-file: 260518_618a48.html
fontsize: 0.85em
---

<style>
.math.display { text-align: left; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
</style>

> **원전**: Pierre Brémaud, *Markov Chains: Gibbs Fields, Monte Carlo Simulation and Queues* (2nd ed., Springer 2020), Chapter 4 "Long-Run Behavior" (pp. 145–165).
> `ref/BOOK_PierreBremaud_MarkovChain.pdf`

---

# Ch.4 Long-Run Behavior — 팟캐스트

> *(스튜디오, 화이트보드 하나. 호스트 H 와 게스트 G 가 마주 앉아 있다.)*

**H.** 안녕하세요, 오늘은 마르코프 체인의 "장기 거동" 을 다뤄봅니다. 지난 시간에 recurrence, ergodicity, stationary distribution 얘기를 했는데 — 오늘 핵심 질문은 하나예요. **"시작점을 잊는가?"**

**G.** 네, 그게 Ch.4 전체를 관통하는 질문이에요. 어떤 초기분포 $\mu$ 에서 출발해도, 충분히 오래 걸으면 체인의 분포가 하나의 정상분포 $\pi$ 로 수렴하느냐. 수렴한다면 얼마나 빨리.

---

## 에피소드 1: Total Variation — "두 분포가 얼마나 다른가"

**H.** 수렴을 말하려면 먼저 "두 분포가 가깝다" 를 재야겠죠.

**G.** 맞아요. 여기서 쓰는 자가 **Total Variation distance** 입니다.

$$d_V(\alpha, \beta) = \sup_{A \subset E} |\alpha(A) - \beta(A)| = \frac{1}{2} \sum_{i \in E} |\alpha(i) - \beta(i)|$$

두 분포가 가장 크게 다를 수 있는 사건 $A$ 를 찾아서 그 확률 차이를 재는 거예요.

**H.** 직관적으로?

**G.** 동전 두 개가 있다고 해요. 하나는 앞면 60%, 하나는 앞면 40%. 그러면 $A = \{\text{앞면}\}$ 에서 차이가 가장 크죠: $d_V = |0.6 - 0.4| / 2 = 0.1$... 아, 아니 이건 두 점이라 그냥 $|0.6 - 0.4| = 0.2$, 그리고 절반이니까 $d_V = 0.1$. 어쨌든 $d_V \in [0, 1]$ 이고, 0 이면 같은 분포, 1 이면 완전히 disjoint 한 support.

**H.** 이걸로 "마르코프 체인의 $n$-step 분포가 $\pi$ 에 가까워지는가" 를 잴 수 있는 거군요.

**G.** 정확해요. $d_V(\mu^\top \mathbf{P}^n, \pi) \to 0$ 인가, 그리고 얼마나 빨리.

---

## 에피소드 2: Coupling — "두 체인을 한 지붕 아래"

**H.** 그런데 TV distance 를 직접 계산하는 건 어렵잖아요?

**G.** 네, 그래서 나온 도구가 **coupling** 입니다. 아이디어가 아름다워요. 두 체인 $\{X'_n\}$, $\{X''_n\}$ 를 **같은 확률 공간** 에 올려놓고, 언젠가 **만나게** 만드는 거예요. 만나는 시점을 $\tau$ 라 하면:

$$n \geq \tau \implies X'_n = X''_n$$

만난 이후로는 같이 걸으니까, 분포 차이가 사라지죠.

**H.** 그래서 핵심 부등식이?

**G.** Brémaud 의 Theorem 4.1.7, **fundamental coupling inequality**:

$$d_V(X'_n, X''_n) \leq \Pr(\tau > n)$$

coupling time $\tau$ 가 빨리 유한해지면 TV distance 가 빨리 0 으로 갑니다. 이 한 줄로 수렴 속도 문제가 "coupling time 의 tail bound" 문제로 변환돼요.

**H.** 증명이 길어요?

**G.** 네 줄이에요. 모든 사건 $A$ 에 대해:

$$\Pr(X'_n \in A) - \Pr(X''_n \in A) = \Pr(X'_n \in A, \tau \leq n) + \Pr(X'_n \in A, \tau > n) - \Pr(X''_n \in A)$$

$\tau \leq n$ 이면 $X'_n = X''_n$ 이니 첫 항은 $\Pr(X''_n \in A, \tau \leq n)$ 과 같아서 상쇄, 남는 건 $\Pr(X'_n \in A, \tau > n) \leq \Pr(\tau > n)$. 끝.

**H.** 네 줄이라니 아름답네요.

---

## 에피소드 3: 정상상태로의 수렴 — Positive Recurrent Case

**H.** 이제 실제 수렴 정리를 봅시다.

**G.** Brémaud Theorem 4.2.1. **ergodic chain** (= irreducible + positive recurrent + aperiodic) 이면:

$$\lim_{n \to \infty} d_V(\mu^\top \mathbf{P}^n, \nu^\top \mathbf{P}^n) = 0$$

**모든** 초기분포 $\mu, \nu$ 에 대해. 특히 $\nu = \pi$ (정상분포) 잡으면:

$$\lim_{n \to \infty} |\mu^\top \mathbf{P}^n - \pi^\top| = 0$$

**어디서 시작하든 $\pi$ 로 간다.**

**H.** 증명 전략은?

**G.** Coupling 입니다. 두 독립 체인 $X^{(1)}_n \sim \mu$, $X^{(2)}_n \sim \nu$ 를 만들고, 둘이 만나면 합치는 **product chain** $Z_n = (X^{(1)}_n, X^{(2)}_n)$ 을 돌려요. Product chain 이 대각선 $\{(i, i)\}$ 에 도달하면 coupling 성공.

핵심: aperiodic + irreducible 이면 product chain 도 irreducible + aperiodic. 정상분포 $\pi(i)\pi(j)$ 가 존재하니 대각선에 양의 확률로 도달 → $\tau < \infty$ a.s. → coupling inequality 로 TV 수렴.

**H.** Aperiodic 이 왜 필요해요?

**G.** Period $d > 1$ 이면 체인이 cyclic class $C_0, C_1, \ldots, C_{d-1}$ 을 돌아요. $C_0$ 에서 출발하면 짝수 시점에만 $C_0$, 홀수에만 $C_1$, ... 이런 식이라 product chain 의 대각선이 $d$ 스텝마다만 접근 가능 → irreducibility 가 깨져요.

**G.** 근데 periodic 도 완전 포기는 아닙니다. Theorem 4.2.3 이 커버: period $d$ 체인은 $\mathbf{P}^{nd}$ 가 각 cyclic class 안에서 수렴. 즉 **$d$ 스텝 건너뛰기로 보면 수렴**.

---

## 에피소드 4: Null Recurrent — "느린 잊기"

**H.** Positive recurrent 가 아니면요?

**G.** Null recurrent (Theorem 4.2.4): $p_{ij}(n) \to 0$ for all $i, j$. 정상분포가 존재하지 않고, 전이확률이 0 으로 깔려요. 수렴이 아니라 **소산**.

**H.** 직관이?

**G.** 1차원 대칭 random walk. 원점 재방문 확률은 1 (recurrent) 이지만 평균 재방문 시간이 $\infty$ (null). 시간이 지나면 walker 가 어디 있는지 점점 더 불확실 — 분포가 넓게 퍼져서 어떤 한 점에 있을 확률이 0 으로.

---

## 에피소드 5: 수렴 속도 — Coupling Rate 와 Perron-Frobenius

**H.** 수렴한다는 건 알겠는데, **얼마나 빨리** 수렴해요?

**G.** 두 가지 도구가 있어요.

### 도구 1: Coupling rate

Theorem 4.3.1: coupling time $\tau$ 에 대해 $\mathbb{E}[\psi(\tau)] < \infty$ 이면

$$|\mu^\top \mathbf{P}^n - \nu^\top \mathbf{P}^n| = o\!\left(\frac{1}{\psi(n)}\right)$$

$\psi$ 가 빨리 커지면 수렴이 빠른 거죠. 유한상태면 $\psi(n) = e^{\alpha n}$ 잡을 수 있어서 **기하적 수렴** (Theorem 4.3.3).

### 도구 2: Perron-Frobenius

유한상태이고 primitive (= irreducible + aperiodic) 면 spectral decomposition:

$$\mathbf{P}^n = \mathbf{1}\pi^\top + O(n^{m_2 - 1}|\lambda_2|^n)$$

$\lambda_2$ = second-largest eigenvalue modulus (**SLEM**). 수렴 속도가 $|\lambda_2|^n$ 에 지배.

**H.** 예를 들면?

**G.** Two-state chain $E = \{1, 2\}$:

$$\mathbf{P} = \begin{pmatrix} 1-\alpha & \alpha \\ \beta & 1-\beta \end{pmatrix}$$

$\lambda_2 = 1 - \alpha - \beta$. 수렴 속도:

$$(\mathbf{P}^n - \mathbf{P}^\infty) = \frac{(1-\alpha-\beta)^n}{\alpha+\beta} \begin{pmatrix} \alpha & -\alpha \\ -\beta & \beta \end{pmatrix}$$

$\alpha + \beta$ 가 클수록 (두 상태 사이 전이가 활발할수록) $|\lambda_2|$ 가 작아서 빨리 수렴. $\alpha = \beta = 0.5$ 면 $\lambda_2 = 0$, 한 스텝에 수렴!

---

## 에피소드 6: Ehrenfest Model — 볼츠만 vs 체르멜로

**H.** 마지막으로 물리 예제 하나 갈게요.

**G.** Ehrenfest model. 상자 두 개 $A$, $B$ 에 입자 $N$ 개. 매 시점 입자 하나를 무작위로 골라 다른 상자로 옮김. 상태 = $A$ 안 입자 수 $X_n \in \{0, 1, \ldots, N\}$.

**H.** 정상분포는?

**G.** $\pi(k) = \binom{N}{k} 2^{-N}$. 이항분포! 양쪽 상자에 반반 있는 게 가장 높은 확률.

**H.** 열역학 비가역성이랑 무슨 상관이에요?

**G.** 볼츠만은 "엔트로피 증가 법칙" 으로 비가역성을 주장했어요. 체르멜로가 반박: "Poincaré recurrence 에 의해 모든 상태를 언젠가 재방문하잖아 — 비가역성이 어딨어?" 

Ehrenfest 가 이 논쟁을 해결했어요. **둘 다 맞지만 시간 스케일이 다릅니다.** 평형 ($N/2$) 에서 벗어난 상태 $L$ 로부터의 평균 재방문 시간은:

$$\frac{1}{2L} \cdot 2^{2L}(1 + O(L))$$

$L = 10^6$ 이면 $10^{-5}$ 초에 평형 도달, 그러나 텅 빈 상태로 재방문은 $\frac{1}{2 \cdot 10^{11}} \times 2^{2 \times 10^6}$ 초 — 우주 나이보다 긴 시간.

**G.** 이게 Brémaud 이 "statistical equilibrium and recurrence are not antagonistic" 이라고 쓴 이유예요. Recurrence 는 맞지만, **관측 가능한 시간 안에는 비가역적으로 보인다**. 마르코프 체인 이론이 그 정량적 설명을 줍니다.

**H.** 그리고 이게 Newton 의 냉각 법칙으로도 연결되죠?

**G.** 네. $N \to \infty$ 극한에서 $\mathbb{E}[X_n - L] = (i - L)(1 - 1/L)^n \to (i - L)e^{-\gamma t}$ — 지수 감쇠. 냉각 법칙이 마르코프 체인의 평균 수렴에서 나옵니다.

---

## HST 와의 연결 (에피소드 보너스)

**H.** 이 Ch.4 가 우리 회사 본업인 HST 랑 어떻게 연결돼요?

**G.** 세 가지.

**첫째**, Theorem 4.2.1 (ergodic 수렴) 이 **HST 에 직접 적용 안 됩니다**. HST 의 $\mathbf{P}(t)$ 가 non-homogeneous (높이 의존) 라서요. 그래서 HST 논문은 Hajnal-Bartlett 의 **weak ergodicity** 로 우회 — "두 초기분포가 서로 잊지만, 특정 $\pi$ 로 수렴하진 않는다."

**둘째**, coupling inequality ($d_V \leq \Pr(\tau > n)$) 는 HST 의 Doeblin minorization 증명 (에이 Theorem A) 에서 핵심 도구. 에이가 쓴 $K^*$-step minorization 이 본질적으로 coupling construction.

**셋째**, Perron-Frobenius 의 SLEM ($|\lambda_2|$) 이 HST_0 (block-free variant) 의 $SD^2/t^2$ 상수를 결정 — 본인이 증명한 `hst0_variance_via_Z.tex` Theorem 2 의 spectral filter $(1+\lambda)/(1-\lambda)$. 유한상태 ergodic chain 이면 SLEM 이 수렴 속도이자 SD 상수의 핵심 파라미터.

**H.** 정리하면: "시작점을 잊는가?" 에 대한 답은 **"잊는다, coupling 덕분에. 속도는 spectral gap 이 결정한다"** — 그리고 HST 는 non-homogeneous 라 이 깔끔한 그림에서 벗어나 있고, 그게 HST 의 도전이자 독특함이다.

**G.** 완벽한 요약이에요.

**H.** Ch.4 마칩니다. 다음 시간에는 Ch.6 의 Fundamental Matrix $\mathbf{Z}$ 로 가보겠습니다. 감사합니다.

**G.** 감사합니다.

---

## 참고

| 절 | Brémaud 정리 | 핵심 |
|---|---|---|
| §4.1 | Lemma 4.1.2, Thm 4.1.3, Thm 4.1.7 | TV distance, maximal coincidence, **coupling inequality** |
| §4.2.1 | Thm 4.2.1 | Ergodic → TV convergence |
| §4.2.2 | Thm 4.2.4 | Null recurrent → $p_{ij}(n) \to 0$ |
| §4.2.3 | Thm 4.2.3 | Periodic ($d > 1$) → $d$-step 수렴 |
| §4.3.1 | Thm 4.3.1, 4.3.3 | Coupling rate, geometric coupling |
| §4.3.2 | Thm 4.3.8 (Perron-Frobenius) | SLEM $|\lambda_2|$ → 기하적 수렴 |
| Ex 4.2.5 | Ehrenfest model | 열역학 비가역성의 마르코프 체인 설명 |
| Ex 4.3.4 | Two-state chain | $\lambda_2 = 1-\alpha-\beta$, 명시적 수렴 |

---

*작성: 연 (연구팀), 2026-05-18*
*원전: Pierre Brémaud, Markov Chains (2nd ed., 2020), Chapter 4.*
