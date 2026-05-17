---
title: 연구 ▷ HST ▷ Three-Regime Convergence 자세히 따라가기
author: 유진
date: 05/15/2026
draft: false
output-file: 260515_8b8636.html
fontsize: 0.85em
---

```{=html}
<style>
.math.display { font-size: 0.9em; text-align: left; }
.callout { font-size: 0.9em; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
</style>
```

> **본 글의 위치**: ABC증명(페이퍼팀)의 `모델A통합증명_ABC.tex` (★★★★★, 2026-05-17) 의 **한국어 해설본**. 원문 섹션 구조 그대로 따라가되, 수식 전개는 한 줄씩 풀어 쓰고 callout 3색(빨강 직관 / 빨강 큰 그림 / 파랑 증명) 짝으로 정리.

Heavy-Snow Transform의 snow distance

$$SD^2_{ij}(t) := \sum_{s=0}^t \bigl(h(v_i, s) - h(v_j, s)\bigr)^2$$

의 점근 거동은 **장기 적립률**

$$\rho_i := \lim_{t \to \infty} \tfrac{1}{t} \sum_{s=1}^t \mathbf{1}\{X_s = v_i\}$$

에 따라 **세 가지 regime**으로 갈린다:

| Regime | 조건 | 스케일링 | 극한 |
|---|---|---|---|
| **A. Balanced** | $\rho_i = 1/n$ for all $i$ | $SD^2_{ij}/t$ | $\to c_{ij}$ (Foster–Lyapunov + Doeblin) |
| **B. Intermediate** | $\rho_i = \rho_j$, but globally unbalanced | $SD^2_{ij}/t^2$ | $\to \sigma^2_{ij}/2$ (diffusive 가정 하) |
| **C. Drift** | $\rho_i \neq \rho_j$ | $SD^2_{ij}/t^3$ | $\to \bar{b}^2(\rho_i - \rho_j)^2/3$ |

`-` **Theorem A (Balanced)** — Foster–Lyapunov drift 와 Doeblin minorization 으로 augmented chain의 positive Harris recurrence 확보. 가장 어려움.
`-` **Theorem B (Intermediate)** — 가중 Cesàro lemma 한 줄. recurrence 기계 불필요.
`-` **Theorem C (Drift)** — 마팅게일 SLLN + 급수 계산. 역시 recurrence 불필요.

> **모델 노트.** 본 글은 **모델 A** (per-step i.i.d. $b'_s \sim \text{Unif}(0,b)$) baseline. 2026-05-17 01:10 대표 결정으로 회사 공식 채택은 모델 B (per-round constant $b'^{(r)}$) 이지만, 본 글은 reference baseline 역할. 모델 B 손질은 후속 포스트로.

---

# §1 세팅과 표기

## 그래프와 파라미터

- $\mathcal{G} = (V, \mathbf{E}, \mathbf{W})$: 연결 가중 그래프, $V = \{v_1, \ldots, v_n\}$, $|V| = n < \infty$
- $b > 0$: snowfall scale (한 step당 증분 상한)
- $T_{\max} \in \mathbb{N}$: 연속 flow 횟수 상한
- $\boldsymbol{\mu}_0$: $V$ 위 확률측도, **full support** ($\mu_{\min} := \min_i \boldsymbol{\mu}_0(v_i) > 0$)
- 보통 $\boldsymbol{\mu}_0 \propto \deg$ (degree-proportional) 또는 uniform

## 상태 변수

각 시점 $t \geq 0$에서:

- $h(v_i, t) \in \mathbb{R}$: 노드 $v_i$의 **눈 높이**
- $X_t \in V$: **walker** 위치
- $Z_t \in \{0, 1, \ldots, T_{\max}\}$: **block timer** (연속 flow 횟수)

초기 조건: $h(v_i, 0) = 0$ for all $i$ (zero signal), $X_0 \sim \boldsymbol{\mu}_0$, $Z_0 := T_{\max}$ (첫 step은 fresh fall).

## Augmented state (증강 상태)

분석 대상 Markov chain:

$$S_t := (\boldsymbol{\delta}_t,\, X_t,\, Z_t) \in \mathcal{X}^* := \mathbb{R}^{n-1} \times V \times \{0, 1, \ldots, T_{\max}\}$$

여기서

$$\boldsymbol{\delta}_t := \bigl(h(v_2, t) - h(v_1, t),\ \ldots,\ h(v_n, t) - h(v_1, t)\bigr) \in \mathbb{R}^{n-1}$$

은 **상대 높이 벡터** ($v_1$ 기준). 절대 높이 $h$는 $t$에 따라 무한히 커지므로 차이만 추적.

`-` $\mathcal{F}_t := \sigma(X_0, b'_1, X_1, \ldots, b'_t, X_t)$ — 자연 filtration
`-` $\mathcal{B}(\mathcal{X}^*)$ — Borel σ-algebra
`-` $P^m(s, B)$ — $S_0 = s$에서 $m$ step 후 $B$ 도달 확률

## 중심화 높이와 Lyapunov 함수

$$\bar h(t) := \tfrac{1}{n}\sum_i h(v_i, t), \qquad \hbar(v_i, t) := h(v_i, t) - \bar h(t), \qquad \Phi(t) := \sum_i \hbar(v_i, t)^2$$

- $\hbar$: **centered height** — 평균을 뺀 상대 높이
- $\Phi$: **Lyapunov 함수** — 노드 간 높이 차이의 quadratic 측도. $\Phi$가 크면 지형이 울퉁불퉁, $\Phi = 0$이면 평평.
- $M(t) := \max_i |\hbar(v_i, t)|$ — (단측) centered height range. $\sum \hbar = 0$이므로 $\max \hbar - \min \hbar \in [M(t), 2M(t)]$.

## 네 가지 후보 알고리즘 (A / B / C / D)

HST에서 walker 규칙 (block / flow)은 정해져 있지만, **step당 증분 $b'$의 추출 방식**은 후보 4종:

| 모델 | 증분 추출 | step당 $b'$ |
|---|---|---|
| **A** (per-step i.i.d.) | 매 step | $b'_s \sim \text{Unif}(0, b)$, i.i.d. across $s$ |
| **B** (per-round constant) | 매 round $r$ | $b'^{(r)} \sim \text{Unif}(0, b)$; round 내 모든 step 동일 |
| **C** (per-round w/ decay) | 매 round $r$ | $b'^{(r)} f(s)$ — round 내 감쇠 함수 $f$ |
| **D** (deterministic) | — | $b'_s = b$ (상수) |

**모멘트**: A, B, C 모두 $\mathbb{E}[b'_s] = b/2$, $\mathbb{E}[(b'_s)^2] = b^2/3$. D는 $b'_s = b$ 상수.

### Algorithm A (본 글의 대상)

매 step $t+1$ ($t \geq 0$):

1. **증분 추첨**: $b'_{t+1} \sim \text{Unif}(0, b)$ — 과거와 독립
2. **walker 선택**:
   - **Block step** ($Z_t = T_{\max}$, timer 만료): $X_{t+1} \sim \boldsymbol{\mu}_0$ (fresh fall)
   - **Flow step**: $X_t$의 이웃 중 $h \leq h(X_t, t)$인 downstream set $\mathcal{D}_t$. 비어 있으면 stall → block 처리. 아니면 $\mathcal{D}_t$ 안에서 degree-weighted 선택.
3. **높이 갱신**: $h(v_i, t+1) = h(v_i, t) + b'_{t+1}\,\mathbf{1}\{v_i = X_{t+1}\}$
4. **timer 갱신**: $Z_{t+1} = 0$ (block step) 또는 $Z_t + 1$ (flow step)

`-` $t_r$: $r$번째 **block-flag time**
`-` $m_r := t_{r+1} - t_r \in [1, T_{\max} + 1]$: 라운드 길이

### Algorithm B/C/D (참고)

- **B**: round 시작 시 $b'^{(r)} \sim \text{Unif}(0, b)$ 1회 추첨, round 내 모든 step에서 같은 $b'^{(r)}$ 사용. 라운드 총 적립 = $m_r \cdot b'^{(r)}$.
- **C**: round당 $b'^{(r)}$ 추첨 + decay 함수 $f$ (지수/선형 등). 라운드 내 증분이 점차 작아짐.
- **D**: $b'_s = b$ 상수. 무작위성 없음. → lattice irreducibility 문제로 채택 불가 (과거 시도).

> 회사 공식 채택 모델: **B** (2026-05-17 01:10). 본 글은 A baseline.

---

# §2 Three-regime 분류

## 정의

::: {.callout-note collapse="false" title="Definition (Regime classification)"}

$\rho_i$가 모든 $i \in V$에 대해 a.s. 존재한다고 가정.

- **Balanced regime**: 모든 $i$에 대해 $\rho_i = 1/n$. $\sum_i \rho_i = 1$이므로 이는 "모든 $i, j$에 대해 $\rho_i = \rho_j$ **and** global balance"와 동치.
- **Drift regime (pair $i, j$)**: $\rho_i \neq \rho_j$.
- **Intermediate regime (pair $i, j$)**: $\rho_i = \rho_j$이지만 **globally unbalanced** (어떤 $k, \ell$에 대해 $\rho_k \neq \rho_\ell$).

:::

regime은 **노드 쌍 $(v_i, v_j)$ 단위**로 매겨진다. 같은 그래프 안에서도 어떤 쌍은 balanced, 어떤 쌍은 intermediate, 또 어떤 쌍은 drift일 수 있다 (Helm 그래프 사례 참조).

## 실험 분류 (소미·유진 41f2b0 인용)

소미 작성 블로그 *Thm A, B의 실제 검증* (41f2b0, 2026-05-15) 의 시뮬레이션 결과를 인용. 모두 **degree-proportional $\boldsymbol{\mu}_0$** 사용 ($\tau \sim 10^6$, $b = 0.5$, Algorithm A).

| 그래프 | Pair $(v_i, v_j)$ | $\hat\rho$ 관계 | Regime |
|---|---|---|---|
| $K_n$ (regular) | any | $\hat\rho_i = 1/n$ | A (balanced) |
| Parity cycle $C_n$ | any | $\hat\rho_i = 1/n$ | A (balanced) |
| Directed cycle $C_n$ | any | $\hat\rho_i = 1/n$ | A (balanced) |
| Wheel $W_n$ | any | $\hat\rho_i \approx 1/n$ | A (balanced) |
| Star $S_n$ | hub–leaf | $\hat\rho_{\mathrm{h}} \neq \hat\rho_{\mathrm{l}}$ | C (drift) |
| Helm $H_k$ | hub–leaf | $\hat\rho_{\mathrm{h}} \neq \hat\rho_{\mathrm{l}}$ | C (drift) |
| Helm $H_k$ | ring–ring | $\hat\rho_{\mathrm{r}} = \hat\rho_{\mathrm{r}}$, ring로 연결됨 | A (balanced within ring) |
| Helm $H_k$ | leaf–leaf | $\hat\rho_{\mathrm{l}} = \hat\rho_{\mathrm{l}}$, globally unbalanced | B (intermediate) |

## Reading guide

::: {.callout-important collapse="true" title="해석: 그래프 구조와 $\boldsymbol{\mu}_0$가 regime을 결정하는 방식"}

**Regular graph + degree-prop $\boldsymbol{\mu}_0$**: degree가 모두 같으므로 $\mu_0 \propto \deg$ = uniform. 모든 노드 fall 균등 → $\hat\rho_i = 1/n$ → 모든 쌍 **balanced**. doubly stochastic kernel (예: directed cycle)도 같은 결론.

**Non-regular graph + degree-prop $\boldsymbol{\mu}_0$**: high-degree 노드가 fall을 더 자주 받음. **Star $S_n$**: hub 차수 $n-1$, leaf 차수 $1$ → $\hat\rho_{\mathrm{h}} \neq \hat\rho_{\mathrm{l}}$ → hub–leaf 쌍 **drift**.

**Helm $H_k$ (3-class hierarchy)**: hub, ring, leaf 세 등급. $\hat\rho_{\mathrm{h}} > \hat\rho_{\mathrm{r}} > \hat\rho_{\mathrm{l}}$.

- **hub–leaf**: $\rho$ 다름 → drift
- **ring–ring**: $\rho$ 같음 + ring sub-graph로 연결 → 그 안에서 balanced dynamics → $SD^2/t$ 수렴 (balanced "within ring")
- **leaf–leaf**: $\rho$ 같음, 그러나 두 leaf는 hub를 거쳐야만 연결. **hub의 linearly growing height가 flow를 통해 noise로 전파, 두 leaf는 직접 연결 X → noise가 독립** → 높이차가 random walk → $|D_t| \sim \sigma \sqrt{t}$ → SD²/t² 수렴. → **intermediate**.

:::

## 완전성 (exhaustiveness)

위 세 regime은 **상호 배반 + 합쳐서 완전**: $\rho_i, \rho_j$가 존재하는 모든 쌍 $(i, j)$가 정확히 한 부류. $t, t^2, t^3$ 외의 다른 power scaling은 발생 안 함. 상세는 [`regime_exhaustiveness.tex`](../paper/260514_guebin/해리스/regime_exhaustiveness.tex) (해리스) 참조.

## $\rho_i$의 존재성

regime 분류는 $\rho_i$의 a.s. 존재를 전제. random-step 변형에서는 $\mu_{\min} > 0$인 모든 연결 그래프에 대해 **무조건 성립** (오큐리의 stochastic approximation + ODE 인자, [`occupation_slln_general.tex`](../paper/260514_guebin/오큐리/occupation_slln_general.tex) 참조). Theorem A는 $\rho_i = 1/n$을 별도 가정하므로 $\rho_i$ 존재성을 자체 입력으로 받음. Theorem B, C에서는 $\rho_i$ 존재를 주어진 입력으로 사용.

---

---

# §3 Theorem A: Balanced regime

Balanced regime ($\rho_i = 1/n$) 에서는 drift가 사라지므로 $SD^2 = \sum \xi(s)^2$이 되고, 이것의 $SD^2/t \to c_{ij}$ 수렴을 보이려면 에르고딕 정리가 필요. 증명 체인:

$$\text{DT (B(S))} \to \text{round drift} \to \text{Doeblin} \to \psi\text{-irred.} \to \text{Harris rec.} \to \text{moment bound} \to \text{LLN}$$

본 글의 모멘트:

$$\bar b := \mathbb{E}[b'_t] = b/2, \qquad \mathbb{E}[(b'_t)^2] = b^2/3, \qquad \text{Var}(b'_t) = b^2/12.$$

`-` $\mathcal{B}_\Phi$ — **Lyapunov 함수 $\Phi$의 sublevel set** (아래첨자 $\Phi$는 "$\Phi$로 정의된 set"):

$$\mathcal{B}_\Phi := \{S \in \mathcal{X}^* : \Phi(S) \leq R\}, \qquad R := n\left(\tfrac{C_{\text{tight}}+1}{\varepsilon}\right)^{\!2}$$

`-` 정상 분포는 $\pi^*$로 표기 (존재성/유일성은 증명 체인의 결과)

## §3.1 Centred height: two preparatory facts

**Fact 1 ($\hbar$의 single-step update).** 시점 $t \geq 0$, step $t+1$에서 증분 $b'_{t+1} \sim \text{Unif}(0,b)$이 walker $X_{t+1}$에 적립될 때:

$$\hbar(v_i,\, t+1) = \begin{cases} \hbar(v_i,\, t) + b'_{t+1}\dfrac{n-1}{n}, & v_i = X_{t+1} \\[4pt] \hbar(v_i,\, t) - \dfrac{b'_{t+1}}{n}, & v_i \neq X_{t+1} \end{cases}$$

::: {.callout-note collapse="true" title="증명: Fact 1 ($\\hbar$의 single-step update)"}

높이 갱신: $h(v_i, t+1) = h(v_i, t) + b'_{t+1}\,\mathbf{1}\{v_i = X_{t+1}\}$.

평균 갱신: $\bar h(t+1) = \tfrac{1}{n}\sum_i h(v_i, t+1) = \bar h(t) + \tfrac{b'_{t+1}}{n}$.

$\hbar(v_i, t+1) := h(v_i, t+1) - \bar h(t+1)$에 대입.

**$v_i = X_{t+1}$:**
$$\hbar(X_{t+1}, t+1) = \bigl[h(X_{t+1}, t) + b'_{t+1}\bigr] - \bigl[\bar h(t) + \tfrac{b'_{t+1}}{n}\bigr] = \hbar(X_{t+1}, t) + b'_{t+1}\,\tfrac{n-1}{n}$$

**$v_i \neq X_{t+1}$:**
$$\hbar(v_i, t+1) = h(v_i, t) - \bigl[\bar h(t) + \tfrac{b'_{t+1}}{n}\bigr] = \hbar(v_i, t) - \tfrac{b'_{t+1}}{n} \qquad \square$$

:::

**Fact 2.** $\sum_i \hbar(v_i, t) = 0$ for all $t$ (평균을 뺀 정의에서 자명).

## §3.2 Single-step $\Phi$ drift

**Lemma ($\Phi$ drift).** 시점 $t \geq 0$, $X_{t+1} \in V$에서:

$$\mathbb{E}\bigl[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t,\, X_{t+1}\bigr] = b\,\hbar(X_{t+1}, t) + \frac{(n-1)b^2}{3n}$$

::: {.callout-important collapse="true" title="해석: $\\Phi$ drift의 의미"}

- 첫째 항 $b\,\hbar(X_{t+1}, t)$: **위치 의존적 drift**. 높은 노드($\hbar > 0$)에 눈이 쌓이면 $\Phi$ 증가 (높이차 더 벌어짐), 낮은 노드($\hbar < 0$)에 쌓이면 $\Phi$ 감소 (균형 회복).
- 둘째 항 $(n-1)b^2/(3n)$: **항상 양수인 noise**. 노드 선택과 무관하게 매 step마다 $\Phi$를 살짝 키움 (분산 효과).

drift = "체계적 변화" + "noise", 그 합의 sign은 $\hbar(X_{t+1})$에 의해 결정.

:::

::: {.callout-note collapse="true" title="증명: $\\Phi$ drift lemma"}

$\Phi(t)$를 $X_{t+1}$ 기준으로 split + Fact 1을 $\Phi(t+1)$에 대입:

$$\Phi(t) = \hbar(X_{t+1}, t)^2 + \sum_{u \neq X_{t+1}} \hbar(u, t)^2$$

$$\Phi(t+1) = \left[\hbar(X_{t+1}, t) + b'_{t+1}\tfrac{n-1}{n}\right]^2 + \sum_{u \neq X_{t+1}}\left[\hbar(u, t) - \tfrac{b'_{t+1}}{n}\right]^2$$

차이를 $\text{Term}_1$ (recipient) + $\text{Term}_2$ (others)로 분해:

$$
\begin{aligned}
\text{Term}_1
&= \left[\hbar(X_{t+1}, t) + b'_{t+1}\tfrac{n-1}{n}\right]^2 - \hbar(X_{t+1}, t)^2 \\
&= 2\hbar(X_{t+1}, t)\cdot\tfrac{(n-1)b'_{t+1}}{n} + \tfrac{(n-1)^2(b'_{t+1})^2}{n^2} && (\because (A+B)^2 - A^2 = 2AB + B^2)
\end{aligned}
$$

$$
\begin{aligned}
\text{Term}_2
&= \sum_{u \neq X_{t+1}} \left[-\tfrac{2b'_{t+1}}{n}\hbar(u, t) + \tfrac{(b'_{t+1})^2}{n^2}\right] && (\because (A+B)^2 - A^2,\ B = -b'/n) \\
&= -\tfrac{2b'_{t+1}}{n}\sum_{u \neq X_{t+1}} \hbar(u, t) + \tfrac{(n-1)(b'_{t+1})^2}{n^2} \\
&= -\tfrac{2b'_{t+1}}{n}\bigl(-\hbar(X_{t+1}, t)\bigr) + \tfrac{(n-1)(b'_{t+1})^2}{n^2} && (\because \text{Fact 2}: \textstyle\sum \hbar = 0) \\
&= \tfrac{2b'_{t+1}}{n}\hbar(X_{t+1}, t) + \tfrac{(n-1)(b'_{t+1})^2}{n^2}
\end{aligned}
$$

합치고 정리 ($\tfrac{n-1}{n} + \tfrac{1}{n} = 1$, $\tfrac{(n-1)^2 + (n-1)}{n^2} = \tfrac{n-1}{n}$):

$$\Phi(t+1) - \Phi(t) = 2b'_{t+1}\,\hbar(X_{t+1}, t) + \tfrac{n-1}{n}(b'_{t+1})^2$$

기대값 ($\mathbb{E}[b'] = b/2$, $\mathbb{E}[(b')^2] = b^2/3$):

$$\mathbb{E}\bigl[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t, X_{t+1}\bigr] = 2 \cdot \tfrac{b}{2} \cdot \hbar(X_{t+1}, t) + \tfrac{n-1}{n} \cdot \tfrac{b^2}{3} = b\,\hbar(X_{t+1}, t) + \tfrac{(n-1)b^2}{3n} \qquad \square$$

:::

## §3.3 Shape decomposition + Round-level $\Phi$ drift

§3.2는 **한 step**의 $\Phi$ 변화. 라운드 전체로 합산하려면 라운드 중간 시점의 $\hbar$를 라운드 시작 시점의 $\hbar$로 표현해야. 그것이 Shape decomposition.

**기호.**

- $t_r$: $r$번째 block-flag time (라운드 시작)
- $m_r := t_{r+1} - t_r \in [1, T_{\max}+1]$: 라운드 길이
- $\tilde X_j := X_{t_r + j}$: 라운드 $r$의 $j$번째 방문 노드
- $\#_u := |\{j : \tilde X_j = u\}|$: 라운드 $r$ 동안 노드 $u$를 방문한 횟수
- $\pi(u, S) := \mathbb{E}[\#_u \mid S_{t_r}]$: 기대 방문 횟수
- $M := \max_i |\hbar(v_i, t_r)|$: 단측 centered height range

**Lemma (Shape decomposition).**

$$\sum_{j=1}^{m_r} \hbar(\tilde X_j,\, t_r + j - 1) = \sum_{u \in V} \#_u\,\hbar(u, t_r) + \Delta\text{Shape}^{(r)}$$

$$|\Delta\text{Shape}^{(r)}| \leq C_{\Delta\text{Shape}} := \tfrac{(T_{\max}+1)\,T_{\max}}{2}\,b$$

좌변: 라운드 중간 시점의 $\hbar$ 합. 우변 첫째 항: 라운드 시작 시점의 $\hbar$로 표현. $\Delta\text{Shape}$: 라운드 도중 높이가 변해서 생기는 오차 (bounded, $M$ 무관).

::: {.callout-note collapse="true" title="증명: Shape decomposition"}

라운드 중간 시점 $\hbar(\tilde X_j, t_r + j - 1)$를 라운드 시작 시점 $\hbar(\tilde X_j, t_r)$ + 보정으로 표현. Fact 1을 $j - 1$번 반복 적용:

$$\hbar(\tilde X_j, t_r + j - 1) = \hbar(\tilde X_j, t_r) + \sum_{k=1}^{j-1} b'_{t_r+k}\!\left[\mathbf{1}\{\tilde X_k = \tilde X_j\} - \tfrac{1}{n}\right]$$

($\tilde X_k$ 자리에 떨어진 step $k$의 눈은: 같은 노드 $\tilde X_j$이면 $+b'(n-1)/n$ = $+b'(1 - 1/n)$, 다른 노드이면 $-b'/n$. 즉 $b'[\mathbf{1}\{\tilde X_k = \tilde X_j\} - 1/n]$ 만큼 변화.)

$j = 1, \ldots, m_r$에 대해 합산:

$$\sum_{j=1}^{m_r} \hbar(\tilde X_j, t_r + j - 1) = \underbrace{\sum_{j=1}^{m_r} \hbar(\tilde X_j, t_r)}_{=:\, \text{Term}_1} + \underbrace{\sum_{j=1}^{m_r}\sum_{k=1}^{j-1} b'_{t_r+k}\!\left[\mathbf{1}\{\tilde X_k = \tilde X_j\} - \tfrac{1}{n}\right]}_{=:\, \text{Term}_2 \,=:\, \Delta\text{Shape}^{(r)}}$$

**Term 1** (visit count로 재인덱싱):

$$\sum_{j=1}^{m_r} \hbar(\tilde X_j, t_r) = \sum_{u \in V} \#_u\,\hbar(u, t_r)$$

**Term 2 = $\Delta\text{Shape}^{(r)}$ bound.** 각 summand 크기:

$$|b'_{t_r+k}[\mathbf{1}\{\cdots\} - 1/n]| \leq b \cdot \tfrac{n-1}{n} \leq b$$

이중합 항 개수: $\sum_{j=2}^{m_r}(j-1) = \tfrac{m_r(m_r-1)}{2} \leq \tfrac{(T_{\max}+1)T_{\max}}{2}$. 따라서:

$$|\Delta\text{Shape}^{(r)}| \leq \tfrac{(T_{\max}+1)T_{\max}}{2}\,b \qquad \square$$

:::

**Lemma (Round-level $\Phi$ drift).** $B(S) := \sum_{u \in V} \pi(u, S)\,\hbar(u, t_r)$ 정의. 라운드 단위 $\Phi$ drift:

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S_{t_r}] = b\,B(S_{t_r}) + C^{(r)}_{\text{rem}}$$

remainder $C^{(r)}_{\text{rem}}$ ($M$ 무관):

$$|C^{(r)}_{\text{rem}}| \leq C_{M\text{-free}} := (T_{\max}+1)\,b^2\!\left[\tfrac{T_{\max}}{2} + \tfrac{n-1}{3n}\right]$$

→ **$M$에 비례하는 유일한 항은 $b\,B(S)$**. 이것이 §3.4 DT 조건의 표적.

::: {.callout-note collapse="true" title="증명: Round-level $\\Phi$ drift"}

$$
\begin{aligned}
\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S_{t_r}]
&= \sum_{j=1}^{m_r} \mathbb{E}\!\left[b\,\hbar(\tilde X_j, t_r{+}j{-}1) + \tfrac{(n-1)b^2}{3n} \,\Big|\, S_{t_r}\right] && (\because \text{§3.2}) \\
&= b\,\mathbb{E}\!\left[\sum_{j=1}^{m_r} \hbar(\tilde X_j, t_r{+}j{-}1) \,\Big|\, S_{t_r}\right] + \tfrac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S_{t_r}] \\
&= b\,\mathbb{E}\!\left[\sum_u \#_u\,\hbar(u, t_r) + \Delta\text{Shape}^{(r)} \,\Big|\, S_{t_r}\right] + \tfrac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S_{t_r}] && (\because \text{Shape decomp.}) \\
&= b\,B(S_{t_r}) + \underbrace{b\,\mathbb{E}[\Delta\text{Shape}^{(r)} \mid S_{t_r}] + \tfrac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S_{t_r}]}_{=:\, C^{(r)}_{\text{rem}}} && (\because \pi(u, S) := \mathbb{E}[\#_u \mid S],\ B := \textstyle\sum \pi\,\hbar)
\end{aligned}
$$

remainder bound:

$$
\begin{aligned}
|C^{(r)}_{\text{rem}}|
&\leq b \cdot C_{\Delta\text{Shape}} + \tfrac{(T_{\max}+1)(n-1)b^2}{3n} && (\because |\Delta\text{Shape}|,\ m_r \text{ bounded}) \\
&= (T_{\max}+1)\,b^2\!\left[\tfrac{T_{\max}}{2} + \tfrac{n-1}{3n}\right] =: C_{M\text{-free}} \qquad \square
\end{aligned}
$$

:::

## §3.4 Drift-tightness (DT) condition

Foster–Lyapunov 정리가 작동하려면 **$M$이 클 때 round drift가 강하게 음수** 라는 조건이 필요. ABC증명의 명명: **Drift-tightness (DT)**.

::: {.callout-note collapse="false" title="Assumption (DT, Drift-tightness $\\mathrm{DT}(\\varepsilon, C_{\\text{tight}})$)"}

$(\mathcal{G}, \boldsymbol{\mu}_0, b, T_{\max})$에 대해 어떤 $\varepsilon \in (0, \infty)$ 와 $C_{\text{tight}} \in [0, \infty)$ 이 존재해, 모든 $S \in \mathcal{X}^*$에 대해

$$\mathbb{E}\bigl[\Phi(t_{r+1}) - \Phi(t_r) \mid S_{t_r} = S\bigr] \leq -\,\varepsilon\,M(S) + C_{\text{tight}} \tag{DT}$$

이름의 유래: 이 부등식이 chain $\{S_t\}$의 tightness (Meyn–Tweedie Thm 14.0.1) 를 끌어내므로 **"drift가 tight를 강제한다"** 는 의미.

:::

### $B(S)$로 환원

§3.3 round-level drift 식 $\mathbb{E}[\Delta\Phi] = b\,B(S) + C^{(r)}_{\text{rem}}$ 에서 $|C^{(r)}_{\text{rem}}| \leq C_{M\text{-free}}$. 따라서 DT는 다음 $B(S)$ 부등식과 동치:

$$\boxed{\ B(S) \leq -\varepsilon_{\text{OB}}\,M(S) + C_{\text{OB}}\quad \forall S \in \mathcal{X}^*\ } \tag{B-form}$$

상수 변환: $\varepsilon = b\,\varepsilon_{\text{OB}}$, $C_{\text{tight}} = b\,C_{\text{OB}} + C_{M\text{-free}}$.

→ 구체 그래프에서 DT를 검증할 때는 (B-form) 이 더 편리. 본 글에서 OB 검증은 모두 (B-form) 기준.

### $B(S)$의 의미

::: {.callout-important collapse="true" title="해석: $B(S)$가 말하는 것"}

$B(S) = \sum_u \pi(u, S)\,\hbar(u, S)$ — walker가 라운드 동안 방문하는 노드들의 round-start centered height를 **방문 횟수로 가중**한 기댓값.

- $B(S) < 0$: walker가 평균보다 **낮은** 노드($\hbar < 0$)를 더 자주 방문 → 거기 눈 쌓여 지형이 **평탄해짐**, $\Phi$ 감소
- $B(S) > 0$: walker가 **높은** 노드 선호 → 불균형 심화, $\Phi$ 증가
- $B(S) = 0$: 가중 평균 상쇄, $\Phi$ 1차 drift 없음

**OB ($B \leq -\varepsilon_{\text{OB}} M$) = "walker가 systematically 낮은 곳으로 흐른다"** 는 정량적 조건. 지형이 울퉁불퉁할수록 ($M$ 클수록) drift가 더 강하게 음수. $\varepsilon_{\text{OB}}$는 그래프의 **자기교정 세기**.

:::

### Balanced 가정 하의 자동 단순화

가정 (A1) Balanced — $\rho_i = 1/n$ — 에서는 $\boldsymbol{\mu}_0 = $ uniform이라 fall step이 각 노드에 균등하게 $1/n$ 기여. 라운드를 (fall 1번) + (flow $m_r - 1$번) 으로 분해:

$$\pi(u, S) = \underbrace{\tfrac{1}{n}}_{\text{fall}} + \underbrace{\pi_{\text{flow}}(u, S)}_{\text{flow part}}, \qquad \pi_{\text{flow}}(u, S) := \mathbb{E}_S\!\left[\sum_{s=1}^{m_r-1} \mathbf{1}\{X_{t_r+s} = u\}\right]$$

$B(S)$에 대입 + Fact 2 ($\sum \hbar = 0$) 사용:

$$B(S) = \tfrac{1}{n}\underbrace{\sum_u \hbar(u, S)}_{=\, 0} + \sum_u \pi_{\text{flow}}(u, S)\,\hbar(u, S) = \sum_u \pi_{\text{flow}}(u, S)\,\hbar(u, S)$$

→ **$B(S)$의 strict negativity는 전적으로 flow step에서 나온다.** Fall은 자동으로 0 기여.

**검증 절차** (balanced regime):

1. $\pi_{\text{flow}}(u, S)$ 계산 — round $S$에서 시작해 노드 $u$를 방문하는 기대 flow 횟수
2. $B(S) = \sum_u \pi_{\text{flow}}(u, S)\,\hbar(u, S)$ 형성
3. $B(S) \leq -\varepsilon_{\text{OB}}\,M(S) + C_{\text{OB}}$ 를 만족하는 $(\varepsilon_{\text{OB}}, C_{\text{OB}})$ 찾기

### Example: $K_3$, $T_{\max} = 5$, $b = 1$

$n = 3$, $\boldsymbol{\mu}_0 = (1/3, 1/3, 1/3)$, $b = 1$, $T_{\max} = 5$. 모든 round-start $S$에 대해 valid한 $(\varepsilon_{\text{OB}}, C_{\text{OB}})$ 유도.

::: {.callout-note collapse="true" title="증명: $K_3$ OB 상수"}

**Well-separated regime**: 모든 인접 rank gap $> b = 1$, 즉 $h_1 - h_2 > 1$, $h_2 - h_3 > 1$.

**Step 1 ($\pi_{\text{flow}}$).** 각 fall (확률 $1/3$):

- Fall $v_1$: $h_1 \to h_1 + 1$. $v_1$의 가장 낮은 이웃 $v_3$. Flow $v_3$ ($+1$). flow 후 $v_3$ 높이 $h_3 + 1$. 그 가장 낮은 이웃 $v_2$, 그런데 $h_2 > h_3 + 1$ (well-sep). **Block**. $m_r = 2$.
- Fall $v_2$: 같은 논리 — flow $v_3$ ($+1$), block. $m_r = 2$.
- Fall $v_3$: $h_3 + 1$, 가장 낮은 이웃 $v_2$인데 $h_2 > h_3 + 1$. **즉시 block**. $m_r = 1$.

$T_{\max} = 5$ 비활성 ($m_r \leq 2$). flow 방문 횟수:

$$\pi_{\text{flow}}(v_3) = \tfrac{2}{3},\quad \pi_{\text{flow}}(v_1) = \pi_{\text{flow}}(v_2) = 0$$

**Step 2 ($B(S)$).** Balanced 자동 단순화:

$$B(S) = \pi_{\text{flow}}(v_3)\,\hbar_3 = \tfrac{2}{3}\,\hbar_3$$

**Step 3 ($M$-linearity).** $\sum \hbar = 0$ + $h_1 = M$ 정의에서 worst case $\hbar = (M, -M/2, -M/2)$ ($\hbar_3 = -M/2$). 따라서:

$$B(S) \leq \tfrac{2}{3} \cdot (-M/2) = -\tfrac{M}{3}$$

→ $\varepsilon_{\text{OB}} = 1/3$.

**확인**: $(h_1, h_2, h_3) = (4, 2, 0)$ → $\bar h = 2$, $\hbar = (2, 0, -2)$, $M = 2$. $B = (2/3)(-2) = -4/3$. $B/M = -2/3$ ≥ $-1/3 \cdot 1$. OB 만족 (margin 있음, 이 state는 worst case 아님).

**Small-$M$ regime**: 어떤 gap $\leq 1$이면 $M \leq 2$. $|B(S)| \leq m_r \cdot M \leq 6 \cdot 2 = 12$. $C_{\text{OB}} = 12$.

**결론**: $\text{OB}(\varepsilon_{\text{OB}}, C_{\text{OB}}) = (1/3,\, 12)$ on $K_3$, $T_{\max} = 5$, $b = 1$, balanced. $\square$

:::

### 통계적 검증 (분석이 어려울 때)

::: {.callout-tip collapse="true" title="보충: DT의 통계적 검증 절차"}

비정규 그래프나 vertex-transitive 아닌 그래프에서는 $\pi_{\text{flow}}$의 분석적 계산이 어렵다. 이 경우 시뮬레이션으로 검증.

**핵심 관찰**: $\mathbb{E}[\hbar(X_t, t-1) \mid \mathcal{F}_{t-1}] = B(S_{t-1})$. 즉 매 step의 적립 위치 $\hbar(X_t, t-1)$ 가 $B$의 unbiased 단일 step 관측.

**절차**:
1. HST를 $\tau$ step 돌리고 매 step $R_t := \hbar(X_t, t-1)/M(t-1)$ 기록 ($M(t-1) > c$ 인 step만)
2. 배치 평균 $K$개 형성, $t$-test → $\hat\kappa := -\bar R$의 CI
3. DT 성립 시: CI가 $0$ 위에 strictly + $\hat\kappa \approx \varepsilon_{\text{OB}}$. 실패 시: CI가 $0$ 포함하거나 아래.

**post-hoc 가능**: 이미 돌린 시뮬의 tail (last $\tau/2$) 만 써도 됨. 추가 실행 불필요.

**권장 hyperparameters**: $\tau \geq 10^6$, $c = \text{median}(M_t)$, $K = 50$, batch size $\geq 10n$.

상세: ABC증명의 `ob_statistical_test.tex` 참조.

:::

## §3.5 Moment bound via iterated Lyapunov function

**Proposition (Moment bound).** $\mathbb{E}_{\pi^*}[\Phi] < \infty$.

Round drift $\mathbb{E}[\Delta\Phi] \leq -\varepsilon M + C_{\text{tight}}$에서 $M \sim \sqrt{\Phi}$ 이므로 drift 차수가 $-\Phi^{1/2}$. Meyn–Tweedie 14.0.1로는 $\mathbb{E}_{\pi^*}[\Phi^r] < \infty$ ($r < 1/2$) 까지만. **해법**: Lyapunov를 $V_{\text{Ly}} := \Phi^2$로 올리고 14.3.7 적용.

::: {.callout-note collapse="true" title="증명: Moment bound"}

$V_{\text{Ly}} := \Phi^2$, $\Delta\Phi_r := \Phi(t_{r+1}) - \Phi(t_r)$로 둠. 우선 전개:

$$V_{\text{Ly}}(t_{r+1}) = (\Phi(t_r) + \Delta\Phi_r)^2 = \Phi(t_r)^2 + 2\Phi(t_r)\Delta\Phi_r + (\Delta\Phi_r)^2$$

$$\Delta V_{\text{Ly}} = 2\Phi(t_r) \cdot \Delta\Phi_r + (\Delta\Phi_r)^2$$

기대값:

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S_{t_r}] = \underbrace{2\Phi \cdot \mathbb{E}[\Delta\Phi_r \mid S_{t_r}]}_{=:\,\text{Term1}} + \underbrace{\mathbb{E}[(\Delta\Phi_r)^2 \mid S_{t_r}]}_{=:\,\text{Term2}}$$

**Fact 3.** $\Phi \geq M^2/2$, $\Phi \leq nM^2$.

$$
\begin{aligned}
\Phi = \sum_i \hbar_i^2
&\geq \hbar_{\max}^2 + \hbar_{\min}^2 \\
&\geq \tfrac{(\hbar_{\max} + |\hbar_{\min}|)^2}{2} = \tfrac{M^2}{2} && (\because (A^2 + B^2) \geq (A+B)^2/2,\ \textstyle\sum \hbar = 0) \\
\Phi
&\leq nM^2 && (\because n\text{개 항 각각 } \leq M^2)
\end{aligned}
$$

**Term1.**

$$
\begin{aligned}
2\Phi \cdot \mathbb{E}[\Delta\Phi_r \mid S_{t_r}]
&= 2\Phi \sum_{j=1}^{m_r} \mathbb{E}\!\left[b\,\hbar(\tilde X_j, t_r{+}j{-}1) + \tfrac{(n-1)b^2}{3n} \,\Big|\, S_{t_r}\right] && (\because \text{§3.2}) \\
&= 2\Phi\left\{ b\,\mathbb{E}\!\left[\sum_{j=1}^{m_r}\hbar(\tilde X_j, t_r{+}j{-}1) \,\Big|\, S_{t_r}\right] + \tfrac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S_{t_r}] \right\} \\
&= 2\Phi\left\{ b\,\mathbb{E}\!\left[\sum_u \#_u\,\hbar(u, t_r) + \Delta\text{Shape}^{(r)} \,\Big|\, S_{t_r}\right] + \tfrac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S_{t_r}] \right\} && (\because \text{Shape decomp.}) \\
&= 2\Phi\left\{ b\,B(S_{t_r}) + C^{(r)}_{\text{rem}} \right\} && (\because \pi(u, S) := \mathbb{E}[\#_u \mid S],\ B := \textstyle\sum \pi\,\hbar) \\
&\leq 2\Phi\left\{ b\,B(S_{t_r}) + C_{M\text{-free}} \right\} && (\because |C^{(r)}_{\text{rem}}| \leq C_{M\text{-free}}) \\
&\leq 2\Phi(-b\,\varepsilon_{\text{OB}} M + b\,C_{\text{OB}} + C_{M\text{-free}}) && (\because \text{DT (B-form, §3.4)}) \\
&= 2\Phi(-\varepsilon M + C_{\text{tight}}) && (\varepsilon := b\,\varepsilon_{\text{OB}},\ C_{\text{tight}} := b\,C_{\text{OB}} + C_{M\text{-free}}) \\
&= -2\varepsilon\,\Phi M + 2C_{\text{tight}}\,\Phi \\
&\leq -\varepsilon M^3 + 2C_{\text{tight}}\,nM^2 && (\because \text{Fact 3}: \Phi \geq M^2/2,\ \Phi \leq nM^2)
\end{aligned}
$$

$$\therefore\ \text{Term1} \leq -\varepsilon M^3 + 2C_{\text{tight}}\,nM^2$$

**Term2.** 한 step의 $\Phi$ 변화 bound:

$$
\begin{aligned}
|\Delta\Phi \text{ in 1 step}|
&\leq 2b(M + T_{\max}b) + \tfrac{(n-1)b^2}{n} && (\because \beta \leq b,\ |\hbar(t_r{+}j)| \leq M + T_{\max}b) \\
|\Delta\Phi_r|
&\leq (T_{\max}{+}1)\!\left[2b(M + T_{\max}b) + \tfrac{(n-1)b^2}{n}\right] =: A_1 M + A_2 && (\because m_r \leq T_{\max}{+}1) \\
\mathbb{E}[(\Delta\Phi_r)^2 \mid S_{t_r}]
&\leq (A_1 M + A_2)^2 \leq 2A_1^2 M^2 + 2A_2^2 =: C_3 M^2 + C_4 && (\because (a+b)^2 \leq 2a^2 + 2b^2)
\end{aligned}
$$

상수: $A_1 := 2b(T_{\max}+1)$, $A_2 := (T_{\max}+1)[2bT_{\max}b + (n-1)b^2/n]$, $C_3 := 2A_1^2$, $C_4 := 2A_2^2$.

$$\therefore\ \text{Term2} \leq C_3 M^2 + C_4$$

**합치기.**

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S_{t_r}] \leq -\varepsilon M^3 + (2C_{\text{tight}}\,n + C_3)\,M^2 + C_4 = -\varepsilon M^3 + C_5 M^2 + C_6$$

상수: $C_5 := 2C_{\text{tight}}\,n + C_3$, $C_6 := C_4$. $M \geq M_0 := 2C_5/\varepsilon$이면 $M^3$ 지배:

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S_{t_r}] \leq -\tfrac{\varepsilon}{2}\,M^3 + C_6$$

$V$로 표현: $\Phi \leq nM^2$ 에서 $M \geq (\Phi/n)^{1/2}$ 이므로 $M^3 \geq \Phi^{3/2}/n^{3/2} = V_{\text{Ly}}^{3/4}/n^{3/2}$.

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S_{t_r}] \leq -\tfrac{\varepsilon}{2n^{3/2}}\,V_{\text{Ly}}^{3/4} + C_6 \cdot \mathbf{1}_{\{M \leq M_0\}}$$

Meyn–Tweedie Thm 14.3.7: $f := \tfrac{\varepsilon}{2n^{3/2}}\,V_{\text{Ly}}^{3/4}$이 $\pi^*$-integrable.

$$\mathbb{E}_{\pi^*}[V_{\text{Ly}}^{3/4}] = \mathbb{E}_{\pi^*}[\Phi^{3/2}] < \infty$$

$\Phi \leq \Phi^{3/2} + 1$ 이므로:

$$\mathbb{E}_{\pi^*}[\Phi] \leq \mathbb{E}_{\pi^*}[\Phi^{3/2}] + 1 < \infty \qquad \square$$

:::

## §3.6 Doeblin minorization

**Lemma (Pushforward density).** $\xi_1, \ldots, \xi_n \stackrel{\text{iid}}{\sim} \text{Unif}(0, b)$, $\Lambda(\xi) := (\xi_2 - \xi_1, \ldots, \xi_n - \xi_1)$. Pushforward 밀도 $\rho_\Lambda(\mathbf{y}) \geq b^{1-n}/2$ on $\mathcal{O} := \{\max_i |y_i| < b/4\}$.

::: {.callout-note collapse="true" title="증명: Pushforward density lemma"}

fibre $\Lambda^{-1}(\mathbf{y}) \cap [0, b]^n$의 길이: $\min(b, b - \max y_i) - \max(0, -\min y_i) > b/2$ on $\mathcal{O}$. 즉 $\Lambda$의 fiber가 길이 $> b/2$인 구간이라 pushforward density $\geq 1/b \cdot 1/b^{n-1} \cdot (b/2)/1 \cdot 1/1 = b^{1-n}/2$. $\square$

:::

**Proposition (Doeblin minorization).** 임의의 bounded set $\mathcal{B}_\Phi \subset \mathcal{X}^*$ 에 대해, 어떤 $\eta > 0$, $K^* \in \mathbb{N}$, $\mathcal{X}^*$ 위의 확률측도 $\nu$가 존재하여

$$P^{K^*}(s, B) \geq \eta\,\nu(B), \qquad \forall\,s \in \mathcal{B}_\Phi,\ \ \forall\,B \in \mathcal{B}(\mathcal{X}^*)$$

::: {.callout-important collapse="true" title="해석: 한 줄 직관 — '공통 운명'"}

Doeblin minorization을 한 줄로 요약하면: **"어디서 출발하든 결국 겹치게 되는 최소한의 공통 운명이 존재한다."**

**비유**: 전국 어디서 출발하든 $K^*$번 이동하면, 최소 $\eta$ (예: 10%) 확률로 모두 같은 목적지 분포 $\nu$ (예: '서울역' 주변)에 떨어지게 만드는 **마법의 규칙**이 숨어있다.

**$P^{K^*}(s, B) \geq \eta\,\nu(B)$가 말하는 것**:

- **독립적인 바닥 (lower bound)**: 출발점 $s \in \mathcal{B}_\Phi$가 아무리 극단적이어도, 도착 분포의 '바닥'을 공통 측도 $\nu$가 받쳐줌. **출발지에 무관한 양의 공통 부분**이 존재.
- **과거 망각 (coupling)**: 이 공통 부분 덕분에 두 다른 궤적이 양의 확률로 **같은 $\nu$ 표본으로 coupling** → 그 시점부터 두 궤적은 같은 분포. 초기 상태의 기억 리셋.
- **결과 (uniform ergodicity)**: coupling 강제로 모든 궤적이 결국 **같은 정상 분포 $\pi^*$로 수렴**. 출발 무관한 균등 에르고딕성.

이게 HST 증명에서 **"유일한 정상 분포 $\pi^*$가 존재하고, 그 분포가 초기 신호 $\mathbf{y}$에 의존하지 않는다"** 를 보이는 핵심 무기.

:::

::: {.callout-important collapse="true" title="해석: 증명의 큰 그림"}

**문제**: 출발 상태 $\boldsymbol{\delta}_0$ 에서 $K^*$ step 후 도달 분포가 모든 target $\boldsymbol{\delta}^* \in \mathcal{B}_\Phi$를 양의 밀도로 cover하는 걸 보여야.

**Randomness 원천**: HST에는 두 가지.

- 이산 randomness: block 직후 Fall에서 노드 선택 $X_s \sim \mu_0$ (i.i.d.)
- **연속 randomness**: 눈 증분 $b'_s \sim \text{Unif}(0, b)$ (i.i.d.) ← Lebesgue 밀도 원천

연속 밀도는 오로지 **block 직후 Fall에서 떨어지는 첫눈** 에서 (flow step의 $b'_s$는 결정론적 노드 이동과 합쳐져 noise처럼 작동).

**전략 4단계**:

1. **충분히 긴 window**: $K = n(T_{\max}{+}1)$ step. round 길이 $\leq T_{\max}{+}1$이므로 이 안에 block이 최소 $n$번.
2. **모든 노드에 첫눈이 한 번씩 떨어지는 사건 $\mathcal{A}_1$**: 그 $n$번의 block 직후 첫눈이 정확히 $v_1, \ldots, v_n$ 순으로. 모든 노드에 fresh $\text{Unif}(0, b)$ 한 번씩 → $n$차원 shift 가능.
3. **Flow step 눈은 작게 ($\mathcal{A}_2$)**: window 안 flow step의 $b'_s \leq \varepsilon$ conditioning. noise만 줄이는 것이고 block 첫눈 분포는 안 건드림 (independence). $L$ window 누적해도 residual $< Lb/8$.
4. **한 window는 작은 ball, $L$ window는 큰 ball**: 한 window 첫눈들이 만드는 shift은 small ball $\mathcal{O} = \{\max|y_i| < b/4\}$ 안. $\mathcal{B}_\Phi$ 전체 ($\|\boldsymbol{\delta}\|_\infty \leq D_\Phi$)를 cover하려면 $L \cdot b/4 > 2D_\Phi$ → $L := \lceil 16 D_\Phi/b + 1 \rceil$.

**그래서 상수**: $K^* = LK$, $\varepsilon = b/(8K)$.

:::

::: {.callout-note collapse="true" title="증명: Doeblin minorization"}

**상수 정의**.

$$K := n(T_{\max}{+}1),\quad \varepsilon := \tfrac{b}{8K},\quad D_\Phi := \sup_{s \in \mathcal{B}_\Phi} \|\boldsymbol{\delta}(s)\|_\infty,\quad L := \lceil 16 D_\Phi / b + 1 \rceil,\quad K^* := LK$$

**한 window의 핵심 사건 두 개**.

- $\mathcal{A}_1$ — 길이 $K$ window 안 첫 $n$번의 block 직후 첫눈이 차례로 $v_1, v_2, \ldots, v_n$에 떨어짐 (모든 노드에 한 번씩 첫눈)
- $\mathcal{A}_2$ — 같은 window의 모든 flow step 눈 증분이 $b'_s \leq \varepsilon$ (noise 억제)

**$P(\mathcal{A}_1) > 0$ 증명**. round 길이 $\leq T_{\max}+1$이므로 $K = n(T_{\max}+1)$ step 안에 block 시점 $\geq n$개 보장. 각 block 직후 첫눈의 노드 선택은 $\mu_0$ i.i.d., $\mathcal{F}_{s-1}$ 독립.

$$P(\mathcal{A}_1) \geq \mu_{\min}^n > 0 \qquad (\mu_{\min} := \min_i \mu_0(v_i))$$

**$P(\mathcal{A}_2) > 0$ + 첫눈 분포 유지**. block/flow 분류는 timer $Z_{s-1}$로 $\mathcal{F}_{s-1}$-measurable, $b'_s \perp \mathcal{F}_{s-1}$. flow step에 $\{b'_s \leq \varepsilon\}$ 조건 걸어도 **block 첫눈의 $b'_s$ marginal은 $\text{Unif}(0, b)$ 그대로**.

$$P(\mathcal{A}_2) \geq (\varepsilon/b)^{n T_{\max}} > 0$$

**Flow step residual bound** ($L$ window 누적).

$$\|\mathbf{R}\|_\infty \leq \underbrace{L n T_{\max}}_{\text{flow step 수}} \cdot \underbrace{\varepsilon}_{= b/(8K)} = L n T_{\max} \cdot \tfrac{b}{8n(T_{\max}+1)} < \tfrac{Lb}{8}$$

**한 window의 첫눈 shift 분포**. $\mathcal{A}_1$ 하 한 window 첫눈 $n$개는 i.i.d. $\text{Unif}(0, b)$. Pushforward density lemma로 $\boldsymbol{\delta}$-shift이 $\mathcal{O} = \{\max_i |y_i| < b/4\}$ 위 밀도 $\geq b^{1-n}/2$.

**$L$ window 합성으로 target 도달**. $L$개 독립 window shift 합은 Minkowski sum으로 $L\mathcal{O}$ 위 밀도 (convolution).

도달해야 할 net shift: $\boldsymbol{\delta}^* - \boldsymbol{\delta}_0 - \mathbf{R}$.

$$\|\boldsymbol{\delta}^* - \boldsymbol{\delta}_0 - \mathbf{R}\|_\infty \leq \underbrace{\|\boldsymbol{\delta}^*\|_\infty}_{\leq D_\Phi} + \underbrace{\|\boldsymbol{\delta}_0\|_\infty}_{\leq D_\Phi} + \underbrace{\|\mathbf{R}\|_\infty}_{< Lb/8} \leq 2D_\Phi + \tfrac{Lb}{8}$$

$L \geq 16D_\Phi/b + 1$이면:

$$2D_\Phi + \tfrac{Lb}{8} < \tfrac{Lb}{8} + \tfrac{Lb}{8} = \tfrac{Lb}{4}$$

→ shift이 $L\mathcal{O}$ 내부 ball 안. block shift 분포로 cover 가능.

**결론**. $\nu := L\mathcal{O}$ 내 작은 ball의 normalized Lebesgue.

$$P^{K^*}(s_0, B) \geq \underbrace{[P(\mathcal{A}_1 \cap \mathcal{A}_2)]^L}_{>\,0} \cdot \underbrace{(b^{1-n}/2)^L}_{\text{밀도 lower bound}^L} \cdot \nu(B) =: \eta\,\nu(B)$$

$\eta > 0$은 $s_0$ 무관. $\square$

:::

## §3.7 Theorem A: statement & proof

::: {.callout-note collapse="false" title="Theorem A (Balanced-regime convergence)"}

다음 가정 하에:
- (Standing) 연결 그래프 $\mathcal{G}$, full-support $\boldsymbol{\mu}_0$ ($\mu_{\min} > 0$);
- **(A1) Balanced**: $\rho_i = 1/n$ for all $i$;
- **(A2) Drift-tightness**: $(\mathcal{G}, \boldsymbol{\mu}_0, b, T_{\max})$가 $\text{DT}(\varepsilon, C_{\text{tight}})$를 만족 ($\varepsilon > 0$)

deterministic 행렬 $C = [c_{ij}]$ ($c_{ij} \geq 0$, 그래프와 $b, T_{\max}$에만 의존)가 존재해

$$\frac{SD^2_{ij}(t)}{t} \xrightarrow{t \to \infty} c_{ij} \quad \text{a.s.}$$

극한 $c_{ij}$는 **초기 신호 $\mathbf{y}$에 무관**.

:::

::: {.callout-important collapse="true" title="해석: Foster-Lyapunov + Doeblin의 역할 분담"}

증명 체인이 왜 두 기계 (Foster drift, Doeblin) 를 **나란히** 쓰는지 — 둘은 **상호보완**.

**비유**: state space를 큰 그릇이라고 보자.

- **바깥** ($\mathcal{B}_\Phi$ 밖, $\Phi$ 큰 곳): 그릇 가장자리. **Foster drift가 안쪽으로 끌어당김** (Lyapunov 강하게 감소 → chain이 가장자리에 오래 못 머묾).
- **안쪽** ($\mathcal{B}_\Phi$ 안, $\Phi$ 작은 곳): 그릇 바닥. drift가 약함. 대신 **Doeblin minorization이 coupling 강제** (출발 무관 공통 분포로 수렴).

**왜 둘 다 필요한가**:

| 가진 것 | 빠진 것 |
|---|---|
| Foster drift 만 | chain이 $\mathcal{B}_\Phi$ 안에 들어와도 분포 결정 불가 → **유일성 보장 안 됨** |
| Doeblin 만 | $\mathcal{B}_\Phi$가 small set이라도 **거기 도달 보장 없음** → 발산 가능 |
| Foster + Doeblin | 가장자리 → 바닥 (Foster) → 모두 coupling (Doeblin) → 유일 $\pi^*$로 수렴 |

**표준 패턴**: Meyn-Tweedie (2009)의 핵심 도구. Thm 11.3.4 (positive Harris), Thm 14.3.7 (moment bound), Thm 17.3.2 (LLN) 모두 "Foster drift + small set" 조합 위에.

증명 체인 순서: **Step 1** (Foster drift) → **Step 2** (Doeblin) → **Step 3** (ψ-irreducibility) → **Step 4** (Skeleton → original) → **Step 5** (Moment bound) → **Step 6** (LLN).

:::

::: {.callout-note collapse="true" title="증명: Theorem A"}

**Step 1 (Round-skeleton Foster drift).** §3.3 round-level $\Phi$ drift + §3.4 DT (B-form):

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S_{t_r}] \leq -\varepsilon\,M(t_r) + C_{\text{tight}}$$

$\Phi \leq nM^2$ 에서 $M \geq \sqrt{\Phi/n}$. 위 우변이 $-\varepsilon\sqrt{\Phi/n} + C_{\text{tight}}$ 이하. $R := n\bigl((C_{\text{tight}} + 1)/\varepsilon\bigr)^2$, $\mathcal{B}_\Phi := \{S : \Phi(S) \leq R\}$로 두면:

$$S_{t_r} \notin \mathcal{B}_\Phi \implies \mathbb{E}[\Delta\Phi \mid S_{t_r}] \leq -1$$

**Step 2 (Doeblin).** §3.6 Doeblin minorization: $\mathcal{B}_\Phi$가 $K^*$-skeleton의 small set.

$$P^{K^*}(s, B) \geq \eta\,\nu(B) \qquad \forall s \in \mathcal{B}_\Phi$$

**Step 3 ($\psi$-irreducibility, 단일 closed Harris class).**

Foster drift (Step 1) 와 optional stopping: $\tau_{\mathcal{B}_\Phi} := \inf\{r : S_{t_r} \in \mathcal{B}_\Phi\}$ 정의. $Y_r := \Phi(t_{r \wedge \tau}) + r \wedge \tau$ 가 supermartingale ($\mathbb{E}[\Delta Y_r] \leq 0$ for $r < \tau$). 따라서 $\mathbb{E}[\tau] \leq \mathbb{E}[Y_0] = \Phi(t_0) < \infty$ → $\tau < \infty$ a.s. 즉 **모든 initial state가 $\mathcal{B}_\Phi$에 a.s. 도달**.

$\mathcal{B}_\Phi$ 안에서 Doeblin smallness로 $\nu$-irreducibility. $\nu$가 absolutely continuous (Lebesgue on a ball) 이므로 $\psi \ll \text{Leb}$. Meyn-Tweedie Prop 5.5.5: $\psi$-irreducible chain + small set → **단일 closed Harris class**. $\pi^*$ unique. $\nu$가 $\boldsymbol{\delta}$-공간에 정의되어 출발 $\boldsymbol{\delta}_0$ (= 초기 신호 $\mathbf{y}$) 와 무관 → **$\pi^*$가 $\mathbf{y}$ 무관**.

**Step 4 (Skeleton → original transfer).** Foster drift + Doeblin은 round-skeleton $\{S_{t_r}\}$에 대해. original chain $\{S_t\}$로 전달:

- $m_r \leq T_{\max} + 1$ — 라운드 길이 bounded
- 따라서 original chain의 $\tau_{\mathcal{B}_\Phi}$ ≤ $(T_{\max}+1) \cdot \tau^{\text{sk}}_{\mathcal{B}_\Phi}$
- $\mathbb{E}[\tau^{\text{sk}}] < \infty$ → $\mathbb{E}[\tau] < \infty$
- Meyn-Tweedie Thm 17.3.2: skeleton의 positive Harris recurrence + LLN이 original chain으로 전달

**Step 5 (Moment bound).** §3.5 Proposition: $\mathbb{E}_{\pi^*}[\Phi] < \infty$. 따라서:

$$g(S) := (\delta_i - \delta_j)^2 \leq 2\Phi \implies \mathbb{E}_{\pi^*}[g] < \infty$$

**Step 6 (LLN).** Positive Harris chain의 LLN (Meyn-Tweedie Thm 17.1.7):

$$\frac{1}{t}\sum_{s=0}^t g(S_s) \xrightarrow{t \to \infty} \mathbb{E}_{\pi^*}[g] =: c_{ij} \quad \text{a.s.}$$

$g(S_s) = (\delta_{i,s} - \delta_{j,s})^2 = (h(v_i, s) - h(v_j, s))^2$ (차이에서 $h(v_1, s)$ 소거). 좌변 = $SD^2_{ij}(t)/t$.

$c_{ij}$는 $\pi^*$의 유일성 (Step 3) 에서 $\mathbf{y}$ 무관. $\square$

:::

---

# §4 Theorem B: Intermediate regime

::: {.callout-note collapse="false" title="Theorem B (Intermediate-regime convergence)"}

다음 가정 하에:
- (Standing) 연결 그래프 $\mathcal{G}$, full-support $\boldsymbol{\mu}_0$;
- **(A3) Diffusive scaling**: $\rho_i, \rho_j$ 가 a.s. 존재하고, 높이차 $D_t := h(v_i, t) - h(v_j, t)$ 가
  $$\frac{D_t^2}{t} \xrightarrow{t \to \infty} \sigma^2_{ij} \in (0, \infty) \quad \text{a.s.}$$
  를 만족 (이는 $\rho_i = \rho_j$를 implies; $\rho_i \neq \rho_j$ 면 $D_t^2/t \to \infty$).

그러면:

$$\frac{SD^2_{ij}(t)}{t^2} \xrightarrow{t \to \infty} \frac{\sigma^2_{ij}}{2} \quad \text{a.s.}$$

:::

::: {.callout-important collapse="true" title="해석: Intermediate regime의 의미"}

조건 $D_t^2/t \to \sigma^2 > 0$ 은 **높이차가 diffusive 거동** ($|D_t| \sim \sigma\sqrt{t}$) 을 한다는 뜻. 비교:

- **Balanced** ($\rho_i = \rho_j$ + global balance): $D_t = O(1)$ stationary. → $\sum D_s^2 \sim t$
- **Intermediate** ($\rho_i = \rho_j$, global imbalance): $|D_t| \sim \sigma\sqrt{t}$ diffusive. → $\sum D_s^2 \sim t^2$ ← 본 정리
- **Drift** ($\rho_i \neq \rho_j$): $D_t \sim \gamma t$ ballistic. → $\sum D_s^2 \sim t^3$

**왜 diffusion이 나오는가**: $\rho_i = \rho_j$이지만 시스템 전체가 unbalanced인 경우. 예: Helm 그래프 leaf–leaf 쌍. hub의 linearly growing height가 flow를 통해 noise로 두 leaf에 전파. 두 leaf는 직접 연결 X → noise가 독립 → $D_t$가 random walk → step variance $\sim (b^2/3) \cdot 2\rho$, $\sigma^2 = 2\rho b^2/3$.

**Theorem B의 특이성**: recurrence 기계 (Foster–Lyapunov, Doeblin) 가 **필요 없다**. 가중 Cesàro lemma 한 줄로 끝. 가정 (A3) 자체가 "recurrence가 깨졌지만 diffusive하게 깨졌다" 는 conditional 정보.

:::

::: {.callout-note collapse="true" title="증명: Theorem B"}

$\rho_i = \rho_j$이므로 drift 항이 사라짐: $D_s = \bar b(\rho_i - \rho_j)s + \xi(s) = \xi(s)$ (즉 $D_s = \xi(s)$).

$a_s := D_s^2/s$로 둠 (for $s \geq 1$). 가정 (A3) 에서 $a_s \to \sigma^2_{ij}$ a.s.

$$\frac{SD^2_{ij}(t)}{t^2} = \frac{D_0^2}{t^2} + \frac{1}{t^2}\sum_{s=1}^t s \cdot a_s$$

첫째 항 $\to 0$. 둘째 항: **가중 Cesàro lemma** (Toeplitz lemma 변형) — 가중치 $w_s = s$, $\sum_{s=1}^t s = t(t+1)/2$. $a_s \to L$이면

$$\frac{\sum_{s=1}^t s\,a_s}{\sum_{s=1}^t s} \to L$$

따라서 $\sum_{s=1}^t s\,a_s \sim L \cdot t(t+1)/2 \sim L\,t^2/2$. $L = \sigma^2_{ij}$ 대입:

$$\frac{SD^2_{ij}(t)}{t^2} \xrightarrow{t \to \infty} \frac{\sigma^2_{ij}}{2} \quad \text{a.s.} \qquad \square$$

:::

::: {.callout-tip collapse="true" title="보충: $\\sigma^2_{ij}$가 양수일 조건"}

가정 (A3) 은 $\sigma^2_{ij}$가 **strictly positive** 라고 가정한다. 만약 $\sigma^2_{ij} = 0$이면 $D_t = o(\sqrt{t})$, 즉 diffusion보다 빠른 mean-reversion → balanced regime ($SD^2/t$ 수렴)으로 갈 가능성.

Helm 그래프 leaf–leaf에서 $\sigma^2 = 2\rho b^2/3 > 0$이 명시적으로 계산됨 ($\rho$ = leaf의 deposition rate). 일반 그래프에서 $\sigma^2$의 양성은 그래프 구조 분석 필요.

:::

---

# §5 Theorem C: Drift regime

::: {.callout-note collapse="false" title="Theorem C (Drift-regime convergence)"}

다음 가정 하에:
- (Standing) 연결 그래프 $\mathcal{G}$, full-support $\boldsymbol{\mu}_0$;
- **(A4) Drift**: $\rho_i, \rho_j$ 가 a.s. 존재하고 $\rho_i \neq \rho_j$.

그러면:

$$\frac{SD^2_{ij}(t)}{t^3} \xrightarrow{t \to \infty} \frac{\bar b^2(\rho_i - \rho_j)^2}{3} \quad \text{a.s.}$$

여기서 $\bar b := \mathbb{E}[b'] = b/2$. 극한은 $\rho_i \neq \rho_j$ 하에서 strictly positive.

:::

::: {.callout-important collapse="true" title="해석: Drift regime의 의미"}

$\rho_i \neq \rho_j$ 이면 노드 $v_i, v_j$의 평균 적립률이 다르다 → 높이차 $D_s = h(v_i, s) - h(v_j, s)$ 가 **선형으로 발산**:

$$D_s \approx \bar b(\rho_i - \rho_j)\,s + o(s)$$

따라서 $SD^2 = \sum D_s^2 \approx \gamma^2 \sum s^2 \sim \gamma^2 t^3/3$ ($\gamma := \bar b(\rho_i - \rho_j)$).

**증명 도구**: 마팅게일 SLLN + 급수 계산. Theorem A의 무거운 recurrence 기계 (Foster–Lyapunov, Doeblin) 불필요 — drift 항이 fluctuation을 압도하므로.

:::

::: {.callout-note collapse="true" title="증명: Theorem C"}

**Step 1 (높이 분해).**

$$h(v_i, t) = h(v_i, 0) + \sum_{s=1}^t b'_s\,\mathbf{1}\{X_s = v_i\}$$

알고리즘 순서상 $X_s$가 $S_{s-1}$에서 결정 (즉 $b'_s$ 추첨 **이전**), $b'_s \perp (\mathcal{F}_{s-1}, X_s)$. 분리:

$$\sum_{s=1}^t b'_s\,\mathbf{1}\{X_s = v_i\} = \bar b \sum_{s=1}^t \mathbf{1}\{X_s = v_i\} + \sum_{s=1}^t (b'_s - \bar b)\,\mathbf{1}\{X_s = v_i\}$$

**첫째 항**: $\rho_i$ 정의에서 $\frac{1}{t}\sum \mathbf{1}\{X_s = v_i\} \to \rho_i$ a.s. → 합은 $\bar b \rho_i\,t + o(t)$.

**둘째 항**: $M_t^{(i)} := \sum_{s=1}^t (b'_s - \bar b)\,\mathbf{1}\{X_s = v_i\}$ 가 $\mathcal{F}_t := \sigma(X_1, b'_1, \ldots, X_t, b'_t)$-마팅게일.

$$
\begin{aligned}
\mathbb{E}[\Delta M_s^{(i)} \mid \mathcal{F}_{s-1}]
&= \mathbb{E}\!\left[\mathbf{1}\{X_s = v_i\}(b'_s - \bar b) \,\Big|\, \mathcal{F}_{s-1}\right] \\
&= \mathbb{E}\!\left[\mathbf{1}\{X_s = v_i\}\,\underbrace{\mathbb{E}[b'_s - \bar b \mid \mathcal{F}_{s-1}, X_s]}_{=\, 0} \,\Big|\, \mathcal{F}_{s-1}\right] && (\because \text{tower property}) \\
&= 0 && (\because b'_s \perp (\mathcal{F}_{s-1}, X_s))
\end{aligned}
$$

증분 bounded ($|\Delta M_s^{(i)}| \leq b/2$) → 마팅게일 SLLN: $M_t^{(i)}/t \to 0$ a.s.

**Combining**: $h(v_i, t) = h(v_i, 0) + \bar b \rho_i\,t + o(t)$ a.s.

**Step 2 (높이차 점근).** $D(s) := h(v_i, s) - h(v_j, s)$, $\gamma := \bar b(\rho_i - \rho_j)$. Step 1에서:

$$D(s) = \gamma\,s + \xi(s), \qquad \xi(s) = o(s) \text{ a.s.}$$

즉 $\xi(s)/s \to 0$ a.s.

**Step 3 ($SD^2$ 전개).**

$$
\begin{aligned}
SD^2_{ij}(t)
&= \sum_{s=0}^t D(s)^2 = \sum_{s=0}^t \bigl(\gamma s + \xi(s)\bigr)^2 \\
&= \gamma^2 \sum_{s=0}^t s^2 + 2\gamma \sum_{s=0}^t s\,\xi(s) + \sum_{s=0}^t \xi(s)^2
\end{aligned}
$$

**주항**: $\sum_{s=0}^t s^2 = t(t+1)(2t+1)/6 = t^3/3 + O(t^2)$.

**교차항 = $o(t^3)$**. $\xi(s)/s \to 0$이므로 $\forall \epsilon > 0,\ \exists S_\epsilon < \infty$ (random a.s.) 에 대해 $|\xi(s)| \leq \epsilon s$ for $s \geq S_\epsilon$.

$$\left|\sum_{s=0}^t s\,\xi(s)\right| \leq \underbrace{\sum_{s=0}^{S_\epsilon} s\,|\xi(s)|}_{=:\, C_\epsilon} + \epsilon \sum_{s=S_\epsilon+1}^t s^2 \leq C_\epsilon + \tfrac{\epsilon\,t^3}{3}$$

따라서 $\limsup_t t^{-3}|\sum s\,\xi(s)| \leq \epsilon/3$ a.s. $\epsilon$ arbitrary → $\sum s\,\xi(s) = o(t^3)$ a.s.

**잔차항 = $o(t^3)$**. 같은 논법: $|\xi(s)|^2 \leq \epsilon^2 s^2$ for $s \geq S_\epsilon$ → $\sum \xi(s)^2 \leq C'_\epsilon + \epsilon^2 t^3/3 = o(t^3)$ a.s.

**결론**:

$$\frac{SD^2_{ij}(t)}{t^3} = \frac{\gamma^2}{3} + o(1) + o(1) \xrightarrow{t \to \infty} \frac{\gamma^2}{3} = \frac{\bar b^2(\rho_i - \rho_j)^2}{3} \qquad \square$$

:::

::: {.callout-tip collapse="true" title="보충: $\\rho_i$의 존재성과 recurrence 불필요성"}

**$\rho_i$ 존재성**: drift regime에서는 a.s. existence가 일반적으로 open problem이지만, 모든 시뮬에서 경험적 확인됨. 충분조건: walker process의 occupation measure에 대한 SLLN — random-step variant에서는 오큐리의 `occupation_slln_general.tex` 가 $\mu_{\min} > 0$ 인 모든 연결 그래프에서 무조건 성립을 증명.

**Recurrence 기계 불필요**: Theorem C는 마팅게일 SLLN + Toeplitz (가중 Cesàro) 만으로 충분. Theorem A의 난점은 정확히 dominant $O(t)$ drift가 없어서 $o(t)$ fluctuation $\xi$ 자체를 recurrence 기계로 통제해야 한다는 점.

:::

---

# §6 Synthesis

## §6.1 통합 분해 관점

세 정리를 하나로 묶는 관점: 높이차를

$$h(v_i, s) - h(v_j, s) = \underbrace{\bar b(\rho_i - \rho_j)\,s}_{\text{drift}} + \underbrace{\xi_{ij}(s)}_{\text{fluctuation}}, \qquad \xi_{ij}(s) = o(s) \text{ a.s.}$$

로 분해하면

$$SD^2_{ij}(t) = \bar b^2(\rho_i - \rho_j)^2 \sum s^2 + 2\bar b(\rho_i - \rho_j)\sum s\,\xi_{ij}(s) + \sum \xi_{ij}(s)^2$$

`-` **Drift regime** ($\rho_i \neq \rho_j$): 첫째 항 $\sim t^3$이 지배. Normalization: $SD^2/t^3$.

`-` **Balanced regime** ($\rho_i = \rho_j$, global balance): drift 소실. $SD^2 = \sum \xi^2$. Harris recurrence가 $\xi$를 stationary ($= O(1)$) 로 만들어 $\sum \xi^2/t \to c_{ij}$. Normalization: $SD^2/t$.

`-` **Intermediate regime** ($\rho_i = \rho_j$, global imbalance): drift 소실, 그러나 $\xi$가 **diffusive** ($\xi^2(s)/s \to \sigma^2 > 0$, global imbalance가 Harris recurrence 깬다). $\sum \xi^2/t^2 \to \sigma^2/2$. Normalization: $SD^2/t^2$.

세 경우 모두 극한이 의미 있는 그래프 구조를 포착:

| Regime | 극한이 반영 | 예시 |
|---|---|---|
| Balanced ($SD^2/t$) | flow/block stationary dynamics | $K_n$, uniform $\boldsymbol{\mu}_0$ |
| Intermediate ($SD^2/t^2$) | equal-$\rho$ 그룹 내 flow noise 확산 | Helm leaf–leaf, deg-prop $\boldsymbol{\mu}_0$ |
| Drift ($SD^2/t^3$) | deposition-rate gap $|\rho_i - \rho_j|$ | Star hub–leaf, deg-prop $\boldsymbol{\mu}_0$ |

## §6.2 $\boldsymbol{\mu}_0$ 선택 가이드

regime은 primarily $\boldsymbol{\mu}_0$가 결정.

**Uniform** $\boldsymbol{\mu}_0 = (1/n, \ldots, 1/n)$. 모든 노드 fall rate $1/n$. flow/block 동력학의 대칭성과 합쳐서 $\rho_i \approx 1/n$. → **balanced regime**, $SD^2/t \to c_{ij}$. finite-time 거동과 fixed-step formulation은 그래프-신호 상호작용을 미세하게 반영 (예: 신호 manifold recovery).

**Degree-proportional** $\boldsymbol{\mu}_0 \propto \deg$. high-degree 노드가 fall을 더 자주 받음: $\mu_0(v_i) = \deg(v_i)/\sum \deg$. 비정규 그래프에서 unequal $\rho_i$ → linear drift. → **drift regime**, $SD^2/t^3 \to \bar b^2(\rho_i - \rho_j)^2/3$. 극한이 degree heterogeneity 포착 (star: hub 분리).

**Intermediate** $\boldsymbol{\mu}_0 = (1 - \alpha)\,\text{unif} + \alpha\,\text{deg}$. small $\alpha$: drift 약함, balanced 거동이 실용 $t$ 영역 지배. large $\alpha$: drift 지배. crossover time $t^* \sim 1/(\alpha\,\Delta\rho)^2$ ($\Delta\rho := \max_{i,j}|\rho_i - \rho_j|$).

**권장**: parent paper의 Algorithm 1은 $\boldsymbol{\mu}_0 \propto \deg$ 사용. 정규 그래프에선 uniform과 같음. 비정규 그래프에서:

- **manifold/신호 복원** 목적 (Example 1–2): uniform $\boldsymbol{\mu}_0$ → balanced regime
- **degree-구조 탐지** (hub identification): degree-prop $\boldsymbol{\mu}_0$ + $SD^2/t^3$ normalization

## §6.3 실험 검증 (유진/소미 시뮬)

**Star $S_{13}$, $y = 0$, $b = 0.5$, random-step (Model A), $\tau = 500{,}000$**:

| $\boldsymbol{\mu}_0$ | $M_{\text{final}}$ | $SD^2_{HL}/t$ | $SD^2_{HL}/t^3$ |
|---|---|---|---|
| uniform | 15 | 수렴 | $\to 0$ |
| degree-prop ($\alpha = 1$) | 16{,}848 | 발산 ($O(t^2)$) | $\to 3.78 \times 10^{-4}$ |

이론 예측과 일치:

$$\frac{\bar b^2(\rho_h - \rho_l)^2}{3} = \frac{0.25^2 \times (0.199 - 0.067)^2}{3} = 3.63 \times 10^{-4} \approx 3.78 \times 10^{-4}\ \checkmark$$

**Helm$_{21}$ (3 regime 공존), degree-prop, $\tau = 10^7$**:

| Pair type | $\rho_i, \rho_j$ | $SD^2/t$ | $SD^2/t^2$ | $SD^2/t^3$ |
|---|---|---|---|---|
| hub–leaf (drift) | $0.054 \neq 0.040$ | 발산 | 발산 | $\to c > 0$ |
| hub–ring (intermediate) | $0.054 \approx 0.054$ | 발산 | $\to 0.49$ | $\to 0$ |
| leaf–leaf (intermediate) | $0.040 = 0.040$ | 발산 | $\to c' > 0$ | $\to 0$ |

leaf–leaf 쌍: $\rho_{L_i} = \rho_{L_j}$ by symmetry이지만 $SD^2/t$ 발산 ($SD^2 \propto t^2$, Theorem B 일치). diffusion rate $\sigma^2$는 hub의 linearly growing height에서 flow noise로 전파.

상세 시뮬: 소미의 41f2b0 블로그 (*Thm A, B의 실제 검증*), 본인의 본 글 (Three-Regime Convergence 자세히 따라가기) 참조.

---

# §7 Appendix: OB on $K_n$ 풀증명

§3.4의 DT 조건은 그래프별로 검증해야 한다. $K_n$ (complete graph) 에서는 **uniform $\boldsymbol{\mu}_0$** 가정 하 closed form으로 증명 가능. 이게 OB가 분석적으로 풀리는 유일한 그래프 클래스 (이외는 §3.4의 통계적 검증).

::: {.callout-note collapse="false" title="Proposition (OB for $K_n$)"}

$K_n$ + $T_{\max} \geq n$ + uniform $\boldsymbol{\mu}_0$ 하에서 OB:

- **General case**: $\varepsilon_{\text{OB}} = 1/(2n^2)$
- **Well-separated sub-case** (모든 인접 rank gap $> b(T_{\max}+1)$): $\varepsilon_{\text{OB}} = (H_n - 1)/n$ (여기서 $H_n := \sum_{k=1}^n 1/k$)

:::

::: {.callout-note collapse="true" title="증명: OB for $K_n$"}

$\delta := b(T_{\max} + 1)$ — 한 round에 노드 높이 최대 변화량.

**Case 1: $M \leq 2n\delta$.** 자명하게 $|B(S)| \leq (T_{\max}+1)\,M \leq 2n(T_{\max}+1)\delta =: C_{\text{triv}}$. $C_{\text{OB}}$에 흡수.

**Case 2: $M > 2n\delta$.** 노드를 $\hbar$ 내림차순 정렬: $v_{(0)}, \ldots, v_{(n-1)}$. round 시작 $X_1 \sim \mu_0 = \text{Unif}(V)$. flow rule: $v_{(j)}$에서 이웃 중 $h \leq h(v_{(j)})$인 곳으로 균등 선택. $K_n$이므로 이웃 = $V \setminus \{v_{(j)}\}$.

### Step 1 (Well-separated, 모든 gap $> \delta$)

인접 rank gap $> \delta$이면 round 내 어떤 perturbation도 ranking을 뒤집지 못함. $v_{(j)}$에서 downstream set은 정확히 $\{v_{(j+1)}, \ldots, v_{(n-1)}\}$, walker가 균등 선택 → transition $1/(n-1-j)$.

귀납으로: $j < i < n-1$ 에서 $P(\text{rank } i \text{ visited} \mid \text{start at } j) = 1/(n-i)$. $i = n-1$ (가장 낮은 rank) 은 항상 도달 → $\pi(v_{(n-1)}) = 1$.

기대 방문 횟수 (초기 draw + 이후 flow):

$$\pi(v_{(i)}) = \tfrac{1}{n} + \sum_{j=0}^{i-1} \tfrac{1}{n} \cdot \tfrac{1}{n-i} = \tfrac{1}{n} + \tfrac{i}{n(n-i)} = \tfrac{1}{n-i} \qquad (0 \leq i \leq n-2)$$

$\pi(v_{(n-1)}) = 1 = 1/(n - (n-1))$. 즉 모든 $i$에 대해 $\pi(v_{(i)}) = 1/(n-i)$.

$B$ 형성:

$$B = \sum_{i=0}^{n-1} \tfrac{\hbar_{(i)}}{n-i}$$

constraint $\sum \hbar = 0$, $\hbar$ 정렬, $\hbar_{(0)} - \hbar_{(n-1)} = M$ 하에서 maximize. weight $1/(n-i)$가 $i$ 따라 증가, $\hbar_{(i)}$가 $i$ 따라 감소 → maximum at one-high profile: $\hbar_{(0)} = (n-1)M/n$, $\hbar_{(i)} = -M/n$ for $i \geq 1$.

$$B \leq \tfrac{(n-1)M/n}{n} - \tfrac{M}{n}\sum_{k=1}^{n-1}\tfrac{1}{k} = \tfrac{(n-1)M}{n^2} - \tfrac{M\,(H_n - 1)}{n}$$

$\tfrac{n-1}{n^2} - \tfrac{H_n - 1}{n} = \tfrac{n-1 - n(H_n - 1)}{n^2} = \tfrac{2n - 1 - nH_n}{n^2}$. 정리하면 $B \leq -(H_n - 1)M/n$ (well-separated $\varepsilon_{\text{OB}}$).

### Step 2 (General, 어떤 gap $\leq \delta$)

rank를 cluster로 partition: 인접 rank $i, i+1$이 같은 cluster ⟺ $\hbar_{(i)} - \hbar_{(i+1)} \leq \delta$. cluster $C_1, \ldots, C_K$ ($K \leq n$). inter-cluster gap $> \delta$ → inter-cluster ordering이 round 내 stable.

**Inter-cluster 동력학**: walker의 cluster 간 이동은 well-sep $K$-super-node 문제와 동일.

**Perturbative contribution**: round 내 총 방문 $m_r \leq T_{\max}+1$. 각 방문 노드의 $\hbar$가 cluster 대표 $h_k^{\min}$와 최대 $n\delta$ 차이.

$$|B_{\text{extra}}| \leq (T_{\max}+1) \cdot n\delta =: C_{\text{pert}}$$

$M$ 무관.

**Clean contribution**: $B_{\text{clean}} := B - B_{\text{extra}} = \sum_k \Pi_k\,h_k^{\min}$ ($\Pi_k$ = cluster $C_k$ 전체 방문 횟수).

::: {.callout-note collapse="true" title="Lemma (Cluster LP)"}

$n$개 노드를 $K \geq 2$ 개 ordered cluster로 (sizes $|C_k|$, heights $h_1 > \cdots > h_K$, $\sum_k |C_k| h_k = 0$, $M_{\text{cl}} := h_1 - h_K$) 분할할 때, $K_n$의 clean contribution:

$$B_{\text{clean}} \leq -M_{\text{cl}}/n^2$$

**증명** ($K = 2$): sizes $[s, n-s]$. 초기 draw 기여 $0$. walker가 $C_1 \to C_2$ 항상 이동 → $\Pi_1 = s/n$, $\Pi_2 = 1$. $h_1 = (n-s)M_{\text{cl}}/n$, $h_2 = -sM_{\text{cl}}/n$ 대입:

$$B_{\text{clean}} = \tfrac{s}{n} \cdot \tfrac{(n-s)M_{\text{cl}}}{n} + 1 \cdot \tfrac{-sM_{\text{cl}}}{n} = -\tfrac{s^2}{n^2}M_{\text{cl}}$$

$s \in \{1, \ldots, n-1\}$ minimize: $s = 1$ → $B_{\text{clean}} \leq -M_{\text{cl}}/n^2$.

**$K \geq 3$**: $K-1$ free variable의 finite LP. **computational 검증**: $n \leq 40$, $K \leq n$ 모든 partition에 대해 $B_{\text{clean}} \leq -M_{\text{cl}}/n^2$, $K = 2$ + $[1, n-1]$이 항상 worst case 도달 — open problem로 일반 $n$의 induction 증명. $\square$

:::

$M_{\text{cl}} \geq M - n\delta \geq M/2$ (Case 2 가정 $M > 2n\delta$).

**Conclusion**:

$$B = B_{\text{clean}} + B_{\text{extra}} \leq -\tfrac{M_{\text{cl}}}{n^2} + C_{\text{pert}} \leq -\tfrac{M/2}{n^2} + C_{\text{pert}} + C_{\text{triv}} =: -\varepsilon_{\text{OB}}\,M + C_{\text{OB}}$$

$\varepsilon_{\text{OB}} := 1/(2n^2) > 0$. $\square$

:::

---

# §8 외부 의존

본 글의 결과는 두 가지 외부 paper-team 산출물에 의존:

`-` **$\rho_i$의 무조건 존재성**: regime 분류는 $\rho_i := \lim_t \tfrac{1}{t}\sum \mathbf{1}\{X_s = v_i\}$의 a.s. 존재를 전제. random-step variant에서 $\mu_{\min} > 0$인 모든 연결 그래프에 대해 **무조건 성립** — 오큐리의 stochastic approximation + piecewise-linear ODE 인자, [`occupation_slln_general.tex`](../paper/260514_guebin/오큐리/occupation_slln_general.tex). Theorem A는 $\rho_i = 1/n$을 자체 가정으로 가지므로 별도 필요 X; Theorem B, C는 $\rho_i$ 존재를 input으로 받음.

`-` **3-regime 완전성**: §2 Definition 의 세 regime이 상호 배반 + 합쳐서 완전 (4번째 power scaling 불필요). 상세: 해리스의 [`regime_exhaustiveness.tex`](../paper/260514_guebin/해리스/regime_exhaustiveness.tex). 본 글은 완전성을 input으로 받고 각 regime 내부 수렴만 증명.

표준 reference:

`-` Choi, G. and Oh, H.-S. (2026). *Heavy-Snow Transform for Analysis of Data on Graphs*. Manuscript.
`-` Meyn, S. and Tweedie, R. L. (2009). *Markov Chains and Stochastic Stability* (2nd ed.). Cambridge UP. — §5.2 (small sets), §11.3 (Foster–Lyapunov), §14.3 (moment bounds), §17 (LLN).

---

# Open items

`-` **OB for general connected graph**: $K_n$만 closed form 증명. 일반 그래프는 §3.4 통계적 검증으로 대체.

`-` **Cluster LP, $K \geq 3$**: $K = 2$ analytic, $K \geq 3$ 은 $n \leq 40$ computational. 일반 $n$의 inductive 증명 open.

`-` **Drift regime에서 $\rho_i$ 존재의 일반 증명**: 오큐리의 occupation_slln_general이 무조건 존재를 보임. 다만 해당 증명의 자체검증 일부 still pending.

`-` **모델 B로의 transfer**: 2026-05-17 01:10 회사 결정으로 공식 모델은 B. 본 글의 모델 A 결과가 모델 B에서 어떻게 바뀌는지는 후속 작업 (Theorem B, C의 leading constants는 보존, Theorem A의 $c_{ij}$는 graph-specific 수치 변동 가능).
