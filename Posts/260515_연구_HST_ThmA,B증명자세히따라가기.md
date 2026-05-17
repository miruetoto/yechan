---
title: 연구 ▷ HST ▷ Three-Regime Convergence 자세히 따라가기
author: 유진
date: 05/15/2026
draft: false
output-file: 260515_48a84c.html
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

# Theorem B 증명

B가 더 쉬우니까 먼저.

### Step 1: 높이 분해

$$h(v_i, t) = h(v_i, 0) + \sum_{s=1}^t b'_s \cdot \mathbf{1}\{X_s = v_i\}$$

$b'_s$와 $X_s$는 독립 (알고리즘 순서: $X_s$ 결정 후 $b'_s$ 추출). 분리하면:

$$\sum b'_s \mathbf{1}\{X_s = v_i\} = \bar{b}\sum \mathbf{1}\{X_s = v_i\} + \underbrace{\sum (b'_s - \bar{b})\mathbf{1}\{X_s = v_i\}}_{=: M_t^{(i)}}$$

첫째 항: $\rho_i$ 정의에 의해 $\sim \bar{b}\rho_i t$.

둘째 항: $M_t^{(i)}$는 $\mathcal{F}_t$-마팅게일이다. 정확히는 자연 filtration $\mathcal{F}_s = \sigma(X_1, b'_1, \ldots, X_s, b'_s)$에서:

$$\mathbb{E}[\Delta M_s^{(i)} \mid \mathcal{F}_{s-1}] = \mathbb{E}\left[\mathbf{1}\{X_s = v_i\} \cdot \underbrace{\mathbb{E}[b'_s - \bar{b} \mid \mathcal{F}_{s-1}, X_s]}_{= 0}\right] = 0$$

tower property를 사용했다. Block step에서 $X_s$는 fresh $\mu_0$-draw이므로 $\mathcal{F}_{s-1}$-measurable이 아니지만, $b'_s$는 $(\mathcal{F}_{s-1}, X_s)$와 독립이므로 내부 조건부 기대값이 0. 증분 bounded ($\leq b/2$)이므로 마팅게일 SLLN에 의해 $M_t^{(i)}/t \to 0$ a.s.

**결론:** $h(v_i, t) = h(v_i, 0) + \bar{b}\rho_i t + o(t)$ a.s.

### Step 2-7: $SD^2$ 전개

$D(s) = h(v_i, s) - h(v_j, s) = \gamma s + \xi(s)$, $\gamma = \bar{b}(\rho_i - \rho_j)$, $\xi(s) = o(s)$.

$$SD^2_{ij}(t) = \sum_{s=0}^t D(s)^2 = \gamma^2 \sum s^2 + 2\gamma \sum s\xi(s) + \sum \xi(s)^2$$

`-` 주항: $\sum s^2 = t^3/3 + O(t^2)$

`-` 교차항: $|\xi(s)| \leq \epsilon s$ for $s \geq S_\epsilon$ 이므로 $|\sum s\xi(s)| \leq C_\epsilon + \epsilon t^3/3 = o(t^3)$

`-` 잔차항: 같은 논법으로 $\sum \xi(s)^2 = o(t^3)$

$$\frac{SD^2_{ij}(t)}{t^3} \to \frac{\gamma^2}{3} = \frac{\bar{b}^2(\rho_i - \rho_j)^2}{3} \qquad \square$$

Theorem B는 Foster-Lyapunov도 Doeblin도 필요 없다. 마팅게일 SLLN + Toeplitz 보조정리만으로 끝난다.

---

# Theorem A 증명

Balanced regime ($\rho_i = 1/n$)에서는 drift가 사라지므로 $SD^2 = \sum \xi(s)^2$이 되고, 이것의 $SD^2/t \to c_{ij}$ 수렴을 보이려면 에르고딕 정리가 필요하다. 증명 체인:

$$\text{OB} \to \text{round drift} \to \text{Doeblin} \to \psi\text{-irred.} \to \text{Harris rec.} \to \text{moment bound} \to \text{LLN}$$

### 가정

`-` 연결 그래프 $\mathcal{G}$, full-support $\mu_0$ ($\mu_{\min} > 0$)

`-` $\rho_i = 1/n$ for all $i$

`-` OB 조건: $\mathcal{B}(S) = \sum_u \pi(u,S)\hbar(u) \leq -\kappa_G M + C_{\text{OB}}$ (아래 설명)

### Augmented state

증명 전체에서 Markov chain으로 다룰 상태는 **augmented state**:

$$S(t) := (\boldsymbol{\delta}_t,\, X_t,\, Z_t) \in \mathcal{X}^* := \mathbb{R}^{n-1} \times V \times \{0, 1, \ldots, T_{\max}\}$$

각 성분:

- $\boldsymbol{\delta}_t := (h(v_2, t) - h(v_1, t),\ \ldots,\ h(v_n, t) - h(v_1, t)) \in \mathbb{R}^{n-1}$ — 노드 1을 기준으로 한 **상대 높이 벡터** (절대 높이 $h$는 $t$에 따라 무한히 커지므로 차이만 추적)
- $X_t \in V$ — 현재 walker 위치
- $Z_t \in \{0, 1, \ldots, T_{\max}\}$ — 연속 flow 카운터 (직전 step부터 몇 번 연속 flow했는지; $T_{\max}$ 도달하면 block 강제)

`-` $\mathcal{B}(\mathcal{X}^*)$ — $\mathcal{X}^*$의 Borel σ-algebra

`-` $P^m(s, B)$ — $S(0) = s$에서 출발해 $m$ step 후 $B \in \mathcal{B}(\mathcal{X}^*)$ 에 도달할 확률

`-` 정상 분포는 $\pi^*$로 표기 (존재성/유일성은 증명 체인의 결과)

`-` $\mathcal{B}_\Phi$ — **Lyapunov 함수 $\Phi$의 sublevel set** (아래첨자 $\Phi$는 "$\Phi$로 정의된 set"이라는 의미):

  $$\mathcal{B}_\Phi := \{S \in \mathcal{X}^* : \Phi(S) \leq R\}, \qquad R := \tfrac{n}{2}\left(\tfrac{C+1}{b\kappa_G}\right)^2$$

  $R$은 **Foster drift가 $\leq -1$이 되는 임계값 밖**의 안쪽 영역을 잘라냄. 즉 $S \notin \mathcal{B}_\Phi$ ⟺ $\Phi(S) > R$ ⟺ "$M$이 충분히 큼" ⟺ $\mathbb{E}[\Delta\Phi \mid S] \leq -1$.

  **왜 별도로 잘라내는가**: Foster drift는 큰 $\Phi$ 영역에서만 강한 음수. **작은 $\Phi$ 영역(= $\mathcal{B}_\Phi$ 안)에서는 drift가 보장 안 됨**. 그 안에서는 다른 기제(Doeblin minorization)로 chain을 컨트롤. → "**큰 $\Phi$에선 Foster, 작은 $\Phi$(= $\mathcal{B}_\Phi$ 안)에선 Doeblin**"의 역할 분담.

### 보조 결과 1: $\Phi$ drift

**정의.** $\hbar(v_i, t) := h(v_i, t) - \bar{h}(t)$ (centered height), $\Phi(t) := \sum_i \hbar(v_i, t)^2$ (높이 분산 $\times n$).

$\Phi$가 크면 노드 간 높이 차이가 크고, $\Phi = 0$이면 모든 노드가 같은 높이.

**Fact 1 ($\hbar$ 변화).** step $t+1$에서 노드 $X_{t+1}$에 $b'_{t+1}$이 쌓이면:

$$\hbar(X_{t+1},\, t+1) = \hbar(X_{t+1},\, t) + b'_{t+1}\frac{n-1}{n}, \qquad \hbar(u,\, t+1) = \hbar(u,\, t) - \frac{b'_{t+1}}{n} \quad (u \neq X_{t+1})$$

::: {.callout-note collapse="true" title="증명: Fact 1 ($\\hbar$ 변화)"}

구하고 싶은 것: step $t+1$ 후 각 노드의 $\hbar$가 어떻게 바뀌는가.

$\hbar = h - \bar{h}$이므로, $h$의 변화와 $\bar{h}$의 변화를 각각 구하자.

$h$의 변화: $X_{t+1}$에만 $b'_{t+1}$이 쌓이므로:

$$h(X_{t+1},\, t+1) = h(X_{t+1},\, t) + b'_{t+1}, \qquad h(u,\, t+1) = h(u,\, t) \quad (u \neq X_{t+1})$$

$\bar{h}$의 변화: 전체 합이 $b'_{t+1}$만큼 늘고 $n$으로 나누므로:

$$\bar{h}(t+1) = \bar{h}(t) + \frac{b'_{t+1}}{n}$$

$\hbar = h - \bar{h}$에 대입:

$$\hbar(X_{t+1},\, t+1) = h(X_{t+1},\, t) + b'_{t+1} - \bar{h}(t) - \frac{b'_{t+1}}{n} = \hbar(X_{t+1},\, t) + b'_{t+1}\left(1 - \frac{1}{n}\right) = \hbar(X_{t+1},\, t) + b'_{t+1}\frac{n-1}{n}$$

$$\hbar(u,\, t+1) = h(u,\, t) - \bar{h}(t) - \frac{b'_{t+1}}{n} = \hbar(u,\, t) - \frac{b'_{t+1}}{n} \qquad (u \neq X_{t+1}) \qquad \square$$

:::

**Fact 2.** $\sum_i \hbar(v_i, t) = 0$ (항상 성립, $\hbar$의 정의에서 평균을 뺐으므로).

**Lemma ($\Phi$ drift).** step $t+1$에서 노드 $X_{t+1}$에 $b'_{t+1} \sim \text{Unif}(0,b)$가 쌓일 때:

$$\mathbb{E}[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t, X_{t+1}] = b\,\hbar(X_{t+1}, t) + \frac{(n-1)b^2}{3n}$$

높은 노드($\hbar > 0$)에 눈이 쌓이면 $\Phi$ 증가, 낮은 노드($\hbar < 0$)에 쌓이면 $\Phi$ 감소. 둘째 항은 항상 양수인 noise.

::: {.callout-note collapse="true" title="증명: $\\Phi$ drift lemma"}

우리가 구하고 싶은 것은 $\mathbb{E}[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t, X_{t+1}]$이다. 이를 위해 먼저 $\Phi(t+1) - \Phi(t)$를 계산하자.

$\Phi(t+1) - \Phi(t)$를 계산하려면 $\Phi(t+1)$과 $\Phi(t)$가 각각 뭔지 알아야 한다.

$\Phi(t)$는 그냥 정의대로:

$$\Phi(t) = \hbar(X_{t+1},\, t)^2 + \sum_{u \neq X_{t+1}} \hbar(u,\, t)^2$$

$\Phi(t+1)$은 Claim의 $\hbar$ 변화를 대입하면:

$$\Phi(t+1) = \left[\hbar(X_{t+1},\, t) + b'_{t+1}\frac{n-1}{n}\right]^2 + \sum_{u \neq X_{t+1}} \left[\hbar(u,\, t) - \frac{b'_{t+1}}{n}\right]^2$$

이제 빼자.

$$\Phi(t+1) - \Phi(t) = \underbrace{\left[\hbar(X_{t+1}, t) + b'_{t+1}\frac{n-1}{n}\right]^2 - \hbar(X_{t+1}, t)^2}_{X_{t+1}\text{ 항: } (A+B)^2 - A^2} + \underbrace{\sum_{u \neq X_{t+1}}\left[\left(\hbar(u,t) - \frac{b'_{t+1}}{n}\right)^2 - \hbar(u,t)^2\right]}_{\text{나머지 항: } (A+B)^2 - A^2}$$

각 항에 $(A+B)^2 - A^2 = 2AB + B^2$을 전개하면:

$$= \underbrace{2\hbar(X_{t+1}, t)\cdot\frac{(n-1)b'_{t+1}}{n} + \frac{(n-1)^2(b'_{t+1})^2}{n^2}}_{X_{t+1}\text{ 항}} + \underbrace{\sum_{u \neq X_{t+1}}\left[-\frac{2b'_{t+1}}{n}\hbar(u,t) + \frac{(b'_{t+1})^2}{n^2}\right]}_{\text{나머지 항}}$$

나머지 항을 정리하자. $\sum_{u \neq X_{t+1}} \hbar(u,t) = -\hbar(X_{t+1},t)$ ($\because \sum_i \hbar_i = 0$)이므로:

$$\text{나머지 항} = \frac{2b'_{t+1}}{n}\hbar(X_{t+1},t) + \frac{(n-1)(b'_{t+1})^2}{n^2}$$

$X_{t+1}$ 항 + 나머지 항:

$$\Phi(t+1) - \Phi(t) = 2\hbar(X_{t+1},t)\cdot\frac{(n-1)b'_{t+1}}{n} + \frac{(n-1)^2(b'_{t+1})^2}{n^2} + \frac{2b'_{t+1}}{n}\hbar(X_{t+1},t) + \frac{(n-1)(b'_{t+1})^2}{n^2}$$

$$= 2b'_{t+1}\,\hbar(X_{t+1},t)\underbrace{\left[\frac{n-1}{n} + \frac{1}{n}\right]}_{=\,1} + (b'_{t+1})^2\underbrace{\left[\frac{(n-1)^2 + (n-1)}{n^2}\right]}_{=\,\frac{n-1}{n}}$$

$$= 2b'_{t+1}\,\hbar(X_{t+1},t) + \frac{n-1}{n}(b'_{t+1})^2$$

이제 원래 목표로 돌아오자. $b'_{t+1} \sim \text{Unif}(0,b)$에 대해 기대값을 취하면 ($\mathbb{E}[b'_{t+1}] = b/2$, $\mathbb{E}[(b'_{t+1})^2] = b^2/3$):

$$\mathbb{E}[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t, X_{t+1}] = 2 \cdot \frac{b}{2} \cdot \hbar(X_{t+1},t) + \frac{n-1}{n} \cdot \frac{b^2}{3}$$

$$= b\,\hbar(X_{t+1},t) + \frac{(n-1)b^2}{3n} \qquad \square$$

:::


### 보조 결과 2: Shape decomposition

보조 결과 1은 **한 스텝**의 $\Phi$ 변화를 구했다. Theorem A 증명에서는 **라운드 전체**의 $\Phi$ 변화가 필요하다. 이를 위해 먼저 기호를 정의한다.

**정의.**

`-` $t_r$: $r$번째 block-flag time (라운드의 시작 시점)

`-` $m_r = t_{r+1} - t_r$: 라운드 $r$의 길이 ($\leq T_{\max}+1$)

`-` $\tilde{X}_j = X_{t_r + j}$: 라운드 $r$의 $j$번째 방문 노드

`-` $\#_u = |\{j : \tilde{X}_j = u\}|$: 라운드 $r$ 동안 노드 $u$를 방문한 횟수

`-` $\pi(u) = \mathbb{E}[\#_u \mid S(t_r)]$: 기대 방문 횟수

`-` $M = \max_i \hbar(v_i, t_r) - \min_i \hbar(v_i, t_r)$: 높이 range

라운드 전체의 $\Phi$ 변화를 구하려면 $\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1)$을 계산해야 한다. 문제는 $\hbar(\tilde{X}_j, t_r + j - 1)$이 **라운드 중간 시점**의 값이라 라운드 시작 값 $\hbar(u, t_r)$과 다르다는 것이다. Shape decomposition은 이 차이가 bounded error임을 말한다.

**Lemma (Shape decomposition).**

$$\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1) = \sum_u \#_u \hbar(u, t_r) + \Theta, \qquad |\Theta| \leq \frac{(T_{\max}+1)T_{\max}}{2}b$$

좌변: 라운드 중간 시점의 $\hbar$ 합. 우변의 $\sum_u$는 모든 노드 $u \in V$에 대한 합이고, $\#_u$는 노드 $u$를 방문한 횟수. 첫째 항: 라운드 시작 시점의 $\hbar$로 계산한 합. $\Theta$: 라운드 중 높이가 변해서 생기는 오차 (bounded).

::: {.callout-note collapse="true" title="증명: Shape decomposition"}

구하고 싶은 것: 좌변 $\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1)$을 라운드 시작 시점의 $\hbar$로 표현하기.

$\hbar(\tilde{X}_j, t_r + j - 1)$은 라운드 $j$번째 방문 노드의 **시점 $t_r + j - 1$에서의** centered height이다. Fact 1을 반복 적용하면, 시점 $t_r$에서의 값 $\hbar(\tilde{X}_j, t_r)$에 라운드 도중의 보정이 더해진다:

$$\hbar(\tilde{X}_j, t_r + j - 1) = \hbar(\tilde{X}_j, t_r) + \sum_{k=1}^{j-1} b'_{t_r+k}\left[\mathbf{1}\{\tilde{X}_k = \tilde{X}_j\} - \frac{1}{n}\right]$$

첫째 항: 라운드 시작 시점의 $\hbar$. 둘째 항: step $k$에서 $\tilde{X}_j$와 같은 노드에 눈이 쌓이면 $(1-1/n)$, 다른 노드면 $-1/n$만큼 $\hbar$가 바뀐 누적 (Fact 1에 의해).

$j = 1, \ldots, m_r$에 대해 합산하자.

$$\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1) = \underbrace{\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r)}_{\text{첫째 항의 합}} + \underbrace{\sum_{j=1}^{m_r}\sum_{k=1}^{j-1} b'_{t_r+k}\left[\mathbf{1}\{\tilde{X}_k = \tilde{X}_j\} - \frac{1}{n}\right]}_{=:\,\Theta}$$

첫째 항의 합: 노드 $u$를 $\#_u$번 방문했으므로:

$$\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r) = \sum_u \#_u\, \hbar(u, t_r)$$

$|\Theta|$의 bound: 이중합에서 각 항은 $|b'_{t_r+k}| \cdot |\mathbf{1}\{\cdots\} - 1/n| \leq b \cdot 1$. 항의 개수는:

$$\sum_{j=1}^{m_r}(j-1) = \frac{m_r(m_r - 1)}{2} \leq \frac{(T_{\max}+1)T_{\max}}{2}$$

따라서:

$$|\Theta| \leq \frac{(T_{\max}+1)T_{\max}}{2}\,b \qquad \square$$

:::


보조 결과 1의 $\Phi$ drift lemma와 이 shape decomposition을 합치면 **라운드 단위 drift**를 얻는다:

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] = b\,\mathcal{B}(S) + (\text{bounded terms})$$

여기서 $\mathcal{B}(S) = \sum_u \pi(u) \hbar(u)$이다. Bounded terms는 $M$과 무관한 상수. 따라서 **$M$에 비례하는 유일한 항은 $\mathcal{B}(S)$**이다.

::: {.callout-note collapse="true" title="증명: 라운드 단위 drift"}

구하고 싶은 것: $\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)]$.

라운드 $r$은 step $t_r + 1$부터 $t_{r+1}$까지 총 $m_r$ step이다. $\Phi$ drift lemma (보조 결과 1)를 각 step에 적용하고 합산하면:

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] = \sum_{j=1}^{m_r} \mathbb{E}\left[b\,\hbar(\tilde{X}_j, t_r + j - 1) + \frac{(n-1)b^2}{3n} \;\middle|\; S(t_r)\right]$$

$$= b\,\mathbb{E}\left[\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1) \;\middle|\; S(t_r)\right] + \frac{(n-1)b^2}{3n}\,\mathbb{E}[m_r \mid S(t_r)]$$

둘째 항: $\mathbb{E}[m_r] \leq T_{\max}+1$이므로 bounded.

첫째 항: shape decomposition (보조 결과 2)을 적용하면:

$$\mathbb{E}\left[\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1) \;\middle|\; S(t_r)\right] = \mathbb{E}\left[\sum_u \#_u\,\hbar(u, t_r) + \Theta \;\middle|\; S(t_r)\right]$$

$$= \sum_u \underbrace{\mathbb{E}[\#_u \mid S(t_r)]}_{=\,\pi(u)}\,\hbar(u, t_r) + \mathbb{E}[\Theta \mid S(t_r)]$$

$$= \mathcal{B}(S) + \mathbb{E}[\Theta \mid S(t_r)]$$

$|\Theta| \leq (T_{\max}+1)T_{\max}b/2$이므로 $\mathbb{E}[\Theta \mid S(t_r)]$도 bounded.

합치면:

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] = b\,\mathcal{B}(S(t_r)) + \underbrace{b\,\mathbb{E}[\Theta \mid S(t_r)] + \frac{(n-1)b^2}{3n}\mathbb{E}[m_r \mid S(t_r)]}_{\text{bounded terms}} \qquad \square$$

:::

### 보조 결과 3: OB (Occupation-Bias) 조건

**OB 조건:** $\mathcal{B}(S) = \sum_u \pi(u)\hbar(u) \leq -\kappa_G M + C_{\text{OB}}$

$K_n$에서 well-separated case: $\kappa_G = (H_n - 1)/n$. Near-tie case: $\kappa_G = 1/(2n^2)$.

::: {.callout-warning collapse="true" title="보충: OB 조건의 의미와 통계적 검증"}
**직관:** 지형이 울퉁불퉁할수록($M$이 클수록) 눈이 낮은 곳에 편향되어 쌓인다는 조건이다. $\kappa_G > 0$은 그래프 $G$의 **자기교정 세기**이고, $C_{\text{OB}}$는 $M$이 작을 때의 여유분이다.

**적용 범위:** OB가 성립하면 Theorem A의 증명 체인이 작동한다. $K_n$에서는 $\kappa_G$를 이론적으로 유도했지만, 일반 그래프에서는 $\kappa_G$와 $C_{\text{OB}}$가 그래프 구조에 따라 달라지며, Star 그래프처럼 아예 OB가 성립하지 않는 경우도 있다(drift regime).

**통계적 검증:** 이론 증명이 없는 그래프에서도 OB를 직접 확인할 수 있다. HST를 $T$ 스텝 돌리면서 매 스텝 $R_t = \hbar(X_t)/M_t$를 기록한 뒤, batch mean의 $t$-test로 $\hat{\kappa} = -\bar{R} > 0$인지 검정하면 된다. $\kappa_G$나 $C_{\text{OB}}$의 구체적 값을 몰라도, 데이터에서 OB 성립 여부와 $\kappa$의 크기를 추정할 수 있다.
:::

### 보조 결과 4: Moment bound

**Proposition (Moment bound).** $\mathbb{E}_{\pi^*}[\Phi] < \infty$.

Round drift $\mathbb{E}[\Delta\Phi] \leq -b\kappa_G M + C$에서 $M \sim \sqrt{\Phi}$이므로 drift 차수가 $-\Phi^{1/2}$. 이것만으로는 $\mathbb{E}_{\pi^*}[\Phi] < \infty$를 얻기에 부족하다 (Meyn-Tweedie로 $r < 1/2$까지만 가능). 해법: Lyapunov를 $V_{\text{Ly}} = \Phi^2$로 올린다.

::: {.callout-note collapse="true" title="증명: Moment bound"}

구하고 싶은 것: $\mathbb{E}_{\pi^*}[\Phi] < \infty$.

이를 위해 $V_{\text{Ly}} = \Phi^2$의 drift를 구한다. 먼저 $\Delta V_{\text{Ly}} = V_{\text{Ly}}(t_{r+1}) - V_{\text{Ly}}(t_r)$를 계산하자.

$$V_{\text{Ly}}(t_{r+1}) = \Phi(t_{r+1})^2 = (\Phi(t_r) + \Delta\Phi_r)^2$$

$$V_{\text{Ly}}(t_r) = \Phi(t_r)^2$$

빼면:

$$\Delta V_{\text{Ly}} = (\Phi + \Delta\Phi_r)^2 - \Phi^2 = 2\Phi \cdot \Delta\Phi_r + (\Delta\Phi_r)^2$$

기대값을 취하면:

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S(t_r)] = \underbrace{2\Phi \cdot \mathbb{E}[\Delta\Phi_r \mid S(t_r)]}_{=:\,\text{Term1}} + \underbrace{\mathbb{E}[(\Delta\Phi_r)^2 \mid S(t_r)]}_{=:\,\text{Term2}}$$

**Fact 3.** $\Phi \geq M^2/2$, $\Phi \leq nM^2$.

$$
\begin{aligned}
\Phi = \sum_i \hbar_i^2
&\geq \hbar_{\max}^2 + \hbar_{\min}^2 \\
&\geq \tfrac{(\hbar_{\max} + |\hbar_{\min}|)^2}{2} = \tfrac{M^2}{2} && (\because (A^2{+}B^2)\geq(A{+}B)^2/2,\ \textstyle\sum_i\hbar_i=0) \\
\Phi
&\leq nM^2 && (\because n\text{개 항 각각 } \leq M^2)
\end{aligned}
$$

**Term1.**

$$
\begin{aligned}
2\Phi \cdot \mathbb{E}[\Delta\Phi_r \mid S(t_r)]
&= 2\Phi \sum_{j=1}^{m_r} \mathbb{E}\!\left[b\,\hbar(\tilde{X}_j, t_r{+}j{-}1) + \tfrac{(n{-}1)b^2}{3n} \,\middle|\, S(t_r)\right] && (\because \text{보조 결과 1}) \\
&= 2\Phi\left\{ b\,\mathbb{E}\!\left[\sum_{j=1}^{m_r}\hbar(\tilde{X}_j, t_r{+}j{-}1) \,\middle|\, S(t_r)\right] + \tfrac{(n{-}1)b^2}{3n}\,\mathbb{E}[m_r \mid S(t_r)] \right\} \\
&= 2\Phi\left\{ b\,\mathbb{E}\!\left[\sum_u \#_u\,\hbar(u, t_r) + \Theta \,\middle|\, S(t_r)\right] + \tfrac{(n{-}1)b^2}{3n}\,\mathbb{E}[m_r \mid S(t_r)] \right\} && (\because \text{보조 결과 2}) \\
&= 2\Phi\left\{ b\,\mathcal{B}(S(t_r)) + b\,\mathbb{E}[\Theta \mid S(t_r)] + \tfrac{(n{-}1)b^2}{3n}\,\mathbb{E}[m_r \mid S(t_r)] \right\} && (\because \pi(u):=\mathbb{E}[\#_u \mid S(t_r)],\ \mathcal{B}(S):=\textstyle\sum_u\pi(u)\hbar(u)) \\
&\leq 2\Phi\left\{ b\,\mathcal{B}(S(t_r)) + C' \right\} && (\because |\Theta|,\ m_r \text{ bounded}) \\
&\leq 2\Phi(-b\kappa_G M + bC_{\text{OB}} + C') && (\because \text{보조 결과 3 (OB)}) \\
&= 2\Phi(-b\kappa_G M + C) && (C := bC_{\text{OB}} + C') \\
&= -2b\kappa_G\,\Phi M + 2C\Phi \\
&\leq -b\kappa_G M^3 + 2CnM^2 && (\because \text{Fact 3})
\end{aligned}
$$

$$\therefore\ \text{Term1} \leq -b\kappa_G M^3 + 2CnM^2$$

**Term2.**

$$
\begin{aligned}
|\Delta\Phi \text{ in 1 step}|
&\leq 2b(M + T_{\max}b) + \tfrac{(n{-}1)b^2}{n} && (\because \beta \leq b,\ |\hbar(t_r{+}j)| \leq M+T_{\max}b) \\
|\Delta\Phi_r|
&\leq (T_{\max}{+}1)\!\left[2b(M+T_{\max}b)+\tfrac{(n{-}1)b^2}{n}\right] =: A_1 M + A_2 && (\because m_r \leq T_{\max}{+}1) \\
\mathbb{E}[(\Delta\Phi_r)^2 \mid S(t_r)]
&\leq (A_1 M+A_2)^2 \leq 2A_1^2 M^2 + 2A_2^2 && (\because (a{+}b)^2 \leq 2a^2 + 2b^2)
\end{aligned}
$$

상수 정리:

$$A_1 := 2b(T_{\max}{+}1),\quad A_2 := (T_{\max}{+}1)\!\left[2bT_{\max}b+\tfrac{(n{-}1)b^2}{n}\right],\quad C_3 := 2A_1^2,\quad C_4 := 2A_2^2$$

$$\therefore\ \text{Term2} \leq C_3 M^2 + C_4$$

**합치기.**

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S(t_r)] \leq -b\kappa_G M^3 + 2CnM^2 + C_3 M^2 + C_4$$

$$= -b\kappa_G M^3 + C_5 M^2 + C_6$$

$M$이 크면 $M^3$이 $M^2$를 압도한다. $M \geq M_0 := 2C_5/(b\kappa_G)$이면:

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S(t_r)] \leq -\frac{b\kappa_G}{2} M^3 + C_6$$

$V$로 표현하자. $\Phi \leq nM^2$에서 $M \geq (\Phi/n)^{1/2}$이므로:

$$M^3 \geq \frac{\Phi^{3/2}}{n^{3/2}} = \frac{V_{\text{Ly}}^{3/4}}{n^{3/2}}$$

대입하면:

$$\mathbb{E}[\Delta V_{\text{Ly}} \mid S(t_r)] \leq -\frac{b\kappa_G}{2n^{3/2}} V_{\text{Ly}}^{3/4} + C_6 \cdot \mathbf{1}_{\{M \leq M_0\}}$$

Meyn-Tweedie Thm 14.3.7: $\mathbb{E}[\Delta V_{\text{Ly}}] \leq -f + C \cdot \mathbf{1}_B$ 형태이면 $\mathbb{E}_{\pi^*}[f] < \infty$. 따라서:

$$\mathbb{E}_{\pi^*}[V_{\text{Ly}}^{3/4}] = \mathbb{E}_{\pi^*}[\Phi^{3/2}] < \infty$$

$\Phi \geq 0$에서 $\Phi \leq \Phi^{3/2} + 1$이므로 ($\Phi \geq 1$이면 $\Phi \leq \Phi^{3/2}$, $\Phi < 1$이면 $\Phi \leq 1$):

$$\mathbb{E}_{\pi^*}[\Phi] \leq \mathbb{E}_{\pi^*}[\Phi^{3/2}] + 1 < \infty \qquad \square$$

:::

### 보조 결과 5: Doeblin minorization

**Proposition (Doeblin minorization).** 임의의 bounded set $\mathcal{B}_\Phi \subset \mathcal{X}^*$ 에 대해, 어떤 $\eta > 0$, $K^* \in \mathbb{N}$, 그리고 $\mathcal{X}^*$ 위의 확률측도 $\nu$가 존재하여 다음이 성립한다:

$$P^{K^*}(s, B) \geq \eta\,\nu(B), \qquad \forall\,s \in \mathcal{B}_\Phi,\ \ \forall\,B \in \mathcal{B}(\mathcal{X}^*)$$

(기호 $\mathcal{X}^*$, $\mathcal{B}(\mathcal{X}^*)$, $P^m$ 은 위 §Augmented state 참조.)

**Pushforward density lemma (보조).** $\xi_1, \ldots, \xi_n \stackrel{\text{iid}}{\sim} \text{Unif}(0,b)$, $\Lambda(\xi) := (\xi_2 - \xi_1, \ldots, \xi_n - \xi_1)$. Pushforward 밀도 $\rho_\Lambda(\mathbf{y}) \geq b^{1-n}/2$ on $\mathcal{O} := \{\max_i |y_i| < b/4\}$.

::: {.callout-important collapse="true" title="해석: 한 줄 직관 — '공통 운명'"}

Doeblin minorization을 한 줄로 요약하면: **"어디서 출발하든 결국 겹치게 되는 최소한의 공통 운명이 존재한다."**

**비유.** 전국 어디서 출발하든 $K^*$번 이동하면, 최소 $\eta$ (예: 10%) 확률로 모두 같은 목적지 분포 $\nu$ (예: '서울역' 주변)에 떨어지게 만드는 **마법의 규칙**이 숨어있다.

**$P^{K^*}(s, B) \geq \eta\,\nu(B)$가 말하는 것:**

- **독립적인 바닥 (lower bound)**: 출발점 $s \in \mathcal{B}_\Phi$가 아무리 극단적이어도, 도착 분포의 '바닥'을 공통 측도 $\nu$가 받쳐줌. **출발지에 무관한 양의 공통 부분**이 존재.
- **과거 망각 (coupling)**: 이 공통 부분 덕분에 두 다른 궤적이 양의 확률로 **같은 $\nu$ 표본으로 coupling** → 그 시점부터 두 궤적은 같은 분포. 초기 상태의 기억 리셋.
- **결과 (uniform ergodicity)**: coupling 강제로 모든 궤적이 결국 **같은 정상 분포 $\pi^*$로 수렴**. 출발 무관한 균등 에르고딕성.

이게 HST 증명에서 **"유일한 정상 분포 $\pi^*$가 존재하고, 그 분포가 초기 신호 $\mathbf{y}$에 의존하지 않는다"** 를 보이는 핵심 무기.

:::

::: {.callout-important collapse="true" title="해석: 증명의 큰 그림"}

**문제.** 출발 상태 $\boldsymbol{\delta}_0$ (높이차 벡터)에서 $K^*$ step 후 도달 분포가 모든 target $\boldsymbol{\delta}^* \in \mathcal{B}_\Phi$를 양의 밀도로 cover하는 걸 보여야 한다.

**Randomness 원천.** HST에는 두 가지가 있음:

- 이산 randomness: block 직후 Fall에서 노드 선택 $X_s \sim \mu_0$ (i.i.d.)
- **연속 randomness**: 눈 증분 $b'_s \sim \text{Unif}(0,b)$ (i.i.d.) ← Lebesgue 밀도 원천

연속 밀도는 오로지 **block 직후 Fall에서 떨어지는 첫눈** 에서 나온다 (flow step의 $b'_s$는 결정론적 노드 이동과 합쳐져 noise처럼 작동).

**전략 4단계.**

1. **충분히 긴 window**: $K = n(T_{\max}{+}1)$ step. 한 round 길이 $\leq T_{\max}{+}1$이므로 이 안에 block(막힘)이 최소 $n$번 발생.
2. **모든 노드에 첫눈이 한 번씩 떨어지는 사건 $\mathcal{A}_1$**: 그 $n$번의 block 직후 첫눈이 정확히 $v_1, v_2, \ldots, v_n$ 순으로 떨어지는 사건. 그러면 모든 노드에 fresh $\text{Unif}(0,b)$ 가 한 번씩 들어가 $n$차원 shift 가능.
3. **Flow step 눈은 작게**: 같은 window 안 flow step의 눈 증분을 $\leq \varepsilon$로 강제 ($\mathcal{A}_2$). 이건 noise만 줄이는 것이고 block 직후 첫눈의 분포는 안 건드림 (independence). $L$ window에 걸쳐 누적해도 residual은 $< Lb/8$.
4. **한 window는 작은 ball, $L$ window는 큰 ball**: 한 window의 첫눈들이 만드는 shift은 small ball $\mathcal{O} = \{\max|y_i| < b/4\}$ 안. $\mathcal{B}_\Phi$ 전체 ($\|\boldsymbol{\delta}\|_\infty \leq D_\Phi$)를 cover하려면 $L \cdot b/4 > 2D_\Phi$ 필요 → $L := \lceil 16 D_\Phi/b + 1 \rceil$.

**그래서 상수 $K^* = LK$, $\varepsilon = b/(8K)$.**

:::

::: {.callout-note collapse="true" title="증명: Doeblin minorization"}

**상수 정의** (위 빨간 callout의 4단계 전략에서 유도).

$$K := n(T_{\max}{+}1),\quad \varepsilon := \tfrac{b}{8K},\quad D_\Phi := \sup_{s \in \mathcal{B}_\Phi} \|\boldsymbol{\delta}(s)\|_\infty,\quad L := \lceil 16 D_\Phi / b + 1 \rceil,\quad K^* := LK$$

**한 window의 핵심 사건 두 개.**

- $\mathcal{A}_1$ — 길이 $K$ window 안 첫 $n$번의 block 직후 첫눈이 차례로 $v_1, v_2, \ldots, v_n$ 에 떨어짐 (모든 노드에 한 번씩 첫눈)
- $\mathcal{A}_2$ — 같은 window의 모든 flow step 눈 증분이 $b'_s \leq \varepsilon$ (noise 억제)

**$P(\mathcal{A}_1) > 0$ 증명.** 한 round 길이가 $\leq T_{\max}{+}1$이므로 $K = n(T_{\max}{+}1)$ step 안에 block 시점이 $\geq n$개 보장됨. 각 block 직후 첫눈의 노드 선택은 $\mu_0$에서 i.i.d.이고 $\mathcal{F}_{s-1}$과 독립.

$$P(\mathcal{A}_1) \geq \mu_{\min}^n > 0 \qquad (\mu_{\min} := \min_i \mu_0(v_i))$$

**$P(\mathcal{A}_2) > 0$ + 첫눈 분포 유지.** 현재 step이 block 직후 첫눈인지 flow step인지는 timer $Z_{s-1}$로 결정되며 $\mathcal{F}_{s-1}$-measurable. $b'_s \perp \mathcal{F}_{s-1}$이므로 flow step에 $\{b'_s \leq \varepsilon\}$ 조건을 걸어도 **block 직후 첫눈의 $b'_s$ marginal은 $\text{Unif}(0,b)$ 그대로**.

$$P(\mathcal{A}_2) \geq (\varepsilon/b)^{n T_{\max}} > 0 \qquad (\because \text{window 내 flow step 수} \leq n T_{\max})$$

**Flow step residual bound** ($L$ window 누적).

$\mathcal{A}_2$ 하에서 $L$ window 동안 flow step 눈 증분 누적은 노드별로:

$$\|\mathbf{R}\|_\infty \leq \underbrace{L n T_{\max}}_{\text{flow step 수}} \cdot \underbrace{\varepsilon}_{=\,b/(8K)} = L n T_{\max} \cdot \tfrac{b}{8 n (T_{\max}{+}1)} < \tfrac{Lb}{8}$$

즉 noise는 $L$ window 합쳐도 $Lb/8$ 이하.

**한 window의 첫눈 shift 분포.** $\mathcal{A}_1$ 하에서 한 window의 첫눈 $n$개는 i.i.d. $\text{Unif}(0,b)$. Pushforward density lemma에 의해, 이것이 만드는 $\boldsymbol{\delta}$-shift는 $\mathcal{O} := \{\max_i |y_i| < b/4\}$ 위에 밀도 $\geq b^{1-n}/2$로 분포.

**$L$ window 합성으로 target $\boldsymbol{\delta}^*$ 도달.** $L$개 독립 window의 shift 합은 Minkowski 합성으로 $L\mathcal{O} = \{\max|y_i| < Lb/4\}$ 위에 밀도 (convolution).

도달해야 할 net shift:

$$\boldsymbol{\delta}^* - \boldsymbol{\delta}_0 - \mathbf{R}$$

크기 bound:

$$\|\boldsymbol{\delta}^* - \boldsymbol{\delta}_0 - \mathbf{R}\|_\infty \leq \underbrace{\|\boldsymbol{\delta}^*\|_\infty}_{\leq D_\Phi} + \underbrace{\|\boldsymbol{\delta}_0\|_\infty}_{\leq D_\Phi} + \underbrace{\|\mathbf{R}\|_\infty}_{< Lb/8} \leq 2D_\Phi + \tfrac{Lb}{8}$$

$L \geq 16D_\Phi/b + 1$이므로:

$$2D_\Phi + \tfrac{Lb}{8} < \tfrac{Lb}{8} + \tfrac{Lb}{8} = \tfrac{Lb}{4}$$

즉 도달해야 할 shift이 $L\mathcal{O}$의 내부 ball 안에 들어감 → block shift 분포로 cover 가능.

**결론.** $\nu := L\mathcal{O}$ 내부 작은 ball 위 normalized Lebesgue. 그러면 임의의 $s_0 \in \mathcal{B}_\Phi$와 Borel $B$에 대해

$$P^{K^*}(s_0, B) \geq \underbrace{[P(\mathcal{A}_1 \cap \mathcal{A}_2)]^L}_{>\,0\,\text{(window 독립)}} \cdot \underbrace{(b^{1-n}/2)^L}_{\text{밀도 lower bound}^L} \cdot \nu(B) =: \eta\,\nu(B)$$

$\eta > 0$은 $s_0$에 무관 ($n, b, T_{\max}, \mu_{\min}, D_\Phi$에만 의존). $\square$

:::


### 증명 체인

::: {.callout-important collapse="true" title="해석: Foster-Lyapunov + Doeblin의 역할 분담"}

증명 체인이 왜 두 기계(Foster drift, Doeblin)를 **나란히** 쓰는지 — 둘은 **상호보완**으로 서로 못 하는 일을 메꿔준다.

**비유**: state space를 큰 그릇이라고 보자.

- **바깥 영역** ($\mathcal{B}_\Phi$ **밖**, $\Phi$ 큰 곳): 그릇 가장자리. 여기 있으면 **Foster drift가 안쪽으로 끌어당김** (Lyapunov가 강하게 감소 → chain이 가장자리에 오래 못 머묾).
- **안쪽 영역** ($\mathcal{B}_\Phi$ **안**, $\Phi$ 작은 곳): 그릇 바닥. 여기서는 drift가 약함. 대신 **Doeblin minorization이 coupling을 강제** (출발 무관 공통 분포로 수렴).

**왜 둘 다 필요한가:**

| 가진 것 | 빠진 것 |
|---|---|
| Foster drift 만 | chain이 $\mathcal{B}_\Phi$ 안에 들어와도 거기서 어떻게 분포가 결정되는지 모름 → 정상 분포 **유일성 보장 안 됨** |
| Doeblin 만 | $\mathcal{B}_\Phi$가 small set이라도 chain이 **거기 도달한다는 보장 없음** → 발산 가능 |
| Foster + Doeblin | 가장자리 → 바닥(Foster) → 모두 coupling(Doeblin) → 유일 $\pi^*$로 수렴 |

**표준 패턴**: Meyn-Tweedie (2009) 책의 핵심 도구. Theorem 11.3.4 (positive Harris recurrence), Theorem 14.3.7 (moment bound), Theorem 17.3.2 (LLN) 모두 "Foster drift + small set" 조합 위에 세워짐.

이래서 증명 체인이 Step 1 (Foster drift) → Step 2 (Doeblin) → 그 다음 Harris recurrence → LLN 순서.

:::

**Step 1 (Round-skeleton drift).** OB $\Rightarrow$ $\mathbb{E}[\Delta\Phi] \leq -b\kappa_G M + C$, $\leq -1$ outside $\mathcal{B}_\Phi$.

**Step 2 (Doeblin).** $\mathcal{B}_\Phi$가 $K^*$-skeleton의 small set.

**Step 3 ($\psi$-irreducibility).** Foster drift로 모든 initial state가 $\mathcal{B}_\Phi$에 도달 (optional stopping: $\mathbb{E}[\tau] \leq \Phi(t_0) < \infty$). Doeblin smallness로 $\nu$-irreducibility. 단일 closed Harris class, $\pi^*$ unique, $\mathbf{y}$-독립.

**Step 4 (Skeleton $\to$ original).** $m_r \leq T_{\max}+1$이므로 hitting time, stationary measure, LLN 모두 round-skeleton에서 original chain으로 전달 (Meyn-Tweedie 17.3.2).

**Step 5 (Moment bound).** $V_{\text{Ly}} = \Phi^2$ iterated Lyapunov로 $\mathbb{E}_{\pi^*}[\Phi] < \infty$.

**Step 6 (LLN).** $g(S) = (\delta_i - \delta_j)^2 \leq 2\Phi$이므로 $\pi^*$-integrable. Positive Harris chain의 LLN (Meyn-Tweedie 17.1.7):

$$\frac{SD^2_{ij}(t)}{t} = \frac{1}{t}\sum_{s=0}^t g(S(s)) \to \mathbb{E}_{\pi^*}[g] =: c_{ij} \qquad \text{a.s.} \qquad \square$$

---

# Open items

`-` OB for $C_n$, 일반 connected graph: $K_n$만 증명됨

`-` Cluster LP $K \geq 3$: $K = 2$ 완전 증명, $K \geq 3$은 $n \leq 40$ computational verification

`-` Drift regime에서 $\rho_i$ 존재: empirical 확인, formal proof open
