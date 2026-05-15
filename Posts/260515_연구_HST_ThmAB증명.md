---
title: 연구 ▷ HST ▷ Thm A, B 증명
author: 신록예찬
date: 05/15/2026
draft: false
---

Snow distance $SD^2_{ij}(t) = \sum_{s=0}^t (h_i(s) - h_j(s))^2$ 의 점근 거동은 장기 적립률 $\rho_i = \lim \frac{1}{t}\sum \mathbf{1}\{X_s = v_i\}$ 에 따라 두 가지로 나뉜다.

`-` **Theorem A (Balanced, $\rho_i = 1/n$):** $SD^2_{ij}(t)/t \to c_{ij}$ a.s. 극한 $c_{ij}$는 그래프와 $b, T_{\max}$에만 의존하고 초기 신호 $\mathbf{y}$에 무관.

`-` **Theorem B (Drift, $\rho_i \neq \rho_j$):** $SD^2_{ij}(t)/t^3 \to \bar{b}^2(\rho_i - \rho_j)^2 / 3$ a.s.

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

### 보조 결과 1: $\Phi$ drift

$\hbar(v_i, t) = h(v_i, t) - \bar{h}(t)$ (centered height), $\Phi(t) = \sum_i \hbar(v_i, t)^2$.

**Lemma ($\Phi$ drift).** 

$$\mathbb{E}[\Phi(t+1) - \Phi(t) \mid \mathcal{F}_t, X_{t+1}] = b\,\hbar(X_{t+1}, t) + \frac{(n-1)b^2}{3n}$$

증명: $\Phi(t+1) - \Phi(t) = 2\beta\hbar(v,t) + (n-1)\beta^2/n$. $\mathbb{E}[\beta] = b/2$, $\mathbb{E}[\beta^2] = b^2/3$을 대입.

### 보조 결과 2: Shape decomposition

Round $r$의 방문 노드 $\tilde{X}_1, \ldots, \tilde{X}_{m_r}$에 대해:

$$\sum_{j=1}^{m_r} \hbar(\tilde{X}_j, t_r + j - 1) = \sum_u N_u \hbar(u, t_r) + \Theta, \qquad |\Theta| \leq \frac{(T_{\max}+1)T_{\max}}{2}b$$

좌변의 $\hbar$는 라운드 중간 시점의 값인데, 라운드 시작 시점의 값 $\hbar(u, t_r)$으로 분해하면 correction $\Theta$가 bounded.

이 두 lemma를 합치면 **라운드 단위 drift**:

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] = b\,\mathcal{B}(S) + (\text{bounded terms})$$

여기서 $\mathcal{B}(S) = \sum_u \pi(u) \hbar(u)$가 **M에 비례하는 유일한 항**이다.

### 보조 결과 3: OB (Occupation-Bias) 조건

**OB가 뭔가?** Walker가 아래로 흘러가므로, 높이가 낮은(음수인) 노드를 더 많이 방문한다. 따라서 $\mathcal{B}(S) = \sum \pi(u)\hbar(u)$는 $M$이 클수록 음수가 된다:

$$\mathcal{B}(S) \leq -\kappa_G M + C_{\text{OB}}$$

이걸 OB 조건이라 부른다. $\kappa_G > 0$이면 $\Phi$ drift가 음이 되어 Lyapunov 조건이 성립.

**$K_n$ (완전 그래프)에서 OB 증명.** 노드를 $\hbar$ 내림차순으로 정렬: $v_{(0)}, \ldots, v_{(n-1)}$.

Well-separated case (인접 rank 간 gap $> \delta = b(T_{\max}+1)$): 라운드 중 순위 역전이 없으므로 rank $j$에서 출발하면 rank $i > j$를 방문할 확률이 정확히 $1/(n-i)$. Expected visit count:

$$\pi(v_{(i)}) = \frac{1}{n} + \frac{i}{n(n-i)} = \frac{1}{n-i} \qquad (0 \leq i \leq n-2)$$

$\pi(v_{(n-1)}) = 1$ (최저 노드는 항상 도달).

Bilinear form: $\mathcal{B} = \sum_{i=0}^{n-1} \frac{\hbar_{(i)}}{n-i}$. **LP** (선형계획법)로 제약 $\sum \hbar = 0$, sorted, range $= M$ 하에서 최대값을 구하면:

$$\mathcal{B} \leq -\frac{H_n - 1}{n} M, \qquad H_n = \sum_{k=1}^n \frac{1}{k}$$

최대점은 one-high profile: $\hbar_{(0)} = (n-1)M/n$, 나머지 $= -M/n$.

Near-tie가 있는 일반 case: 노드를 cluster로 묶고, perturbative contribution $\leq (T_{\max}+1)n\delta$ ($M$과 무관)을 분리한 뒤, clean contribution에 $K=2$ cluster LP ($\kappa \geq 1/n^2$)를 적용. 최종 $\kappa_G = 1/(2n^2)$.

### 보조 결과 4: Moment bound

Round drift $\mathbb{E}[\Delta\Phi] \leq -b\kappa_G M + C$는 $-\Phi^{1/2}$ 차수라서 Meyn-Tweedie로는 $\mathbb{E}_{\pi^*}[\Phi^r] < \infty$ ($r < 1/2$)만 나온다. LLN에는 $r = 1$이 필요.

**해법: Lyapunov를 $V = \Phi^2$로 올린다.**

$$\mathbb{E}[\Delta V \mid S] = 2\Phi \cdot \mathbb{E}[\Delta\Phi] + \mathbb{E}[(\Delta\Phi)^2]$$

`-` 첫째 항: $2\Phi \cdot (-b\kappa_G M) \leq -b\kappa_G M^3$ (pigeonhole: $\Phi \geq M^2/2$)

`-` 둘째 항: $(\Delta\Phi)^2 \leq (A_1 M + A_2)^2 = O(M^2)$

$M^3$ vs $M^2$: cubic이 이기므로 $\mathbb{E}[\Delta V] \leq -cV^{3/4} + C\mathbf{1}_{\mathcal{B}'}$.

Meyn-Tweedie 14.3.7 적용: $\mathbb{E}_{\pi^*}[\Phi^{3/2}] < \infty \Rightarrow \mathbb{E}_{\pi^*}[\Phi] < \infty$.

### 보조 결과 5: Doeblin minorization

Bounded set $\mathcal{B}_\Phi$가 small set임을 보인다. $K = n(T_{\max}+1)$ steps의 window를 잡으면:

1. $\geq n$ 번의 block-flag time이 보장됨 (각 round $\leq T_{\max}+1$ steps)
2. Non-block steps에서 $b'_s \leq \varepsilon$로 conditioning (block/non-block 분류는 pre-step timer $Z_{s-1}$로 결정, $b'_s$와 독립이므로 block increments는 여전히 $\text{Unif}(0,b)$)
3. Block increments의 pushforward가 Lebesgue density를 줌 (fibre length $> b/2$)
4. $L$ windows 합성으로 target coverage

### 증명 체인

**Step 1 (Round-skeleton drift).** OB $\Rightarrow$ $\mathbb{E}[\Delta\Phi] \leq -b\kappa_G M + C$, $\leq -1$ outside $\mathcal{B}_\Phi$.

**Step 2 (Doeblin).** $\mathcal{B}_\Phi$가 $K^*$-skeleton의 small set.

**Step 3 ($\psi$-irreducibility).** Foster drift로 모든 initial state가 $\mathcal{B}_\Phi$에 도달 (optional stopping: $\mathbb{E}[\tau] \leq \Phi(t_0) < \infty$). Doeblin smallness로 $\nu$-irreducibility. 단일 closed Harris class, $\pi^*$ unique, $\mathbf{y}$-독립.

**Step 4 (Skeleton $\to$ original).** $m_r \leq T_{\max}+1$이므로 hitting time, stationary measure, LLN 모두 round-skeleton에서 original chain으로 전달 (Meyn-Tweedie 17.3.2).

**Step 5 (Moment bound).** $V = \Phi^2$ iterated Lyapunov로 $\mathbb{E}_{\pi^*}[\Phi] < \infty$.

**Step 6 (LLN).** $g(S) = (\delta_i - \delta_j)^2 \leq 2\Phi$이므로 $\pi^*$-integrable. Positive Harris chain의 LLN (Meyn-Tweedie 17.1.7):

$$\frac{SD^2_{ij}(t)}{t} = \frac{1}{t}\sum_{s=0}^t g(S(s)) \to \mathbb{E}_{\pi^*}[g] =: c_{ij} \qquad \text{a.s.} \qquad \square$$

---

# Open items

`-` OB for $C_n$, 일반 connected graph: $K_n$만 증명됨

`-` Cluster LP $K \geq 3$: $K = 2$ 완전 증명, $K \geq 3$은 $n \leq 40$ computational verification

`-` Drift regime에서 $\rho_i$ 존재: empirical 확인, formal proof open
