---
title: 연구 ▷ HST ▷ Thm A, B 증명 자세히 따라가기
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

**Proposition (Doeblin minorization).** 임의의 bounded set $\mathcal{B}_\Phi$에 대해 $\eta > 0$, $K^* < \infty$, 확률측도 $\nu$가 존재하여 $P^{K^*}(s,B) \geq \eta\,\nu(B)$.

::: {.callout-important collapse="true" title="해석: Doeblin minorization이 뭘 말하는가"}

$\mathcal{B}_\Phi$ 안의 아무 상태에서 출발해도, $K^*$ step 후에는 **어디든 양의 확률로 갈 수 있다.**

HST의 상태는 연속 ($\mathbb{R}^{n-1}$의 높이 벡터)이라서, 특정 상태에서 특정 상태로 갈 확률이 0일 수 있다. 하지만 block event에서 i.i.d. $\text{Unif}(0,b)$가 들어오므로 연속 density가 생긴다. 이를 이용하면 bounded 영역 안에서는 어디든 양의 확률로 도달 가능.

Moment bound (보조 결과 4)에서 $\{M \leq M_0\}$가 petite set이라고 쓴 근거가 이것이다.

:::


### 증명 체인

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
