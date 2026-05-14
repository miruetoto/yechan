---
title: 연구 ▷ HST ▷ Snow 거리는 발산할 수도 있다
author: 신록예찬
date: 05/15/2026
draft: false
output-file: 260515_c0531a.html
---

Snow distance $SD^2_{ij}(t) = \sum_{s=0}^t (h_i(s) - h_j(s))^2$ 가 항상 수렴하는 건 아니다. 초기분포 $\mu_0$에 따라 두 가지 regime이 존재한다.

각 노드의 장기 적립률을 $\rho_i = \lim_{t\to\infty} \frac{1}{t}\sum_{s=1}^t \mathbf{1}\{X_s = v_i\}$ 라 하자. 모든 노드에 눈이 균등하게 쌓이면 ($\rho_i = 1/n$) **balanced**, 아니면 **drift**.

`-` **Balanced** ($\rho_i = 1/n$): $SD^2/t$ 가 유한 상수로 수렴. 극한은 그래프-신호 상호작용을 반영.

`-` **Drift** ($\rho_i \neq \rho_j$): $SD^2/t$ 는 발산. 대신 $SD^2/t^3$ 이 수렴한다:

$$\frac{SD^2_{ij}(t)}{t^3} \;\to\; \frac{\bar{b}^2(\rho_i - \rho_j)^2}{3}$$

**Drift regime 증명.** 노드 $v_i$의 누적 높이는

$$h(v_i, t) = h(v_i, 0) + \sum_{s=1}^t b'_s \cdot \mathbf{1}\{X_s = v_i\}$$

여기서 $b'_s \sim \mathrm{Unif}(0, b)$는 매 스텝의 눈량이고 $X_s$는 눈이 쌓이는 노드이다. $b'_s$와 $X_s$는 독립이므로 합을 분리할 수 있다:

$$\sum_{s=1}^t b'_s \cdot \mathbf{1}\{X_s = v_i\} = \bar{b}\sum_{s=1}^t \mathbf{1}\{X_s = v_i\} + \sum_{s=1}^t (b'_s - \bar{b})\cdot\mathbf{1}\{X_s = v_i\}$$

첫째 항: $\rho_i$의 정의에 의해 $\frac{1}{t}\sum \mathbf{1}\{X_s = v_i\} \to \rho_i$ (a.s.)이므로 $\bar{b}\rho_i \cdot t$ 로 성장.

둘째 항은 마팅게일 강대수법칙으로 처리한다. $M_t^{(i)} := \sum_{s=1}^t (b'_s - \bar{b}) \cdot \mathbf{1}\{X_s = v_i\}$ 라 놓으면, $b'_s - \bar{b}$는 평균 0이고 $\mathcal{F}_{s-1}$ (시점 $s-1$까지의 정보)과 독립이므로

$$\mathbb{E}[M_t^{(i)} - M_{t-1}^{(i)} \mid \mathcal{F}_{t-1}] = \mathbb{E}[(b'_t - \bar{b}) \mid \mathcal{F}_{t-1}] \cdot \mathbf{1}\{X_t = v_i\} = 0$$

즉 $\{M_t^{(i)}\}$는 마팅게일이다. 마팅게일 강대수법칙은 **각 증분의 분산이 bounded이면 $M_t/t \to 0$ a.s.** 를 말한다. 정확히는:

> $\{M_t\}$가 마팅게일이고 $\sum_{t=1}^\infty \frac{\mathrm{Var}(M_t - M_{t-1})}{t^2} < \infty$ 이면 $\frac{M_t}{t} \to 0$ a.s.

우리 경우 $\mathrm{Var}(M_t^{(i)} - M_{t-1}^{(i)}) = \mathrm{Var}((b'_t - \bar{b}) \cdot \mathbf{1}\{X_t = v_i\}) \leq \mathrm{Var}(b'_t) = b^2/12$ 로 상수에 bounded. $\sum_{t=1}^\infty \frac{b^2/12}{t^2} < \infty$ ($p$-급수, $p=2$) 이므로 조건이 만족되어 $M_t^{(i)}/t \to 0$ a.s., 즉 둘째 항은 $o(t)$.

따라서 두 노드의 높이 차이는

$$h(v_i, t) - h(v_j, t) = \bar{b}(\rho_i - \rho_j)\cdot t + o(t)$$

이걸 제곱해서 더하면:

$$SD^2_{ij}(t) = \sum_{s=0}^t \left[\bar{b}(\rho_i - \rho_j)\cdot s + o(s)\right]^2 = \bar{b}^2(\rho_i - \rho_j)^2 \sum_{s=0}^t s^2 + \text{(작은 항들)}$$

주항에서 $\sum_{s=0}^t s^2 = \frac{t(t+1)(2t+1)}{6} \sim \frac{t^3}{3}$ 이고, 나머지 교차항 $\sum s \cdot o(s)$과 잔차항 $\sum o(s)^2$은 모두 $o(t^3)$이다. $t^3$으로 나누면:

$$\frac{SD^2_{ij}(t)}{t^3} \;\to\; \frac{\bar{b}^2(\rho_i - \rho_j)^2}{3} \qquad \square$$

---

regime을 결정하는 건 $\mu_0$이다. 정규 그래프에서는 degree-proportional이든 uniform이든 $\rho_i = 1/n$이라 balanced. 비정규 그래프에서 degree-proportional $\mu_0$를 쓰면 높은 차수 노드에 눈이 더 자주 쌓여서 drift에 빠진다.

아래 Star 그래프로 검증해보았다.

![](attachments/260515_c0531a_05.png)

hub의 degree는 12, leaf는 1이다. degree-proportional $\mu_0$를 쓰면 $\mu_0(\text{hub}) = 12/24 = 0.5$로 hub에 reset이 집중된다. 실측 적립률은 $\rho_{\text{hub}} = 0.33$, $\rho_{\text{leaf}} = 0.056$ ($\mu_0$와 정확히 같지는 않다 — block/flow dynamics가 재분배하기 때문). $\rho$가 균등하지 않으므로 drift regime. $SD^2/t^3$이 이론값(점선)으로 수렴하는 과정:

![](attachments/260515_c0531a_04.gif)

결국 두 regime 모두 $\tau \to \infty$에서 초기 신호 $\mathbf{y}$는 씻겨나간다. 차이는 **어떤 그래프 정보가 남느냐**: balanced는 flow/block dynamics의 세밀한 구조, drift는 단순히 적립률 격차.
