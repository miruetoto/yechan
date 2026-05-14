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

증명은 간단하다. $\rho_i \neq \rho_j$이면 높이 차이가 선형으로 벌어진다: $h_i(t) - h_j(t) \approx \bar{b}(\rho_i - \rho_j) \cdot t$. 이걸 제곱해서 더하면 $\sum s^2 \sim t^3/3$. 끝. 강대수법칙 하나로 된다.

regime을 결정하는 건 $\mu_0$이다. 정규 그래프에서는 degree-proportional이든 uniform이든 $\rho_i = 1/n$이라 balanced. 비정규 그래프에서 degree-proportional $\mu_0$를 쓰면 높은 차수 노드에 눈이 더 자주 쌓여서 drift에 빠진다.

Star $S_{13}$에서 degree-prop $\mu_0$로 검증. $SD^2/t^3$이 이론값(점선)으로 수렴하는 과정:

![](attachments/260515_c0531a_04.gif)

결국 두 regime 모두 $\tau \to \infty$에서 초기 신호 $\mathbf{y}$는 씻겨나간다. 차이는 **어떤 그래프 정보가 남느냐**: balanced는 flow/block dynamics의 세밀한 구조, drift는 단순히 적립률 격차.