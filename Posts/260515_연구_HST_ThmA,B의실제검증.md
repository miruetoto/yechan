---
title: 연구 ▷ HST ▷ Thm A, B의 실제 검증
author: 소미
date: 05/15/2026
draft: false
output-file: 260515_41f2b0.html
---

각 노드의 장기 적립률을 $\rho_i = \lim_{t\to\infty} \frac{1}{t}\sum_{s=1}^t \mathbf{1}\{X_s = v_i\}$ 라 하자. $\rho$가 균등한지 여부에 따라 SD의 스케일링이 달라진다.

`-` **Thm A — Balanced** ($\rho_i = 1/n$): $SD^2/t \to c$ (유한 상수). 극한은 그래프-신호 상호작용을 반영.

`-` **Thm B — Drift** ($\rho_i \neq \rho_j$): $SD^2/t \to \infty$. 대신 $\frac{SD^2_{ij}(t)}{t^3} \to \frac{\bar{b}^2(\rho_i - \rho_j)^2}{3}$.

본 글은 두 부분으로 나뉜다. **A 그래프 구조 검증** ($f=0$, 24개 예제) 에서는 그래프 모양·엣지 방향이 regime 분류에 미치는 영향을, **B 신호 변형 검증** ($f \neq 0$, 3개 예제) 에서는 초기 신호의 효과를 본다.

## A. 그래프 구조 검증 ($f = 0$)

[「그래프 도메인에서의 거리」](./260516_연구_HST_그래프도메인에서의거리.html) 에서 정의한 24개 테스트 그래프 (7 모양 × 변형) 에 대해 $f = 0$ 으로 두고 HST 시뮬을 수행했다. 예제 네이밍은 거리 블로그와 동일 (1-1 ~ 7-2). 시뮬 설정: $b = 0.05$, $T_{\max} = 20$, $\tau = 10^7$, seed=42, **모델 B random-step variant**, degree-prop $\mu_0$.

`-` 모양 1 — Star ($n=13$): 1-1 Inward / 1-2 Outward / 1-3 Bidir / 1-4 Gaussian
`-` 모양 2 — Wheel ($n=11$): 4변형
`-` 모양 3 — Helm ($n=21$): 4변형
`-` 모양 4 — Double Helm ($n=31$): 4변형
`-` 모양 5 — Path ($n=60$): 5-1 Directed / 5-2 Bidir / 5-3 Gaussian
`-` 모양 6 — Ring ($n=60$): 3변형
`-` 모양 7 — Cylinder ($n=60$): 7-1 Inter-dir / 7-2 Inter-bidir

각 예제마다 3패널 애니메이션: (그래프 시그널 3D | $\sqrt{SD^2/t}$ 3D MDS embedding | $SD^2/t$ log-log 수렴 곡선).

### A.1 Star Inward (1-1)

hub 단방향 inward. degree-prop $\mu_0$ 가 hub 에 집중되지만 inward 방향 때문에 flow 가 hub 에서 멈춤. $\hat c \approx 0.31$ — **약발산** (Thm A 경계).

![](attachments/260515_41f2b0_10.gif)

### A.2 Star Outward (1-2)

hub outward. flow 가 leaf 로 빠지면서 $\rho_\text{hub} \gg \rho_\text{leaf}$. $\hat c \approx 6.4 \times 10^7$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_11.gif)

### A.3 Star Bidir (1-3)

양방향. degree-prop $\mu_0$ 가 hub 에 0.5 집중. $\hat c \approx 1.6 \times 10^7$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_12.gif)

### A.4 Star Gaussian (1-4)

Gaussian kernel 가중치. hub-leaf 격차 완화. $\hat c \approx 0.011$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_13.gif)

### A.5 Wheel Inward (2-1)

Wheel = Star + outer cycle. outer cycle 덕분에 hub 우회 가능. $\hat c \approx 0.016$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_14.gif)

### A.6 Wheel Outward (2-2)

hub outward 이지만 outer cycle 이 격차 완화. $\hat c \approx 0.021$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_15.gif)

### A.7 Wheel Bidir (2-3)

양방향. $\hat c \approx 0.093$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_16.gif)

### A.8 Wheel Gaussian (2-4)

Gaussian. $\hat c \approx 0.011$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_17.gif)

### A.9 Helm Inward (3-1)

Helm = Wheel + pendant leaves. hub-ring-pendant 3계층 격차. inward 방향. $\hat c \approx 7.6 \times 10^6$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_18.gif)

### A.10 Helm Outward (3-2)

hub outward. $\hat c \approx 5.4 \times 10^6$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_19.gif)

### A.11 Helm Bidir (3-3)

양방향. $\hat c \approx 2.5 \times 10^6$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_20.gif)

### A.12 Helm Gaussian (3-4)

Gaussian kernel. hub-pendant 큰 격차는 Gaussian 으로도 완화 부족. $\hat c \approx 0.40$ — **약발산**.

![](attachments/260515_41f2b0_21.gif)

### A.13 D-Helm Inward (4-1)

Double Helm = Helm + outer ring. hub-inner-pendant-outer 4계층. $\hat c \approx 1.9 \times 10^5$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_22.gif)

### A.14 D-Helm Outward (4-2)

hub outward. $\hat c \approx 2.0 \times 10^5$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_23.gif)

### A.15 D-Helm Bidir (4-3)

양방향. 격차 완화되지만 여전히 큼. $\hat c \approx 1.2 \times 10^4$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_24.gif)

### A.16 D-Helm Gaussian (4-4)

Gaussian kernel. 4계층 격차가 너무 커 Gaussian 도 완화 못함. $\hat c \approx 2.4 \times 10^5$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_25.gif)

### A.17 Dir Path (5-1)

단방향 path. 양 끝 (degree 1) vs 내부 (degree 2). $\hat c \approx 768$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_26.gif)

### A.18 Bidir Path (5-2)

양방향. 양 끝 lighter. $\hat c \approx 494$ — **drift (Thm B)**.

![](attachments/260515_41f2b0_27.gif)

### A.19 Gauss Path (5-3)

Gaussian kernel. 양 끝 격차 완화. $\hat c \approx 1.4$ — **약발산**.

![](attachments/260515_41f2b0_28.gif)

### A.20 Dir Ring (6-1)

단방향 ring (shift matrix). doubly stochastic → $\rho_i = 1/n$. $\hat c \approx 0.010$ — **balanced (Thm A)**, 비대칭에서도 성립.

![](attachments/260515_41f2b0_29.gif)

### A.21 Bidir Ring (6-2)

양방향 ring. 정규. $\hat c \approx 0.18$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_30.gif)

### A.22 Gauss Ring (6-3)

Gaussian kernel ring. $\hat c \approx 0.014$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_31.gif)

### A.23 Cylinder Inter-dir (7-1)

5 ring × 12 node, ring 사이 단방향 결합. ring별 적립률 격차. $\hat c \approx 8.0$ — **약발산**.

![](attachments/260515_41f2b0_32.gif)

### A.24 Cylinder Inter-bidir (7-2)

5 ring × 12 node, ring 사이 양방향. 정규. $\hat c \approx 0.029$ — **balanced (Thm A)**.

![](attachments/260515_41f2b0_33.gif)

### A 정리

`-` **Wheel 의 4 변형 모두 balanced**: outer cycle 이 차수 격차를 완화해서 $\rho_i \approx 1/n$.

`-` **Star/Helm/D-Helm 의 directed/bidir 변형은 모두 drift**: hub 가 reset 을 집중적으로 받아 $\rho_\text{hub} \gg \rho_\text{leaf}$.

`-` **6-1 Directed Ring**: doubly stochastic 이라 비대칭임에도 $\rho_i = 1/n$ → balanced. 「비대칭 그래프에서도 Thm A」 라는 비자명한 결론.

`-` **Gaussian kernel** 은 정규 그래프(Wheel, Ring)와 작은 격차(Star)에서 잘 작동하지만, hub-pendant 격차가 큰 Helm/D-Helm 에서는 완화 부족.

`-` **Cylinder 7-1 (inter-dir)** 은 ring 사이 단방향 결합으로 ring별 적립률 격차 → 약발산. **7-2 (inter-bidir)** 은 정규 → balanced.

## B. 신호 변형 검증 ($f \neq 0$)

§A의 그래프 중 일부 위에 초기 신호를 얹어 Thm A 가 신호-그래프 상호작용을 어떻게 반영하는지 본다. 모두 balanced regime 그래프 위에서.

### B.1 Parity Cycle $C_{60}$ — $f_i = (-1)^i$

§A 6-2 Bidir Ring 과 동일한 양방향 cycle. 신호는 Nyquist 주파수. 이웃끼리 신호가 반대이므로 block 이 자주 발생한다.

![](attachments/260515_41f2b0_04.gif)

### B.2 Directed Cycle $C_{60}$ — $f = \pm 1$

§A 6-1 Dir Ring 위에 반원 경계 $\pm 1$ 신호. $f \neq 0$ 이면 매 시점 $(h_i - h_j)^2$ 에 초기 신호 차이 $(f_i - f_j)^2$ 가 상수항으로 누적되어 $SD^2/t$ 가 큰 값에서 수렴 (scaling 자체는 Thm A 그대로).

![](attachments/260515_41f2b0_07.gif)

### B.3 Outlier Cylinder

§A 7-2 Cyl Inter-bidir 와 비슷한 cylinder 구조 (Gaussian kernel, 위=-3, 아래=+3) 에서 +3 그룹 정중앙 한 노드만 -3 으로 flip. 정규 → balanced. 국소적 outlier 가 SD 임베딩에서 어떻게 분리되는지.

![](attachments/260515_41f2b0_08.gif)

---

결국 두 regime 모두 $\tau \to \infty$ 에서 초기 신호 $\mathbf{y}$ 는 씻겨나간다. 차이는 **어떤 그래프 정보가 남느냐**: balanced 는 flow/block dynamics 의 세밀한 구조, drift 는 단순히 적립률 격차.

---

## Appendix: random-step variant 모델 정의 민감도

본문의 실제 검증은 random-step variant의 **모델 A** (매 step별 $b'_s \sim \text{Unif}(0,b)$ 추첨, $\mathbb{E}[b'] = b/2$) 가정 위에서 수행되었다. 그런데 "한 라운드 동안 $b'$이 step 수만큼인지, 라운드당 1개인지"에 대한 정의가 코드와 직관 사이에 모호할 수 있다. 후보 모델을 다음과 같이 정리한다.

| 모델 | $b'$ 추첨 시점 | 라운드 내 적립 |
|---|---|---|
| **A** (코드/이론) | 매 step | $b'_s$ 매번 새로 추첨 |
| **B** (round-const) | 라운드 시작 시 1회 | 라운드 동안 같은 $b'$ |
| **C** (감쇠) | 라운드 시작 시 1회 | step $k$에서 $b' \cdot r^k$ ($r<1$) |

추가예제 6 (Directed Cycle $C_{60}$, $f=0$, $b=0.05$, $T_{\max}=20$, $\tau=50\text{M}$, seed=42) 으로 네 가지 모델을 동일 조건에서 비교했다.

![](attachments/260515_41f2b0_09.png)

| 모델 | 정의 | $SD^2/t$ (pair 0,30) | $SD^2/t^2$ | 비고 |
|---|---|---|---|---|
| A | per-step iid | **0.0135** | $2.7\times 10^{-10}$ | 본문 기준 |
| B | round-const | **0.0155** | $3.1\times 10^{-10}$ | +15% vs A |
| C, $r=0.9$ | 점진 감쇠 | **0.0109** | $2.2\times 10^{-10}$ | $-19\%$ vs A |
| C, $r=0.5$ | 급격 감쇠 | **0.1214** | $2.4\times 10^{-9}$ | 큰 폭 변동 |

### 관찰

1. **Thm A 스케일링은 모델 정의 무관하게 성립**한다. 네 모델 모두 $SD^2/t$는 일정값으로 수렴하고 $SD^2/t^2$는 $0$으로 감소한다. 따라서 추가예제 6의 정성적 결론 ("비대칭 그래프에서도 Thm A 성립") 은 모델 선택에 robust 하다.

2. **$c$ 값의 절댓값은 모델에 민감**하다. A 대비 B는 $+15\%$, C($r=0.9$)는 $-19\%$, C($r=0.5$)는 $\sim 9$배 증폭. 모델 C 의 급감쇠($r=0.5$)는 매 라운드의 첫 step에 적립이 집중되어 노드별 높이 격차가 크게 벌어지는 효과 때문이다.

3. 따라서 본문의 $\hat c \approx 0.0136$ 같은 **수치**는 모델 A 가정에 종속적이지만, "수렴한다"는 **정성**은 모델 정의의 자유도에 영향받지 않는다.

### 비고

이 부록은 모델 A/B/C 가운데 어느 것이 "맞는" 정의인지 결정하지 않는다. 그 판단은 random-step variant 의 물리적 직관, 기존 증명들의 호환성, 그리고 다른 그래프 (Helm, Wheel, Star) 에서의 동일 비교가 필요하다.
