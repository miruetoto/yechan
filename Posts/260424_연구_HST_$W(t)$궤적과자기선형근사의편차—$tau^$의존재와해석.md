---
title: 연구 ▷ HST ▷ $W(t)$ 궤적과 자기 선형근사의 편차 — $\tau^*$의 존재와 해석
author: 신록예찬
date: 04/24/2026
draft: false
output-file: 260424_09cf43.html
---

# 0. 동기: HST는 자기 자신의 선형근사에서 얼마나 벗어나 있는가

[보충아이디어 정리](260208_21543c.html)의 Direction A.7에서 제기한 문제를 본격적으로 전개한다.

Heavy-snow weight 행렬 $W^{\mathrm{HST}}(t)$는 $t = 0$에서 순수 유클리드 정보, $t \to \infty$에서 순수 그래프 정보(Harris 극한 $W^{G}$)를 반영한다. 중간 $t$에서는 두 정보가 섞여 있는데, 이 섞임이 **자명하게 선형**이라면 HST의 고유 기여는 없을 것이다. 본 포스트는 다음을 묻는다:

> **HST의 $W(t)$는 두 endpoint $W^E$와 $W^G$를 잇는 자명한 선형결합과 얼마나 다른가?**

본 포스트의 전략: 같은 두 endpoint를 잇는 **두 궤적** — HST path와 선형결합 path — 을 구성하고, Frobenius 편차 $\Delta(t)$를 분석한다. 편차가 최대가 되는 $\tau^*$는 HST의 **비선형 섞임이 최대로 표출되는 시점**이며, 이 시점의 성질이 "HST가 단순 가중합이 아닌 이유"의 정량적 특성화가 된다.

**본 포스트는 Internal 관점.** $W^G$를 **HST 자신의 Harris 극한** $\lim_t W^{\mathrm{HST}}(t)$로 잡는다. 이 선택은 개념적으로 깨끗하다: HST 궤적이 **자기 자신의 선형근사** 대비 얼마나 휘어있는가를 묻는 self-consistency 체크. 두 endpoint가 동일하므로 $\Delta(t)$가 양 끝에서 0으로 수렴하여 내부 최댓값이 깨끗하게 존재.

**구성.** §1에서 notation. §2에서 $\Delta(t)$와 $\tau^*$의 기본 성질. §3에서 $\tau^*$ 존재 정리와 해석. §4에서 matrix path의 위상 — 둘러싸는 면적을 HST의 대역 비선형성 지표로. §5는 ring·bi-cluster·MCU·random graph 등에서의 실험 계획. §6은 열린 문제와 결론.

> **External 버전.** HST를 **외재적** classical graph metric (diffusion distance, effective resistance 등)과 비교하는 관점은 별도 포스트 [HST와 외재 Graph Distance의 편차](260424_ext.html)에서 다룬다. 본 포스트는 순수 Internal 분석.

---

# 1. Setup

## 1.1 공통 기호

$\mathcal{G} = (\mathcal{V}, \mathcal{E}, \mathbf{W})$: 연결 가중 그래프, $n = |\mathcal{V}|$.
$\mathbf{y} \in \mathbb{R}^n$: 노드 신호.
$b > 0, T_{\max} \geq 1$: HST 파라미터.
$\{X_t\}$: Algorithm 1이 생성하는 walker 위치.
$h(v, t)$: 시각 $t$에서의 snow height.

논문 Definition 2.2에 따른 snow weight 행렬:

$$W^{\mathrm{HST}}_{ij}(t) := \begin{cases} \exp\!\left(-\dfrac{SD_{ij}^2(t)}{2\theta^2 t}\right) & t > 0,\; i \neq j \\ \exp\!\left(-\dfrac{(y(v_i)-y(v_j))^2}{2\theta^2}\right) & t = 0,\; i \neq j \\ 0 & i = j \end{cases}$$

## 1.2 두 극한: $W^E$와 $W^G$

**유클리드 극한 (출발점).** $t = 0$에서 snow distance가 $|y(v_i) - y(v_j)|$로 퇴화:

$$W^E_{ij} := W^{\mathrm{HST}}_{ij}(0) = \exp\!\left(-\frac{(y(v_i) - y(v_j))^2}{2\theta^2}\right), \qquad i \neq j$$

**그래프 극한 (도착점) — Harris 극한.** [*260407 Main Theorem*](260407_eecd78.html)의 결과:

$$W^G_{ij} \;:=\; \lim_{t \to \infty} W^{\mathrm{HST}}_{ij}(t) \;=\; \exp\!\left(-\frac{c_{ij}}{2\theta^2}\right)$$

여기서 $c_{ij} = \lim_t \overline{SD}^2_{ij}(t) = \mathbb{E}_\pi[(\delta_i - \delta_j)^2]$. drift 조건 $(\mathrm{DC}): \Delta_{\mathrm{drift}} < 0$ 하에서 존재하며 a.s. 수렴. $W^G$는 $\mathbf{y}$의 거시적 스케일에 독립이고 $\mathbf{y} \bmod b$에만 국소적으로 종속.

**요점.** $W^E$와 $W^G$는 **모두 HST 자신이 결정**하는 양. HST 궤적 $W^{\mathrm{HST}}(t)$는 이 두 점을 $t: 0 \to \infty$에 걸쳐 잇는다.

## 1.3 Linear combination path (자기 선형근사)

같은 두 endpoint를 잇는 **자명한 mixing**:

$$W^{\mathrm{lin}}(\alpha) := (1-\alpha)\, W^E + \alpha\, W^G, \qquad \alpha \in [0, 1]$$

- $\alpha = 0$: $W^E$ (HST의 $t = 0$과 일치)
- $\alpha = 1$: $W^G$ (HST의 $t = \infty$와 일치)

이 path는 HST의 두 endpoint를 직선으로 잇는, **HST 자신의 "선형근사 buddy"**이다. 자기가 만든 start/end를 자기가 직선으로 이으면 원래 궤적과 비교할 수 있다.

## 1.4 Reparametrization $\alpha(t)$

두 path의 편차를 비교하려면 HST의 "진행도" $t \in [0, \infty]$를 선형 path의 "진행도" $\alpha \in [0, 1]$로 매칭시켜야 한다.

**Endpoint-calibrated 선택 (본 포스트의 기본):**

$$\alpha(t) \;:=\; \frac{\|W^{\mathrm{HST}}(t) - W^E\|_F}{\|W^G - W^E\|_F}$$

이 선택은 HST가 "이미 이동한 거리"를 비율로 측정하여 같은 진행도에서 비교할 수 있게 한다. 다음 성질을 만족:

- $\alpha(0) = 0$, $\alpha(t) \to 1$ as $t \to \infty$
- 단조증가 (HST의 snow distance 단조성에서 유도)
- $[0, \infty] \to [0, 1]$ 연속 bijection

**대안 (비교용):**
- **단순 시간 기반**: $\alpha(t) = t/(t + t_0^*)$, 여기서 $t_0^*$는 §A.2 (보충아이디어)의 위상적 전환 시각.
- **지수형**: $\alpha(t) = 1 - e^{-t/\tau_{\mathrm{mix}}}$, 혼합 시간 $\tau_{\mathrm{mix}}$.

본 포스트의 기본 분석은 endpoint-calibrated을 사용한다.

---

# 2. $\Delta(t)$와 $\tau^*$의 기본 성질

## 2.1 정의

$$\Delta(t) \;:=\; \big\| W^{\mathrm{HST}}(t) - W^{\mathrm{lin}}(\alpha(t)) \big\|_F$$

$$\tau^* \;:=\; \arg\max_{t \geq 0} \Delta(t)$$

(argmax가 유일하지 않을 경우 최소의 $t$로 관습.)

**해석.** $\Delta(t)$는 HST가 자기 자신의 직선 근사에서 시각 $t$에 얼마나 떨어져 있는지. $\tau^*$는 이 편차가 최대가 되는 "가장 휘어있는" 시점.

## 2.2 경계값

**$t = 0$에서:** $W^{\mathrm{HST}}(0) = W^E$이고 $\alpha(0) = 0$이므로 $W^{\mathrm{lin}}(\alpha(0)) = W^E$. 따라서

$$\Delta(0) = 0$$

**$t \to \infty$에서:** Harris 극한 하에서 $W^{\mathrm{HST}}(t) \to W^G$이고 $\alpha(t) \to 1$이므로 $W^{\mathrm{lin}}(\alpha(t)) \to W^G$. 따라서

$$\Delta(t) \xrightarrow{t \to \infty} 0$$

**깔끔한 대칭.** 양 끝에서 정확히 0으로 수렴하므로, $\Delta$가 비자명하면 반드시 내부 최댓값이 존재. 이것이 Internal 관점의 이론적 이점.

## 2.3 연속성

$SD^2_{ij}(t)$는 $t$의 결정론적 선형 증분($+b$)으로만 변하므로, 그 지수변환인 $W^{\mathrm{HST}}_{ij}(t)$도 이산 $t$를 적절히 보간하면 연속. $\alpha(t)$ 역시 연속. 따라서:

**Lemma 2.1.** $\Delta(t)$는 $[0, \infty]$에서 연속이다. $\square$

## 2.4 비자명성

**Proposition 2.2 (비자명 조건).** 다음이 동치:

(a) 모든 $t \in [0, \infty]$에서 $\Delta(t) = 0$, 즉 HST path가 정확히 직선.

(b) $W^{\mathrm{HST}}(t) = (1 - \alpha(t)) W^E + \alpha(t) W^G$ for all $t$.

**언제 (a)가 성립하는가?** 매우 예외적. 예: 모든 snow distance 성분 $SD^2_{ij}(t)$가 시간에 대해 "지수 함수 평균 하에서 정확히 선형으로 진행"하는 경우. 이는 일반적인 그래프·신호 조합에서는 깨진다 (HST의 flow-block 동역학이 본질적으로 비선형). 따라서 일반적 상황에서 (a)는 성립하지 않고, $\Delta(t)$는 비자명하다.

**엄밀한 비자명성 증명은 열린 문제** (§6 O2).

---

# 3. $\tau^*$의 존재 정리와 해석

## 3.1 존재 정리

**Theorem 3.1 ($\tau^*$의 존재).** 연결 그래프 $\mathcal{G}$가 drift 조건 (DC)를 만족하고, $W^E \neq W^G$ (즉 $\mathbf{y}$의 위상 효과가 그래프 극한을 바꾼다)이며, $W^{\mathrm{HST}}$ path가 직선이 아니다 (Proposition 2.2의 (a)가 거짓)라 하자. 그러면:

(i) $\Delta(t)$는 $[0, \infty]$에서 연속이고 $\Delta(0) = \Delta(\infty) = 0$.

(ii) $\max_t \Delta(t) > 0$이고 내부 최댓값 $\tau^* \in (0, \infty)$가 존재한다.

(iii) $\Delta(\tau^*) > 0$이 **HST의 내재적 비선형성 지표**이다.

**증명 스케치.** (i)은 §2.2, §2.3. (ii)는 비퇴화 가정 + 경계 0 + Weierstrass 최댓값 정리. (iii)은 정의. $\square$

## 3.2 $\tau^*$의 해석

$\tau^*$에서의 의미를 여러 각도로:

**물리적 그림 (눈 쌓임의 feedback).** HST는 "눈이 쌓이며 지형이 변하고, 변한 지형이 다시 눈의 흐름을 바꾸는" feedback loop를 가진다. 이 feedback이 없다면 $W(t)$는 단순 보간. Feedback이 가장 강하게 작용하는 시점이 $\tau^*$다.

**기계학습적 그림 (feature engineering).** HST의 "비선형 mixing"이 새 feature를 가장 많이 창출하는 시점. 논문의 cliff detection (Fig 5), hierarchical clustering separation (§5) 같은 HST 고유 효과가 $\tau^*$ 근처에 집중될 것이라는 예측.

**통계적 그림 (effective degrees of freedom).** 선형 보간은 두 extreme matrix 사이의 1차원 직선. HST는 $\tau^*$ 부근에서 **직선에서 가장 멀리 떨어진** matrix를 만들어 내므로, 분석 대상이 1차원 저차원 구조가 아닌 더 높은 **effective degrees of freedom**을 가진다고 해석 가능.

## 3.3 그래프 class별 의존성 (추측)

**추측 3.2.** drift 조건 (DC) 하의 그래프에 대해, $\tau^*$는 다음 세 지표의 조합에 대략 의존:

1. **그래프 spectral gap** $\gamma(\mathcal{G})$: 작을수록 mixing 느림 → $\tau^*$ 큼
2. **$\mathbf{y}$의 variation** $\mathrm{Var}(\mathbf{y})$: 클수록 유클리드 정보 강함 → Phase I 길어짐 → $\tau^*$ 큼
3. **파라미터 $b$**: 작을수록 slow dynamics → $\tau^*$ 큼

**후보 scaling:**

$$\tau^* \;\sim\; \frac{\mathrm{Var}(\mathbf{y})^{1/2}}{b \cdot \gamma(\mathcal{G})} \quad \text{(추측)}$$

exponent 정확성은 실험과 엄밀 증명으로 검증해야 할 대상.

**graph class별 정성적 예측:**

| 그래프 class | $\gamma$ | 예상 $\tau^*$ | 예상 $\max \Delta$ |
|:---|:---:|:---:|:---:|
| 완전 그래프 $K_n$ | 큼 | 작음 | 작음 (구조 거의 없음) |
| Ring / Cycle | 중간 | 중간 | 중간 |
| $k$-NN / random geometric | 중간 | 중간 | **큼** |
| Bottleneck graph (bi-cluster) | 작음 | **큼** | **큼** (cliff 효과) |
| Expander | 큼 | 작음 | 중간 |

Bottleneck 그래프에서 $\tau^*$가 크고 편차도 큼 — HST의 강점이 가장 돋보이는 class라는 예측.

---

# 4. Matrix path의 위상과 면적

## 4.1 두 궤적, 하나의 공간

$n \times n$ symmetric·nonneg 행렬 공간 $\mathcal{S}_+^n$ 위에서:

- $t \mapsto W^{\mathrm{HST}}(t)$: 확률적 궤적 (각 sample path별)
- $\alpha \mapsto W^{\mathrm{lin}}(\alpha)$: 결정론적 직선

둘은 같은 두 점 $W^E, W^G$를 잇는다. 두 path의 합집합이 $\mathcal{S}_+^n$ 내 **닫힌 곡선(closed loop)**을 형성.

## 4.2 Enclosed area — HST 비선형성의 대역 지표

**정의 (Path-integrated deviation):**

$$\mathcal{A}_{\mathrm{HST}} \;:=\; \int_0^\infty \Delta(t)\, |\alpha'(t)|\, dt$$

또는 등가적으로

$$\mathcal{A}_{\mathrm{HST}} \;=\; \int_0^1 \big\|W^{\mathrm{HST}}(t(\alpha)) - W^{\mathrm{lin}}(\alpha)\big\|_F\, d\alpha$$

이것은 두 궤적이 둘러싸는 영역의 "Frobenius-면적"이다 (반경형 측정).

**의미.**

- $\mathcal{A}_{\mathrm{HST}} = 0$ ⟺ HST path = 선형결합 path (HST가 자명한 mixing만 함)
- $\mathcal{A}_{\mathrm{HST}}$가 클수록 HST가 선형 null 대비 많이 휘어진 궤적을 따름 → HST의 **대역적 비선형성**

**$\tau^*$와 $\mathcal{A}$의 상보성.**

- $\tau^*$: **최대 순간 편차**의 시점
- $\mathcal{A}$: **누적 편차**의 총량

같은 $\mathcal{A}$를 주는 그래프 둘이라도 $\tau^*$의 위치가 다를 수 있다. 전자는 "얼마나 휘었나", 후자는 "언제 가장 휘었나".

**그래프 class 지표로.**

$$\mathcal{A}_{\mathrm{HST}}(\mathcal{G}, \mathbf{y}) \;\text{: HST가 } (\mathcal{G}, \mathbf{y}) \text{에서 얼마나 의미 있는가의 대역 지표}$$

이 지표를 실험적으로 계산하여 그래프 class별 HST의 "가치"를 순서화할 수 있다.

## 4.3 Homotopy 관점

두 궤적은 $\mathcal{S}_+^n$에서 같은 endpoint를 잇지만 **다른 homotopy class에 속할 수 있다**. 예컨대 HST path가 특정 "singular locus" (예: 특정 고유값이 0이 되는 표면)를 **우회**한다면 선형결합 path와 homotopically 다른 경로.

실용적으로는 $\mathcal{S}_+^n$이 convex이므로 global homotopy는 trivial이지만, 의미 있는 **codimension-1 subvariety**를 singular로 정의하면 두 path의 winding 관점이 생긴다. 예:

- Laplacian의 두 번째 고유값 $\lambda_2 = 0$이 되는 표면 (연결성이 깨지는 표면)
- 특정 고유벡터의 부호가 바뀌는 표면 (hierarchical clustering의 birth/death 표면)

이런 singular locus를 path가 지나가느냐 여부가 topological invariant 후보.

(상세 전개는 §6 O4 참조.)

---

# 5. 실험 계획

각 실험은 공통으로 (a) $\Delta(t)$의 시각화, (b) $\tau^*$의 측정, (c) $\mathcal{A}_{\mathrm{HST}}$의 계산을 수행.

## 5.1 Ring graph (논문 Example 2)

**세팅.** $n = 60$, ring 위에 Gaussian kernel weight. $\mathbf{y} = \mathbf{0}$ (순수 그래프 case) vs $\mathbf{y} = \pm 3$ (modified) 두 가지.

**측정.**
- $\Delta(t)$ vs $t$ 플롯
- $\tau^*$ 위치와 $\Delta(\tau^*)$ 값
- $\mathcal{A}_{\mathrm{HST}}$ 값
- **비교 지표**: $\mathbf{y} = \mathbf{0}$과 $\mathbf{y} = \pm 3$ 사이에서 $\tau^*, \mathcal{A}$가 어떻게 달라지는가? Phase I 유무가 보임.

**예상.** $\mathbf{y} = \mathbf{0}$에서는 Phase I이 없어 $\tau^*$가 작고 $\mathcal{A}$도 작음. $\mathbf{y} = \pm 3$에서는 두 클러스터가 섞이는 transition이 $\tau^*$ 근처.

## 5.2 Bi-cluster with cliff (논문 Fig 5)

**세팅.** Ring 위에 $y$가 불연속적으로 바뀌는 지점 (cliff) 포함. 논문의 cliff detection 결과가 나타나는 세팅 재현.

**측정.**
- $\Delta(t)$ peak이 cliff 감지 발생 시점과 일치하는가?
- $\tau^*$와 논문의 "적절한 $t$" 추천 범위 비교

**예상.** $\tau^*$ 근처가 cliff 감지의 "sweet spot" — HST 실전 활용의 이론적 근거.

## 5.3 MCU 데이터 (논문 §6)

**세팅.** 논문의 MCU 영화 network. 23 nodes, asymmetric $\mathbf{W}$, $\mathbf{y}$는 box office.

**측정.**
- $\Delta(t)$ vs $t$
- 논문의 권장 $t = 10^6$, $b = 1$, $T_{\max} = 500$이 $\tau^*$ 근처인가?
- MCU의 hierarchical clustering (multi-hero / core movies / others) 분리능이 $\tau^*$에서 최대인지 확인

**예상.** 논문이 경험적으로 선택한 $t$ 값이 $\tau^*$와 대략 일치하면 HST의 실용적 guideline을 이론화한 것.

## 5.4 Random graph family

**세팅.** 다음 random graph family에서 파라미터을 varying시키며 실험:
- Erdős–Rényi $G(n, p)$ with $p$ 변화
- Stochastic block model with community size/strength 변화
- Random geometric graph with radius 변화

**측정.**
- 각 class에서 $\tau^*$와 $\mathcal{A}$의 empirical distribution
- 그래프의 spectral gap $\gamma$와 $\tau^*$의 scaling 관계 — 추측 3.2 검증

**예상.** spectral gap이 작을수록 $\tau^*$가 크게 나오는 scaling. Stochastic block model에서 bottleneck이 강할수록 $\mathcal{A}$가 큼.

## 5.5 Bottleneck graph (극단 case)

**세팅.** 비대칭 바벨, 두 클러스터 + 좁은 bridge. Drift 조건이 거의 깨지는 경계 영역.

**측정.**
- (DC) 만족 경계에서의 $\tau^*$ behavior
- (DC) 깨지면 $\tau^*$가 어떻게 되는가 — $\infty$로 발산? 유한이지만 비정형?

**예상.** 이 세팅이 HST 이론의 경계를 탐색하는 실험.

## 5.6 파라미터 sensitivity

모든 실험에서 공통으로:

- $b$의 영향: $b \to 0$에서 $\tau^*$가 어떻게 scaling?
- $T_{\max}$의 영향: 큰 $T_{\max}$에서 $\tau^*$ 안정성
- $\theta$의 영향: snow weight의 bandwidth — $\tau^*$ location이 $\theta$에 의존?

---

# 6. 열린 문제와 결론

## 6.1 열린 문제

**O1. $\tau^*$의 closed-form expression.** 추측 3.2의 $\tau^* \sim \mathrm{Var}(\mathbf{y})^{1/2}/(b\gamma)$의 엄밀화. spectral gap, $\mathbf{y}$ variation, $b$의 exponent 정확히 결정.

**O2. 비자명성의 엄밀 판정.** $W^{\mathrm{HST}}$ path가 직선일 조건 (Proposition 2.2의 (a))을 정확히 특성화. 어떤 $(\mathcal{G}, \mathbf{y})$에서 HST가 자명한 mixer가 되는가? (직관적으로 거의 항상 비자명이지만 엄밀 증명 필요.)

**O3. $\mathcal{A}_{\mathrm{HST}}$의 해석적 공식.** 누적 면적이 $(\mathcal{G}, \mathbf{y}, b, T_{\max})$에서 어떻게 결정되는가. Spectral quantity로 표현 가능한가?

**O4. Singular locus와 homotopy class.** §4.3의 topological invariant 엄밀 정의와 계산. HST path가 피하는 singular locus의 특성화.

**O5. 다른 norm.** Frobenius 대신 operator norm, $\ell_\infty$, Bregman divergence 등에서 $\tau^*$가 어떻게 달라지는가. 특히 spectral-analysis 관점에서 operator norm이 더 자연스러울 수 있음.

**O6. Computational tractability.** $\Delta(t)$와 $\tau^*$를 계산하려면 $W^{\mathrm{HST}}(t)$ trajectory 전체가 필요. Batch-parallel ([보충아이디어 B](260208_21543c.html))로 효율화 가능한가?

**O7. $\tau^*$ 근처에서 HST 효과의 집중.** "cliff detection이 $\tau^*$ 근처에서 극대화"라는 §3.2의 가설의 실험적·이론적 검증.

## 6.2 결론

본 포스트는 HST의 "자기 자신의 선형근사 대비 비자명한 mixing"을 정량화하기 위해 **선형결합 path와의 편차** $\Delta(t)$와 최대 편차 시각 $\tau^*$를 도입했다.

핵심 기여:

1. **형식화 (Internal)**: $W^G$를 HST의 Harris 극한으로 잡아, 두 endpoint가 HST 내재적으로 결정되도록 깔끔하게 구성
2. **존재 정리 (Theorem 3.1)**: drift 조건 하에서 내부 최댓값 $\tau^* \in (0, \infty)$가 존재함을 간결하게 증명
3. **해석**: $\tau^*$는 HST의 feedback loop가 최대로 작용하는 시점. Cliff detection, hierarchical clustering 등 HST 실전 효과의 **이론적 스위트 스팟** 후보
4. **대역 지표**: $\mathcal{A}_{\mathrm{HST}}$ (path-integrated deviation)가 그래프별 HST 가치의 global measure
5. **실험 계획**: ring, bi-cluster, MCU, random graph family에서의 구체적 시뮬레이션 로드맵

[보충아이디어 A.7](260208_21543c.html)에서 "매우 중요"로 표시된 이 주제가 본 포스트에서 본격 전개되었다. 다음 단계:

- 실험 (§5)을 실제 코드로 구현하여 $\tau^*$·$\mathcal{A}$의 empirical behavior 확인
- 추측 3.2의 scaling을 random graph family에서 검증
- cliff/bottleneck 사례에서 $\tau^*$가 실용 guideline과 일치하는지 검증

이 결과들이 축적되면, HST의 "$t$ 선택 이론"이라는 실용 문제가 이론적 지지대를 얻는다.

**External 관점과의 보완.** 본 포스트는 HST 자신의 비선형성(self-consistency)에 집중. HST를 **classical graph metric (diffusion, effective resistance 등)**과 비교하여 "HST가 기존 접근과 얼마나 다른가"를 묻는 External 관점은 별도 포스트 [HST와 외재 Graph Distance의 편차](260424_ext.html)에서 다룬다. 두 관점은 상보적이다 — Internal은 HST 내부 구조, External은 HST의 타자 대비 독창성.

---

# 부록: 보충아이디어 문서와의 연결

본 포스트는 다음과 크로스-참조된다:

- [보충아이디어 A.7](260208_21543c.html): 본 포스트의 주제 제기 (짧은 개요)
- [보충아이디어 A.5](260208_21543c.html): Phase I/II 점근 프레임워크 — $t_0^*$가 재매개화 후보에 등장
- [보충아이디어 A.6](260208_21543c.html): Gromov-Hausdorff 수렴 — 극한 거리공간 $\mathcal{M}_\infty$
- [보충아이디어 D](260208_21543c.html): Harris 에르고드 수렴 — $W^G$ 존재의 근거
- [260407 Main Theorem](260407_eecd78.html): $W^G$의 a.s. 수렴 엄밀 증명
