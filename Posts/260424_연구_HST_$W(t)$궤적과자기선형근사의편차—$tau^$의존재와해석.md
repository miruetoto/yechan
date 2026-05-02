---
title: 연구 ▷ HST ▷ $W(t)$ 궤적과 자기 선형근사의 편차 — $\tau^*$의 존재와 해석
author: claude
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

## 5.7 ex2 첫 empirical 결과 — 신호 dominance regime의 발견

§5.1 의 ring (cylinder $y=\pm 3$) 세팅에서 실제 시뮬레이션 ([260424_예제2_편차분석.py](../260424_예제2_편차분석.py), [260425_예제2_거리비교.py](../260425_예제2_거리비교.py), [260425_예제2_선형결합실험.py](../260425_예제2_선형결합실험.py)) 을 돌려본 결과, **본 포스트의 이론 frame 자체는 잘 작동하지만, ex2 가 frame 의 강점을 보이기엔 부적합**하다는 점이 드러났다.

### 5.7.1 측정값 요약 ($\tau = 5 \times 10^5$, $b = 0.05$)

세 metric 으로 잰 $\tau^*$ 와 $\mathcal{A}_{\mathrm{HST}}$:

| metric | $\tau^*_{\mathrm{int}}$ | $\max \Delta^{\mathrm{int}}$ | $\mathcal{A}^{\mathrm{int}}_{\mathrm{HST}}$ |
|:---|---:|---:|---:|
| Frobenius (raw) | 598 | 75.18 | 46.15 |
| Procrustes (shape) | 17 286 | 0.005 | 0.003 |
| PD bottleneck ($H_1$) | 1 640 | 0.007 | 0.005 |

해석:
- **$\tau^*$ 자체는 존재한다 — Theorem 3.1 의 예측이 empirical 로 확인.**
- **그러나 max $\Delta$ 가 작다.** 특히 Procrustes 와 PD 는 거의 0 에 가까움. 즉 HST path 가 linear path 위에서 *별로 벗어나지 않음*. §3.1 에서 가정한 "비퇴화 가정" 이 ex2 에서 **양적으로는 거의 퇴화에 가까운** 상태.
- **세 metric 의 $\tau^*$ 가 한 자릿수~두 자릿수 차이로 다름.** "어느 metric 으로 정의된 $\tau^*$ 가 의미 있는가" 라는 질문이 새로 등장.

### 5.7.2 결정적 발견 — Linear-combo path 가 HST 와 거의 같음

별도 실험: $D_{\mathrm{lin}}(\alpha) = (1-\alpha)\hat D^E + \alpha \hat D^G$ 에 대해 $\alpha \in [0, 1]$ sweep, cMDS-3D 임베딩 후 HST 임베딩과 **Procrustes disparity** 를 계산.

| $D^G$ 후보 | 최적 $\alpha^*$ | min Procrustes disparity to HST |
|:---|---:|---:|
| ring angular dist² | **0.02** | **0.0029** |
| shortest path² | **0.02** | **0.0029** |
| effective resistance | **0.04** | **0.0031** |

세 후보 모두 $\alpha^* \approx 0$ — **거의 pure Euclidean ($\hat D^E$ 만)** 에서 HST shape 에 사실상 도달. ex2 에서는 HST cylinder 의 cMDS 임베딩이 **신호-only Euclidean 의 cMDS 임베딩과 좌표 수준에서 구분되지 않음**.

### 5.7.3 왜 그런가 — 신호 dominance via eigenvalue gap

cMDS 의 분산 분해:
$$\text{Var}(\text{embedding}) \;=\; \sum_k \lambda_k(D_{\text{centered}})$$

Cylinder 신호 ($y = \pm 3$) 의 squared scale = $36$, ring radius (=1) 의 squared scale = $\mathcal{O}(1)$. → $D^{\mathrm{HST}}$ 의 첫 eigenvalue $\lambda_1$ (signal split 축) 이 $\lambda_2, \lambda_3$ (ring 축) 을 **압도** ($\sim 36 : 1$).

Procrustes 는 standardize (Frob 정규화) 후 비교하므로 $\sum \lambda_k$ 으로 나눔 → dominant axis 만 정렬되면 disparity 는 대략

$$\text{disparity} \;\approx\; \frac{\lambda_2 + \lambda_3 + \cdots}{\sum_k \lambda_k} \;\ll\; 1$$

**즉 $\Delta(t)$ 가 작은 것은 HST 의 비선형성이 작은 것이 아니라, signal-dominated regime 의 자연스러운 귀결**. Linear path 와의 편차가 magnitude 적으로 미미.

### 5.7.4 그러나 Procrustes 가 잡지 못하는 차원 — HST 의 *남는* 가치

위 disparity 0.003 은 좌표-수준의 닮음만 측정. Procrustes 가 무시하는 차원에서는 HST 와 pure Eucl 이 여전히 다름:

- **위상 ($H_1$)**: HST 의 cylinder 안 ring 은 closed loop → $H_1 \neq 0$. Pure Eucl 은 두 cluster point → $H_1 = 0$. **위상적으로 다름**.
- **Cluster 내부 정렬**: HST 는 같은 cluster 내 노드를 ring angle 순으로 매끄럽게 배치, Eucl 은 noise 수준 (~$0.04$) 로 무작위.
- **Multi-scale 제어**: HST 는 $\tau$ 로 임베딩 차원성·복잡도 조절. Eucl 은 단일 scale.

즉 ex2 에서 **HST 의 진정한 contribution 은 "small but qualitatively essential" 한 ring perturbation**. Δ-frame (Frobenius) 으로는 놓치고, **PD H1** 또는 **within-cluster ordering correlation** 같은 metric 으로만 보임.

### 5.7.5 Regime A/B/C 분류

본 결과를 일반화하면 신호 vs 그래프 스케일 비율로 3 개 regime:

| regime | 조건 | HST 의 행동 | $\tau^*$ |
|:---|:---|:---|:---|
| **A** (신호 dominant) | $\|\mathbf y\|^2 \gg$ (graph scale) | $W^{\mathrm{HST}}(t) \approx W^E$ 부근에서 거의 안 움직임. Linear interpolation 잘 맞음. | 매우 작음, $\Delta$ 작음 |
| **B** (그래프 dominant) | $\|\mathbf y\|^2 \ll$ (graph scale) | Pure $\mathbf y = 0$ 이면 §8.1 에서 본 commute-time 류로 환원. Linear path 도 그래프 metric 에 수렴. | 의미 있는 $\tau^*$, 그러나 HST 만의 가치 약함 |
| **C** (균형 또는 이질결합) | $\|\mathbf y\|^2 \sim$ graph scale, 또는 신호와 구조가 *반대 방향* | 신호와 그래프가 진정으로 nonlinearly 결합. Linear path 로 도달 불가능한 shape. | 의미 있고, $\Delta$ 도 큼 — **Theorem 3.1 의 비퇴화 가정이 진짜 의미를 가짐** |

ex2 ($y = \pm 3$, ring r=1) 는 **regime A** 에 속한다. §3.3 의 추측 $\tau^* \sim \mathrm{Var}(\mathbf y)^{1/2} / (b\gamma)$ 는 신호가 *상대적으로* 클수록 $\tau^*$ 가 *커진다* 고 예측하는데, regime A 에서는 *오히려 max $\Delta$ 가 작아져* HST path 가 trivial 해짐. 추측은 부분적으로만 옳고, **regime classification 이 더 근본적인 변수**임을 시사.

### 5.7.6 ex2 의 demo 적격성 재고

본 포스트 (Δ-frame) 와 [HST 와 외재 거리의 편차](260424_ext.html) 의 두 frame 모두에서 **ex2 ($y=\pm 3$) 는 HST 의 강점이 잘 드러나지 않는 example** 임이 확인됨. 향후 frame 의 power 를 보이려면:

- regime C 인 example 로 이동: **예제5 (cliff), 예제6 (parity), 예제7 (twin rings), 예제8 (directed cycle)** ([260424 추가예제](260424_extra.html))
- 또는 ex2 자체에서 신호 amplitude 를 $\pm 0.5 \sim \pm 1$ 로 낮춰 regime A → C 천이 관찰

→ §6.1 에 새 열린 문제로 추가.

### 5.7.7 별도 정리

본 5.7 의 발견은 [예제2 재검토 — HST 는 Euclidean 과 얼마나 다른가?](260425_ex2review.html) 에 자세한 narrative 로 정리. 본 포스트의 frame 검증 관점에서는 위 요약으로 충분.

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

**O8. Regime classification.** §5.7.5 에서 도입한 A/B/C 분류 (신호 dominance vs 그래프 dominance vs 균형) 의 정량적 경계. 어떤 invariant — 예: dominance index $\rho := \lambda_1(D^{\mathrm{HST}}) / \sum_k \lambda_k$ — 가 regime 을 잡는가? $\rho$ 의 임계값은 무엇인가? Regime C 에서 추측 3.2 의 scaling 이 비로소 검증되는가.

**O9. Norm 별 $\tau^*$ 의 일치/불일치.** §5.7.1 에서 Frob / Procrustes / PD 의 $\tau^*$ 가 한두 자릿수 차이 — 어느 norm 의 $\tau^*$ 가 "이론적으로 자연스러운가". 후보 비교: Frob (코드량 small), Procrustes (shape, regime A 에서 거의 0), PD (topology, signal-graph 결합 잡음). $\tau^*$ 들이 *언제* 일치하는가의 sufficient condition.

**O10. ex2 외 예제들에서 frame 의 거동.** §5.7.6 의 후속 — 예제5–8 (cliff, parity, twin rings, directed) 에서 동일 sweep. Regime C 가 실현되는 example 의 분포. 본 frame 의 *적용 영역* 을 정의하는 문제.

**O11. 신호 amplitude scan.** ex2 에서 $\mathbf y$ 의 amplitude 를 $0.1 \to 3$ 으로 변화시키며 $\tau^*, \mathcal{A}, \alpha^*$ (linear-combo 최적), Procrustes-min disparity 의 변화 추적. **Regime A → C 천이가 어디서 일어나는가** 의 구체적 수치 baseline.

## 6.2 결론

본 포스트는 HST의 "자기 자신의 선형근사 대비 비자명한 mixing"을 정량화하기 위해 **선형결합 path와의 편차** $\Delta(t)$와 최대 편차 시각 $\tau^*$를 도입했다.

핵심 기여:

1. **형식화 (Internal)**: $W^G$를 HST의 Harris 극한으로 잡아, 두 endpoint가 HST 내재적으로 결정되도록 깔끔하게 구성
2. **존재 정리 (Theorem 3.1)**: drift 조건 하에서 내부 최댓값 $\tau^* \in (0, \infty)$가 존재함을 간결하게 증명
3. **해석**: $\tau^*$는 HST의 feedback loop가 최대로 작용하는 시점. Cliff detection, hierarchical clustering 등 HST 실전 효과의 **이론적 스위트 스팟** 후보
4. **대역 지표**: $\mathcal{A}_{\mathrm{HST}}$ (path-integrated deviation)가 그래프별 HST 가치의 global measure
5. **실험 계획**: ring, bi-cluster, MCU, random graph family에서의 구체적 시뮬레이션 로드맵

[보충아이디어 A.7](260208_21543c.html)에서 "매우 중요"로 표시된 이 주제가 본 포스트에서 본격 전개되었다. 다음 단계:

- 실험 (§5)을 실제 코드로 구현하여 $\tau^*$·$\mathcal{A}$의 empirical behavior 확인 → §5.7 (ex2 첫 결과) 에서 부분 수행
- 추측 3.2의 scaling을 random graph family에서 검증
- cliff/bottleneck 사례에서 $\tau^*$가 실용 guideline과 일치하는지 검증

이 결과들이 축적되면, HST의 "$t$ 선택 이론"이라는 실용 문제가 이론적 지지대를 얻는다.

**6.3 ex2 첫 empirical 결과의 함의 (§5.7 요약 — 본 frame 자체에 대한 메시지).**

§5.7 의 ex2 시뮬레이션에서 다음이 드러났다:

1. **$\tau^*$ 자체는 잘 정의되고 측정된다** — Theorem 3.1 의 존재성 주장은 empirical 로 확인.
2. **그러나 ex2 ($y=\pm 3$) 에서 $\Delta$ 가 양적으로 미미** — Procrustes 0.003, PD 0.007. **신호 dominance regime A** (§5.7.5) 의 자연스러운 귀결로, HST path 가 linear path 위에서 거의 안 벗어남.
3. **Linear-combo path ($\alpha^* \approx 0$, 거의 pure Eucl) 가 HST shape 를 사실상 재현** — *좌표 수준* 에선. 위상·cluster 내 정렬 등 Procrustes 가 잡지 않는 차원에서만 HST 가 unique.

**즉 본 frame ($\Delta$, $\tau^*$, $\mathcal{A}_{\mathrm{HST}}$) 의 *predictive power* 는 regime C (신호와 그래프가 이질적으로 결합) 에서만 충분히 발휘**된다. ex2 (regime A) 는 frame 의 demo 로 부적절. 향후 frame 의 의미를 보이려면 cliff/parity/twin/directed 같은 regime C example 로 옮겨가야 함 (O10).

**또 한 가지 — Frame 의 *metric* 선택의 중요성.** §5.7.1 의 Frob / Procrustes / PD 의 $\tau^*$ 가 일치하지 않음 (O9). 본 포스트의 형식화는 Frobenius 를 묵시적 default 로 채택했지만, regime A 에서 Frobenius $\Delta$ 는 dominant signal axis 의 magnitude 를 반영하고, Procrustes/PD 는 더 미세한 shape/topology 를 본다. **각 norm 이 잡는 "비선형성" 의 종류가 다르다** — 이는 §6.1 O5 (다른 norm) 의 실증적 시작점.

**External 관점과의 보완.** 본 포스트는 HST 자신의 비선형성(self-consistency)에 집중. HST를 **classical graph metric (diffusion, effective resistance 등)**과 비교하여 "HST가 기존 접근과 얼마나 다른가"를 묻는 External 관점은 별도 포스트 [HST와 외재 Graph Distance의 편차](260424_ext.html)에서 다룬다. 두 관점은 상보적이다 — Internal은 HST 내부 구조, External은 HST의 타자 대비 독창성. 흥미롭게도 ex2 의 regime A 발견은 두 frame 모두에 영향: *둘 다 ex2 에선 약하게 나옴*.

**연결 포스트:**
- [예제2 재검토 — HST 는 Euclidean 과 얼마나 다른가?](260425_ex2review.html) — regime A 발견의 단독 narrative
- [그래프 도메인에서의 거리](260425_distancesurvey.html) — HST 가 위치 잡혀야 할 family 지형
- [예제2의 diffusion distance 정의 오류](260425_diffusionbug.html) — 본 검토의 출발점

---

# 부록: 보충아이디어 문서와의 연결

본 포스트는 다음과 크로스-참조된다:

- [보충아이디어 A.7](260208_21543c.html): 본 포스트의 주제 제기 (짧은 개요)
- [보충아이디어 A.5](260208_21543c.html): Phase I/II 점근 프레임워크 — $t_0^*$가 재매개화 후보에 등장
- [보충아이디어 A.6](260208_21543c.html): Gromov-Hausdorff 수렴 — 극한 거리공간 $\mathcal{M}_\infty$
- [보충아이디어 D](260208_21543c.html): Harris 에르고드 수렴 — $W^G$ 존재의 근거
- [260407 Main Theorem](260407_eecd78.html): $W^G$의 a.s. 수렴 엄밀 증명
