---
title: 연구 ▷ HST ▷ HST와 외재 Graph Distance의 편차
author: 클로드
date: 04/24/2026
draft: false
output-file: 260424_6ae2dc.html
---

# 0. 동기: HST는 classical graph metric과 어떻게 다른가

[Internal 버전 포스트](260424_a5d784.html)에서는 HST가 **자기 자신의 선형근사** 대비 얼마나 휘어있는지를 분석했다. 두 endpoint가 모두 HST 내재적이므로 $\Delta(t)$가 양 끝에서 0으로 수렴하여 내부 최댓값 $\tau^{*,\mathrm{int}}$가 깔끔하게 존재했다.

본 포스트는 **External 관점**을 다룬다: HST를 $\mathbf{y}$를 무시하는 기존 graph metric — diffusion distance, effective resistance, shortest path 등 — 과 **직접 경쟁**시킨다. 이 관점의 가치:

> **"HST는 classical graph theory가 보지 못하는 것을 본다"**는 정성적 주장을 **정량적 지표**로 환원.

Internal 관점이 "HST 내부 구조의 비자명성"을 묻는다면, External 관점은 "HST의 타자 대비 독창성"을 묻는다. 둘은 독립적 질문이고, 외재적 null과의 비교는 HST 연구의 external validity에 직접 기여한다.

**개념적 대가.** External 버전은 Internal보다 이론이 까다롭다:

- 두 endpoint가 **다르다** — HST는 Harris 극한 $W^{G,\mathrm{int}}$로, classical은 $W^{G,\mathrm{ext}}$로 끝남
- 따라서 $\Delta^{\mathrm{ext}}(t)$는 $t \to \infty$에서 0으로 가지 **않고** 양의 상수 $\Delta^{\mathrm{ext}}_\infty$로 수렴
- 내부 최댓값 $\tau^{*,\mathrm{ext}}$가 유한할 수도, 무한할 수도 있음 — 그래프 class별로 정성적 다름

이 복잡성이 오히려 External 버전을 "HST 우위의 경계 탐색"으로서 흥미롭게 만든다.

**구성.** §1 외재적 $W^G$의 세 후보 (diffusion, effective resistance, shortest path). §2 $\Delta^{\mathrm{ext}}(t)$의 경계 거동 분석. §3 $\tau^{*,\mathrm{ext}}$의 존재 조건 — transient peak vs 점근 편차. §4 두 graph metric 사이의 hierarchy — HST가 어떤 classical metric과 가장 멀리 떨어져 있나. §5 실험 계획. §6 열린 문제.

---

# 1. 외재적 $W^G$의 세 후보

## 1.1 기본 조건

외재적 $W^G_{\mathrm{ext}}$는 다음을 만족해야 한다:

- (C1) **$\mathbf{y}$-무관**: $y$에 의존하지 않고 $\mathcal{G}$의 구조만으로 결정
- (C2) **대칭적**: $W^{G,\mathrm{ext}}_{ij} = W^{G,\mathrm{ext}}_{ji}$
- (C3) **양정치성**: Gaussian kernel로 환산 가능, 즉 거리 $d_{ij}$로부터 $W^{G,\mathrm{ext}}_{ij} = \exp(-d_{ij}^2/2\theta^2)$

HST의 $W^{G,\mathrm{int}} = \lim_t W^{\mathrm{HST}}(t)$는 $\mathbf{y} \bmod b$에 국소 종속이므로 (C1)을 엄격히 위반. 따라서 Internal 버전의 그래프 극한은 진정한 외재적 metric이 아니다. External 버전은 (C1)을 만족하는 엄밀한 null을 필요로 한다.

## 1.2 후보 A: Diffusion Distance

**정의.** 그래프 Laplacian $\mathbf{L} = \mathbf{D} - \mathbf{W}$의 고유분해 $\mathbf{L} = \sum_k \lambda_k \boldsymbol{\psi}_k \boldsymbol{\psi}_k^\top$. Diffusion distance는

$$DD^2_{ij}(s) \;:=\; \sum_{k \geq 1} e^{-2\lambda_k s} (\psi_k(i) - \psi_k(j))^2$$

여기서 $s > 0$은 diffusion scale parameter.

**$W$ 행렬로 변환.**

$$W^{G,\mathrm{diff}}_{ij}(s) \;:=\; \exp\!\left(-\frac{DD^2_{ij}(s)}{2\theta^2}\right)$$

**특징.**

- **Random walk 기반** — 주요 random walk based 방법의 대표
- **Scale parameter $s$** 선택이 필요 (HST의 $t$와 유사한 knob)
- **Closed form** — 고유분해로 직접 계산
- **$\mathbf{y}$ 완전 무관** — 순수 $\mathcal{G}$-결정

**HST와의 관계.** 논문 §3의 배경 설명에서 diffusion distance가 직접 비교 대상으로 언급. HST는 "이걸 넘어서 $\mathbf{y}$도 반영한다"는 주장을 내세움. External $\Delta$는 이 주장의 정량화.

## 1.3 후보 B: Effective Resistance

**정의.** 그래프를 전기회로로 보고, edge $(i,j)$의 저항을 $1/w_{ij}$로 설정. 두 노드 사이의 effective resistance:

$$R_{\mathrm{eff}}(i, j) \;:=\; (\boldsymbol{e}_i - \boldsymbol{e}_j)^\top \mathbf{L}^\dagger (\boldsymbol{e}_i - \boldsymbol{e}_j)$$

여기서 $\mathbf{L}^\dagger$는 Moore-Penrose pseudoinverse.

**$W$ 행렬로 변환.**

$$W^{G,\mathrm{res}}_{ij} \;:=\; \exp\!\left(-\frac{R_{\mathrm{eff}}(i, j)}{2\theta^2}\right)$$

**특징.**

- **Parameter-free** — scale parameter 없음
- **Commute time 관련** — $R_{\mathrm{eff}}(i,j) = C_{ij}/(2|E|)$ where $C_{ij}$ is expected commute time
- **Scale invariance** 부족 — edge weight의 absolute scale에 민감

**HST와의 관계.** HST의 $c_{ij} = \mathbb{E}_\pi[(\delta_i - \delta_j)^2]$가 **commute time과 관련이 있다**는 추측이 [보충아이디어 H의 열린 문제 1](260208_21543c.html)에서 제기. Effective resistance를 null로 쓰면 이 관계를 정량화할 수 있다.

## 1.4 후보 C: Shortest Path

**정의.**

$$d_{\mathrm{sp}}(i, j) \;:=\; \min_{\text{path } i \to j} \sum_{\text{edges in path}} \frac{1}{w_{ij}}$$

**$W$ 행렬로 변환.**

$$W^{G,\mathrm{sp}}_{ij} \;:=\; \exp\!\left(-\frac{d_{\mathrm{sp}}(i, j)^2}{2\theta^2}\right)$$

**특징.**

- **Parameter-free**
- **결정론적** (random walk 기반 아님)
- **Metric** — 삼각부등식 엄격

**HST와의 관계.** 가장 "결정론적인 graph metric" — random walk 기반의 HST와 가장 구조가 다름. External $\Delta$에서 가장 큰 편차를 기대.

## 1.5 세 후보의 비교

| | Diffusion | Effective Resistance | Shortest Path |
|:---|:---:|:---:|:---:|
| Parameter | $s$ (scale) | 없음 | 없음 |
| 기반 | Random walk | Random walk (commute) | 결정론적 |
| Closed form | ✓ (고유분해) | ✓ (pseudoinverse) | ✓ (Dijkstra) |
| Metric? | Pseudo | Pseudo | True metric |
| HST와 관련 | 논문에서 직접 비교 | $c_{ij}$ 추측과 연결 | 가장 멀 것으로 예상 |

**본 포스트의 기본 선택.** 가장 자연스럽게 HST와 비교되는 **diffusion distance**를 주로 사용. 다른 두 후보는 §4의 hierarchy 비교 및 §5의 실험에서 활용.

---

# 2. $\Delta^{\mathrm{ext}}(t)$의 경계 거동

## 2.1 정의

$W^{G,\mathrm{ext}}$를 위의 외재적 후보 중 하나로 고정. Linear combination path:

$$W^{\mathrm{lin,ext}}(\alpha) := (1-\alpha)\, W^E + \alpha\, W^{G,\mathrm{ext}}, \qquad \alpha \in [0, 1]$$

재매개화 $\alpha(t): [0, \infty] \to [0, 1]$은 Internal 버전과 동일 (단순 시간 기반, 지수형, endpoint-calibrated 중 선택). 본 포스트의 기본은 **단순 시간 기반** $\alpha(t) = t/(t + t_0^*)$ — endpoint가 다르므로 endpoint-calibrated이 의미가 다르다 (Internal과 달리 HST 종점이 아닌 classical 종점을 목표로 calibrate).

편차:

$$\Delta^{\mathrm{ext}}(t) \;:=\; \big\| W^{\mathrm{HST}}(t) - W^{\mathrm{lin,ext}}(\alpha(t)) \big\|_F$$

$$\tau^{*,\mathrm{ext}} \;:=\; \arg\max_{t \geq 0} \Delta^{\mathrm{ext}}(t)$$

## 2.2 경계값: $t = 0$

$W^{\mathrm{HST}}(0) = W^E$이고 $\alpha(0) = 0$이므로 $W^{\mathrm{lin,ext}}(\alpha(0)) = W^E$. 따라서

$$\Delta^{\mathrm{ext}}(0) = 0$$

Internal과 동일 (출발점 $W^E$는 HST와 classical 모두 공유).

## 2.3 경계값: $t \to \infty$ — **핵심 차이**

HST는 Harris 극한 $W^{G,\mathrm{int}}$로 수렴:

$$W^{\mathrm{HST}}(t) \xrightarrow{t \to \infty} W^{G,\mathrm{int}}$$

반면 linear path는 classical $W^{G,\mathrm{ext}}$로 수렴:

$$W^{\mathrm{lin,ext}}(\alpha(t)) \xrightarrow{t \to \infty} W^{G,\mathrm{ext}}$$

일반적으로 $W^{G,\mathrm{int}} \neq W^{G,\mathrm{ext}}$이므로:

$$\Delta^{\mathrm{ext}}(t) \xrightarrow{t \to \infty} \Delta^{\mathrm{ext}}_\infty \;:=\; \|W^{G,\mathrm{int}} - W^{G,\mathrm{ext}}\|_F \;\geq\; 0$$

**$\Delta^{\mathrm{ext}}_\infty$의 의미.** HST의 Harris 극한이 classical graph metric과 얼마나 다른가의 정량 지표. 이 자체가 HST의 **점근적 독창성**을 측정한다.

- $\Delta^{\mathrm{ext}}_\infty = 0$: HST가 점근적으로 classical metric과 수렴 (고유성 없음)
- $\Delta^{\mathrm{ext}}_\infty > 0$: HST가 어떤 classical metric과도 본질적으로 다름

**추측.** 일반적인 $(\mathcal{G}, \mathbf{y})$에서 $\Delta^{\mathrm{ext}}_\infty > 0$. 즉 Harris 극한은 classical graph metric과 결코 같지 않다. 증명은 $W^{G,\mathrm{int}}$의 $\mathbf{y} \bmod b$ 종속성으로부터 — classical metric은 $\mathbf{y}$에 전혀 의존하지 않으므로, $\mathbf{y}$가 격자 경계를 움직이면 $W^{G,\mathrm{int}}$가 변하지만 $W^{G,\mathrm{ext}}$는 그대로. 이 불일치가 $\Delta^{\mathrm{ext}}_\infty$의 소멸을 막는다.

## 2.4 경계 거동의 표

| | $t = 0$ | $t \to \infty$ |
|:---|:---:|:---:|
| Internal $\Delta^{\mathrm{int}}(t)$ | 0 | 0 |
| External $\Delta^{\mathrm{ext}}(t)$ | 0 | $\Delta^{\mathrm{ext}}_\infty \geq 0$ |

Internal은 닫힌 곡선을 만들고, External은 열린 곡선. 이 차이가 $\tau^*$ 분석의 근본 차이를 만든다.

---

# 3. $\tau^{*,\mathrm{ext}}$의 존재 조건

## 3.1 존재의 분기

$\Delta^{\mathrm{ext}}(t)$는 $[0, \infty]$에서 연속이고 $\Delta^{\mathrm{ext}}(0) = 0$, $\Delta^{\mathrm{ext}}(\infty) = \Delta^{\mathrm{ext}}_\infty$. 따라서 $\max_t \Delta^{\mathrm{ext}}(t)$는 $[0, \infty]$ compact한 의미에서 존재하지만, **내부 최대**인지 점근 극한의 경계인지 구분 필요:

**Case A: 단조증가.** $\Delta^{\mathrm{ext}}$가 $t$에 대해 단조비감소이면

$$\tau^{*,\mathrm{ext}} = \infty, \qquad \max \Delta^{\mathrm{ext}} = \Delta^{\mathrm{ext}}_\infty$$

이 경우 HST와 classical metric의 차이는 "시간이 갈수록 계속 커지며 점근에 수렴". **end-to-end 비선형성**.

**Case B: 비단조 — Transient Peak.** 어떤 유한 $t^{\mathrm{peak}} \in (0, \infty)$에서 $\Delta^{\mathrm{ext}}(t^{\mathrm{peak}}) > \Delta^{\mathrm{ext}}_\infty$이면

$$\tau^{*,\mathrm{ext}} = t^{\mathrm{peak}} \in (0, \infty), \qquad \max \Delta^{\mathrm{ext}} > \Delta^{\mathrm{ext}}_\infty$$

이 경우 HST의 transient 거동이 점근보다 더 크게 classical과 편차를 보임. **중간 거동의 고유성**.

**두 case의 해석 차이.** Case A는 "HST가 점진적으로 classical metric과 떨어져 최종적으로 $\Delta^{\mathrm{ext}}_\infty$만큼 다름". Case B는 "transient 중에 뭔가 독특한 일이 벌어지며, 이후 점근은 상대적으로 덜 극적". Case B가 HST 이론상 더 흥미 — HST 고유의 **중간 단계 효과**가 존재.

## 3.2 Theorem 3.1 (경계 거동과 $\tau^{*,\mathrm{ext}}$ 분기)

**Theorem 3.1.** drift 조건 (DC)와 $W^{G,\mathrm{int}} \neq W^{G,\mathrm{ext}}$ 하에서:

(i) $\Delta^{\mathrm{ext}}(t)$는 $[0, \infty]$에서 연속이고 $\Delta^{\mathrm{ext}}(0) = 0$, $\Delta^{\mathrm{ext}}(\infty) = \Delta^{\mathrm{ext}}_\infty > 0$.

(ii) Case A ($\tau^{*,\mathrm{ext}} = \infty$) vs Case B (유한 $\tau^{*,\mathrm{ext}}$)의 분기는 다음과 동치:
- Case A ⟺ $\sup_{t < \infty} \Delta^{\mathrm{ext}}(t) \leq \Delta^{\mathrm{ext}}_\infty$
- Case B ⟺ $\sup_{t < \infty} \Delta^{\mathrm{ext}}(t) > \Delta^{\mathrm{ext}}_\infty$

(iii) Case B이면 유한 $t^{\mathrm{peak}} \in (0, \infty)$이 존재하여 $\Delta^{\mathrm{ext}}(t^{\mathrm{peak}}) > \Delta^{\mathrm{ext}}_\infty$.

**증명 스케치.** (i)은 §2.3. (ii)는 sup의 정의에서 자명. (iii)은 연속 함수가 $(0, \infty)$에서 점근값을 상회하면 최댓값이 내부에 있음. $\square$

## 3.3 어떤 그래프가 Case B인가? (추측)

**추측 3.2.** HST가 classical metric과 transient에서 더 다르게 움직이는 그래프는 다음 특성을 갖는다:

1. **Bottleneck 구조**: 두 클러스터 + 좁은 bridge. HST는 $\mathbf{y}$ 정보로 bridge를 "발견"하지만 classical은 graph 구조만 봐서 bridge 특수성을 놓침. Transient 동안 HST가 이 차이를 극대화.

2. **Cliff / Discontinuity**: $\mathbf{y}$가 급격히 변하는 지점. HST는 cliff 근처에서 snow 축적 패턴이 변하며 classical과 크게 차이. 논문 Figure 5의 상황.

3. **Hierarchical structure**: 다층 클러스터. HST의 $t$ 변화에 따라 서로 다른 계층이 순차적으로 드러나면서 classical과 transient 차이 발생.

**Case A가 예상되는 그래프:**

1. **Homogeneous graph**: Complete graph, regular graph, expander — $\mathbf{y}$에 무관한 구조에서 HST 고유성이 최종 극한에서 드러나고 transient에서는 monotone.

2. **Smooth $\mathbf{y}$**: $\mathbf{y}$가 그래프 위에서 매끄럽게 변하면 cliff/bottleneck 효과 없음.

## 3.4 Case B의 직관: 왜 transient peak이 나타나는가

HST의 동역학을 보면, 초기에는 $\mathbf{y}$의 구배가 snow의 흐름을 강하게 지시한다 (Phase I). 중간기에는 이 구배가 snow 축적으로 평탄화되어 "원래 $\mathbf{y}$"의 효과는 잔류만 남는다 (Phase I→II 전환). 후기에는 그래프 구조만 남는다 (Phase II).

Classical graph metric은 **$\mathbf{y}$의 효과가 전혀 없으므로 Phase II의 근사**. 따라서:

- 초기 ($t = 0$): HST $\approx$ 유클리드, classical $\approx$ 그래프, 둘 다 같은 "자명한" 경로 출발 → $\Delta = 0$
- 중간 ($t \sim \tau^*$): HST가 **$\mathbf{y}$의 non-trivial mixing 중** — classical은 여전히 $\mathbf{y}$ 무시하여 linear mixing만. 편차 최대.
- 후기 ($t \to \infty$): HST는 $W^{G,\mathrm{int}}$에 수렴, classical은 $W^{G,\mathrm{ext}}$. 편차는 점근 상수.

Transient peak이 존재하는 조건은 "중간 단계에서 HST가 $\mathbf{y}$를 적극 활용하는 정도"가 "점근에서의 차이 $\Delta^{\mathrm{ext}}_\infty$"보다 큰 경우 — bottleneck·cliff·hierarchical 그래프에서 이 조건이 충족된다는 예측.

---

# 4. Classical Null들 사이의 Hierarchy

## 4.1 HST와 어떤 classical metric의 편차가 가장 큰가

세 후보 중 어떤 것을 null로 쓰느냐에 따라 $\Delta^{\mathrm{ext}}_\infty$와 $\tau^{*,\mathrm{ext}}$가 달라진다. 이를 비교하면 HST가 "어떤 classical 접근과 가장 다른가"의 hierarchy가 나온다.

**정의.**

$$\Delta^{\mathrm{diff}}_\infty := \|W^{G,\mathrm{int}} - W^{G,\mathrm{diff}}\|_F$$

$$\Delta^{\mathrm{res}}_\infty := \|W^{G,\mathrm{int}} - W^{G,\mathrm{res}}\|_F$$

$$\Delta^{\mathrm{sp}}_\infty := \|W^{G,\mathrm{int}} - W^{G,\mathrm{sp}}\|_F$$

## 4.2 예상 Hierarchy

**추측 4.1 (Null의 유사도 순서).** 일반적 연결 그래프 $\mathcal{G}$에 대해:

$$\Delta^{\mathrm{diff}}_\infty \;\leq\; \Delta^{\mathrm{res}}_\infty \;\leq\; \Delta^{\mathrm{sp}}_\infty$$

**근거.**

- **Diffusion이 HST와 가장 가까움**: 둘 다 random walk 기반. Diffusion은 homogeneous walk의 극한, HST는 $\mathbf{y}$-biased non-homogeneous walk. 구조적으로 가장 유사.
- **Effective Resistance는 중간**: 역시 random walk 기반 (commute time 관련)이지만 Laplacian의 pseudoinverse에 의존하여 spectral 성질이 diffusion과 다름.
- **Shortest Path가 가장 멂**: 결정론적, random walk 없음. HST의 확률적 성격과 가장 멀리 떨어짐.

## 4.3 Transient Peak 여부의 hierarchy

위와 유사한 추측:

**추측 4.2 (Transient Peak의 likeliness).** HST가 transient peak ($\tau^{*,\mathrm{ext}} < \infty$)을 보일 likeliness:

$$\mathrm{Diffusion} \;>\; \mathrm{Effective\;Resistance} \;>\; \mathrm{Shortest\;Path}$$

즉 HST와 가까운 classical metric일수록 transient에서 편차가 두드러질 수 있다. 멀리 떨어진 shortest path는 처음부터 끝까지 꾸준히 멀어 Case A 경향.

이 추측은 실험적으로 검증 가능하며, 긍정된다면 "HST와 classical random walk 접근" 사이의 fine-grained difference를 잡을 수 있는 도구가 된다.

## 4.4 추가 null 후보: Commute Time Distance 직접

$R_{\mathrm{eff}}$ 대신 commute time $C_{ij} = \mathbb{E}[\tau_{i \to j} + \tau_{j \to i}]$을 그대로 쓸 수 있다:

$$W^{G,\mathrm{ct}}_{ij} := \exp(-C_{ij}/2\theta^2)$$

$C_{ij}$가 HST의 $c_{ij}$와 **직접적으로 관련**된다는 추측 (보충아이디어 H 열린 문제 1) 때문에, 이 null과의 비교는 특히 해석적으로 의미 있다. $\Delta^{\mathrm{ct}}_\infty$가 작다는 결과는 "HST의 Harris 극한 $\approx$ commute time distance"를 시사하고, 크다는 결과는 HST가 commute time 이상을 한다는 증거.

---

# 5. 실험 계획

모든 실험에서 (a) $\Delta^{\mathrm{ext}}(t)$ 시각화, (b) $\tau^{*,\mathrm{ext}}$ 측정, (c) $\Delta^{\mathrm{ext}}_\infty$ 계산, (d) Case A/B 판정.

## 5.1 세 null에 대한 비교

각 graph에 대해 세 null (diffusion, resistance, shortest path)에서 모두 계산하여 §4의 hierarchy 추측 검증.

## 5.2 Ring Graph (논문 Example 2)

**세팅.** $n = 60$, ring, $\mathbf{y} = \mathbf{0}$ vs $\mathbf{y} = \pm 3$.

**측정.**
- 세 null에 대한 $\Delta^{\mathrm{ext}}(t)$
- Case A인지 Case B인지
- $\mathbf{y} = \mathbf{0}$ (유클리드 없음)에서 $\Delta^{\mathrm{ext}}_\infty$이 작을 것으로 예상 (HST $\approx$ classical)
- $\mathbf{y} = \pm 3$에서 $\Delta^{\mathrm{ext}}_\infty$이 커질 것 — $\mathbf{y}$ 효과가 HST의 점근값을 classical에서 멀어지게

## 5.3 Bi-cluster with Cliff (논문 Fig 5)

**세팅.** Cliff가 있는 bi-cluster. 가장 **Case B가 예상되는** 세팅.

**측정.**
- $\Delta^{\mathrm{ext}}(t)$의 transient peak 위치
- Diffusion null에서 peak이 cliff 감지 시점과 일치하는지
- **HST 우위의 empirical 증거**: peak에서 HST가 cliff를 잡고, 그 시점에 diffusion은 못 잡는 것을 보이는 그림

## 5.4 MCU 데이터 (논문 §6)

**세팅.** MCU 영화 network, 23 nodes.

**측정.**
- Diffusion distance와 HST의 MCU 분석 결과 비교
- 두 접근의 주요 component (multi-hero/core/others) 분리능이 다른지
- $\Delta^{\mathrm{ext}}(t)$의 peak이 논문 권장 $t = 10^6$ 근처인지

## 5.5 Expander vs Bottleneck

**세팅.**
- Expander: Ramanujan graph, random $d$-regular
- Bottleneck: Bi-cluster with narrow bridge

**측정.**
- Expander에서 Case A 예상 검증
- Bottleneck에서 Case B 예상 검증
- 추측 3.2 (어떤 그래프가 Case B) 실험적 확인

## 5.6 Commute Time과의 Hierarchy 검증

추측 4.1, 4.2를 여러 random graph family에서 검증:

$$\Delta^{\mathrm{diff}}_\infty, \Delta^{\mathrm{res}}_\infty, \Delta^{\mathrm{sp}}_\infty, \Delta^{\mathrm{ct}}_\infty$$

의 순서가 예측과 일치하는지 empirical distribution으로 확인.

---

# 6. 열린 문제와 결론

## 6.1 열린 문제

**O1. $\Delta^{\mathrm{ext}}_\infty > 0$의 엄밀 증명.** §2.3의 추측 — Harris 극한 $W^{G,\mathrm{int}}$이 어떤 classical $W^{G,\mathrm{ext}}$와도 본질적으로 다르다는 주장. $\mathbf{y} \bmod b$ 의존성 논변을 formal하게.

**O2. $\tau^{*,\mathrm{ext}}$의 Case A/B 판정 조건.** 어떤 $(\mathcal{G}, \mathbf{y}, \mathrm{null})$에서 transient peak이 발생하는지의 엄밀 판정.

**O3. 추측 3.2의 검증.** Bottleneck, cliff, hierarchical 그래프에서 Case B가 예상되는 논리의 엄밀 formulation.

**O4. 추측 4.1, 4.2의 hierarchy.** 세 null 간 편차 크기의 순서가 일반적으로 성립하는가.

**O5. HST와 commute time의 관계 — $c_{ij} \overset{?}{\approx} C_{ij}$.** [보충아이디어 H 열린 문제 1](260208_21543c.html)의 핵심 추측. $\Delta^{\mathrm{ct}}_\infty$가 작다면 이 관계의 empirical 증거이고, 이론적으로는 Foster-Lyapunov 하의 정상분포 $\pi$가 commute time의 structure와 어떻게 관련되는지 탐구 필요.

**O6. Multiple null의 결합.** 세 null과의 편차 벡터 $(\Delta^{\mathrm{diff}}, \Delta^{\mathrm{res}}, \Delta^{\mathrm{sp}})(t)$를 한꺼번에 분석하는 multi-dimensional $\Delta$ 구성. HST가 "어떤 classical 접근도 개별적으로는 잡지 못하는" 정보를 제공하는지의 joint 증거.

**O7. 외재 null로서의 graph kernel.** Weisfeiler-Lehman kernel, graphlet kernel 등 더 현대적인 graph kernel을 null로 쓰는 extension.

## 6.2 결론

본 포스트는 HST의 "**classical graph theory 대비 독창성**"을 정량화하기 위해 외재적 null과의 편차 $\Delta^{\mathrm{ext}}(t)$를 도입했다.

핵심 기여:

1. **Null의 체계화 (§1)**: Diffusion distance, effective resistance, shortest path 세 classical graph metric을 공통 프레임으로 제시, HST와의 비교 대상으로서의 특성 분석.

2. **경계 거동의 이해 (§2)**: Internal과 달리 External $\Delta^{\mathrm{ext}}(\infty) \geq 0$이라는 핵심 차이 identification. $\Delta^{\mathrm{ext}}_\infty$ 자체가 HST의 점근적 독창성 지표.

3. **Case A/B 분기 (§3)**: $\tau^{*,\mathrm{ext}}$가 유한(transient peak)인지 무한(end-to-end growth)인지의 분기를 명시. 각 case의 해석 대조.

4. **Hierarchy 추측 (§4)**: 세 null 간 $\Delta$ 크기의 순서 추측. HST와 classical random walk 방법의 fine-grained difference 탐색 도구.

5. **실험 계획 (§5)**: Cliff, expander, bottleneck, MCU에서 각 null에 대한 체계적 실험.

**Internal 관점과의 보완.** 본 External 관점은 [Internal 포스트](260424_a5d784.html)의 "HST의 내재적 비선형성"과 상보적. 둘을 합치면 HST의 특성에 대한 완전한 그림:

- **Internal $\tau^{*,\mathrm{int}}$**: HST가 자기 자신의 직선 근사에서 얼마나 휘어있나 → 내부 구조의 풍부함
- **External $\tau^{*,\mathrm{ext}}$**: HST가 classical 접근에서 얼마나 멀어져 있나 → 타자 대비 독창성
- **두 관점의 일치점**: 아마 동일한 현상(bottleneck, cliff)에서 둘 다 큰 값. 이 일치가 HST 기여의 **robustness**를 보장.

다음 단계:

- §5의 실험을 코드로 구현 (특히 diffusion distance와 HST의 concrete 비교)
- 추측 3.2, 4.1, 4.2를 random graph family에서 검증
- O5 ($c_{ij} \approx C_{ij}$)의 empirical 확인이 되면 HST와 classical commute time theory 사이의 직접 연결이 열림

---

# 부록: 관련 포스트

- [Internal 버전 포스트](260424_a5d784.html): HST가 자기 자신의 선형근사에서 얼마나 휘어있나
- [보충아이디어 A.7](260208_21543c.html): 본 포스트의 주제 제기
- [보충아이디어 G](260208_21543c.html): $W(t)$에서 입력 복원 — 본 포스트의 External 관점이 정보이론적 관점과 연결
- [보충아이디어 H 열린 문제 1](260208_21543c.html): $c_{ij}$와 commute time의 관계 — O5와 직접 연결
- [260407 Main Theorem](260407_eecd78.html): $W^{G,\mathrm{int}}$ 존재의 근거
