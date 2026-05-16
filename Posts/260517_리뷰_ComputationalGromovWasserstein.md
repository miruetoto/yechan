---
title: 리뷰 ▷ Computational Gromov-Wasserstein distance
author: 연
date: 05/17/2026
draft: false
output-file: 260517_c4f8a2.html
fontsize: 0.85em
---

<style>
.math.display { text-align: left; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
</style>

<https://doi.org/10.1016/j.procs.2022.01.086>

> Zheng, L., Xiao, Y., Niu, L. (2022). *A brief survey on Computational Gromov-Wasserstein distance*. Procedia Computer Science **199**, 697–702.

### 초록 (옮김)

그래프는 점(node) 자체뿐 아니라 점 사이의 관계까지 담는 일반화된 측도 공간(metric measure space)이므로, 그래프 두 개 사이의 거리(distance)는 일반적인 분포 거리로는 잘 잡히지 않는다. 최근 **Gromov-Wasserstein (GW) discrepancy** 가 그래프 위에서의 의사-거리(pseudometric)로 제안되어 이론적 근거를 갖추었지만, 실제 응용은 아직 적다. 본 논문은 GW 의 기본 아이디어·응용·발전을 정리하고, 이차계획(quadratic programming) 형태인 GW 를 **반복 알고리즘**으로 푸는 방법을 다룬다.

---

## 1. 왜 또 거리인가

::: {.callout-important collapse="false" title="한 줄 직관"}
**두 분포 간 거리**(Wasserstein)와 **두 그래프 간 거리**(Gromov-Wasserstein)는 다르다. Wasserstein 은 노드값을 정렬해 옮기지만, GW 는 노드 사이의 **관계 구조**(거리 행렬) 까지 같이 옮긴다.
:::

머신러닝에서 두 분포의 차이를 잴 때 보통 KL divergence, JS divergence, Total Variation, $L^2$ 같은 양을 쓴다. 하지만 데이터가 **그래프** 라면 다음 두 가지가 동시에 필요하다.

- 노드(점)들이 어떤 분포를 따르는가 (질량, 빈도)
- 노드 사이 관계가 어떻게 연결돼 있는가 (인접성, 거리)

전통적 Wasserstein 은 첫째만 다룬다. 둘째까지 같이 다루는 것이 **Gromov-Wasserstein**.

본 회사([HST 주식회사](#))의 본업인 *Heavy-Snow Transform (HST)* 가 같은 동기(노드 값 + 노드 구조 동시 반영)에서 출발한다는 점에서, GW 는 비교 대상이자 보완 대상이다. 단, HST 는 **한 그래프 안에서 노드 쌍 거리**(intra-graph)를, GW 는 **그래프 두 개 사이 전체 거리**(inter-graph)를 잰다는 점이 결정적인 차이다.

---

## 2. Optimal Transport 와 Wasserstein 거리

Optimal Transport (OT) 는 Monge (1781) 가 흙 운반 문제로 제기했고, Kantorovich 가 선형계획(LP) 으로 완화했다.

::: {.callout-note collapse="true" title="증명 식: 이산 Wasserstein"}

운반 문제는 다음과 같다. $n$ 개 창고에 원료 $a_i$ 단위, $m$ 개 공장에 수요 $b_j$ 단위. $i$ 창고에서 $j$ 공장으로 한 단위 옮기는 비용이 $C_{ij}$. 운반 계획 $T \in \mathbb{R}^{n \times m}_+$ 의 총비용을 최소화:

$$
\begin{aligned}
L_C(p, q) &= \min_{T \in \Pi(p,q)} \langle C, T \rangle = \min_T \sum_{i,j} C_{ij} T_{ij} && (\because \text{운반 총비용}) \\
\Pi(p, q) &= \{ T \in \mathbb{R}^{n \times m}_+ : T \mathbf{1}_m = p,\ T^\top \mathbf{1}_n = q \} && (\because \text{주변분포 보존})
\end{aligned}
$$

비용이 $C = D^k$ (거리의 $k$ 제곱) 이면

$$
W_k(p, q) = L_{D^k}(p, q)^{1/k}
$$

는 대칭·양·삼각부등식 모두 만족하는 **k-Wasserstein 거리** 가 된다.

:::

---

## 3. Gromov-Wasserstein discrepancy

핵심 차이점은 **두 공간이 서로 다른 차원·메트릭** 일 때도 거리를 잴 수 있게 만든 것이다.

두 metric measure space $(X, d_X, \mu_X)$, $(Y, d_Y, \mu_Y)$ 에 대해

$$
GW(\mu_X, \mu_Y) = \inf_{\pi \in \Pi(\mu_X, \mu_Y)} \iint_{X \times Y} \iint_{X \times Y} L(x, y, x', y')\, d\pi(x,y)\, d\pi(x',y')
$$

여기서 손실 함수 $L(x, y, x', y') = f\big(d_X(x, x'),\ d_Y(y, y')\big)$. 즉 **두 공간에서 짝지어진 점쌍의 내부거리가 얼마나 다른가** 를 잰다.

::: {.callout-important collapse="false" title="해석"}
Wasserstein 이 "각 점을 어디로 옮길까?" 라면, GW 는 "각 **점쌍** 이 다른 공간의 어떤 점쌍에 대응되는가?" 를 묻는다. **점 자체의 좌표는 의미가 없고, 점 사이의 관계만 본다.** 그래서 차원이 달라도 비교 가능.
:::

### 이산 형태

그래프를 dissimilarity matrix $C \in \mathbb{R}^{N \times N}$ 와 분포 $p \in \Sigma_N$ 의 쌍 $(C, p)$ 로 본다.

$$
GW(C, \tilde C, p, q) = \min_{T \in \Pi(p, q)} \varepsilon_{C, \tilde C}(T),
\quad
\varepsilon_{C, \tilde C}(T) = \sum_{i, j, k, l} L(C_{ik}, \tilde C_{jl})\, T_{ij}\, T_{kl}
$$

4-way 텐서 $\mathcal{L}(C, \tilde C) = \big(L(C_{ij}, \tilde C_{kl})\big)$ 와 텐서-행렬 곱

$$
(\mathcal{L}(C, \tilde C) \otimes T)_{i,j} = \sum_{k,l} L(C_{ik}, \tilde C_{jl})\, T_{kl}
$$

로 쓰면 간결히

$$
\varepsilon_{C, \tilde C}(T) = \langle \mathcal{L}(C, \tilde C) \otimes T,\ T \rangle.
$$

$T$ 에 대한 **이차계획(quadratic programming, QP)**. 일반적으로 NP-hard.

---

## 4. 응용 — graph matching, partitioning

| 작업 | 의미 |
|---|---|
| graph matching | 두 그래프 사이의 노드 대응 (예: 단어 임베딩 정렬, cross-lingual 매칭) |
| graph partitioning | 노드를 군집으로 묶기 (barycenter graph 학습) |
| multi-graph matching/partitioning | 여러 그래프 동시 정렬 |

Alvarez-Melis et al. (2018) 은 GW 로 단어 임베딩을 OT 문제로 풀어, 하이퍼파라미터 거의 없이도 cross-lingual word alignment 에서 SOTA 급 성능을 보였다. Xu Hongteng et al. (2019) 은 GW 학습 (Gromov-Wasserstein Learning, GWL) 으로 대규모 그래프 정렬·분할의 통합 패러다임을 만들었다.

---

## 5. 알고리즘 — Sinkhorn 과 Proximal

QP 라서 직접 푸는 게 어렵다. 두 가지 우회로:

### 5.1 엔트로피 정규화 (Mémoli, Peyré 등)

$$
GW_\epsilon(C, \tilde C, p, q) = \min_{T \in \Pi(p, q)} \varepsilon_{C, \tilde C}(T) - \epsilon H(T),
\quad H(T) = -\sum_{i,j} T_{ij} \log T_{ij}
$$

투영 경사법(projected gradient descent) 으로 다음 갱신을 반복:

$$
T^{(k+1)} \leftarrow \arg\min_{T \in \Pi(p,q)} \langle \mathcal{L}(C, \tilde C) \otimes T^{(k)},\ T \rangle.
$$

각 스텝은 비용 행렬 $\mathcal{L}(C, \tilde C) \otimes T^{(k)}$ 의 **OT 부문제(sub-OT)**. **Sinkhorn projection** 으로 빠르게 풀린다. 단, 전역 수렴 보장은 없다.

### 5.2 Proximal Point (Xu et al.)

KL divergence 를 proximal term 으로 추가:

$$
T^{(k+1)} = \arg\min_{T \in \Pi(p, q)} \langle L(C, \tilde C, T^{(k)}),\ T \rangle + \gamma\, \mathrm{KL}(T \,\Vert\, T^{(k)})
$$

복잡한 비볼록 문제를 **일련의 볼록 부문제** 로 분해. Sinkhorn-Knopp 로 거의 선형 시간에 수렴. 안정점에 **전역 수렴 보장**.

::: {.callout-tip collapse="true" title="보충: 두 방법의 차이"}
- 엔트로피 정규화: 단순하고 빠르지만 정규화 강도 $\epsilon$ 에 민감, 수렴 보장 없음.
- Proximal: 매 스텝 KL 부담을 묶어 비볼록을 볼록 부문제로 풀어냄. 안정성 좋고 큰 그래프에 적합.
:::

---

## 6. Fused Gromov-Wasserstein (FGW)

GW 는 **구조(edge relation)** 만 본다. 그러나 그래프가 **노드 속성(feature)** 도 가질 때는 둘 다 써야 한다. 예시: "boy" vs "girl" 과 "football" vs "basketball" 은 관계 구조가 비슷하지만 의미는 완전히 다르다.

Vayer Titouan et al. (2019) 의 **Fused GW**:

$$
FGW(p, q) = \min_{T \in \Pi(p, q)} \big\langle (1 - \alpha) M + \alpha\, \mathcal{L}(C, \tilde C) \otimes T,\ T \big\rangle
$$

여기서 $M_{ij}$ 는 한 그래프의 $i$ 번 점에서 다른 그래프의 $j$ 번 점으로 옮기는 **노드 특징 비용**. $\alpha \in [0, 1]$ 가 구조/특징 가중치를 조절.

::: {.callout-important collapse="false" title="HST 와의 연결"}
**HST 가 한 그래프 안에서 Euclidean (노드값) + Graph (구조) 를 동시에 반영하는 거리** 라면, **FGW 는 두 그래프 사이에서 노드 특징 + 관계 구조를 동시에 반영하는 거리** 이다. 둘 다 "두 도메인 정보를 동시에 다룬다" 는 motivation 을 공유하지만, **scope** (intra-graph vs inter-graph)와 **연산 방식** (random walk 누적 vs OT 최소화)이 다르다.
:::

조건부 경사법(conditional gradient) 으로 비볼록 문제이지만 국소 stationary point 로 수렴.

---

## 7. 회사 본업과의 비교 정리

| 항목 | HST (회사 본업) | Gromov-Wasserstein (이 논문) |
|---|---|---|
| **scope** | 한 그래프 안 두 노드 거리 (intra-graph) | 두 그래프 사이 전체 거리 (inter-graph) |
| **두 도메인 통합** | Euclidean ($y$) + Graph 구조 ($\mathcal{E}$) | 노드 특징 ($M$) + 거리 행렬 ($C$) (Fused GW) |
| **연산** | 눈 누적 random walk, $SD^2(t)$ | OT 최적화, $\min_T \langle \mathcal{L} \otimes T, T \rangle$ |
| **정렬/매칭** | 노드 임베딩, 클러스터링 | 노드 대응, 그래프 정렬·분할 |
| **수렴 보장** | weak ergodicity (Theorem A.5) | proximal: 전역 수렴, entropy: 보장 X |
| **계산량** | $O(t)$ 시뮬 (대규모는 $\mathbf{P}^t$ 직접 못 구함) | 한 스텝 Sinkhorn, NP-hard 본문제 |

두 방법은 **경쟁 관계가 아니라 보완 관계**. HST 가 한 그래프 안의 풍부한 다중스케일 거리를 만들면, 그 거리 행렬을 **GW 의 $C$** 로 넣어 두 그래프를 정렬·비교하는 응용이 자연스럽다.

---

## 8. 정리

- GW 는 **그래프-그래프 거리**의 사실상 표준 도구로 자리잡고 있다.
- Sinkhorn 기반 엔트로피 정규화는 빠르지만 수렴 보장이 약하고, proximal 방법은 더 안정적·확장적.
- Fused GW 는 노드 속성까지 통합, HST 의 motivation 과 자연스럽게 만난다.
- 한국·국내 그래프 연구에서 HST 와 GW 를 **결합한 응용 사례** 는 아직 보이지 않는다. 회사 차원에서 잠재적 방향 하나.

다음 단계 (개인 노트):
- HST 의 $SD(t)$ 행렬을 GW 의 $C$ 로 넣은 **HST-GW 합성** 시도 검토
- 본 논문이 인용한 원전 (Mémoli, Peyré, Vayer Titouan, Xu Hongteng) 정독

---

*작성: 연 (연구팀), 2026-05-17 01:05*
