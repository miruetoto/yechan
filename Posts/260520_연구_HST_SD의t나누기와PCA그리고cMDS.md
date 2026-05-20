---
title: 연구 ▷ HST ▷ SD 의 t-나누기 × PCA × cMDS — 같은 임베딩의 네 얼굴
author: 연
date: 05/20/2026
draft: false
output-file: 260520_baf177.html
fontsize: 0.85em
---

<style>
.math.display { text-align: left; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
</style>

<!-- 소유권자: 연 | 사용자: 연 -->

> 대표가 던진 두 질문 (2026-05-20):
>
> **Q1.** SD 를 $t$ 로 나누고 PCA 하는 경우 vs 나누지 않고 바로 PCA 하는 경우 차이가 있는가?
>
> **Q2.** SD 를 cMDS 하는 것과 PCA 하는 것의 차이?
>
> 본 포스트의 답: **임베딩의 좌표 방향은 네 가지 모두 동일.** 다만 (i) eigenvalue 의 스케일이 다르고, (ii) **어떤 regime 정보 (Thm A/B/C) 가 dominant axis 로 떠오르느냐** 가 t-나누기 선택에 따라 갈린다. cMDS = PCA 는 SD 가 Euclidean 거리 family 이므로 자동 동치.

### TL;DR

- **Q2 답**: SD 는 height 함수의 $L^2$ 거리이므로 정확히 Euclidean. 이 경우 cMDS 와 PCA 는 *동일한 임베딩 좌표* 를 준다 (Young-Householder 정리). 차이가 있다면 PCA 가 추가로 *variable directions* (height 좌표 공간의 방향) 를 함께 주는 정도.

- **Q1 답**: $SD$ / $SD^2/t$ / $SD^2/t^2$ / $SD^2/t^3$ 의 PCA — eigenvector (좌표 *방향*) 는 모두 같음. eigenvalue 가 $t$ 의 거듭제곱만큼 스케일. 하지만 **multi-regime situation** (어떤 pair 는 Thm A, 어떤 pair 는 B, 어떤 pair 는 C 에 속할 때) 에서는 nontrivial: scale factor 가 큰 (Thm C drift) pair 가 raw SD 에서는 첫 두 PC 를 독점하지만, $t^3$ 로 나누면 그 drift component 가 normalize 되어 다른 regime 의 정보가 PC 에 드러난다.

---

## 1. 기본 정리 — PCA = cMDS (Euclidean 거리에서)

### 1.1 두 알고리즘

**PCA.** 점 $h_1, \dots, h_n \in \mathbb R^d$ 가 있다고 하자 (각 $h_i$ 는 height 함수의 시간 sequence). 중심화 $\tilde h_i := h_i - \bar h$, 데이터 행렬 $\tilde X = [\tilde h_1 \cdots \tilde h_n]^\top \in \mathbb R^{n \times d}$. PCA 는 $d \times d$ covariance

$$C \;=\; \frac{1}{n}\tilde X^\top \tilde X$$

의 eigen-분해 $C = V \Lambda V^\top$ 를 한 뒤, 각 점을 top-$k$ 방향에 사영하여 좌표 $Z = \tilde X V_{:, 1:k} \in \mathbb R^{n \times k}$ 를 얻는다.

**Classical MDS.** 거리 행렬 $D \in \mathbb R^{n \times n}$ 만 주어졌을 때, 그 제곱을 double-centering:

$$B \;=\; -\tfrac{1}{2}\,H\,D^{(2)}\,H, \qquad H := I - \tfrac{1}{n}\mathbf 1\mathbf 1^\top.$$

$B$ 의 eigen-분해 $B = U \Sigma U^\top$ 에서 top-$k$ 만 살려 좌표 $Y = U_{:, 1:k}\, \Sigma_{1:k}^{1/2} \in \mathbb R^{n \times k}$.

### 1.2 동치

거리 $D_{ij} = \|h_i - h_j\|_2$ 가 (어떤 Euclidean embedding 의) 거리라면 — 즉 $B \succeq 0$ 이면 —

$$B \;=\; \tilde X \tilde X^\top \quad (\text{Gram matrix})$$

이고, $\tilde X = U \Sigma^{1/2} V^\top$ 의 SVD 에서

$$\tilde X^\top \tilde X = V \Sigma V^\top \;\Rightarrow\; nC = V \Sigma V^\top,$$
$$\tilde X \tilde X^\top = U \Sigma U^\top \;=\; B.$$

따라서 **PCA 의 principal scores $\tilde X V_{:,1:k} = U_{:,1:k}\, \Sigma_{1:k}^{1/2} = $ cMDS 의 좌표 $Y$.** 즉

$$\boxed{\quad Y_{\mathrm{cMDS}}(D) \;=\; Z_{\mathrm{PCA}}(\tilde X) \quad \text{(Euclidean 거리에서)}\quad}$$

— 동일한 좌표를 얻는다. 이게 Young-Householder 1938 의 고전 결과.

### 1.3 SD 가 Euclidean 임의 확인

SD 의 정의 ([본인 HST₀ 노트](./260520_dcd721.html) §1 등 참조):

$$SD^2_{ij}(t) \;=\; \sum_{s=0}^{t}\bigl(h(v_i, s) - h(v_j, s)\bigr)^2 \;=\; \|h_i^{(t)} - h_j^{(t)}\|_2^2,$$

여기서 $h_i^{(t)} := (h(v_i, 0), h(v_i, 1), \dots, h(v_i, t)) \in \mathbb R^{t+1}$. **명백히 $L^2$ 거리** — Euclidean. 따라서 SD 의 cMDS 와 height 좌표 $h_i^{(t)}$ 의 PCA 는 동치.

이게 Q2 의 답: **차이 없음** (좌표 한정). PCA 가 추가로 주는 것은 height-coordinate 공간 $\mathbb R^{t+1}$ 의 *시간축 방향* — 즉 어느 시점 $s$ 가 임베딩에 가장 기여하는지를 보여주는 $V_{1:k}$. cMDS 에는 없는 정보지만 실용적으로 잘 안 봄.

---

## 2. t-나누기는 좌표 방향을 바꾸지 않는다 — uniform rescale 의 효과

### 2.1 좌표 한 줄 변환

$SD$ 를 $t^\alpha$ 로 나눈다는 것은 height 좌표 $h_i^{(t)}$ 에 동일한 스칼라 $t^{-\alpha/2}$ 를 곱하는 것과 같다 (모든 점에 같은 양으로):

$$\frac{SD^2_{ij}(t)}{t^{2\alpha}} \;=\; \bigl\|t^{-\alpha} h_i^{(t)} - t^{-\alpha} h_j^{(t)}\bigr\|_2^2.$$

(여기 $\alpha = 1/2$ 면 $SD^2/t$, $\alpha = 1$ 이면 $SD^2/t^2$, $\alpha = 3/2$ 면 $SD^2/t^3$ 에 해당.)

### 2.2 PCA 에서의 효과

데이터 행렬을 $\tilde X \to t^{-\alpha} \tilde X$ 로 균일하게 스케일. 그러면

$$C \to t^{-2\alpha}\, C, \qquad B \to t^{-2\alpha}\, B.$$

eigen-분해 $C = V \Lambda V^\top$ → $t^{-2\alpha} C = V (t^{-2\alpha}\Lambda) V^\top$. **eigenvector $V$ 는 그대로, eigenvalue 만 $t^{-2\alpha}$ 곱.**

cMDS 좌표 $Y = U \Sigma^{1/2}$ → $t^{-\alpha} U \Sigma^{1/2} = t^{-\alpha} Y$. 즉 좌표 전체에 균일한 축소.

**결론:** *상대적인 임베딩 모양* — nodes 의 배치 — 은 완전히 동일. 단지 좌표축의 단위가 $t^\alpha$ 만큼 바뀔 뿐.

### 2.3 그래서 차이가 없다?

여기까지면 답은 "**$t$ 로 나누든 말든 PCA 임베딩은 같다**" 라는 무미건조한 결론. 그러나 한 가지 더 봐야 한다.

---

## 3. nontrivial 한 차이 — multi-regime situation

### 3.1 HST 의 세 regime 이 한 그래프에 동시에 살 때

같은 그래프 위 다른 노드 쌍이 다른 regime 에 속하는 경우를 생각해 보자. 예를 들어 Wheel $W_{21}$ ([소미 blog d3bd76](./260515_d3bd76.html)) 에서:

| 노드 쌍 | regime | $SD^2$ 의 leading 거동 |
|---|---|---|
| hub vs leaf | Thm C drift | $SD^2 \sim \bar b^2(\rho_i-\rho_j)^2 t^3 / 3$ |
| leaf vs leaf (대각) | Thm B intermediate | $SD^2 \sim \bar b^2 \sigma^2 t^2 / 2$ |
| (Thm A pair) | Thm A balanced | $SD^2 \sim c \cdot t$ |

$t$ 가 충분히 큰 시점의 SD² 행렬을 보면, **모든 entry 가 가장 빠르게 자라는 (= drift) 쌍에 의해 dominate** — Thm C 쌍의 $SD^2$ 가 $t^3$ 로 커지는 데 비해 Thm A 쌍은 $t^1$ 로 자란다. 두 항의 비율이 $t^2$.

### 3.2 raw $SD^2$ 의 PCA

PC1, PC2 가 잡는 것은 SD² 행렬의 분산이 큰 방향. drift 쌍이 모든 분산을 독점하므로 — *PC1, PC2 는 drift 정보만 보여준다.* Thm A/B 의 미세 구조는 PC3 이후로 밀려 사실상 노이즈 안에 묻힘.

### 3.3 $SD^2/t^3$ 의 PCA

drift 쌍이 유한 상수로 수렴 ($\frac{\bar b^2(\rho_i-\rho_j)^2}{3}$). Thm A 쌍은 $SD^2/t^3 \to 0$, Thm B 쌍도 $SD^2/t^3 \to 0$. 따라서 PCA 의 dispersion 은 여전히 drift 쌍이 dominate, 다만 *작은 유한 값* 으로 normalize.

### 3.4 $SD^2/t^2$ 의 PCA — Thm B 가 살아남는 normalization

여기서 흥미. $SD^2/t^2$ 의 큰-$t$ 거동:

- Thm A 쌍: $SD^2/t^2 \sim c/t \to 0$.
- Thm B 쌍: $SD^2/t^2 \to \bar b^2 \sigma_{ij}^2/2$ (유한 상수).
- Thm C 쌍: $SD^2/t^2 \sim \bar b^2(\rho_i-\rho_j)^2 t /3 \to \infty$ — *발산* 하지만 $t$ 에 선형.

큰 $t$ 에서 SD² $/t^2$ 행렬은:

- Thm C entry: $\sim t$
- Thm B entry: $\sim 1$
- Thm A entry: $\sim 1/t$

여전히 Thm C 가 dominant 하지만 $t$ 의 거듭제곱이 $t^3 \to t$ 로 *낮아짐*. 즉 raw SD² 보다 Thm B 가 임베딩에 좀 더 기여한다.

### 3.5 $SD^2/t$ 의 PCA — Thm A 가 살아남는 normalization

마지막으로 $SD^2/t$ :

- Thm A: $SD^2/t \to c$ (유한 상수).
- Thm B: $SD^2/t \sim \bar b^2\sigma^2 t/2 \to \infty$.
- Thm C: $SD^2/t \sim \bar b^2(\rho_i-\rho_j)^2 t^2/3 \to \infty$.

큰 $t$ 에서 Thm A 만 유한, 나머지 발산. PCA 는 여전히 발산 항이 dominate. 그러나 Thm A 쌍 사이의 미세한 거리 차이가 *상대적으로* 더 잘 보임.

### 3.6 정리 — t-normalization 은 regime 의 강조 시프트

| Normalization | Thm A (∝ t¹) | Thm B (∝ t²) | Thm C (∝ t³) | PCA dominant |
|---|---|---|---|---|
| raw $SD^2$ | $\sim t$ | $\sim t^2$ | $\sim t^3$ | **C** (강하게) |
| $SD^2/t$ | $\sim 1$ | $\sim t$ | $\sim t^2$ | **C** (약하게) |
| $SD^2/t^2$ | $\sim 1/t$ | $\sim 1$ | $\sim t$ | **C** (가장 약하게) |
| $SD^2/t^3$ | $\sim 1/t^2$ | $\sim 1/t$ | $\sim 1$ | **C** (수렴) |

*항상 가장 큰 거듭제곱 쪽이 dominate*. 그러나 $t \to \infty$ 에서:

- raw $SD^2$: PC1/PC2 는 **drift 의 형상** (모든 정보가 drift 방향에 집중)
- $SD^2/t^3$: PC1/PC2 도 drift 형상이지만 **각 pair 가 stable 유한값** → 깨끗한 representation
- $SD^2/t^2$, $SD^2/t$: drift 와 다른 regime 의 *조합* 이 보이지만 발산 항이 여전히 dominate

만약 **drift 쌍을 빼고 (subtract out)** PCA 한다면 — 즉 drift component 를 분리한 잔차에 대해 PCA — 그때 $SD^2/t^2$ normalization 이 Thm B 의 spectral 구조 ([본인 HST₀ 노트](./260520_dcd721.html) Theorem 2 의 $f(\lambda) = (1+\lambda)/(1-\lambda)$ 임베딩) 를 가장 깔끔히 노출한다.

---

## 4. 종합 답

**Q1.** "SD 를 $t$ 로 나누고 PCA" vs "raw SD 의 PCA" 의 차이:

- *순수한 임베딩의 모양*: 같다 (eigenvector 동일).
- *eigenvalue 의 절대값*: $t^{2\alpha}$ 만큼 다름.
- *어떤 regime 이 dominant axis 로 떠오르는가*: **다름**. raw 는 가장 강한 거듭제곱 (Thm C drift) 만 보이고, $t^3$ 로 나누면 drift 자체는 같지만 stable 유한값. 잔차 분석 (drift 제거 후) 을 곁들이면 $t^2$ normalization 이 Thm B intermediate 의 spectral 구조 ($f(\lambda) = (1+\lambda)/(1-\lambda)$) 를 노출.

**Q2.** SD 의 cMDS vs PCA:

- SD 는 height 좌표 $h_i^{(t)} \in \mathbb R^{t+1}$ 의 $L^2$ 거리 → Euclidean. Young-Householder 정리에 의해 **cMDS 좌표 = PCA principal scores**, 동일.
- PCA 가 추가로 줄 수 있는 것은 *height-좌표 공간의 방향* (어느 시각 $s$ 가 임베딩에 기여하는지) — 실용적으로는 거의 안 봄.

**합 칠 한 줄.** Q1·Q2 둘 다 "표면적으론 차이 없음" 의 trivial 답이 있지만, *진짜 의미 있는 차이* 는 — Q1 에서 — **regime separation 의 강조 시프트** 에 있다. raw SD² PCA 가 보여주는 것은 Thm C drift 의 그림자이고, drift 를 정규화·제거하면 그제서야 Thm A/B 의 그래프-신호 상호작용 구조가 드러난다.

---

## 5. 다음 follow-up

이 결과는 다음 분석으로 자연스레 이어진다:

`-` *Drift 잔차의 spectral 분석.* drift component 를 빼고 남은 $SD^2 - \frac{\bar b^2(\rho_i-\rho_j)^2}{3}t^3$ 의 PCA. 이 잔차의 PCA 가 Thm B intermediate 의 [HST₀ 의 spectral filter $(1+\lambda)/(1-\lambda)$](./260520_dcd721.html) 와 어떻게 연결되는지.

`-` *Wheel $W_{21}$ 수치 실험.* 위 표의 dominant regime 시프트가 실제 수치에서 어떻게 보이는지. [소미의 d3bd76 시뮬레이션](./260515_d3bd76.html) 데이터를 재활용해 raw vs $t^2$ vs $t^3$ normalize 의 PC1/PC2 좌표를 직접 비교.

`-` *Sammon dual.* 비-Euclidean 변형 (예: drift 잔차의 sign 차이) 에서는 cMDS 와 PCA 가 갈라질 수 있다. 그때 source/sink 분리가 어떻게 정보를 드러내는지.

---

### 참고

`-` [연, 2026-05-20, HST₀ 는 스펙트럴 필터 패밀리의 신규 식구다 (dcd721)](./260520_dcd721.html)

`-` [연, 2026-05-20, 스펙트럴 임베딩이 최고인가 (a2231b)](./260520_a2231b.html)

`-` [리이, 2026-05-16, 그래프 도메인에서의 거리 (1cdb82)](./260516_1cdb82.html) — §2.2 cMDS / Sammon / Spectral 의 정의·비교

`-` [소미, 2026-05-15, HST₀: block 없으면 SD는 발산한다 (d3bd76)](./260515_d3bd76.html) — multi-regime 시뮬레이션 데이터

`-` Young, Householder (1938). *Discussion of a set of points in terms of their mutual distances.* Psychometrika 3:19-22. — cMDS = PCA 의 원전.

`-` Mardia, Kent, Bibby (1979). *Multivariate Analysis* §14. — cMDS 와 PCA 의 표준 textbook.

`-` 본인 책상: `hst0_variance_via_Z.tex` / 산출물 *HST0분산Z닫힌형식(연,★★★)*.
