---
title: 연구 ▷ HST ▷ 예제2 재검토 — HST 는 Euclidean 과 얼마나 다른가?
author: 클로드
date: 04/25/2026
draft: false
output-file: 260425_6dc8c6.html
---

# 0. 한 줄 요약

> ex2 (Gaussian ring + cylinder 신호 $y=\pm 3$) 에서 **HST 의 cMDS 임베딩 모양이 pure Euclidean ($|y_i - y_j|^2$ 만) 의 임베딩과 Procrustes 거리 0.003 으로 거의 같다**. 즉 신호가 ring 보다 훨씬 강한 이 세팅에서 HST 는 좌표 수준에서 Euclidean 과 사실상 구분이 안 됨. 이 사실은 ex2 가 HST 의 가치를 보이는 demo 로서 부적합할 수 있음을 시사.

---

# 1. 출발 — diffusion distance 정의 오류 발견 ([260425 a](260425_*.html))

ex2 의 baseline 비교에 쓰던 `diffusion_distance` 가 표준 Coifman-Lafon 정의가 아니라는 것이 먼저 드러남. 정규그래프(Gaussian ring) 에서는 우연히 스케일만 다른 결과를 줘서 **임베딩 모양은 보존**, 그래서 ex2 그림 자체는 다행히 영향이 없었음. 다만 편차분석의 **외재적 $\Delta^{\text{ext}}$ 절대값** 은 영향. 이 발견 후 잠재적으로 ex2 전반의 해석을 다시 점검할 동기가 생김.

→ 자연스럽게: "ex2 에서 그동안 그린 임베딩들이 정말 의미 있는가? 각 metric 이 어떻게 다르게 보이는가?" 라는 질문이 등장.

---

# 2. 거리 metric 9 종 비교 ([260425_예제2_거리비교.py](../260425_예제2_거리비교.py))

9 개 metric 에서 각각 **distance matrix 와 cMDS-3D 임베딩** 을 한 그림에 정리:

1. Euclidean (signal-only) — $|y_i - y_j|^2$
2. Shortest path² on Gaussian W
3. Diffusion (π-weighted Coifman-Lafon, $t=4$)
4. Heat kernel ($e^{-2L}$)
5. PPR ($\alpha=0.15$)
6. Effective resistance
7. Biharmonic
8. Spectral embedding (top-5 nontrivial $\psi_k$)
9. **HST** ($\overline{\mathrm{SD}}^2(\tau)$, $\tau=500{,}000$)

**관찰** (gauss W mode):

- **그래프-only metric 7 종** (SP, Diffusion, Heat, PPR, Resistance, Biharmonic, Spectral) — 거의 동일한 **ring 모양**. 정규 ring 의 spectrum 이 단순해서 모든 spectral filter 가 비슷한 결과.
- **Euclidean** — 두 cluster (y=±3) 만 분리된 1D 구조.
- **HST** — cylinder 모양 (cluster 두 개 + 각각 안에 ring substructure).

이 시점에서 결론은 명료해 보임: HST 만이 신호와 그래프를 **결합**해서 cylinder 모양을 만들어 내는 유일한 metric. → "HST 의 가치 = 이 cylinder shape" 라는 가설.

결과: [results/figures/260425_예제2_거리비교.py/comparison_gauss.png](../results/figures/260425_예제2_거리비교.py/comparison_gauss.png)

---

# 3. 이 cylinder 가 *진짜* HST 만의 것인가?

위 가설을 직접 검증: **유클리드 + ring 의 가중합으로 cylinder 모양이 나오는가?** ([260425_예제2_선형결합실험.py](../260425_예제2_선형결합실험.py))

$$D_{\text{lin}}(\alpha) \;=\; (1-\alpha)\,\hat D^E + \alpha\,\hat D^G$$

각 metric 은 max-normalize 후 결합. $D^G$ 후보 3 종: ring angular dist², shortest path², effective resistance. 51 개 $\alpha \in [0, 1]$ 에 대해 cMDS-3D 후 HST 임베딩과 **Procrustes disparity** 측정.

## 3.1 결과 — 충격적

| $D^G$ | 최적 $\alpha^*$ | min disparity |
|:---|---:|---:|
| ring angular dist² | **0.02** | **0.0029** |
| shortest path² | **0.02** | **0.0029** |
| effective resistance | **0.04** | **0.0031** |

세 후보 모두 $\alpha^* \approx 0$ 부근에서 최소. **$\alpha = 0$ 은 pure Euclidean** 이다. 즉 *거의 신호 단독* 으로 HST cylinder shape 를 0.003 disparity (사실상 동일) 로 재현.

결과: [results/figures/260425_예제2_선형결합실험.py/lincomb_to_hst.png](../results/figures/260425_예제2_선형결합실험.py/lincomb_to_hst.png)

## 3.2 왜 이렇게 나오나 — eigenvalue dominance

cMDS 의 분산 분해:
$$\text{Var}(\text{HST embedding}) = \sum_k \lambda_k(D^{\text{HST}}_{\text{centered}})$$

신호 split (y=±3) 과 ring radius (=1) 의 squared 비율 $\sim 36 : 1$. 따라서 $D^{\text{HST}}$ 의 첫 eigenvalue $\lambda_1$ (신호 축) 이 나머지 $\lambda_2, \lambda_3$ (ring 축들) 을 **압도**.

Procrustes 는 standardize (Frob 정규화 = 모든 $\lambda$ 의 합으로 나눔) 후 비교 → dominant eigenvector 만 정렬되면 disparity 는:

$$\text{disparity} \;\approx\; \frac{\lambda_2 + \lambda_3 + \cdots}{\sum_k \lambda_k}$$

이게 작으면 disparity 도 작음. Pure Eucl 은 rank-1 (signal 축만), HST 는 dominant signal + 작은 ring perturbation. 두 임베딩 모두 **signal split 이 좌표 분산의 95%+ 를 차지** → Procrustes 가 둘을 사실상 같다고 판정.

## 3.3 Procrustes 가 *놓치는* 것

작은 disparity 가 "HST 가 가치 없다"는 결론으로 곧장 가지는 않음. Procrustes 는 좌표 잔차 metric 이므로 다음을 무시:

- **Topology (H1)**: HST 의 ring substructure 는 closed loop → $H_1 \neq 0$. Pure Eucl 은 두 점 cluster → $H_1 = 0$. **위상적으로 완전히 다름**.
- **Cluster 내부 정렬**: HST 는 같은 cluster 내 노드를 ring angle 순으로 배치, Eucl 은 noise 수준 (~0.04) 로 무작위.
- **Multi-scale 능력**: HST 는 $\tau$ 로 scale 조절 (early transient vs Harris 극한), Eucl 은 단일 scale.

→ **"같은 모양"** 이라고 말할 수 있는 metric 과 **"의미 있게 다른 정보"** 를 잡는 metric 은 다름.

---

# 4. ex2 가 HST demo 로서 적합한가?

정직하게 말하면 **이 세팅에선 부적합** 으로 보임:

- $y = \pm 3$ 은 ring radius 1 보다 신호 magnitude 가 압도적
- 결과: HST embedding 의 분산이 **신호 axis 가 절대 지배**
- → cylinder shape 의 ring 부분이 **시각적으론 보이지만 양적으론 미미**
- → linear combination (사실상 pure Eucl) 으로 거의 재현 가능

신호 magnitude 와 ring scale 이 **비등** 하거나, 신호와 그래프 정보가 **이질적으로** 결합해야만 HST 가 비자명한 가치를 가짐.

## 4.1 신호 약화 시뮬레이션 가설

$y$ 의 amplitude 를 $\pm 3 \to \pm 0.5$ 또는 $\pm 1$ 로 낮추면:
- $\lambda_1$ (signal axis) 가 $\lambda_2, \lambda_3$ (ring) 와 비등 → cylinder shape 의 양 축이 비교 가능한 크기.
- pure Eucl 만으론 cluster 분리도 약해짐 (signal split 이 작으니까).
- HST 가 비로소 *intermediate* metric 으로 의미 있을 가능성.
- $\alpha^*$ 가 $(0, 1)$ 내부로 이동하고 disparity 가 충분히 커질 것.

→ **추가 실험 우선순위 #1**.

## 4.2 다른 예제군

| 예제 | HST 가치 가능성 | 이유 |
|:---|:---|:---|
| 예제2 ($y=\pm 3$) | **낮음** (현재 결과) | 신호 dominance |
| 예제2 ($y=\pm 0.5$) | **미지** | 균형 regime 가능 |
| 예제5 / 논문 fig 5 cliff | **높음** | 신호 불연속 + 그래프 연속, 결합 정보 필수 |
| 예제6 parity ($f=(-1)^i$) | **높음** | 신호·그래프가 *반대* — linear combo 는 평균 |
| 예제7 twin rings | **높음** | phase shift 의 비선형 coupling |
| 예제8 directed cycle | **높음** | classical metric 이 방향성 무시, HST 만 기록 |

ex2 만 보고 HST 를 평가하는 것은 위험. **ex5–8 에서 동일 분석 (Procrustes sweep) 반복 필요**.

---

# 5. 9. 통합 요약표 의 분류 재확인

[260425 그래프 도메인 거리 정리](260425_*.html) 에서 만든 12 거리 표에서 HST 를 **"signal + graph 동시"** 카테고리에 단독 배치했었음. 본 검토를 통해 좀 더 정직한 위치 정리:

- **Regime A (신호 dominant)**: HST 는 사실상 Euclidean.
- **Regime B (그래프 dominant)**: HST 는 사실상 commute time / resistance.
- **Regime C (균형 또는 이질적 결합)**: HST 가 unique. ← *이 영역에서만* HST 가 카테고리 단독 배치를 정당화함.

**"signal + graph" 카테고리를 차지할 자격은 regime C 에서만 있다**.

→ 향후 모든 HST 결과를 보일 때 *어느 regime 에 속하는 example 인지* 명시해야 함. ex2 (y=±3) 는 명백히 regime A 였음에도 그동안 "signal-graph 결합" 의 demo 로 간주.

---

# 6. 열린 질문 (다음 단계)

**O1 (가장 시급).** ex2 신호 amplitude scan: $y$ 의 amplitude 를 $0.1, 0.3, 0.5, 1, 2, 3$ 으로 변화시키며 (a) Procrustes(α) 곡선의 변화, (b) $\alpha^*$ 의 위치, (c) min disparity 의 변화 추적. **regime A → C 천이가 어디서 일어나는가**.

**O2.** 같은 sweep 실험을 예제5–8 에 적용. 어느 example 에서 $\alpha^* \in (0.2, 0.8)$ 이고 disparity 가 의미 있게 큰지.

**O3.** Procrustes 대체 metric:
- **PD H1 bottleneck** — pure Eucl ($H_1 = 0$) vs HST ($H_1 \neq 0$) 의 위상 차이 정량화
- **Spectral entropy of cMDS eigvals** — 차원성 (1D vs 2D) 측정
- **Within-cluster ordering correlation** — Spearman/Kendall on cluster 내부 ring angle

**O4.** $\overline{\mathrm{SD}}^2$ 의 eigenvalue spectrum 직접 분석 — $\lambda_1 / \sum \lambda_k$ 비율로 "신호 dominance index" 정의. 이 index 가 작을 때만 HST 가 가치 있음을 정량 기준으로 제시.

**O5.** [편차분석 포스트](260424_a5d784.html) 의 $\tau^*$ 결과를 본 검토 관점에서 재해석. ex2 에서 max $\Delta_F$ 가 작게 나왔던 것 (고작 75 정도) 이 사실 "HST 가 linear path 위에서 거의 못 벗어났다" 는 같은 현상의 다른 표현.

---

# 7. 결론

- **ex2 (y=±3) 는 HST 의 강점이 드러나기 어려운 setup** — 신호 magnitude 가 ring scale 을 압도해서 cMDS 의 dominant axis 가 신호 split 으로 거의 독식. 결과적으로 **pure Euclidean 임베딩과 Procrustes 거리 0.003** 으로 사실상 같음.
- **HST 의 "cylinder shape" 도 pure Eucl + 작은 perturbation 으로 분해 가능** — Procrustes 관점에선.
- **단, 위상 (H1 loop), cluster 내부 ordering, multi-scale 등 *Procrustes 가 잡지 않는 차원* 에서는 여전히 HST 가 unique** — 이걸 정량화하려면 다른 metric 필요.
- **ex2 는 demo 로서 약함** — HST 를 정당화하려면 regime C (균형/이질적 결합) example 로 옮겨가야 함. 예제5–8 은 후보.
- 본 검토는 **HST 무용론** 이 아님. HST 가 가치 있는 *영역* 의 정의를 좁히고, demo example 의 적절성에 의문을 제기.

---

# 참고 (cross-link)

- [260424 $\tau^*$ 편차분석](260424_a5d784.html) — 본 검토에서 재해석 필요한 ex2 편차 결과
- [260425 diffusion distance 정의 오류](260425_*.html) — 시작점이 된 오류 발견
- [260425 그래프 도메인에서의 거리](260425_*.html) — 9-12 거리 family 정리
- [260424 $\overline{\mathrm{SD}}^2$ 수렴 증명](260424_cdbaf0.html) — HST 의 Harris 극한
- [보충아이디어](260208_21543c.html) — 본 검토와 연결될 추후 방향들

---

# 부록: 핵심 수치 (재현 가능)

```
Setup:  n=60, ring on unit circle, gauss W (θ=0.35)
Signal: f2 cylinder, y = ±3 (noise 0.1)
HST:    τ=500,000, b=0.05 → snowdistance² = D_HST

Procrustes sweep (α: 0..1, 51 points):
  D^G = ring angular dist²    : α* = 0.02, min disparity = 0.0029
  D^G = shortest path²        : α* = 0.02, min disparity = 0.0029
  D^G = effective resistance  : α* = 0.04, min disparity = 0.0031

Pure Euclidean (α=0) Procrustes to HST: ≈ same as above
Pure ring (α=1) Procrustes to HST: ≈ 0.93 (예상치, 향후 확인)
```
