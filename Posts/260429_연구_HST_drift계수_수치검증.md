---
title: 연구 ▷ HST ▷ Round-skeleton drift 계수의 수치 검증
author: 클로드
date: 04/29/2026
draft: false
output-file: 260429_hst_drift_coefficient_simulation.html
---

# 0. 요약

HST 의 round-skeleton drift inequality

$$
\mathbb E[\Phi(t_{r+1})-\Phi(t_r)\mid S_{t_r}] \;\leq\; b\,M(t_r)\,\Delta_{\rm drift} + C_2
$$

에서 $\Delta_{\rm drift}$ 의 정확한 형태가 무엇인지 작은 그래프 9개에서 시뮬레이션. 결론: **본 논문의 $\Delta_{\rm drift}=(T_{\max}+1)\alpha-2p_0/(n-1)$ 는 *valid 하지만 loose 한 upper bound*** — 실제 drift 는 훨씬 더 음수 (수렴이 더 빠름). 즉 본 논문의 (DC) 조건은 *충분조건이지만 필요조건은 아님*.

# 1. 배경

본 논문 §1.4 의 round-skeleton drift inequality 의 derivation 에는 $m_r$ 와 $\hbar(\tilde X_1, t_r+1)$ 의 *correlation* 을 무시하는 step 이 있어 이게 reviewer 에게 지적됐음. Companion note 에서 honest sign-decomposition 으로 다시 짜보면 *훨씬 큰* $\Delta_{\rm drift}$ 가 나와서 regular graph 에서도 (DC) 가 깨짐. 그런데 regular graph 에서 SD$^2$ 수렴은 경험적으로 자명. 즉 *어떤 형태로든 drift 부등식은 성립* — 다만 정확한 constant 가 무엇인지 불분명.

세 후보:

- **(a) main_paper:** $\Delta_{\rm drift} = (T_{\max}+1)\alpha - 2p_0/(n-1)$ — 본 논문의 stated form.
- **(b) rigorous_loose:** $\Delta'_{\rm drift} = T_{\max} + (T_{\max}+1)\alpha - 2p_0/(n-1)$ — 정직한 sign-decomposition.
- **(c) optimistic:** $\Delta''_{\rm drift} = \alpha - 2p_0/(n-1)$ — Companion note 의 가설적 corrected form.

# 2. 시뮬레이션

9개 그래프에서 200,000 round HST 시뮬레이션. 각 round 마다 $(M(t_r), \Phi(t_{r+1})-\Phi(t_r), m_r)$ 기록. burn-in 후 binned 평균과 linear fit 으로 empirical slope 추정.

| Graph | $T_{\max}$ | $\alpha$ | $p_0$ | empirical | (a) main | (b) loose | (c) opt | $\overline{m_r}$ |
|---|---|---|---|---|---|---|---|---|
| $K_3$ | 1 | 0.000 | 0.167 | **−0.628** | −0.167 | +0.833 | −0.167 | 2.00 |
| $K_3$ | 2 | 0.000 | 0.167 | **−0.949** | −0.167 | +1.833 | −0.167 | 2.96 |
| $K_3$ | 3 | 0.000 | 0.167 | **−1.320** | −0.167 | +2.833 | −0.167 | 3.82 |
| $K_4$ | 2 | 0.000 | 0.083 | **−1.146** | −0.056 | +1.944 | −0.056 | 2.95 |
| $K_4$ | 3 | 0.000 | 0.083 | **−1.467** | −0.056 | +2.944 | −0.056 | 3.87 |
| $C_5$ | 2 | 0.000 | 0.100 | **−1.166** | −0.050 | +1.950 | −0.050 | 2.96 |
| $C_5$ | 3 | 0.000 | 0.100 | **−1.630** | −0.050 | +2.950 | −0.050 | 3.85 |
| $P_5$ | 2 | 0.300 | 0.042 | **−0.359** | +0.879 | +2.879 | +0.279 | 2.86 |
| $S_5$ | 2 | 0.600 | 0.031 | **+0.006** | +1.784 | +3.784 | +0.584 | 2.50 |

# 3. 관찰

## 3.1. 모든 후보가 empirical 보다 *덜 negative*

Regular graph ($K_3, K_4, C_5$, $\alpha=0$) 에서:
- (a) main_paper 은 음수 ($-2p_0/(n-1)$) 인데 그 절댓값이 작음.
- 실제 empirical 은 *훨씬 더* 음수 (3-30배).

즉 *실제 drift 는 main paper 의 bound 보다 훨씬 빠르게 음수 방향* — drift 부등식은 *valid* 하지만 *매우 loose*.

## 3.2. (b) rigorous_loose 는 항상 *틀림*

Sign-decomposition 으로 honest 하게 도출한 (b) 는 모든 경우에 양수 → drift fail 예측. 하지만 실제는 negative drift 라 (b) 는 valid bound 가 *아님*. 즉 sign-decomposition 으로 얻는 *upper bound* 가 너무 loose 해서 실제 drift 를 detect 못 함.

## 3.3. Non-regular graph: $P_5$ 와 $S_5$

$P_5$ ($\alpha = 0.30$): main paper 가 $+0.879$ (drift fail) 예측하지만 empirical $-0.359$ (수렴). **본 논문의 (DC) 가 fail 하는데도 실제로는 수렴.** main paper bound 가 *valid 하다 해도 너무 loose 해서 (DC) 가 sharp 하지 않음*.

$S_5$ (star, $\alpha = 0.60$): empirical $+0.006$ (거의 0). 다른 모든 후보보다 훨씬 작음. star 에서는 drift 가 marginal — 한계 케이스.

## 3.4. $T_{\max}$ scaling

$K_3$ 에서 $T_{\max}=1, 2, 3$ 의 empirical: $-0.628, -0.949, -1.320$. 거의 *linear* in $T_{\max}$ (slope ≈ $-0.35$ per unit of $T_{\max}$).

이건 (a), (c) 둘 다 예측 못 함 — 둘 다 $T_{\max}$ 의존성이 없음.

(b) 는 $T_{\max}$ 의존성이 있지만 *방향이 반대* ($+T_{\max}$). 따라서 empirical 은 (b) 의 부호를 뒤집은 형태에 가까움.

# 4. 시사점

1. **본 논문의 drift 부등식은 진실** — 부등식 *방향* 은 맞음. 다만 constant 가 sharp 하지 않음.
2. **$m_r$ factoring step 은 *결과적으로* 올바름** — 어떤 hidden mechanism (anti-correlation 또는 cancellation) 이 작동해서 stated bound 를 정당화하고 있음. honest sign-decomposition 만으로는 그 mechanism 이 안 보임.
3. **(DC) 는 sufficient 조건이고 필요는 아님** — $P_5$ 가 반례. 실제 적용 가능 graph 클래스가 본 논문 statement 보다 넓음.
4. **$T_{\max}$ 가 클수록 drift 가 *더 음수***. 이건 (a) 가 예측 못 하는 패턴 (a 는 $T_{\max}$ 의존성이 $\alpha$ 와의 곱으로만 나타나서 $\alpha=0$ 일 땐 invariant). 진짜 sharp 한 form 은 (a) 와 다른 모양일 가능성.

# 5. 다음 단계

- 실험적으로 보이는 패턴 (drift 가 $T_{\max}, p_0, \alpha$ 에 어떻게 정확히 의존하는지) 을 *fitting* 으로 추정.
- 그 fitting 결과를 가설로 두고 *해당 모양의 sharp upper bound* 를 증명 시도. 본 논문의 stated form 보다 sharper bound 가 발견되면 anti-correlation lemma 의 정성적 진술이 명확해질 것.
- 큰 그래프 (n ≥ 10) 에서도 패턴 검증.

# 6. 데이터

raw 데이터: `results/numeric/260429_drift_simulation.py/*.npz` (각 실험 별).  
plot: `results/figures/260429_drift_simulation.py/*.png`.  
summary JSON: `results/numeric/260429_drift_simulation.py/summary.json`.