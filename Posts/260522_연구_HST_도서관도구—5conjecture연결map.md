---
title: 연구 ▷ HST ▷ 도서관 도구 — 5 conjecture 연결 map
author: 연
date: 05/22/2026
draft: false
output-file: 260522_db44b2.html
fontsize: 0.85em
---
<!-- 소유권자: 연 | 사용자: 연 -->

```{=html}
<style>
.math.display { font-size: 0.9em; text-align: left; }
mjx-container[display="true"] { text-align: left !important; margin-left: 0 !important; }
.katex-display { text-align: left !important; }
.katex-display > .katex { text-align: left !important; }
.callout { font-size: 0.9em; }
</style>
```

회사 도서관에 비치된 책 (Brémaud, Meyn--Tweedie, Durrett, Djurić--Richard) 안에서, HST 본업과 직접 연결되는 정리·보조정리·도구를 한 페이지로 정리한다. 본인의 *cancellation tier framework* 의 다섯 conjecture (Tier 1.5 cross-over / Tier 2 closed form / Tier 3 emptiness / Hodge $k$-cochain / filter--tier correspondence) 가 어느 도구를 요구하는지까지 함께 표시한다. 작성 동기는 단순한데, 본인 노트들이 다 적힌 뒤 보니 "이건 결국 책의 어디 정리를 HST 에 옮긴 것일 뿐" 인 부분과 "책엔 비슷한 것조차 없고 새로 만들어야 하는 부분" 이 섞여 있었다는 점이다. 둘을 같은 표 안에서 분리하면, 다음 한 달 동안 어디를 정독하고 어디를 새로 쓸지가 보인다.

# A. 이미 본인 `.tex` 에 박힌 도구 (확정)

지금까지의 HST 노트 (`hst0_variance_via_Z.tex`, `sd_t3_upper_bound.tex`, `cancellation_tier_framework.tex`) 안에 이미 인용되어 있는 도구들. 이들은 *책에서 가져와 HST 표기로 옮긴* 것이지 새로 만든 것이 아니다.

| 도구 | 형태 | 출처 | 본업 연결 |
|---|---|---|---|
| **Brémaud (6.31)** | $\sigma^2 = \langle f, f\rangle_\pi - \langle f, \Pi f\rangle_\pi + 2\langle f, (\mathbf{Z}-\mathbf{I})f\rangle_\pi$ | Brémaud Ch.6.3.3 | Tier 1 (intermediate) leading constant |
| **Brémaud (6.34)** | $\sigma^2 = \sum_k \frac{1+\lambda_k}{1-\lambda_k}\pi^2(v_k(i)-v_k(j))^2$ | Brémaud Ch.6.3.3 | Tier 1 spectral filter — HST$_0$ 가 신규 식구로 등재 |
| **Toeplitz lemma** | $a_s \to L \,\Rightarrow\, \frac{\sum w_{ts} a_s}{\sum w_{ts}} \to L$ | Knopp §43, 1928 | Thm A/B/C 의 $\SD^2/t^k$ 점근 상수 정리 |
| **Cesàro means** | $a_s \to L \,\Rightarrow\, \bar a_t \to L$ | Brémaud Ch.4 | Toeplitz 의 특수 케이스, $\rho_i(t) \to \rho_i$ |
| **Fundamental matrix $\mathbf{Z}$** | $(\mathbf{I}-\mathbf{P}+\mathbf{1}\boldsymbol\pi^\top)^{-1}$ | Brémaud Ch.6.3 | (6.31)·(6.34) 의 backbone, Conj 2 $\mathbf{Z}^{\mathrm{HST}}$ 의 출발점 |

이 다섯 개가 Tier 1 (Thm B) 의 닫힌 형식을 받쳐주는 *전체* 도구다. Brémaud 한 권만으로 HST$_0$ 의 intermediate regime 이 끝난다는 것이 본인 `hst0_variance_via_Z.tex` 의 결론이었다.

# B. 본인이 아는데 아직 안 박은 도구 (5 conjecture 에 필요)

본인이 책에서 본 적은 있고 어디에 쓸지 머릿속에 있지만, 아직 HST 노트로 옮기지 않은 도구들. 각 항목 옆에 어느 conjecture 와 연결되는지 표시한다. **B5 가 cancellation framework 의 단기 핵심**: Brémaud (6.31) 은 *time-homogeneous* $\mathbf{P}$ 를 요구하는데 full HST 는 $\mathbf{P}(s)$ 가 height-history 에 의존하는 비균질 체인이고, mixing CLT 가 그 간격을 메우는 정본 도구다.

| 도구 | 출처 | 어떤 conjecture 에 필요 |
|---|---|---|
| **Foster's theorem** (Foster--Lyapunov drift) | Brémaud Ch.5.3, Meyn--Tweedie Ch.11 | Conj 2 Tier 2 closed form (에이가 이미 통합증명 Lemma 3 에서 사용) |
| **Doeblin minorization** | Brémaud Ch.5.4, Durrett Ch.5.8.3 | Conj 2·3 — exponential mixing 의 정량 표현 |
| **Geometric drift / $V$-uniform ergodicity** | Meyn--Tweedie Ch.15 | Conj 3 Tier 3 emptiness 의 martingale lower bound |
| **Harris recurrence** (general state space) | Meyn--Tweedie Ch.9, Durrett Ch.5.8 | Conj 2 non-homogeneous cocycle 의 base |
| **CLT for stationary sequences with mixing** | Durrett Ch.8.3 | Conj 2 — Brémaud (6.31) 의 time-비균질 일반화 |
| **Donsker's theorem** | Durrett Ch.8.1 | Conj 1 Tier 1.5 cross-over 의 functional CLT 형식 |
| **Birkhoff's ergodic theorem** | Durrett Ch.6.2, Brémaud Ch.4 | $\rho_i$ 존재성 (오큐리 Occupation SLLN 의 backbone) |
| **Kingman's subadditive ergodic theorem** | Durrett Ch.6.4 | Conj 2 non-homogeneous cocycle limit |
| **Donsker--Varadhan rate function** | Durrett Ch.2.7 (기초), Dembo--Zeitouni (정본) | Conj 2 approach (ii) — empirical measure LDP 의 variational form |

# C. 본인이 잘 안 본 도구 (도서관에 있으나 정독 X)

본인이 graph distance 표를 만들 때 표 안에는 넣었지만, HST 와의 연결까지는 짚어보지 않은 도구들. **C5 (spectral filter family) 가 cancellation framework Conj 5 의 직접 데이터**: 책 안에 Diffusion·Heat·PPR·Resistance·Biharmonic 의 필터 형식이 다 정리되어 있고, 각 멤버가 $\lambda \to 0$, $\lambda \to 1^-$ 에서 어떻게 행동하는지를 분석하면 어느 tier 에 속하는지가 결정된다.

| 도구 | 출처 | 본업 관련 |
|---|---|---|
| **Itô's formula** | Durrett Ch.7.6 | $\xi(s)$ 의 SDE 연속시간 한계 (HST 의 scaling limit) |
| **Feynman--Kac formula** | Durrett Ch.9.4 | $\SD^2$ 의 generator-level 표현 (open) |
| **Multidim. Brownian motion + PDE** | Durrett Ch.9 | $h(\cdot, t)$ 의 continuum-graph 한계 |
| **Effective resistance / Commute time** | Brémaud Ch.5.6, Djurić--Richard Ch.4 | Resistance distance — graph distance 비교 base |
| **Spectral filter family** (Diffusion / Heat / PPR / Resistance / Biharmonic) | Djurić--Richard Ch.3·4 | Conj 5 filter--tier correspondence 의 직접 데이터 |
| **Hodge Laplacian $\Delta_k$ on simplicial complexes** | 도서관 외 (Lim 2020, Schaub 2020) | Conj 4 $k$-cochain extension |
| **Renewal theory** | Durrett Ch.2.6 | block 사이 간격 분포 — Tier 1 → Tier 2 downshift 의 정량 |
| **Local limit theorems** | Durrett Ch.3.5 | $\xi(s)$ 의 local density — Tier 2 $c_{ij}$ 의 Gaussian form 가능성 |

# D. 우선순위 (대표 보고용)

5 conjecture 중 어느 것을 단기·중기·장기 로 풀 수 있는지를 도구 가용성 기준으로 정리한다. 도구가 도서관에 *이미 있고* 본인이 *아는* 경우가 빠르고, 도구가 외부 ref (Lim 2020, Dembo--Zeitouni 등) 인 경우는 그 자료를 먼저 도서관에 추가해야 한다.

- **즉시 활용 (단기)**: **B5 — Durrett Ch.8.3 mixing CLT** 가 Tier 1 식의 time-비균질 일반화의 가장 짧은 path. Brémaud (6.31) 이 *iid · time-homogeneous* 만 다루는데, full HST 는 $\mathbf{P}(s)$ 비균질이므로 mixing 조건으로 일반화한 CLT 가 정본 도구. 이걸 정독하면 Conj 2 의 (i) cocycle 접근의 첫 한 페이지가 나온다.
- **단기 활용**: **B1·B2 — Foster + Doeblin** 은 에이가 이미 본업 통합증명에 박았다. 본인 노트로 끌어오기만 하면 Conj 2 ·3 의 base 가 채워진다.
- **중기 활용**: **B8 — Kingman subadditive ergodic theorem** 이 non-homogeneous cocycle limit 의 정본 도구. Conj 2 의 cocycle 정의 자체를 정당화한다.
- **장기 활용**: **B9 — Donsker--Varadhan rate function** 은 Conj 2 의 가장 야심찬 형식 (variational $c_{ij}$). Dembo--Zeitouni 정본을 도서관에 추가해야 한다.
- **C5 spectral filter family** 는 Conj 5 검증에 직접 필요한 데이터. Djurić--Richard Ch.3·4 정독이 다음 본업으로 가장 적절하다 — 이론 + 시뮬레이션 양쪽 가능.

# E. 한 줄 결론

본인 cancellation framework 의 5 conjecture 중 **Conj 1·2 는 도서관 도구로 풀 수 있는 path 가 보이고, Conj 4·5 는 도구 추가 (Hodge ref) 가 선행되어야 한다**. Conj 3 (Tier 3 emptiness) 은 martingale lower bound + spectral gap 의 표준 조합이라 본인 짐작으로는 한 페이지짜리 정리로 끝날 가능성이 있다.

# 참조

[1] P. Brémaud, *Markov Chains: Gibbs Fields, Monte Carlo Simulation and Queues* (2nd ed.), Springer, 2020 — Ch.4 (Cesàro), Ch.5 (Foster, Doeblin), Ch.6.3 (fundamental matrix Z), Ch.6.3.3 (eq. 6.31, 6.34).

[2] S. P. Meyn and R. L. Tweedie, *Markov Chains and Stochastic Stability* (2nd ed.), Cambridge University Press, 2009 — Ch.9 (Harris recurrence), Ch.11 (Foster--Lyapunov), Ch.15 ($V$-uniform ergodicity).

[3] R. Durrett, *Probability: Theory and Examples* (5th ed.), 2019 — Ch.2.4 (SLLN), Ch.2.6 (Renewal), Ch.2.7 (LDP basics), Ch.3.5 (Local limit), Ch.5.8 (general state space), Ch.6.2 (Birkhoff), Ch.6.4 (subadditive), Ch.7.6 (Itô), Ch.8.1 (Donsker), Ch.8.3 (mixing CLT), Ch.9.4 (Feynman--Kac).

[4] P. M. Djurić and C. Richard (eds.), *Cooperative and Graph Signal Processing*, Academic Press, 2018 — Ch.3·4 (spectral graph theory + graph distances).

[5] A. Dembo and O. Zeitouni, *Large Deviations Techniques and Applications* (2nd ed.), Springer, 1998 — 도서관 추가 필요 (Conj 2 approach (ii)).
