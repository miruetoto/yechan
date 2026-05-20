---
title: "연구 ▷ HST ▷ Kite 실험 (K_{2,20} + RIGHT)"
author: 소미
date: 05/20/2026
draft: false
output-file: 260520_e4f2a8.html
toc: false
page-layout: full
---
<!-- 소유권자: 소미 | 사용자: 소미 -->

## §1. 그래프 정의

K_{2,20} + RIGHT 양식. 양방향 그래프.

- **LEFT** (주황 □): 회색 노드 20개를 거쳐 ANCHOR 와 연결
- **MID** (회색 ○, 20개): LEFT 와 ANCHOR 사이 병렬 path 의 internal node
- **ANCHOR** (빨강 ○): hub. LEFT 와 RIGHT 양쪽 연결
- **RIGHT** (보라 □): ANCHOR 와 직접 연결 (internal 없음)

| 노드 | degree | 비고 |
|---|---|---|
| LEFT | 20 | hub (회색 20개 거침) |
| MID (×20) | 2 | LEFT, ANCHOR 양쪽 연결 |
| ANCHOR | 21 | LEFT 의 모든 회색 + RIGHT 1개 |
| RIGHT | 1 | ANCHOR 와 직접 |

|V|=23, |E|=41.

![](attachments/260520_e4f2a8_graph.png)

## §2. 시뮬 설정

- **signal** $f = 0$ (모든 노드)
- **b** = 0.05, $T_{\max}$ = 20
- **$\mu_0$** = uniform (degree-proportional)
- **τ** = $10^6$, **N_CKPT** = 30 (geomspace $1 \sim 10^6$, 초기 t 포함)
- random-step variant (b'_s ~ Unif(0, b) i.i.d.)

## §3. 결과

### 2x3 layout

- (0,0) **graph** — 본 그래프 시각화
- (0,1) **수렴 곡선** — SD²/t (Thm A, 보라), SD²/t² (intermediate, 회색), SD²/t³ (Thm B, 초록)
- (0,2) — [비움]
- (1,0) **spectral embedding** — 정규화 Laplacian 의 fiedler 등 eigvec 사용
- (1,1) **HST embedding** — √(SD²/t) 의 2D MDS
- (1,2) **resistance embedding** — pseudoinverse of $L$ 의 resistance distance

![](attachments/260520_e4f2a8_anim.gif)

### 분석

**Regime**: degree heterogeneity 매우 큼 (LEFT/ANCHOR deg ≈ 20, MID deg = 2, RIGHT deg = 1) → balanced 가정 위배 → **Thm B (drift) regime** 예상.

- **SD²/t**: 발산 (우상향)
- **SD²/t²**: 발산 (우상향)
- **SD²/t³**: 수렴 (안정) → **t³ scaling** 확인. mean off-diagonal SD² 가 t³ 으로 증가.

**HST embedding**: LEFT/ANCHOR 가 매우 가깝게 (둘 다 hub, 거의 같은 distance profile to MID). RIGHT 는 ANCHOR 근방. MID 20개는 LEFT-ANCHOR 사이 분포.

**Spectral**: hub 노드 (LEFT/ANCHOR) 가 양 끝에 분리 + MID 중앙 cluster. RIGHT 는 ANCHOR 근방. graph 의 구조 (path/hub) 잘 반영.

**Resistance**: Spectral 과 유사 양식 — resistance distance 가 spectral distance 와 closely related.

### Comparing distances

세 distance (HST / spectral / resistance) 모두 본 graph 의 **hub-MID-leaf 구조**를 다른 방식으로 표현. HST 는 시간 진화 양식 (random walk 위의 height accumulation) 이라 graph 외 dynamic 정보 포함.

## §4. 결론

K_{2,20} + RIGHT 양식은 degree heterogeneity 가 dominant 한 그래프로, HST 의 SD² 가 **t³ scaling** (Thm B drift regime) 으로 발산. embedding 양식에서 hub 와 leaf 의 분리가 명확. spectral / resistance 와의 비교에서 HST 가 dynamic (random walk) 정보를 추가 활용함을 확인.
