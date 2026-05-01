---
title: 리뷰 ▷ TDL Deep Dive 2편 등장인물
author: 클로드
date: 04/30/2026
draft: false
output-file: 260430_327a7a.html
---

# 🎙️ TDL Deep Dive — 시리즈 2편 "이 동네 사람들 — 누가 누군지 정리해드림"

> ⚠️ **면책**: 본 글은 **엔터테인먼트성 내러티브 분석**이다. "두목", "동맹", "본거지" 같은 표현은 모두 **수집된 364편 데이터셋 안에서 본 패턴**에 대한 편집적 라벨이지, 학계의 공식 분류가 아니다. 단, **개별 논문의 저자·소속·arXiv ID·venue·요지**는 PDF 본문에서 직접 확인한 사실. 인용 번호 [n]은 글 끝 부록 표와 일치.

> 364편의 논문을 깠더니 같은 이름이 자꾸 나온다.
> 이 분야의 등장인물 정리편.
> 진행: 이 선배 (이론파), 박 후배 (응용파)

**📑 시리즈 인덱스**
- [시리즈 1편 — 5막 드라마](260430_리뷰_TDLDeepDive_1편5막드라마.md)
- **시리즈 2편 — 이 동네 사람들** ← 여기
- [시리즈 3편 — 수학 무기고](260430_리뷰_TDLDeepDive_3편수학무기고.md)

---

**박 후배**: 지난 편 보고 댓글이 하나 왔어요. "진영은 알겠는데, 그래서 누가 그 진영의 두목인지 모르겠다."

**이 선배**: 좋은 지적이에요. 오늘은 **이 분야의 캐릭터 시트**를 만들어드릴게요. 등장인물 9팀.

---

## 👑 1. Bianconi 학파 — "Dirac operator를 들고 온 사람들"

**이 선배**: 본거지는 **Queen Mary University of London**. 두목은 **Ginestra Bianconi**.

**박 후배**: 이 사람 원래 통계물리 출신 아니에요?

**이 선배**: 맞아요. 복잡계 네트워크의 거장. 그런데 최근 **Dirac operator를 graph/hypergraph에 갖다 붙이는 작업**으로 폭발 중. 본 데이터셋 대표작이 directed hypergraph + Dirac operator로 ABIDE 자폐 데이터를 풀어낸 논문 [1]. 톤이 거의 **Nature Communications**급.

같은 학파가 Sapienza의 Barbarossa와 합작으로 Dirac spectral transform을 정식화한 논문도 있어요 [2]. **영국-이탈리아 라인**이 지금 TSP 최강 동맹.

---

## 🇮🇹 2. Sapienza 학파 — "Hodge decomposition으로 뇌를 분해하는 사람들"

**이 선배**: **Sapienza University of Rome**. 두목은 **Sergio Barbarossa**, 부두목은 **Stefania Sardellitti**.

**박 후배**: Barbarossa는 신호처리 학계의 거장이죠.

**이 선배**: 그렇죠. 이 사람들이 **TSP(Topological Signal Processing)** framework를 만들었어요. 단순 그래프 위 신호처리(GSP)를 simplicial complex / cell complex로 확장한 것.

대표작은 **HCP 데이터셋 100명의 dMRI+fMRI에 Hodge decomposition**을 돌린 논문 [3]. 뇌 신호를 **발산-회전-조화** 세 부분으로 분해. 같은 학파가 **3D Point Cloud**에도 같은 도구를 적용 [4]. **Discrete Exterior Calculus** 까지 끌어옴.

**박 후배**: 이 사람들 진짜 multitasking이네요.

**이 선배**: 동시에 4-5개 라인을 굴립니다. 거의 **모든 핵심 논문에 IEEE 저널 톤**.

---

## 🇬🇧 3. Cambridge 학파 — "Sheaf와 simplicial attention"

**이 선배**: **University of Cambridge**. 두목은 **Pietro Liò**.

**박 후배**: Liò는 진짜 어디에나 있는 사람이에요.

**이 선배**: 양적으로 압도적. 우리 데이터셋만 봐도 sheaf, simplicial, hypergraph 거의 모든 클러스터에 끼어 있어요.

대표작 셋:
- Sapienza와 합작한 **Polynomial Neural Sheaf Diffusion** [5] — Sheaf Laplacian의 K차 다항식.
- **N-simplicial Attention** [6] — Attention의 simplicial 일반화 + RoPE.
- 항저우-홍콩 그룹과 합작한 AAAI-26 Sheaflet HNN [7], 그리고 같은 컨퍼런스에 또 한 편 PEF-HNN [8].

**박 후배**: 마지막 두 편은 홍콩-항저우-싱가포르-케임브리지 합작이네요.

**이 선배**: 네. **이게 최근 1년의 가장 강력한 hypergraph 동맹**입니다.

---

## 🦁 4. Oxford 학파 — "Bronstein의 GFM 진영"

**이 선배**: **University of Oxford**. 두목은 **Michael Bronstein**, 보좌는 **Xiaowen Dong**.

**박 후배**: Bronstein은 Geometric Deep Learning 책의 그 사람이죠.

**이 선배**: 맞아요. 그런데 최근 행보가 흥미로워요. **Graph Foundation Model이 안 된다는 걸 본인이 직접 증명함** [9]. 결론: 고정 backbone GFM은 task-dependent architectural attribute에 robust하지 않다.

**박 후배**: 다음 단계는 architecture-adaptive GFM이라는 거죠.

**이 선배**: 이게 시드 토픽이에요. 6개월 안에 후속 논문 떼거지로 나올 겁니다.

---

## ⚖️ 5. EPFL / Aachen 학파 — "스펙트럴-메시지패싱 종전 협정파"

**이 선배**: 두 진영을 통합한 **종전 협정문**을 쓴 사람들 [10]. **Frossard, Vandergheynst, Schaub, Morris, Wolf, Levie** — 양 학파 거장 다 모았어요.

**박 후배**: 진짜 어벤져스네요.

**이 선배**: 이 논문 한 편으로 **GSP 학파와 MPNN 학파가 같은 framework로 묶입니다.** 박사 1년차한테 던져주기 좋은 핵심 레퍼런스.

---

## 🇰🇷 6. 한국 분파 — "Yonsei 라인"

**이 선배**: 한국에서도 한 팀이 들어와 있어요.

대표작 [11]: **"Are Graph Transformers Necessary?"** — Graph Transformer 없이 MPNN + fractal node로 long-range 처리.

**박 후배**: 이거 도발적인 제목이네요.

**이 선배**: 의도된 도발이죠. **"Graph Transformer는 필요 없다"** 가 본문 주장이에요. 또 다른 한국 그룹은 고전 TDA 쪽 [12] — RCLA 알고리즘으로 **Poisson noise 모델 하 stability theorem** 증명.

---

## 🧠 7. Padova 학파 — "이론신경과학의 부활"

**이 선배**: **University of Padova (이탈리아)**. 두목은 **Samir Suweis**.

대표작 [13]: heterogeneous DMFT(동적 평균장 이론)로 **MICrONS connectome(쥐 시각피질 1mm³ 전체 시냅스)** 검증.

**박 후배**: q-bio 카테고리네요. 이거 PRL/PNAS 노리는 톤?

**이 선배**: PRX급. 이론신경과학에서 굉장히 무거운 논문이에요.

---

## 🤖 8. Google DeepMind — "트랜스포머 위상 비판자"

**이 선배**: 외부 도전자.

대표작 [14]: **Mike Mozer 팀**의 "The Topological Trouble With Transformers". 요지는 "Feedforward transformer는 state tracking이 위상학적으로 불가능, RNN/SSM 회귀해야 한다."

**박 후배**: Mozer는 옛날 RNN 시절 거장이죠.

**이 선배**: **연어 회귀**. 옛날 RNN 학파가 위상학으로 무장하고 돌아온 셈이에요.

---

## 🔬 9. Harvard Med / GRaM 학파 — "GPT-2의 곡률을 잰 사람"

**이 선배**: 단신 외인구단.

대표작 [15]: **Asif Khan**의 GPT-2 / GPT-2-medium에 **Ollivier-Ricci curvature** 측정. 깊은 레이어로 갈수록 contractive support, H1 lifetime 감소, H0 skeleton 유지.

**박 후배**: 이거 mechanistic interpretability하고 합쳐지면 무서워지겠네요.

**이 선배**: 이미 합쳐지는 중입니다.

---

## 🗺️ 캐릭터 시트 한 장 요약

| 학파 | 본거지 | 두목 | 무기 | 주요 적/동맹 | 대표 인용 |
|---|---|---|---|---|---|
| Bianconi | QMUL (런던) | Bianconi | Dirac operator, directed hypergraph | 동맹: Sapienza | [1, 2] |
| Sapienza | 로마 | Barbarossa, Sardellitti | Hodge, DEC, TSP | 동맹: Bianconi | [3, 4] |
| Cambridge | 케임브리지 | Pietro Liò | Sheaf, simplicial attention | 동맹: 항저우/홍콩/싱가포르 | [5, 6, 7, 8] |
| Oxford | 옥스퍼드 | Bronstein | GFM 비판 | 적: 자기 자신(GFM 옛날 주장) | [9] |
| EPFL/Aachen | 로잔/아헨 | Frossard, Vandergheynst, Schaub, Morris | 통합 framework | 중재자 | [10] |
| Yonsei (KR) | 서울 | Jeongwhan Choi, Noseong Park | MPNN + fractal nodes | 적: Graph Transformer | [11] |
| 경북대 추정 (KR) | — | Choi, Yang 등 | Persistent Homology stability | 동맹: TDA 클래식 | [12] |
| Padova | 파도바 | Suweis | DMFT + connectome | 동맹: 이론신경과학 | [13] |
| DeepMind | Mountain View | Mike Mozer | 위상학적 비판 | 적: Transformer 제국 | [14] |
| Harvard Med | 보스턴 | Asif Khan | Ollivier-Ricci curvature | 동맹: mechanistic interp | [15] |

> "두목"·"동맹"·"적" 등은 모두 편집적 라벨. 데이터셋 안 패턴 기반.

---

**박 후배**: 이거 보니까 **이탈리아 라인이 진짜 무섭네요.** Bianconi-Sapienza-Padova가 한 줄로 서면.

**이 선배**: 정확해요. **이탈리아가 TDL의 실리콘밸리**예요. 다음 편(시리즈 3편)에서 그들이 휘두르는 **수학적 무기들**을 정리해드릴게요.

---

## 📚 부록 — 인용된 논문 15편

| # | 논문 | 저자 (소속) | venue / arXiv | 날짜 |
|---|---|---|---|---|
| 1 | Topological cluster sync via Dirac spectral programming on directed hypergraphs | Yupeng Guo, Ahmed A. A. Zaid, Xueming Liu (HUST), Ginestra Bianconi (QMUL) | arXiv:2512.14729 [physics.soc-ph] | 2025-12-10 |
| 2 | Learning Dirac Spectral Transforms for Topological Signals | Leonardo Di Nino, Tiziana Cattai, Sergio Barbarossa (Sapienza), Ginestra Bianconi (QMUL), Paolo Di Lorenzo | arXiv:2602.14590 [eess.SP] | 2026-02-16 |
| 3 | Multimodal Higher-Order Brain Networks: A TSP Perspective | Breno C. Bispo, Stefania Sardellitti, Juliano B. Lima, Fernando A. N. Santos | arXiv:2603.29903 [q-bio.NC], IEEE 톤 | 2026-03-31 |
| 4 | Topological Signal Processing for 3D Point Cloud Data | Tiziana Cattai, Stefania Sardellitti (Mercatorum), Sergio Barbarossa (Sapienza), Stefania Colonnese | arXiv:2602.19636 [eess.SP] | 2026-02-23 |
| 5 | Polynomial Neural Sheaf Diffusion | Alessio Borgi (Sapienza+Cambridge), Fabrizio Silvestri, Pietro Liò | arXiv:2512.00242 [cs.LG] | 2025-11-28 |
| 6 | How Smoothing is N-simplicial Attention? | Alexandre Dussolle (Cambridge / ENPC), Pietro Liò (Cambridge) | arXiv:2512.15600 [cs.LG] | 2025-12-17 |
| 7 | High-Pass Matters: Sheaflet-Based Design for HNN | Ming Li (Zhejiang Normal), Yujie Fang, Dongrui Shen, Han Feng, Xiaosheng Zhuang (CityU HK), Kelin Xia (NTU), Pietro Liò (Cambridge) | **AAAI-26** | 2026 |
| 8 | Permutation Equivariant Framelet-based HNN (PEF-HNN) | Ming Li, Yi Wang, Chengling Gao, Lu Bai (Beijing Normal), Yujie Fang, Xiaosheng Zhuang, Pietro Liò | **AAAI-26** | 2026 |
| 9 | Can Graph Foundation Models Generalize Over Architecture? | Benjamin Gutteridge, Michael Bronstein, Xiaowen Dong (Oxford + AITHYRA) | **GRaM @ ICLR 2026** Proceedings | 2026 |
| 10 | Position: Message-passing and Spectral GNNs are Two Sides of the Same Coin | A. Vasileiou (RWTH), J. Cervino (Penn), P. Frossard (EPFL), C. Kanatsoulis (Penn), C. Morris (RWTH), M. Schaub (RWTH), P. Vandergheynst (EPFL), Z. Wang (Florida), G. Wolf (Montreal), R. Levie (Technion) | arXiv:2602.10031 [cs.LG], Position paper | 2026-02-10 |
| 11 | Are Graph Transformers Necessary? Efficient Long-Range Message Passing with Fractal Nodes in MPNNs | Jeongwhan Choi, Seungjun Park, Sumin Park, Sung-Bae Cho, Noseong Park (Yonsei + KAIST) | preprint | 2026 |
| 12 | Denoising Data Reduction for TDA (RCLA) | Seonmi Choi, Semin Oh, Jeong Rye Park, Seung Yeop Yang | arXiv:2603.29248 [cs.CG] | 2026-03-31 |
| 13 | Topological Origin of Diversity of Timescales in RNN Circuits | Marco Zenari, Luca Taffarello, Luca Mazzucato (Oregon), Amos Maritan, Samir Suweis (Padova) | arXiv:2603.04149 [q-bio.NC] | 2026-03-04 |
| 14 | The Topological Trouble With Transformers | Michael C. Mozer, Shoaib A. Siddiqui, Rosanne Liu (Google DeepMind) | arXiv:2604.17121 [cs.LG] | 2026-04-28 |
| 15 | The Geometrical and Topological Signatures of Transformers | Asif Khan (Harvard Med) | **GRaM @ ICLR 2026** Tiny Paper | 2026 |

---

> 📁 출처: `data/TDL/downloaded_papers/` 364편
> 📅 분석일: 2026-04-30
