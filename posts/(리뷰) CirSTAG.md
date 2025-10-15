---
title: (리뷰) CirSTAG
author: 신록예찬
date: 06/29/2025 
format: html
---

`-` CirSTAG은 **최신 집적 회로(ICs)의 안정성(Stability)을 분석하기 위한 혁신적인 프레임워크**를 제시한다. 회로 설계에서 안정성 분석은 매우 중요하지만, 기존 방법들은 수많은 시뮬레이션을 필요로 하여 시간과 비용이 많이 소모되는 단점이 있다. 

`-` CirSTAG은 이러한 문제점을 해결하기 위해 **사전 학습된 그래프 신경망(GNN) 모델**을 활용한다. GNN은 그래프 형태의 데이터를 효율적으로 처리하는 딥러닝 모델로, 회로를 그래프로 표현하여 GNN이 회로 동작을 모방하도록 학습될 수 있다.

`-` CirSTAG의 핵심 아이디어는 **GNN의 입출력 그래프 기반 매니폴드(Manifold) 사이의 매핑 왜곡(Mapping Distortion)을 정량화**하는 것이다. 

`-` 입력 매니폴드에서 서로 가까웠던 두 노드(데이터 샘플)가 GNN을 거쳐 출력 매니폴드에서 갑자기 멀리 떨어져 버린다면, 이는 GNN이 해당 부분에서 불안정하다는 것을 나타내며, 이는 곧 회로 자체도 해당 부분에서 민감하다는 것을 의미한다.

## CirSTAG의 세 가지 핵심 단계

`-` CirSTAG은 이 분석 과정을 세 가지 주요 단계로 나누어 수행한다:

### **1단계: 입력/출력 임베딩 행렬 구성 (Phase 1: Input/Output Embedding Matrices Construction)**
   
`-` **입력 임베딩 행렬 ($U_M$)**: 원본 회로 그래프(고차원)의 구조적 특성을 유지하면서 차원을 축소하기 위해 **스펙트럼 그래프 임베딩** 방식을 사용한다. 이는 Laplacian Eigenmap 알고리즘에서 영감을 받았다. 주어진 그래프의 인접 행렬 $A$와 차수 행렬 $D$를 사용하여 정규화된 그래프 라플라시안 행렬 $L_{norm}$을 계산한다:
        $$L_{norm} = I - D^{-\frac{1}{2}}AD^{-\frac{1}{2}}$$
        여기서 $I$는 항등 행렬이다. $L_{norm}$의 상위 $M$개의 가장 작은 고유값 $\tilde{\lambda}_1, \dots, \tilde{\lambda}_M$과 해당 고유 벡터 $\tilde{u}_1, \dots, \tilde{u}_M$을 사용하여 입력 스펙트럼 임베딩 행렬 $U_M$을 구성한다:
        $$U_M = \left[ \sqrt{|1-\tilde{\lambda}_1|}\tilde{u}_1, \dots, \sqrt{|1-\tilde{\lambda}_M|}\tilde{u}_M \right]$$
        이 $U_M$은 입력 그래프의 필수 구조적 특성을 보존하면서 차원을 효과적으로 축소한다. 이 임베딩 행렬 $U_M$은 다음 단계에서 입력 그래프 기반 매니폴드를 구성하는 데 사용된다.

`-` **출력 임베딩 행렬 ($Y$)**: 사전 학습된 GNN 모델의 출력은 일반적으로 저차원 노드 임베딩이므로, 이 **GNN이 생성한 노드 임베딩을 직접 출력 임베딩 행렬로 사용**한다. 이 출력 임베딩 행렬도 다음 단계에서 그래프 토폴로지 학습을 통해 출력 그래프 기반 매니폴드로 변환된다.

### **2단계: PGM(확률 그래프 모델)을 통한 그래프 기반 매니폴드 구성 (Phase 2: Constructing Graph-Based Manifolds via PGMs)**


이 단계에서는 1단계에서 얻은 입력($U_M$) 및 출력($Y$) 임베딩 행렬을 사용하여 **저차원 그래프 기반 매니폴드를 효율적으로 구성**한다. 이 과정은 PGM 기반의 그래프 토폴로지 학습(SGL)에서 영감을 받았지만, 확장성을 높이기 위해 효율적인 스펙트럼 희소화(spectral sparsification) 기법을 도입한다.

매니폴드 구성을 위한 **목적 함수($F(\Theta)$)**는 다음의 볼록 최적화 문제로 정식화된다:
        $$\max_{\Theta} F(\Theta) = \log \det(\Theta) - \frac{1}{M} \text{Tr}(X^\top \Theta X)$$
        여기서 $X \in \mathbb{R}^{N \times M}$는 입력 데이터 행렬(1단계에서 구성된 $U_M$ 또는 GNN 출력 $Y$에 해당)이며, $\Theta$는 정밀도 행렬($\Sigma^{-1}$)이다. $\Theta$는 유효한 그래프 라플라시안 행렬 $L$과 $\Theta = L + \frac{1}{\sigma^2}I$ 관계를 가진다.

이 목적 함수는 두 가지 주요 목표를 달성하도록 설계되었다:

- **$F_1 = \log \det(\Theta)$**: **그래프 위에서 신호의 부드러움(Signal Smoothness over the Graph)**을 장려한다.
- **$-\frac{1}{M} F_2 = -\frac{1}{M} \text{Tr}(X^\top \Theta X)$**: **추정된 그래프 토폴로지의 희소성(Sparsity of the Estimated Graph Topology)**을 장려한다. 이는 데이터를 나타내는 노드 간의 불필요한 연결을 줄여 매니폴드를 단순화한다.

이 함수를 최대화하기 위해, 각 엣지 $p,q$에 대한 편미분은 유효 저항(Effective Resistance) $R_{eff}^{p,q}$와 데이터 벡터 간의 ℓ2 거리 $D_{data}^{p,q}$로 표현된다:
        $$\frac{\partial F_1}{\partial w_{p,q}} = R_{eff}^{p,q}$$
        $$\frac{\partial F_2}{\partial w_{p,q}} = D_{data}^{p,q} = \|X^\top e_{p,q}\|_2^2$$
        여기서 $e_{p,q} = e_p - e_q$이며, $e_p$는 $p$번째 항목에 1을 가진 표준 기저 벡터이다. 이 그래디언트를 기반으로 **스펙트럼 왜곡 메트릭 $\eta_{p,q}$**를 사용하여 중요하지 않은 엣지를 가지치기(prune)한다:
        $$\eta_{p,q} = \frac{R_{eff}^{p,q}}{D_{data}^{p,q}} = w_{p,q}R_{eff}^{p,q}$$
        $\eta_{p,q}$ 값이 작은 엣지(즉, 낮은 유효 저항과 큰 데이터 거리)는 가지치기되어 희소한 그래프 구조를 만든다. 초기 조밀 그래프는 k-최근접 이웃(kNN) 알고리즘을 사용하여 효율적으로 구성된다. PGM은 스펙트럼 그래프 희소화 기법인 짧은 사이클 분해(short-cycle decomposition)를 통해 정제된다.

### **3단계: 매니폴드에서의 안정성 분석 (Phase 3: Stability Analysis on the Manifolds)**

`-` GNN 모델의 안정성을 정량화하기 위해 **거리 매핑 왜곡(Distance Mapping Distortion, DMD)** 메트릭을 사용한다. DMD는 입력 매니폴드 $G_X$의 두 데이터 샘플 $p, q$ 사이의 거리 $d_X(p,q)$에 대한 출력 매니폴드 $G_Y$의 거리 $d_Y(p,q)$의 비율로 정의된다:
        $$\delta_F(p,q) \stackrel{\text{def}}{=} \frac{d_Y(p,q)}{d_X(p,q)}$$
        $d_X$와 $d_Y$는 지오데식 거리 대신 효율적으로 계산할 수 있고 전역 그래프 구조를 포착하는 **유효 저항 거리(effective resistance distance)**를 채택한다.


`-` CirSTAG은 입력 $X$가 GNN 모델 $F$에 의해 출력 $Y$로 매핑될 때, 입력 및 출력 매니폴드의 라플라시안 행렬 $L_X$ 및 $L_Y$를 활용하여 노드 안정성을 평가한다. 가장 큰 $s$개의 일반화된 고유값 $\zeta_1 \ge \zeta_2 \ge \dots \ge \zeta_s$와 해당 고유 벡터 $v_1, v_2, \dots, v_s$를 포함하는 가중 고유 부분 공간 행렬 $V_s$를 계산한다:
        $$V_s = \left[ \frac{v_1}{\sqrt{\zeta_1}}, \dots, \frac{v_s}{\sqrt{\zeta_s}} \right]$$
        이 $V_s$를 사용하여 입력 그래프 $G_X$를 $s$차원 벡터로 임베딩한다.
- 엣지 $(p,q) \in E_X$의 안정성은 엣지 양 끝 노드의 임베딩 거리로 정량화된다:
        $$\|V_s^\top e_{p,q}\|_2^2$$
        여기서 $e_{p,q} = e_p - e_q$이며, $e_p$는 표준 기저 벡터이다.
- 최종적으로, 각 노드 $p$의 **안정성 점수($\text{Score}_F(p)$)**는 $p$의 이웃 노드 $q_i$에 대한 DMD 값의 평균을 기반으로 추정된다:
        $$\text{Score}_F(p) \stackrel{\text{def}}{=} \frac{1}{|N_X(p)|} \sum_{q_i \in N_X(p)} \left( \|V_s^\top e_{p,q_i}\|_2^2 \right) \propto \frac{1}{|N_X(p)|} \sum_{q_i \in N_X(p)} \left( \delta_F(p, q_i) \right)^3$$
        여기서 $N_X(p)$는 $G_X$에서 노드 $p$의 이웃 집합을 나타낸다. 이 노드 안정성 점수는 매니폴드 설정에서 GNN 모델의 **'로컬 립시츠 상수(Local Lipschitz constant)'와 유사**한 역할을 한다. 이 점수를 통해 **회로 성능에 큰 영향을 미칠 수 있는 가장 중요하거나 민감한(sensitive) 회로 요소들을 식별**할 수 있게 된다.

## CirSTAG의 장점 및 적용 분야

`-` CirSTAG은 **거의 선형적인 실행 시간 복잡도(near-linear runtime complexity)**를 가지며, 그 데이터 중심적인 특성 덕분에 **다양한 GNN 아키텍처와 호환**된다. 실제 회로 설계에 대한 경험적 평가를 통해, CirSTAG이 다양한 파라미터 변화(노드 특징 및 그래프 토폴로지 교란) 하에서 각 회로 요소의 안정성을 정확하게 추정할 수 있음을 입증한다. 이는 대규모 집적 회로 설계의 안정성을 평가하기 위한 확장 가능한 방법론을 제공한다.

`-` CirSTAG의 응용 분야로는 **회로 타이밍 예측** 및 **회로 기능 역공학(reverse engineering)**과 같은 다양한 EDA(Electronic Design Automation) 문제가 있다. 실험 결과, 차원 축소 단계는 불안정성 순위를 정확하게 예측하는 데 중요함을 보여준다. 또한, 불안정한 노드에 대한 위상 교란(topology perturbation)은 안정적인 노드에 비해 GNN 성능 저하(F1 매크로 점수)가 더 급격하게 나타남을 관찰한다.

