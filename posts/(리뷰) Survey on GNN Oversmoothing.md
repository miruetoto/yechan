---
title: (리뷰) Survey on GNN Oversmoothing
author: 신록예찬
date: 07/22/2025
draft: false
---

<https://arxiv.org/abs/2303.10993>

```         
@article{rusch2023survey,
  title={A survey on oversmoothing in graph neural networks},
  author={Rusch, T Konstantin and Bronstein, Michael M and Mishra, Siddhartha},
  journal={arXiv preprint arXiv:2303.10993},
  year={2023}
}
```

# 빠른요약

-   **서론: GNN의 깊이 한계와 오버-스무딩 현상**

    -   그래프 신경망(GNN)은 관계형 및 상호작용 데이터 학습에 강력한 도구로 부상함.
    -   대부분의 GNN은 컴퓨터 비전 분야의 CNN과 달리 얕은 층을 가지며, 깊은 GNN의 성능을 저해하는 문제점 중 하나가 **오버-스무딩**임.
    -   오버-스무딩은 GNN의 깊이가 증가함에 따라 **노드 특성(node features)이 유사해지는 현상**을 의미하며, 모든 노드 특성이 동일한 상수 값으로 기하급수적으로 수렴하는 것을 말함.
    -   과도한 스무딩은 정보가 없는 한계로 수렴하게 만들며, 특히 노드 레이블이 이웃 노드 레이블과 다른 이종 그래프(heterophilic graphs) 처리 능력을 심각하게 저해할 수 있음.

-   **오버-스무딩의 정의 (정의 1)**

    -   오버-스무딩은 적절한 노드-유사도 측정($\mu$)이 **층 수(**$n$)에 따라 기하급수적으로 $0$에 수렴하는 것으로 공리적으로 정의됨. 즉, $\mu(X_n)$이 $C_1e^{-C_2n}$ 형태로 수렴함.
    -   노드-유사도 측정 $\mu$는 다음 두 가지 조건을 만족해야 함:
        1.  모든 노드 $i$에 대해 $X_i$가 상수 $c$인 경우에만 $\mu(X) = 0$이 성립해야 함 (노드 특성이 일정한 값으로 수렴할 때 오버-스무딩이 발생한다는 개념을 공식화함).
        2.  삼각 부등식 또는 부분 가법성($\mu(X + Y) \le \mu(X) + \mu(Y)$)을 만족해야 함.
    -   이 정의는 기존의 다양한 오버-스무딩 접근 방식들을 통합하며, **기하급수적 수렴**이 오버-스무딩의 중요한 특징임을 강조함.

-   **오버-스무딩 측정 방법**

    -   **디리클레 에너지(Dirichlet energy)**: $$E(X_n) = \frac{1}{v} \sum_{i \in V} \sum_{j \in N_i} |X_n^i - X_n^j|^2_2$$
        -   $\mu(X_n) = \sqrt{E(X_n)}$은 정의 1의 모든 조건을 만족하는 유효한 노드-유사도 측정임.
        -   $L_p$-norm ($p > 1$) 기반의 디리클레 에너지도 가능함.
    -   **평균 절대 거리(Mean Average Distance, MAD)**: $$\mu(X_n) = \frac{1}{v} \sum_{i \in V} \sum_{j \in N_i} \left(1- \frac{X_n^i \cdot X_n^j}{|X_n^i ||X_n^j |}\right)$$
        -   MAD는 정의 1의 조건 1과 2를 만족하지 않아 유효한 노드-유사도 측정이 아님. 특히 스칼라 값의 경우 항상 $0$이 되어 오버-스무딩의 충분 조건이 되지 못함.
        -   그러나 다차원($m > 1$)의 경우 층 수가 증가함에 따라 기하급수적으로 $0$에 수렴하는 경향을 보여 조건 3을 만족할 수 있음.
        -   **디리클레 에너지가 더 안정적이고 모든 조건을 만족하므로 선호됨**.
    -   **실험적 평가**: GCN, GAT, GraphSAGE와 같은 표준 GNN 아키텍처는 Texas, Cora, Cornell5 등 다양한 실제 그래프 데이터셋에서 모두 **층별 측정값(디리클레 에너지 및 MAD)이 층 수가 증가함에 따라 기하급수적으로** $0$에 수렴함을 보임.

-   **오버-스무딩 완화 방법**

    -   **1. 정규화 및 규제 (Normalization and Regularization)**:
        -   명시적 규제: Energetic GNNs (EGNNs)는 디리클레 에너지 편차에 페널티를 부과.
        -   암시적 규제: DropEdge, Graph DropConnect (GDC), PairNorm, Differentiable Group Normalization (DGN), NodeNorm 등이 있음.
            -   PairNorm은 다음 수식을 통해 노드 특성 $X$를 정규화하여 쌍별 거리(pairwise distances)를 각 층에서 일정하게 유지함: $$\hat{X}_i = X_i - \frac{1}{v} \sum_{j=1}^v X_j$$ $$X_i = s\frac{\hat{X}_i}{\sqrt{\frac{1}{v} \sum_{j=1}^v |\hat{X}_j|^2_2}}$$ 여기서 $s > 0$는 하이퍼파라미터임.
    -   **2. GNN 동역학 변경 (Change of GNN dynamics)**:
        -   GCN의 확산과 유사한 동역학을 대체함.
        -   Graph-Coupled Oscillator Network (GraphCON), PDE-GCN, Allen-Cahn Message Passing (ACMP), Gradient Flow Framework (GRAFF), Gradient Gating (G2) 등이 있음.
            -   GraphCON은 비선형 진동자(non-linear oscillators)를 사용하여 확산과 유사한 GCN의 동역학을 대체함: $$Y_n = Y_{n-1} + \Delta t[\sigma(F_{\theta_n}(X_{n-1},G)) - \gamma X_{n-1} - \alpha Y_{n-1}]$$ $$X_n = X_{n-1} + \Delta t Y_n$$ 여기서 $Y_n$은 보조 노드 특성($\Delta t > 0$는 시간 단계)임.
            -   Gradient Gating (G2)는 그래프 기울기(graph-gradient)를 활용하는 학습 가능한 노드별 조기 종료 메커니즘을 게이팅 함수를 통해 구현함: $$\hat{\tau}_n = \sigma(\hat{F}_{\hat{\theta}_n}(X_{n-1},G))$$ $$\tau_n^{ik} = \tanh\left(\sum_{j \in N_i} |\hat{\tau}_n^{jk} - \hat{\tau}_n^{ik}|^p\right)$$ $$X_n = (1-\tau_n) \odot X_{n-1} + \tau_n \odot \sigma(F_{\theta_n}(X_{n-1},G))$$ 여기서 $p \ge 0$임.
    -   **3. 잔차 연결 (Residual connections)**:
        -   ResNet에서 영감을 받아 깊은 GNN에 잔차 연결을 추가함.
        -   Res-GCN, GCNII, Jumping Knowledge Networks (JKNets), Deep Adaptive Graph Neural Networks (DAGNNs) 등이 있음.
            -   Res-GCN은 다음과 같은 잔차 연결을 포함함: $$X_n = X_{n-1} + F_{\theta_n}(X_{n-1},G)$$
            -   GCNII는 모든 층에 초기 노드 특성 $X_0$의 스케일된 잔차 연결을 추가함: $$X_n = \sigma \left[ \left( (1-\alpha_n)\hat{D}^{- \frac{1}{2}} \hat{A} \hat{D}^{- \frac{1}{2}} X_{n-1} + \alpha_n X_0 \right) \left( (1-\beta_n)I + \beta_n W^n \right) \right]$$ 여기서 $\alpha_n, \beta_n \in$는 하이퍼파라미터임.
    -   **완화 방법의 실험적 평가**:
        -   DropEdge-GCN과 Res-GCN은 여전히 오버-스무딩을 겪음 (디리클레 에너지가 기하급수적으로 $0$에 수렴).
        -   G2-GCN, GraphCON-GCN, PairNorm, GCNII는 층별 디리클레 에너지를 거의 일정하게 유지하여 오버-스무딩을 완화함.

-   **표현력 희생의 위험: 완화는 필수적이지만 충분하지 않음**

    -   오버-스무딩을 완화하는 것이 깊은 GNN 구축에 **필수적(necessary)**이지만, 이것만으로는 충분하지 않음.
    -   단순히 편향(bias) 벡터를 추가한 GCN이나 PairNorm은 디리클레 에너지를 일정하게 유지하여 오버-스무딩을 완화하지만, 층 수가 증가할수록 성능(정확도)이 크게 감소함.
        -   예를 들어, 층별 디리클레 에너지를 약 $1$로 일정하게 유지하는 GCN에 편향을 추가한 모델과 PairNorm 모델은 층 수가 증가할 때 성능이 급격히 저하됨.
    -   **오버-스무딩 완화만을 목표로 GNN의 표현력(expressive power)을 희생하는 것이 주요 함정임**.
    -   G2와 같이 오버-스무딩을 완화하면서도 층 수가 증가함에 따라 성능이 약간 증가하는 방법도 있음.

-   **연속시간 GNN으로의 확장**

    -   그래프 표현 학습의 새로운 분야로, GNN을 **깊이에 따라 연속적인 모델**로 다룸.
    -   메시지 패싱 전파를 상미분 방정식(ODEs) 또는 편미분 방정식(PDEs)으로 모델링된 그래프 동적 시스템으로 정식화함.
        -   이는 다음과 같은 미분 방정식으로 모델링될 수 있음: $$X'(t) = \sigma(F_{\theta}(X(t),G))$$ 여기서 $X(t)$는 시간 $t \ge 0$에서의 노드 특성을 의미함.
    -   **연속시간 오버-스무딩 (정의 2)**: 노드 유사성 측정 $\mu$가 **시간(**$t$)에 따라 기하급수적으로 $0$에 수렴하는 것으로 정의됨 (즉, $\mu(X(t))$가 $C_1e^{-C_2t}$ 형태로 수렴).
    -   Graph Neural Ordinary Differential Equations (GDEs), Continuous Graph Neural Networks (CGNN), Graph-Coupled Oscillator Networks (GraphCON) 등이 연속시간 GNN의 예시임.

-   **결론**

    -   깊은 GNN은 장거리 상호작용이 있는 관계형 데이터, 특히 이종 그래프 데이터 학습에 필수적임.
    -   **오버-스무딩은 깊은 GNN 구축의 핵심 과제**임.
    -   논문은 오버-스무딩에 대한 **공리적 정의**를 제시하고, 오버-스무딩 완화가 깊은 GNN 성능 향상을 위한 **필수 조건**이지만 GNN의 **표현력 손실 없이** 이루어져야 함을 강조함.
    -   이 정의는 **연속시간 GNN** 분야로 확장될 수 있음.