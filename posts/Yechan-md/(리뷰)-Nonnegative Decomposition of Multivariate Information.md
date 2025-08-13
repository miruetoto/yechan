---
title: (리뷰) Nonnegative Decomposition of Multivariate Information
author: 신록예찬
date: 07/31/2025
draft: false
bibliography: ref.bib
---



### 1. 중복성(Redundancy)의 새로운 정의: $I_{min}$

PID는 **중복성($I_{min}$)에 대한 새로운 정의에서 시작합니다. 이는 주어진 변수($S$)의 각 가능한 결과에 대해 _어떤 소스라도 제공하는 최소한의 정보_를 모든 가능한 결과에 대해 평균한 값입니다.

- **특정 정보(Specific Information):** 먼저, 소스 $A$가 변수 $S$의 특정 결과 $s$에 대해 제공하는 정보인 **특정 정보 $I(S=s;A)$**를 정의합니다. 이는 $S=s$의 놀람(surprise)에 대한 평균 감소를 정량화합니다. $$I(S=s;A) = \sum_a p(a|s) \left[ \log \frac{1}{p(s)} - \log \frac{1}{p(s|a)} \right] \quad \text{[12, Eq. 2]}$$ 여기서 $p(s)$는 $S=s$의 확률, $p(a|s)$는 $S=s$일 때 $A=a$의 조건부 확률입니다. $I(S=s;A)$는 항상 음수가 아닙니다.
    
- **중복성 $I_{min}$의 공식:** $k$개의 소스($A_1, A_2, \ldots, A_k$)에 대한 중복성 $I_{min}$은 $S$의 각 결과에 대해 어떤 소스라도 제공하는 최소 정보의 기댓값으로 정의됩니다. $$\mathbf{I_{min}(S; {A_1,A_2, \ldots,A_k}) = \sum_s p(s) \min_{A_i} I(S=s;A_i)} \quad \text{[13, Eq. 3]}$$ 이 정의는 모든 소스에 공통된 정보(어떤 소스라도 제공하는 최소 정보)를 포착하며, 소스들이 $S$의 다른 결과에 대해 정보를 제공할 수 있다는 점을 고려합니다.
    
- **$I_{min}$의 특성:**
    
    - $I_{min}$은 **음수가 아닙니다**.
    - $I_{min}$은 모든 $A_i$에 대해 $I(S;A_i)$보다 작거나 같으며, 모든 소스가 $S$에 대해 정확히 동일한 정보를 제공하는 경우에만 같아집니다.
    - 주어진 소스 $A$에 대한 중복 정보량은 "자기 중복성($I_{min}(S;{A}) = I(S;A)$)"에 의해 최대화됩니다.
- **중복성 격자(Redundancy Lattice):** $I_{min}$의 도메인은 $R$의 비어있지 않은 부분집합들의 모든 집합 중에서, 어떤 소스도 다른 소스의 상위 집합이 아닌 집합들의 컬렉션인 **$A(R)$**로 줄일 수 있습니다. $$A(R) = {\alpha \in P_1(P_1(R)) : \forall A_i,A_j \in \alpha, A_i \not\subset A_j} \quad \text{[16, Eq. 4]}$$ $A(R)$의 요소에 대한 부분 순서($\triangleleft$)는 중복성 격자를 생성하며, 이는 중복성의 구조에 대한 통찰력을 제공합니다.
    

### 2. 부분 정보 함수(Partial Information Function): $\Pi_R$

$I_{min}$은 중복성을 누적 정보 함수로 포착하는 반면, **부분 정보 함수($\Pi_R$)**는 각 특정 소스 컬렉션이 **고유하게 기여하는 부분 정보**를 측정합니다. 이 부분 정보들이 전체 섀넌 정보를 분해하는 원자들을 형성합니다.

- **정의:** $\Pi_R$은 $I_{min}$의 뫼비우스(Möbius) 역함수로 정의됩니다. $$I_{min}(S;\alpha) = \sum_{\beta \triangleleft \alpha} \Pi_R(S;\beta) \quad \text{[21, Eq. 6]}$$
    
- **재귀적 계산:** $\Pi_R$은 다음과 같이 재귀적으로 계산할 수 있습니다. $$\mathbf{\Pi_R(S;\alpha) = I_{min}(S;\alpha) - \sum_{\beta \prec \alpha} \Pi_R(S;\beta)} \quad \text{[21, Eq. 7]}$$ 여기서 $\beta \prec \alpha$는 중복성 격자에서 $\alpha$보다 낮은 모든 $\beta$를 의미합니다. 즉, $\Pi_R(S;\alpha)$는 $\alpha$의 소스들에 의해 중복적으로 제공되지만, 더 간단한 소스 컬렉션(즉, 중복성 격자에서 $\alpha$보다 낮은 $\beta$)에 의해 제공되지 않는 정보를 정량화합니다.
    
- **닫힌 형태(Closed Form):** $\Pi_R$은 또한 닫힌 형태로 표현될 수 있습니다. $$\mathbf{\Pi_R(S;\alpha) = I_{min}(S;\alpha) - \sum_s p(s) \max_{\beta \in \alpha^-} \min_{B \in \beta} I(S=s;B)} \quad \text{[22, Eq. 8]}$$ 여기서 $\alpha^-$는 중복성 격자에서 $\alpha$ 바로 아래에 있는 노드들을 나타냅니다.
    
- **$\Pi_R$의 특성:**
    
    - $\Pi_R$은 **음수가 아닙니다**. 이는 정보량으로서 명확하게 해석될 수 있음을 의미합니다.

### 3. 상호 정보(Mutual Information)의 분해

부분 정보(PI) 항들은 상호 정보(Mutual Information)를 합산하여 분해할 수 있게 합니다. $$I(S;A) = I_{min}(S;{A}) = \sum_{\beta \triangleleft {A}} \Pi_R(S;\beta) \quad \text{[22, Eq. 9]}$$

- **3변수 시스템($R = {R_1, R_2}$)의 예시:**
    - $S$와 $R_1$ 사이의 상호 정보: $$I(S;R_1) = \Pi_R(S;{1}) + \Pi_R(S;{1}{2}) \quad \text{[23, Eq. 10]}$$
    - $S$와 $R_1, R_2$ 사이의 상호 정보: $$\mathbf{I(S;R_1, R_2) = \Pi_R(S;{1}) + \Pi_R(S;{2}) + \Pi_R(S;{1}{2}) + \Pi_R(S;{12})} \quad \text{[23, Eq. 11]}$$ 이 방정식들의 각 항은 다음과 같이 해석됩니다:
    - $\mathbf{\Pi_R(S;{1}{2}) = I_{min}(S;{1}{2})}:$ $R_1$과 $R_2$ 사이의 **중복 정보(redundancy)**를 나타냅니다. 이는 $R_1$과 $R_2$가 동일하거나 겹치는 정보를 제공하는 경우입니다.
    - $\mathbf{\Pi_R(S;{1}) = I(S;R_1) - I_{min}(S;{1}{2})}:$ $R_1$의 **고유 정보(unique information)**를 나타냅니다. 이는 $R_2$가 제공하지 않는 $R_1$만의 정보입니다. $\Pi_R(S;{2})$도 마찬가지입니다.
    - $\mathbf{\Pi_R(S;{12})}:$ $R_1$과 $R_2$의 조합에 의해 제공되는 **시너지 정보(synergy)**를 나타냅니다. 이는 개별적으로는 얻을 수 없지만 함께 있을 때만 얻을 수 있는 정보입니다.

PI 원자(PI-atom)는 중복성과 시너지의 다양한 조합으로 해석됩니다. 고유 정보는 중복성 또는 시너지의 퇴화된 형태로 간주될 수 있습니다.

이러한 수식적 정의와 해석을 통해 PID는 다변수 시스템 내에서 정보가 어떻게 분포되고 상호 작용하는지에 대한 명확하고 음수가 아닌 설명을 제공합니다.