---
title: (리뷰) 첸이통, IEEE Signal Process. Lett
author: 신록예찬
date: 04/09/2025
draft: false
bibliography: 2025-04-09-ref.bib
---


`-` $N$개의 노드(=에이전트)가 독립적인 손실함수를 가지고 있음. 각각의 손실함수는 모두 L변수 함수임. 따라서 우리가 최적화해야 할 파라메터의 수는 $N \times L$ 임 

- ${\cal N} = \{1,2,\dots,N\}$ 은 노드들의 집합, 혹은 에이전트의 집합을 의미
- 임의의 에이전트 $k \in {\cal N}$ 는 손실함수 $J_k({\boldsymbol w}_k)$ 를 가짐. 여기에서 ${\boldsymbol w}_k$는 길이가 $L$인 벡터임
- ${\boldsymbol w}_k^\star$는 $J_k({\boldsymbol w}_k)$의 unique minimizer 임. (당연히 ${\boldsymbol w}_k^\star$ 역시 lenght L vector) 


`-` 최적화해야할 파라메터를 모은 NxL 매트릭스를 ${\bf W}$라고 정의하자. 

