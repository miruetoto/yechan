# HST 논문 기호 정리 (hst(guebinOHS)_ver2.pdf 기준)

## 그래프 기본

| 기호 | 의미 | 비고 |
|------|------|------|
| $\mathcal{G} = (V, \mathbf{E}, \mathbf{W})$ | 가중 그래프 | |
| $V = \{1, 2, \ldots, n\}$ | 노드(정점) 집합 | $n = \|V\|$ |
| $\mathcal{E}$ | 간선 집합 | $\{(v_i, v_j) : v_i \leftrightarrow v_j\}$ |
| $\mathbf{E}$ | $n \times n$ 인접 행렬 | $E_{ij} = 1$ if $(v_i, v_j) \in \mathcal{E}$ |
| $\mathbf{W}$ | $n \times n$ 가중치 행렬 | $0 \leq W_{ij} \leq 1$ |
| $v_i$ | $i$번째 노드 | |
| $N_i$ | $v_i$의 이웃 집합 | $\{j : (v_i, v_j) \in \mathcal{E}\}$ |
| $\deg(v_i)$ | $v_i$의 차수 | $\|N_i\|$ |
| $n$ | 노드 수 | $\|V\|$ |

## 초기 데이터

| 기호 | 의미 | 비고 |
|------|------|------|
| $y(v)$ | 노드 $v$에서 관측된 신호값 | |
| $\mathbf{y} = (y(v_1), \ldots, y(v_n))^\top$ | 신호 벡터 | |

## Heavy-Snow Process (Algorithm 1)

| 기호 | 의미 | 비고 |
|------|------|------|
| $h(v, t)$ | 시각 $t$에서 노드 $v$의 높이 | $h(v, 0) := y(v)$ |
| $b$ | 업데이트 양 (눈의 양) | $b > 0$ |
| $T_{max}$ | 최대 연속 flow 횟수 | |
| $\boldsymbol{\mu}_0$ | 초기 낙하 분포 | $\mu_0(v_i) \propto \deg(v_i)$ |
| $X_t$ | 시각 $t$에서 눈이 존재하는 노드 인덱스 | |
| $Z_t$ | 연속 flow 카운터 | $Z_t \in \{0, 1, \ldots, T_{max}\}$ |
| $\mathcal{D}_t$ | 하류 영역 (downstream area) | $y(v_i) \geq y(v_j)$이면 flow 발생 |

## 전이 확률

| 기호 | 의미 | 비고 |
|------|------|------|
| $p_{ij}$ | $\mu_0$에 의해 $v_j$가 선택될 확률 | |
| $p_{v_i \to v_j}$ | $v_i$에서 $v_j$로 flow할 확률 | |
| $p_{v_i \to \cdot}$ | block 사건의 확률 | |

## Snow Distance

| 기호 | 의미 | 비고 |
|------|------|------|
| $SD_{ij}(t)$ | 노드 $v_i$, $v_j$ 사이의 snow distance | $\sqrt{\sum_{s=0}^{t}(h(v_i,s)-h(v_j,s))^2}$ |
| $\mathbf{SD}(t)$ | Snow distance 행렬 | |
| $\mathbf{W}(t)$ | Snow weight 행렬 | |

## 기타 파라미터

| 기호 | 의미 | 비고 |
|------|------|------|
| $\theta$ | 연결 강도 조절 파라미터 | Gaussian kernel에서 사용 |
| $\sigma^2$ | Gaussian kernel 분산 | 예: $\sigma^2 = 0.35$ |
| $\varepsilon_i$ | 노이즈 | $\varepsilon_i \sim N(0, 0.1^2)$ |
| $\|\cdot\|_2$ | 유클리드 노름 | |

## 주의: 증명 문서와의 기호 충돌

논문에서 이미 사용 중인 기호이므로, 증명 문서에서 **다른 의미로 사용하면 안 되는** 기호:

| 논문 기호 | 논문에서의 의미 | 증명에서의 충돌 위험 |
|----------|---------------|------------------|
| $V$ | 노드 집합 | Lyapunov 함수로 사용 중 (**충돌!**) |
| $\mathbf{E}$ | 인접 행렬 | 상태 공간 $E$와 혼동 가능 |
| $\mathbf{W}$ | 가중치 행렬 | |
| $\mathcal{D}_t$ | 하류 영역 | |
| $X_t$ | 눈 위치 | 증명에서도 동일 의미로 사용 (OK) |
| $Z_t$ | flow 카운터 | 증명에서도 동일 의미로 사용 (OK) |
| $h(v,t)$ | 높이 | 증명에서도 동일 의미로 사용 (OK) |
| $b$ | 업데이트 양 | 증명에서도 동일 의미로 사용 (OK) |
| $n$ | 노드 수 | 증명에서도 동일 의미로 사용 (OK) |
| $\boldsymbol{\mu}_0$ | 초기 분포 | 증명에서도 동일 의미로 사용 (OK) |
| $\varepsilon_i$ | 노이즈 | Foster-Lyapunov의 $\epsilon$과 혼동 가능 |
