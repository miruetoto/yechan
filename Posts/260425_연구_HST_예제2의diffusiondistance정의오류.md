---
title: 연구 ▷ HST ▷ 예제2의 diffusion distance 정의 오류 — 왜 그림은 괜찮았는가
author: 신록예찬
date: 04/25/2026
draft: false
---

# 0. 한 줄 요약

> 지금까지 예제2(Gaussian ring)에서 "diffusion distance" 라고 그린 것은 **표준 Coifman-Lafon 정의가 아니다**. 코드가 `W` 를 **스칼라 한 개**로 나누는 잘못된 정규화를 쓰고, π-가중도 빠져 있다. **그럼에도 예제2의 그림은 쓸 만하다** — 링 위 Gaussian kernel은 **정규(regular) 그래프** 라서 스칼라 정규화가 우연히 올바른 row-stochastic 정규화와 같아지고, π-가중 누락은 **일정한 스케일 차이**만 만들어 PCA 임베딩 모양은 동일하다.

단, 비정규 그래프 (MCU 네트워크, 실데이터, ER 랜덤, Bi-cluster 등) 에서는 재계산이 필요하다.

---

# 1. 버그 진단

`pyhst/distances.py` 의 기존 구현:

```python
def diffusion_distance(W, t=4):
    col_sum = W.round(3).sum(axis=0)
    denom = col_sum.flat[0]               # ← 0번째 column sum 하나만!
    P = W / denom                         # ← W 전체를 스칼라 하나로 나눔
    Pt = np.linalg.matrix_power(P, t)
    return l2distance(Pt)
```

세 가지 문제:

1. **$P$ 가 row-stochastic 이 아님.** 표준은 $P = D^{-1} W$ 에서 $D = \operatorname{diag}(\text{row}_i(W))$. 즉 각 노드의 out-degree 로 **행마다 따로** 나눠야 한다. 현재 구현은 **column sum 의 첫 성분** 하나로 전체를 나눠 전이행렬이 안 됨.

2. **$\pi$-가중 누락.** Coifman-Lafon 정의:
   $$d_t^2(i,j) \;=\; \sum_k \frac{\bigl(P^t_{ik} - P^t_{jk}\bigr)^2}{\pi_k}$$
   여기서 $\pi$ 는 $P$ 의 stationary distribution. 코드는 분모 없이 단순 L2:
   $$\tilde d_t^2(i,j) \;=\; \sum_k \bigl(P^t_{ik} - P^t_{jk}\bigr)^2$$

3. `W.round(3)` — 수치 아티팩트. 의미 불명.

---

# 2. 예제2가 괜찮았던 이유 — 정규성(regularity)

예제2의 W 는 링 위 Gaussian kernel (각 노드 좌표 거리로 계산):

$$W_{ij} \;=\; \exp\!\left(-\frac{\|\mathbf{x}_i - \mathbf{x}_j\|^2}{2\theta^2}\right), \quad \mathbf{x}_i = (\cos\theta_i, \sin\theta_i)$$

링의 **대칭성** 덕에 모든 노드 $i$ 에서 $\sum_j W_{ij}$ 가 동일 (=정규 그래프). 수치 확인:

```
Gaussian ring (n=60):   row_sum min=11.453  max=11.453  std=0.000
```

이 경우 `W / col_sum[0]` 은 **우연히** `W / row_sum[i]` 와 정확히 같다 (모든 row_sum 이 11.453으로 동일하므로).

따라서:

- **스칼라 정규화 문제 → 발생 안 함.** OLD $P$ = 표준 $P$.
- **$\pi$-가중 누락 → 스케일 차이만 발생.** 정규 그래프에서는 $\pi \equiv 1/n$ 이므로
  $$d_t^2(i,j)_{\text{표준}} \;=\; n \cdot \tilde d_t^2(i,j)$$
  즉 표준 정의 $= n \times$ 기존 계산값.

수치 확인:

```
Gaussian ring, t=4:
  OLD       : min=0.0000  max=0.0609      ← 지금까지 써온 값
  NEW (π=1) : min=0.0000  max=0.0608      ← OLD와 거의 동일 (0.09% diff)
  NEW (π-w) : min=0.0000  max=3.6509      ← 정확히 60배 ≈ n × OLD
```

**PCA 임베딩 관점에서 스케일 factor 는 무관** — 거리 행렬에 상수곱을 해도 classical MDS 임베딩의 모양(회전·반사 무시)은 동일. 따라서 **figure 5 의 diffusion-distance 임베딩 plot 은 올바르다**.

---

# 3. 그럼에도 고쳐야 하는 이유

## 3.1 Frobenius 스케일에 의존하는 지표

편차분석 ([260424 $\tau^*$ 포스트](260424_a5d784.html)) 의 외재적 버전:

$$\Delta^{\text{ext}}(t) \;=\; \|D(t) - D^{\text{lin,ext}}(\alpha(t))\|_F, \qquad D^{G,\text{diff}} = (\text{diffusion distance}^2)$$

여기서는 **상수 스케일 차이가 결과에 직접 반영**. $\alpha$-재매개화가 endpoint-calibrated 라 비율은 보존되지만, **$A^{\text{ext}}_{\text{HST}}$ 의 절대값**, **$\Delta^{\text{ext}}_\infty$**, 그리고 내재 vs 외재 비교는 스케일에 따라 해석이 달라진다.

예제2 결과 (TAU=500k, gauss, 세 metric):

| metric | OLD diff (π=1) | NEW diff (π-weighted) |
|:---|:---|:---|
| Frob ext τ* | 598, max=75.2 | 500000, max=87.5 (단조증가) |
| Proc ext τ* | 500000, max=0.037 | 500000, max=0.896 |

$\pi$-가중 버전에서는 external $\Delta$ 가 **끝까지 단조증가** — Harris 극한 $D^{G,\text{int}}$ 이 π-가중 diffusion distance 로부터 매우 멀리 떨어져 있다. 이것이 "HST가 classical graph metric 과 근본적으로 다르다"의 *정량적* 증거이자, 스케일에 따라 완전히 다른 서사가 나옴.

## 3.2 비정규 그래프에서의 치명성

MCU 네트워크, stochastic block model, ER 랜덤, 실데이터 등은 **정규 아님**. 이 경우:

- OLD 의 $P$ 는 row-stochastic 이 아니며, 그 거듭제곱 $P^t$ 는 전이확률 해석이 무너짐.
- `col_sum[0]` 이 우연히 큰/작은 특이값이면 전체 `P` 의 스케일이 왜곡.
- `P^t` 의 entry 들이 확률이 아니라 임의의 양수 → 이후 L2 거리가 "diffusion" 이라는 이름값을 잃음.

논문의 **MCU figure** 를 포함해 비정규 graph 에서 그린 모든 diffusion-distance 결과는 재계산해야 한다.

---

# 4. 올바른 정의 (수정 후)

`pyhst/distances.py`:

```python
def diffusion_distance(W, t=4, weighted=True):
    """Coifman-Lafon diffusion distance².

    P = D_out^{-1} W,  π = stationary(P),
    d²(i,j) = Σ_k (P^t[i,k] - P^t[j,k])² / π_k      (weighted=True)
            = Σ_k (P^t[i,k] - P^t[j,k])²             (weighted=False)
    """
    W = np.asarray(W, dtype=float)
    row_sum = W.sum(axis=1)
    mask = row_sum > 0
    P = np.zeros_like(W)
    P[mask] = W[mask] / row_sum[mask, None]         # row-stochastic

    Pt = np.linalg.matrix_power(P, t)

    if not weighted:
        return l2distance(Pt)

    # π: left eigenvector of P (power iteration)
    n = P.shape[0]
    pi = np.ones(n) / n
    for _ in range(2000):
        pi_new = pi @ P
        pi_new = pi_new / pi_new.sum()
        if np.max(np.abs(pi_new - pi)) < 1e-12:
            pi = pi_new; break
        pi = pi_new

    pi_safe = np.where(pi > 1e-14, pi, 1.0)
    diff2 = (Pt[:, None, :] - Pt[None, :, :]) ** 2
    return np.sum(diff2 / pi_safe[None, None, :], axis=-1)
```

`weighted=False` 로 호출하면 **기존 OLD와 정확히 동일한 값** (정규 그래프 기준) 이 나오므로, 과거 실험 결과 재현에 쓸 수 있다.

---

# 5. 실질적 가이드

| 상황 | 권장 |
|:---|:---|
| 예제2 (Gaussian ring) 의 기존 figure 5 임베딩 plot | **그대로 사용 OK.** 스케일만 다르고 모양은 동일. |
| 편차분석 (Δ(t), τ*_ext, A_ext 등) 절대값 해석 | 수정된 정의로 재실행 필요. α-재매개화 덕분에 비율·순위는 비슷하지만 절대 스케일 바뀜. |
| MCU, ER, SBM, 실데이터의 diffusion-distance 결과 | **재계산 필수.** 정규 아님 → OLD 는 의미 없는 값이었음. |
| 논문 figure 중 비정규 graph 관련 diffusion-distance | **전수 재검토 필요.** |

---

# 6. 열린 문제

**O9. diffusion scale $t$ 선택.** 현재 $t=4$ 를 관례로 썼지만, 정규화가 제대로 된 지금은 $t$ 선택이 결과에 더 명확히 반영된다. Random walk mixing time 과의 관계, 그래프 spectral gap 과의 scaling 등이 재조사 대상.

**O10. $\pi$-가중 vs 비가중.** 편차분석에서 $D^{G,\text{diff}}$ 를 π-가중으로 쓸지 비가중으로 쓸지는 "비교 기준" 의 성격 문제. 비가중은 PCA 임베딩 모양 비교에 더 자연스럽고, π-가중은 확률적 의미가 더 명료. 두 버전으로 병렬 실험 필요.

---

# 7. 결론

- **예제2 기존 그림**: 정의가 틀렸지만 정규성 덕에 모양은 맞음 → 그대로 써도 무방.
- **편차분석 절대값**: 재실행 필요.
- **비정규 그래프 결과 전반**: 전수 재검증.
- `pyhst/distances.py` 는 표준 Coifman-Lafon 정의로 수정. 구 behavior 는 `weighted=False` 로 복원 가능.

이 오류는 오래된 노트북 (`old/2022-OHS-HST/...`) 에서 비롯된 것으로 추정되며, 당시에는 ex2 같은 정규 그래프에서만 검증돼서 미발견 상태로 남았던 것으로 보인다.
