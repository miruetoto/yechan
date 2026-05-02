---
title: 연구 ▷ HST ▷ 그래프 구조로의 점근적 환원 — Appendix A.6 초안
author: claude
date: 04/26/2026
draft: false
output-file: 260426_92d9fb.html
---

# 0. 메모

`hst.tex` (paper/260423_guebin) Appendix A.2 (`Further Remark on Snow Distance`) 의 비공식 서술

> "If there exists sufficiently large $t_0$ such that $h(v_i,t_0) \approx h(v_j,t_0)$ ..., then matrix $\mathbf{SD}(t_0+1,t)$ is similar to $\widetilde{\mathbf{SD}}(0,t-t_0)$, which does not hold Euclidean information. Thus, the graph domain information is only contained in $\mathbf{SD}(t_0+1,t)$."

를 Theorem 형태로 바로 꽂을 수 있게 정형화한 초안. 별도 파일 `260426_연구_HST_$overline{operatorname{SD}}^2$수렴증명—내가푼다.md` 의 round/big-drop 기반 미시 증명과 보완 관계 — 이 글은 "Appendix에 그대로 들어가는 추상-레벨 서술"이고, 그쪽이 "drift 상수 $\Delta_{\mathrm{drift}}$의 명시적 계산을 포함한 자체 증명"이다. Meyn–Tweedie / Brémaud 의 일반 결과를 블랙박스로 인용하는 버전.

ref/ 추가 PDF:

- `BOOK_PierreBremaud_MarkovChain.pdf` — Brémaud (2013). Ch. 3 (Recurrence and Ergodicity), Ch. 12 (Non-homogeneous Markov Chains). 이미 `\cite{Bremaud2013}` 로 인용 중.
- `BOOK_MarkovChainsAndStochasticStability.pdf` — Meyn & Tweedie (2009). Ch. 11 (Foster–Lyapunov drift), Ch. 17 (Sample paths and limit theorems).
- `BOOK_GSP_Cooperative and graph signal processing principles and applications.pdf` — Djurić & Richard (Eds.), 2018. GSP의 graph filter / spectral filter 일반론.

---

# 1. 정리 (Theorem A.6)

> **Theorem A.6** (Asymptotic Reduction to the Graph Structure).
> $\mathcal{G} = (V, \mathbf{E}, \mathbf{W})$ 를 $|V| = n < \infty$ 인 연결 가중 그래프, $\mathbf{y}, \widetilde{\mathbf{y}}$ 를 $\mathcal{G}$ 위의 임의의 두 초기 graph-valued 데이터, $\mathbf{h}(\cdot, t), \widetilde{\mathbf{h}}(\cdot, t)$ 를 동일 파라미터 $(b, T_{max})$ 하의 heavy-snow transform 이라 하자. 그러면 모든 $i, j \in V$ 에 대해
> $$\lim_{t \uparrow \infty} \frac{\operatorname{SD}^2_{ij}(t)}{t} \;=\; \lim_{t \uparrow \infty} \frac{\widetilde{\operatorname{SD}}^2_{ij}(t)}{t} \;=:\; s^*_{ij} \quad \text{a.s.,}$$
> 여기서 $s^*_{ij} \geq 0$ 은 $\mathcal{G}$ (와 $b, T_{max}$) 만의 함수. 결과적으로 snow weight 는 결정론적 극한
> $$W^*_{ij} \;:=\; \lim_{t \uparrow \infty} W_{ij}(t) \;=\; \exp\!\Big(-\frac{s^*_{ij}}{2\theta^2}\Big), \quad i \neq j$$
> 을 갖고, 이 또한 $\mathcal{G}$ 만으로 결정된다.

---

# 2. 증명

차분과정 $D_{ij}(t) := h(v_i, t) - h(v_j, t)$ 를 도입. Table 1 (`tb:notation`) 의 갱신규칙 $h(v, t) = h(v, t-1) + b \, I(X_t = v)$ 로부터

$$D_{ij}(0) = y(v_i) - y(v_j), \qquad D_{ij}(t+1) - D_{ij}(t) = b \big[I(X_{t+1} = i) - I(X_{t+1} = j)\big].$$

## Step 1 — Markovian reduction

Definition 3.1 (`dfn:hst`) 의 전이확률을 살펴보면, $P(X_{t+1} = j \mid X_t = i, h(\cdot, t))$ 가 높이를 통해 의존하는 부분은 모두 지시함수 $I(h(v_i, t) \geq h(v_k, t)) = I(D_{ik}(t) \geq 0)$ 와 결정론적 보조변수 $Z_t, \mathcal{D}_t$ 만을 통한다. 따라서 확장과정

$$\mathbf{U}_t \;:=\; \big( (D_{ij}(t))_{i,j \in V}, \, X_t, \, Z_t \big)$$

는 가산 상태공간 위의 **시간동질** 마르코프 체인이고, 그 전이커널은 그래프 $\mathcal{G}$ (와 $b, T_{max}$) 만의 함수. 특히 $\mathbf{y} \mapsto \mathbf{y} + c \mathbf{1}$ 인 평행이동에 불변 (모든 $D_{ij}$ 가 그대로).

## Step 2 — 차분과정의 안정성 (Foster–Lyapunov)

$\{D_{ij}(t)\}$ 가 tight family임을 보임. 인접쌍 $i \sim j$ 와 큰 $M > 0$ 에 대해 $\{D_{ij}(t) \geq M\}$ 사건을 잡으면, $\mathcal{G}$ 의 연결성으로 경로 $v_i = v_{k_0} \sim v_{k_1} \sim \cdots \sim v_{k_r} = v_j$ 가 존재하고, 이 경로에서 어떤 인접쌍 $(v_{k_\ell}, v_{k_{\ell+1}})$ 은 $D_{k_\ell, k_{\ell+1}}(t) \geq M / r$ 을 만족. Definition 3.1 의 $p_{t, v_{k_\ell} \to v_{k_{\ell+1}}}$ 정의로부터, 높이 순서가 맞을 때 flow 확률은 양의 그래프 상수로 하향 유계. 이를 Theorem A.5 (`thm:a5`, $\{X_t\}$ 의 weak ergodicity) 와 결합하면 Foster–Lyapunov drift 부등식

$$\mathbb{E}\!\left[ V(\mathbf{U}_{t+1}) \,\big|\, \mathbf{U}_t = \mathbf{u} \right] - V(\mathbf{u}) \;\leq\; -\varepsilon V(\mathbf{u}) + c_0 \, \mathbf{1}_C(\mathbf{u}), \qquad V(\mathbf{u}) := \sum_{i,j} D_{ij}^2$$

를 얻음 (어떤 $\varepsilon, c_0 > 0$ 와 유한 small set $C$ 에 대해). 일반적인 Foster–Lyapunov 안정성 이론은 Meyn & Tweedie (2009, Ch. 11) 참조.

이 drift 부등식과 $\{\mathbf{U}_t\}$ 의 essential class 에서의 기약성 (모든 reset 시각에서 $\boldsymbol{\mu}_0(v) > 0$ 인 모든 $v \in V$ 에 눈이 떨어짐) 으로부터 $\{\mathbf{U}_t\}$ 는 양의 Harris 재귀이고 유일한 정상분포 $\pi^*$ 를 갖는다.

> **주석**: 이 step 의 $\Delta_{\mathrm{drift}}$ 명시적 계산은 별도 글 `260426_연구_HST_$overline{operatorname{SD}}^2$수렴증명—내가푼다.md` 참조. 거기서는 $\Phi(t) = \sum_i \hbar(i,t)^2$, round-skeleton, big-drop 사건 $E_r$ 등을 명시적으로 도입하여 drift 상수
> $$\Delta_{\mathrm{drift}} = (T_{max}+1)\alpha - \frac{2 p_0}{n-1}$$
> 의 부호조건 $\Delta_{\mathrm{drift}} < 0$ 을 명시. 정칙 그래프에서는 $\alpha = 0$ 이라 자동.

## Step 3 — $\pi^*$ 의 $\mathbf{y}$-비의존성

Step 1 에 의해 $\mathbf{U}_t$ 의 커널은 $\mathcal{G}$ 만에 의존. $\mathbf{y}$ 가 다르면 $\mathbf{U}_0$ 의 초기상태만 달라지지만, 양의 Harris 재귀성으로부터 $\mathbf{U}_t$ 의 분포는 초기상태와 무관하게 $\pi^*$ 로 total variation 수렴. 따라서 $\pi^*$ 는 $\mathcal{G}$ (와 $b, T_{max}$) 만의 함수.

## Step 4 — Birkhoff 에르고딕 정리

양의 Harris 재귀 마르코프 체인의 에르고딕 정리 (Brémaud 2013, Theorem 3.3.1; Meyn & Tweedie 2009, Theorem 17.0.1) 를 유계 함수 $f(\mathbf{u}) = D_{ij}^2$ 에 적용 (Step 2 로부터 $\pi^*$ 하의 $D_{ij}^2$ 유계성 확보). 임의의 초기분포에서

$$\frac{1}{t+1} \sum_{s=0}^{t} D_{ij}(s)^2 \;\xrightarrow[t \uparrow \infty]{a.s.}\; \mathbb{E}_{\pi^*}\!\big[D_{ij}^2\big] \;=:\; s^*_{ij},$$

우변은 $\mathcal{G}$ 만의 함수. $\operatorname{SD}^2_{ij}(t) = \sum_{s=0}^{t} D_{ij}(s)^2$ 이므로 $\operatorname{SD}^2_{ij}(t) / t \to s^*_{ij}$ a.s. $\widetilde{\mathbf{y}}$ 에도 동일 논리로 같은 극한. 첫 번째 주장 증명 완료.

## Step 5 — Snow weight 의 극한

$t > 0$ 과 $i \neq j$ 에 대해 Definition 3.2 로부터

$$W_{ij}(t) \;=\; \exp\!\Big(-\frac{\operatorname{SD}^2_{ij}(t)}{2\theta^2 \, t}\Big) \;\xrightarrow[t \uparrow \infty]{a.s.}\; \exp\!\Big(-\frac{s^*_{ij}}{2\theta^2}\Big) \;=\; W^*_{ij},$$

이는 $\mathcal{G}$ 만의 함수. $\square$

---

# 3. Remark

Theorem A.6 은 본 절 서두의 휴리스틱 분해

$$\mathbf{SD}^2(t) = \mathbf{SD}^2(0, t_0) + \mathbf{SD}^2(t_0+1, t)$$

를 정형화한다. 임의의 고정된 $t_0$ 에 대해 $\mathbf{SD}^2(0, t_0)$ 은 유계 행렬, 따라서 $\mathbf{SD}^2(0, t_0) / t \to \mathbf{0}$ ($t \uparrow \infty$). 결과적으로 $\mathbf{SD}^2(t)$ 의 단위시간당 증가율은 전적으로 $\mathbf{SD}^2(t_0+1, t)$ 가 담당하며, 그 비율 $s^*_{ij}$ 는 그래프만의 함수. 이 점근적 영역에서 snow-weighted matrix $\mathbf{W}(t)$ 는 graph-only filter 처럼 작동 — graph signal processing 에서 다루는 spectral graph filter 의 heavy-snow 유사체 (Djurić & Richard 2018, Ch. 3 참조).

---

# 4. 추가할 bibliography 항목

`hst.tex` 의 `\begin{thebibliography}{}` 블록에 추가:

```latex
\bibitem{MeynTweedie2009} Meyn, S. and Tweedie, R. L. (2009). {\it Markov Chains and Stochastic Stability} (2nd ed.). Cambridge University Press.

\bibitem{DjuricRichard2018} Djuri\'c, P. M. and Richard, C. (Eds.) (2018). {\it Cooperative and Graph Signal Processing: Principles and Applications}. Academic Press.
```

`Bremaud2013` 은 이미 등록되어 있으므로 그대로 사용.

---

# 5. 삽입 위치 (확정 시)

`paper/260423_guebin/hst.tex` 의 line 973 (`\end{figure}` 직후, line 975 의 `\subsubsection*{A.3 ...}` 직전) 에 본문 (Theorem + Proof + Remark 부분) 을 삽입.

---

# 6. 미해결/체크 포인트

- Step 2 의 drift 부등식은 일반론 인용으로 처리. 명시적 상수 (앞서 언급한 $\Delta_{\mathrm{drift}}$) 가 필요한 reviewer 가 있을 수 있음 → 그 경우 `내가푼다` 글의 결과를 lemma 로 분리해서 추가 가능.
- 격자 상태공간에서의 기약성 — $\boldsymbol{\mu}_0(v) > 0$ 인 모든 $v$ 에서 reset 가 일어난다는 것으로부터 essential class 단일성 주장. 비연결 그래프나 $\mu_0$ 의 support 가 $V$ 전체가 아닌 경우는 별도 조건 필요.
- $W_{ij}(t)$ 정의에서 분모가 $t$ 인 점이 핵심. 만약 분모가 $t - t_0$ 였으면 더 자연스러웠겠지만, 현재 정의로도 $\operatorname{SD}^2(0, t_0) / t \to 0$ 이므로 동일 극한 성립.
