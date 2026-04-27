---
title: 연구 ▷ HST ▷ $\overline{\operatorname{SD}}^2$ 수렴 증명 — 내가푼다
author: 신록예찬
date: 04/26/2026
draft: false
output-file: 260426_c05541.html
---
---

**Theorem.** 연결 가중 그래프 $\mathcal{G}$와 파라미터 $(b, T_{max})$가 주어질 때, 다음 drift 조건이 성립한다고 하자:

$$\Delta_{\mathrm{drift}} := (T_{max}{+}1)\alpha - \frac{2p_0}{n-1} < 0 \tag{DC}$$

여기서 $\alpha = \|\boldsymbol{\mu}_0 - \mathbf{u}\|_1$ (차수 비균등도), $\mu_{min} = \min_i \mu_0(v_i)$, $p_0 =\mu_{min}P_{RW,min}$ (1라운드 big-drop확률하한)이다. DC가 만족되면 상수 행렬 $C = [c_{ij}]_{n \times n}$, $c_{ij} = c_{ij}(\mathcal{G}, b, T_{max}, \mathbf{y} \bmod b) \geq 0$이 존재하여, 임의의 초기 신호 $\mathbf{y}$에 대해:

$$\overline{\mathbf{SD}}^2(t) \xrightarrow{t \to \infty} C \quad \text{a.s.}$$

이때 한계값 $c_{ij}$는 초기 신호 $\mathbf{y}$의 거시적인 절대 스케일에는 완전히 독립적이며, 오직 이산 업데이트 규칙의 격자 경계를 결정짓는 $b$ 단위의 미세한 위상차(fractional part modulo $b$)에만 국소적으로 종속된다.

---

이걸 증명하기 위해 가장 도전적인 부분은 아래를 보이는것임. 

Proposition (Height Recurrence). 높이 범위 $M(t)$에 대해, **drift 조건**

$$\Delta_{\mathrm{drift}} := (T_{max}{+}1)\alpha - \frac{2p_0}{n-1} \;<\; 0 \tag{DC}$$

이 성립하면 ($\mathcal{G}$가 정칙이면 $\alpha = 0$이므로 자동 성립), $(\mathcal{G}, b, T_{max})$에만 의존하는 상수 $C > 0$이 존재하여, **높이차 과정 $\{S(t)\}$는 유한 집합 $\mathcal{B}_\Phi$로 유한 기대 시간 내에 반복 복귀하는 양의 재귀 마르코프 체인**이다. 특히 $M(t)$는 $\{M \leq C\}$로 반복 복귀.

> 질문: $\mathcal{B}_\Phi$ 은 뭐지?? 아무렇게나 잡으면 되나?? 유한기대시간내에 반복복귀하도록 하는 $\mathcal{B}_\Phi$가 존재한다는 것이지? $\{S(t)\}$가 $\mathcal{B}_\Phi$에 반복복귀한다는거랑 $M(t)$가 $\{M\leq C\}$로 반복복귀한다는건 같은 뜻아닌가? 다른거야?? 그리고 $M(t)$는 랜덤아니야? 따라서 $C$라고 쓰면안되고 $O_p(1)$ 이런식으로 바운드를 잡아야하지 않어??

> **답**:
>
> **(a) $\mathcal{B}_\Phi$ 의 정체 — 아무렇게나 잡는 게 아님.** Foster-Lyapunov 기준은 "**어떤** 유한 집합 $\mathcal{B}$에 대해 그 바깥에서 음의 drift가 성립하면 양의 재귀" 라는 정리. 즉 $\mathcal{B}$는 "존재성"으로 주어지는 것이지 임의 선택이 아님. 본 증명에서는 구체적으로
> $$\mathcal{B}_\Phi := \{S : \Phi(\boldsymbol{\delta}) \leq C_1\}, \qquad C_1 := n\big(C_3/(b|\Delta_{\mathrm{drift}}|)\big)^2$$
> 로 잡는다 ($C_1$은 그래프-결정 상수 — 아래 질문2/3에서 보임). 이때 격자 $\boldsymbol{\delta}(0) + b\mathbb{Z}^{n-1}$ 과 사블레벨 $\{\Phi \leq C_1\}$ 의 교집합은 유한이라 $\mathcal{B}_\Phi$ 가 유한 집합이 됨 (Foster의 전제 충족).
>
> **(b) "$\{S\}$ 반복복귀" vs "$M$ 반복복귀" — 동치는 아니지만 양방향 implication.** $\Phi(t) = \sum_i \hbar(i,t)^2$ 이고 $\Phi \geq M^2/2$ (왜냐: 두 극단점 $\hbar_{\max}, \hbar_{\min}$ 이 $\hbar_{\max}^2 + \hbar_{\min}^2 \geq (\hbar_{\max} - \hbar_{\min})^2/2 = M^2/2$). 따라서:
> - $\Phi \leq C_1 \Longrightarrow M \leq \sqrt{2C_1}$. ($S \in \mathcal{B}_\Phi \Rightarrow M \leq C := \sqrt{2C_1}$)
> - 역도 정성적으로 성립: $\sum_i \hbar^2 \leq nM^2$ 이므로 $M \leq C \Longrightarrow \Phi \leq nC^2$, 즉 약간 큰 사블레벨에 들어감.
>
> 두 표현 모두 "유한 집합으로 무한히 자주 복귀" 라는 같은 정성적 사실을 다른 좌표(높이 분산 vs. 높이 범위)로 본 것. Proposition 본문은 후자 쪽 ("$M(t) \leq C$ 반복복귀") 만 명시한 것.
>
> **(c) $C$ vs $O_p(1)$ — $C$ 가 맞고 $O_p(1)$ 이 아님.** $M(t)$ 자체는 확률변수. 하지만 명제는 "**결정적 상수** $C$ 가 존재해서 $M(t) \leq C$ 인 사건이 a.s. 무한히 자주 일어난다" 라는 거. 즉 random path 하나하나가 무작위지만, 모든 경로가 같은 결정적 박스 $\{M \leq C\}$ 안을 무한 자주 방문함 (양의 Harris 재귀의 정의). $O_p(1)$ 은 단일 시각 $t$에서의 분포 tightness 개념이라 약함 — 본 증명은 그보다 강한 a.s. 양의 재귀.

Proposition을 증명하기 위해서는 Foster-Lyapunov정리를 써야함. 그런데 이걸 위해서는 적당한 함수 $\Phi(t)$를 잡아서 그 함수가 $t$가 커질때 점점 작아진다는 느낌을 줘야함. 그러니까 $\Phi(t+1) - \Phi(t)$를 음수로 바운드하면 좋음. 그런데 그게 쉽진않음. 그래서 우리는 적당히 round의 개념을 도입하여 1라운드 이후 $\Phi$의 값이 줄어듦을 보일것임. 

$$\mathbb{E}[\Phi(t_{r+1})-\Phi(t_r) | S(t_r)] \leq 0$$

그런데 식을 전개하다보면 아래가 보여진다. 

$$\mathbb{E}[\Phi(t_{r+1})-\Phi(t_r) | S(t_r)] \leq bM(t_r)\left[(T_{max}+1)\alpha - \frac{2p_0}{n-1}\right] +C_3$$

대괄호속에있는게 $\Delta_{\text{drift}}$인데 이것이 음수라고 가정했으므로 $bM(t_r)$가 충분히 크면 음수가 나온다. 여기가 헷갈리는 부분인데 

1. $M(t_r)$이 충분히 크면 $\Phi(t_r)$ 가 천천히 줄어든다는 논리전개가 가능
2. $M(t_r)$가 충분히 크지 않으면 원래 잘 수렴하므로 상관없다(? )

이런식임. 사실 이런식의 논리전개가 가능한지 모르겠음. 딱 1과 2의 경계가 되는 $M(t)$가 있었으면 좋겠지? 

> 질문1: 이거 이런논리전개가 가능해? 2의 경우 충분히 크지 않아도 잘 수렴하는지 체크해야하지않어 
> 질문2: 경계가 되는 $M(t)$는?

> **답**:
>
> **(a) "충분히 크지 않으면 어떻게 되나?" — Foster-Lyapunov가 우아하게 처리.** 핵심: Foster-Lyapunov 정리는 **refuge $\mathcal{B}_\Phi$ 바깥에서만** 음의 drift 를 요구함. refuge 안에서는 drift 가 양수든 음수든 상관없음 (단, 한 스텝 변화량이 유한이면 됨). 직관:
> - $M(t_r)$ 큼 ($\Rightarrow$ refuge 밖): 음의 drift, $\Phi$ 가 평균적으로 줄어듦 → 결국 refuge 들어감.
> - $M(t_r)$ 작음 ($\Rightarrow$ refuge 안): drift 가 양수일 수 있어도 무한 발산 못함. 왜냐하면 (i) 한 스텝에 $|\Phi$ 변화량$|$ 유한 ($\leq O(b M_{\text{round}}) + O(b^2)$), (ii) 어떤 양수 변화로 refuge 를 빠져나가면 그 다음부터 다시 음의 drift 가 작용해 끌어당김. 즉 "refuge 밖으로 나가도 다시 끌려들어오는" 메커니즘이 자동.
>
> 그래서 "2의 경우" 를 별도로 체크할 필요 없음 — Foster-Lyapunov 가 한 묶음으로 처리. (경로 단위 a.s. 양의 재귀가 정리의 결론.)
>
> **(b) 경계 $M^\star$ 의 명시적 식.** Drift 부등식
> $$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \leq bM(t_r)\Delta_{\mathrm{drift}} + C_3$$
> 의 우변이 $\leq -1$ 이 되려면
> $$bM(t_r)\Delta_{\mathrm{drift}} + C_3 \leq -1 \;\iff\; M(t_r) \;\geq\; \frac{C_3 + 1}{b\,|\Delta_{\mathrm{drift}}|} \;=:\; M^\star.$$
> 따라서 refuge $\mathcal{B}_\Phi := \{\Phi \leq C_1\}$ 에서 $C_1$ 을 잡을 때 $\Phi \leq C_1 \Rightarrow M \leq M^\star$ 이 되도록, $\Phi \leq nM^2$ 관계에서 $C_1 := n(M^\star)^2$ 정도 잡으면 됨. 정확히는
> $$M^\star := \frac{C_3 + 1}{b\,|\Delta_{\mathrm{drift}}|}, \qquad C_1 := n(M^\star)^2.$$
> 이 두 값이 모두 $(\mathcal{G}, b, T_{max})$ 에만 의존하는 결정적 상수.

암튼 이 식까지 나오면 나머지는 나름대로 자명해보임. 이 식을 유도하는 과정이 문제가 있는지 살펴보자. 증명을 위해서 아래를 계산한다. 

`-` Lemma1: $\Phi(t) - \Phi(t-1) = 2b\,\hbar(\ell, t) + \frac{(n-1)b^2}{n}$ 

`-`  ($\dagger$): $\mathbb{E}[\Phi(t+1) - \Phi(t) \mid S(t)] = 2b \cdot \mathbb{E}\big[\hbar(X_{t+1}, t+1) \mid S(t)\big] + \frac{(n-1)b^2}{n}$

`-` 라운드별 기여 상한
$$\sum_{j=1}^{m_r}\hbar(\tilde{X}_{j}, t_r+j)<m_r \hbar(\tilde{X}_1,t_r+1)+C_{round}$$
한번 정리하면 아래와 같이 쓸 수 있다. 

$$\sum_{j=1}^{m_r}\left(\hbar(\tilde{X}_{j}, t_r+j) - \hbar(\tilde{X}_1,t_r+1)\right)<C_{round}$$
`-` 그 라운드에 만약에 큰 낙차가 포함되어있다?? 그러면 아래와 같이 됨. 

$$\sum_{j=1}^{m_r}\left(\hbar(\tilde{X}_{j}, t_r+j)-\hbar(\tilde{X}_{1}, t_r+1)\right)\leq -(m_r-1)\frac{M(t_r)}{n-1} +C_{round}$$
이때 $m_r \geq 2$ 이므로 

$$\sum_{j=1}^{m_r}\left(\hbar(\tilde{X}_{j}, t_r+j)-\hbar(\tilde{X}_{1}, t_r+1)\right)\leq -\frac{M(t_r)}{n-1} +C_{round}$$

`-` 이제 아래에 관심을 가지자. 

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \;=\; 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\hbar(\tilde X_j,\, t_r{+}j)\,\Big|\,S(t_r)\bigg] - m_r\,\tfrac{(n-1)b^2}{n} $$
좀더 정리하면 아래와 같이 쓸 수 있다. 

$=\; 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r{+}j) -\hbar(\tilde{X}_1,t_r+1)\big)\Big| S(t_r)\bigg] +2bm_r\mathbb{E}\left[\hbar(\tilde{X}_1,t_r+1) \mid S(t_r) \right]-  m_r\,\tfrac{(n-1)b^2}{n}$
이때 $\mathbb{E}\left[\hbar(\tilde{X}_1, t_r+1) \mid S(t_r)\right]=\sum_{i=1}^{n}\mu_0(i)\hbar(i,t_r)+\frac{b(n-1)}{n}\leq \frac{\alpha}{2}M(t_r)+\frac{b(n-1)}{n}$를 쓰면 

$\leq 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r{+}j) -\hbar(\tilde{X}_1,t_r+1)\big)\Big| S(t_r)\bigg] +bm_r\alpha M(t_r)+m_r \frac{b^2(n-1)}{n}$

이 계산된다. 그런데 남은 부분의 합을 정리해야 한다.

**(II) 큰 낙차 항의 기댓값 상한.** 합의 $j=1$ 항은 $\hbar(\tilde X_1, t_r+1) - \hbar(\tilde X_1, t_r+1) = 0$ 이므로 실제로는

$$\mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)\big) \,\Big|\, S(t_r)\bigg] = \mathbb{E}\bigg[\sum_{j=2}^{m_r}\big(\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)\big) \,\Big|\, S(t_r)\bigg]$$

big-drop 사건 $E_r$ ($\tilde X_1 = i^\star$, $\tilde X_2 = j^\star$, 확률 $\Pr(E_r \mid S(t_r)) \geq p_0$) 의 발생 여부로 케이스 분해:

- $E_r$ 발생 시 ($\spadesuit$): $\sum_{j=2}^{m_r} (\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)) \leq -(m_r-1)\frac{M(t_r)}{n-1} + C_{round} \leq -\frac{M(t_r)}{n-1} + C_{round}$  ($m_r \geq 2$)
- $E_r^c$: B.3 일반 bound $\sum (\hbar - \hbar_1) \leq C_{round}$

기댓값 결합 ($\Pr(E_r) \geq p_0 \Rightarrow$ 음의 항이 $\geq -p_0 M/(n-1)$, 부등호 방향 주의):

$$\mathbb{E}\bigg[\sum_{j=2}^{m_r}\big(\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)\big) \,\Big|\, S(t_r)\bigg] \leq -\frac{p_0\,M(t_r)}{n-1} + C_{round}$$

양변에 $2b$ 곱하면:

$$2b\,\mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)\big) \,\Big|\, S(t_r)\bigg] \;\leq\; -\frac{2b\,p_0\,M(t_r)}{n-1} + 2b\,C_{round}$$

**최종 결합.** 앞서 정리한

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \;\leq\; 2b\,\mathbb{E}\bigg[\sum_{j=1}^{m_r}(\hbar(\tilde X_j, t_r+j) - \hbar(\tilde X_1, t_r+1)) \,\Big|\, S(t_r)\bigg] + b\,m_r\,\alpha\,M(t_r) + m_r\frac{b^2(n-1)}{n}$$

에 위 (II) bound 와 $m_r \leq T_{max}+1$ 을 대입:

$$\begin{aligned}
\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] &\leq -\frac{2b\,p_0\,M(t_r)}{n-1} + 2b\,C_{round} + b(T_{max}+1)\,\alpha\,M(t_r) + (T_{max}+1)\frac{b^2(n-1)}{n} \\
&= b\,M(t_r)\underbrace{\Big[(T_{max}+1)\alpha - \frac{2p_0}{n-1}\Big]}_{=\,\Delta_{\mathrm{drift}}} + \underbrace{\Big[2b\,C_{round} + (T_{max}+1)\frac{b^2(n-1)}{n}\Big]}_{=:\,C_3}
\end{aligned}$$

따라서 정리하면

$$\boxed{\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \;\leq\; bM(t_r)\,\Delta_{\mathrm{drift}} + C_3} \tag{$\star\star$}$$

여기서 $C_3 := 2bC_{round} + (T_{max}+1)b^2(n-1)/n$ 은 $M(t_r)$ 에 무관하고 $(\mathcal{G}, b, T_{max})$ 에만 의존하는 결정적 상수.

---

결국 DC아래에서 위의 $M(t_r)$이 충분히크면 음이되고, (충분히 크지 않으면??) Foster-Lyapunov정리를 쓸 수 있다. 

> 질문: 아까 했던거랑 비슷한데 충분히 크지 않으면 어떻게 되지?

> **답** (앞 답변 (a) 의 재진술 + 정리):
>
> $M(t_r) < M^\star$ ($:= (C_3+1)/(b|\Delta_{\mathrm{drift}}|)$) 인 영역에서 drift 가 양수가 될 수 있는데, 이건 **걱정할 필요가 없는 영역**임. 이유:
>
> 1. **그 영역 자체가 refuge.** $\mathcal{B}_\Phi := \{\Phi \leq C_1 := n(M^\star)^2\}$ 으로 잡으면 $M < M^\star \Rightarrow \Phi \leq nM^2 < n(M^\star)^2 = C_1 \Rightarrow S \in \mathcal{B}_\Phi$. 즉 **"$M$ 작음" = "이미 refuge 안"**.
>
> 2. **Foster-Lyapunov 의 형식적 요구.** 정리는 다음 형태의 drift 만 요구:
>    $$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \;\leq\; -1 + \beta\,\mathbf{1}_{\mathcal{B}_\Phi}(S(t_r)), \qquad \beta < \infty.$$
>    - **밖** ($S \notin \mathcal{B}_\Phi$, 즉 $M \geq M^\star$): 우변 $= -1$. ($\star\star$) 의 $bM\Delta_{\mathrm{drift}} + C_3 \leq -1$ 이 $M \geq M^\star$ 정의로 자동.
>    - **안** ($S \in \mathcal{B}_\Phi$): 우변 $= -1 + \beta$. drift 가 양수든 음수든 한 스텝 변화량 유한 (예: $|2bM \cdot \alpha + b^2| \leq 2bM^\star + b^2 =: \beta + 1$) 이라 $\beta < \infty$ 만족. 즉 안에선 drift 부호 무관, 변화량만 유한이면 OK.
>
> 3. **결과.** Foster-Lyapunov → $\mathbb{E}[\tau_{\mathcal{B}_\Phi}] < \infty$ (refuge 도달 기대시간 유한) → 양의 재귀. 경로 단위로 $S(t)$ 가 무한히 자주 $\mathcal{B}_\Phi$ 에 복귀, 등가적으로 $M(t) \leq \sqrt{2C_1} =: C$ 에 무한 자주 복귀 (a.s.).
>
> 한 줄 요약: **"$M$ 큼"**은 음의 drift 로 끌어내리고, **"$M$ 작음"**은 이미 refuge 안. 둘 사이 경계 $M^\star$ 가 결정적 상수로 잡혀서 case 분리가 깔끔.