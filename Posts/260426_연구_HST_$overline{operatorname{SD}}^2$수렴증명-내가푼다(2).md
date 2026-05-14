---
title: 연구 ▷ HST ▷ $\overline{\operatorname{SD}}^2$ 수렴 증명 - 내가푼다 (2)
author: 신록예찬
date: 04/26/2026
draft: false
output-file: 260426_03a7dc.html
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

> 여기에서  $C_1 := n\big(C_3/(b|\Delta_{\mathrm{drift}}|)\big)^2$

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

`-` 이제 아래에 관심을 가지자. 

$$\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] \;=\; 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\hbar(\tilde X_j,\, t_r{+}j)\,\Big|\,S(t_r)\bigg] - m_r\,\tfrac{(n-1)b^2}{n} $$
좀더 정리하면 아래와 같이 쓸 수 있다. 

$=\; 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r{+}j) -\hbar(\tilde{X}_1,t_r+1)\big)\Big| S(t_r)\bigg] +2bm_r\mathbb{E}\left[\hbar(\tilde{X}_1,t_r+1) \mid S(t_r) \right]-  m_r\,\tfrac{(n-1)b^2}{n}$
이때 $\mathbb{E}\left[\hbar(\tilde{X}_1, t_r+1) \mid S(t_r)\right]=\sum_{i=1}^{n}\mu_0(i)\hbar(i,t_r)+\frac{b(n-1)}{n}\leq \frac{\alpha}{2}M(t_r)+\frac{b(n-1)}{n}$를 쓰면 

$\leq 2b\, \mathbb{E}\bigg[\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r{+}j) -\hbar(\tilde{X}_1,t_r+1)\big)\Big| S(t_r)\bigg] +bm_r\alpha M(t_r)+m_r \frac{b^2(n-1)}{n}$

이 계산된다.  그런데 $\sum_{j=1}^{m_r}\big(\hbar(\tilde X_j, t_r{+}j) -\hbar(\tilde{X}_1,t_r+1)\big)=\bigstar$ 에 대하여 아래의 식을 유도할 수 있다. 

- 라운드내 낙차가 없는 경우: $\bigstar<C_{round}$
- 라운드내 낙차가 있는 경우:$\bigstar \leq -\frac{M(t_r)}{n-1} +C_{round}$
 
이를 이용하여 정리하면 

$$\;\leq\; 2bC_{round} - \frac{2b p_0 M(t_r)}{n-1} +bm_r\alpha M(t_r)+m_r \frac{b^2(n-1)}{n}$$
정리하면

$$=bM(t_r)\left[m_r \alpha -\frac{2p_0}{n-1}\right]+\left[2bC_{round}+m_r\frac{b^2(n-1)}{n} \right]$$

$m_r \leq T_{max}+1$로 바운드치면 

$$\leq bM(t_r)\left[(T_{max}+1) \alpha -\frac{2p_0}{n-1}\right]+\left[2bC_{round}+(T_{max}+1)\frac{b^2(n-1)}{n} \right]$$
따라서 

$$\begin{aligned}
\mathbb{E}[\Phi(t_{r+1}) - \Phi(t_r) \mid S(t_r)] &\leq -\frac{2b\,p_0\,M(t_r)}{n-1} + 2b\,C_{round} + b(T_{max}+1)\,\alpha\,M(t_r) + (T_{max}+1)\frac{b^2(n-1)}{n} \\
&= b\,M(t_r)\underbrace{\Big[(T_{max}+1)\alpha - \frac{2p_0}{n-1}\Big]}_{=\,\Delta_{\mathrm{drift}}} + \underbrace{\Big[2b\,C_{round} + (T_{max}+1)\frac{b^2(n-1)}{n}\Big]}_{=:\,C_3}
\end{aligned}$$

여기서 $C_3 := 2bC_{round} + (T_{max}+1)b^2(n-1)/n$ 은 $M(t_r)$ 에 무관하고 $(\mathcal{G}, b, T_{max})$ 에만 의존하는 결정적 상수.

---

결국 DC아래에서 위의 $M(t_r)$이 충분히크면 음이되고, (충분히 크지 않으면??) Foster-Lyapunov정리를 쓸 수 있다. 

