---
title: 연구 ▷ HST ▷ $\overline{\operatorname{SD}}^2$ 수렴 증명 — 내가푼다
author: 신록예찬
date: 04/26/2026
draft: true
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

Proposition을 증명하기 위해서는 Foster-Lyapunov정리를 써야함. 그런데 이걸 위해서는 적당한 함수 $\Phi(t)$를 잡아서 그 함수가 $t$가 커질때 점점 작아진다는 느낌을 줘야함. 그러니까 $\Phi(t+1) - \Phi(t)$를 음수로 바운드하면 좋음. 그런데 그게 쉽진않음. 그래서 우리는 적당히 round의 개념을 도입하여 1라운드 이후 $\Phi$의 값이 줄어듦을 보일것임. 

$$\mathbb{E}[\Phi(\tau_{n+1})-\Phi(\tau_n) | S(\tau_n)] \leq 0$$ 
그런데 식을 전개하다보면 아래가 보여진다. 

$$\mathbb{E}[\Phi(t+T)-\Phi(t) | S(t)] \leq bM(t)\left[(k+1)(T_{max}+1)\alpha - \frac{2q_0}{n-1}\right] +C_3$$

대괄호속에있는게 $\Delta_{\text{drift}}$인데 이것이 음수라고 가정했으므로 $bM(t)$는 음수값이다. 여기가 헷갈리는 부분인데 

1. $M(t)$가 충분히 크면 $\Phi(t)$ 가 천천히 줄어든다는 논리전개가 가능
2. $M(t)$가 충분히 크지 않으면 원래 잘 수렴한다? 

이런식임. 사실 이런식의 논리전개가 가능한지 모르겠음. 딱 1과 2의 경계가 되는 $M(t)$가 있었으면 좋겠지? 

> 질문1: 이거 이런논리전개가 가능해? 
> 질문2: 경계가 되는 $M(t)$는?

암튼 이 식까지 나오면 나머지는 나름대로 자명해보임. 이 식을 유도하는 과정이 문제가 있는지 살펴보자. 증명을 위해서 아래를 계산한다. 

`-` Lemma1: $\Phi(t+1) - \Phi(t) = 2b\,\hbar(\ell, t) + \frac{(n-1)b^2}{n}$ 

`-`  ($\dagger$): $\mathbb{E}[\Phi(t+1) - \Phi(t) \mid S(t)] = 2b \cdot \mathbb{E}\big[\hbar(X_{t+1},\, t) \mid S(t)\big] + \frac{(n-1)b^2}{n}$

`-` $(\dagger\dagger)$: $\mathbb{E}[\Phi(t{+}T) - \Phi(t) \mid S(t)]  = 2b \sum_{s=0}^{T-1} \mathbb{E}\big[\hbar(X_{t{+}s{+}1},\, t{+}s) \mid S(t)\big] + T \cdot \tfrac{(n-1)b^2}{n}$

`-` 결국 아래를 계산해서 싹더하는게 포인트다. 

- $\mathbb{E}\big[\hbar(X_{t+1},\, t) \mid S(t)\big]$ = $t$시점에서 눈이 평균높이보다 더 낮은곳으로 쌓일것으로 기대되면 음수, 그렇지 않으면 양수
- $\mathbb{E}\big[\hbar(X_{t+2},\, t+1) \mid S(t)\big]$ = $t+1$시점에서 눈이 평균높이보다 더 낮은곳으로 쌓일것으로 기대되면 음수, 그렇지 않으면 양수
- $\mathbb{E}\big[\hbar(X_{t+3},\, t+2) \mid S(t)\big]$ = $t+2$시점에서 눈이 평균높이보다 더 낮은곳으로 쌓일것으로 기대되면 음수, 그렇지 않으면 양수
- ... 
- $\mathbb{E}\big[\hbar(X_{t+T},\, t+T-1) \mid S(t)\big]$ = $t+T-1$시점에서 눈이 평균높이보다 더 낮은곳으로 쌓일것으로 기대되면 음수, 그렇지 않으면 양수

`-` 이걸 다 따지면 어려우니까 $\{t,t+1,\dots, t+T-1\}$의 시간을 라운드별로 쪼개자. 최소한 이 시간구간에 $k$개의 완전한 라운드가 포함됨을 보장할 수 있다. 여기에서 $k$는 모든 노드가 평균1번정도는 "눈떨어짐"을 경험할 수 있는 수 이다. (아무리 연결이 적은 노드라도 $k$라운드중 1번은 눈 떨어짐을 평균적으로 경험할것이라 기대함) 이 숫자는 big-drop을 강제발생시키기 위한 장치인데, big-drop이 발생하기 위해서는 여러가지 가능성이 있지만 별다른 증명과정없이 자명하게 알수있는건 눈이 우연히 $i^\star$로 떨어져서 우연히 $j^\star$로 흘러간 경우이다. 실제로는 이 case말고도 big-drop이 발생할 수 있으므로, 논문에서 설정한 $T$는 매우 매우 보수적인 셋팅이다. 

`-` big-drop에 왜케 집착하느냐? big-drop이 일어나면 $\Phi(t)$가 드라마틱하게 변하기 때문이다. 직관적인 설명으로는 "수치예시: $\Phi$가 제곱합이라서 생기는 차이" 를 참고하자. 

`-` 빅드랍이 없는 라운드 ($t_r$이 시작점이라고 하자): 아래를 각각 계산해서 

- $\hbar(X_{t_r}, t_r)=\hbar(\tilde{X}_{0}, t_r)$
- $\hbar(X_{t_r+1}, t_r+1)=\hbar(\tilde{X}_{1}, t_r+1)$
- ... 
- $\hbar(X_{t_r+m_r-1}, t_r+m_r-1)=\hbar(\tilde{X}_{m_r-1}, t_r+m_r-1)$

> 질문: 근데 인덱싱이 안맞는거 아님?? $\mathbb{E}\big[\hbar(X_{t+1},\, t) \mid S(t)\big]$ 이거면 $\hbar(\star, \ast)$ 꼴에서 $\star$와 $\ast$는 한시점차이나는듯한데? 

싹 더하면? 

$$\sum_{j=0}^{m_r-1}\hbar(\tilde{X}_{j}, t_r+j)<m_r \hbar(\tilde{X}_0,t_r)+C_{round}$$
한번 정리하면 아래와 같이 쓸 수 있다. 

$$\sum_{j=0}^{m_r-1}\left(\hbar(\tilde{X}_{j}, t_r+j) - \hbar(\tilde{X}_0,t_r)\right)<C_{round}$$


그런데, 그 라운드에 만약에 큰 낙차가 포함되어있다?? 그러면 아래와 같이 됨. 

$$\sum_{j=0}^{m_r-1}\left(\hbar(\tilde{X}_{j}, t_r+j)-\hbar(\tilde{X}_{0}, t_r)\right)\leq -(m_r-1)\frac{M(t_r)}{n-1} +C_{round}$$
이때 $m_r \geq 2$ 이므로 

$$\sum_{j=0}^{m_r-1}\left(\hbar(\tilde{X}_{j}, t_r+j)-\hbar(\tilde{X}_{0}, t_r)\right)\leq -\frac{M(t_r)}{n-1} +C_{round}$$

`-` 이제 아래에 관심을 가지자. 

$$2b \sum_{s=0}^{T-1} \mathbb{E}\big[\hbar(X_{t{+}s{+}1},\, t{+}s) \mid S(t)\big] + T \cdot \tfrac{(n-1)b^2}{n}$$
위의식은 아래와 같이 바꿀 수 있다. 

$$2b \sum_{r=1}^{R}\sum_{j=0}^{m_r-1} \mathbb{E}\big[\hbar(X_{t{+}s{+}1},\, t{+}s) \mid S(t)\big] + T \cdot \tfrac{(n-1)b^2}{n}$$
여기에서 $R$은 $T$에 들어있는 라운드의 수이다. (확률변수) 그런데 $T$시간동안 최소 $k$개의 라운드가 있음은 무조건 보장되므로 


`-` T시간 동안 최소한 완전한 $k$개의 라운드가 있다. 이 $k$개의 라운드에서 몇개의 라운드에서나 big-drop이 일어날까? 그 숫자는 $\sum_{r=1}^{k}{\bf 1}_{E_r}$로 표현할 수 있을 것이다. 그러면 이제 아래에 관심을 가져보자. 

$$\sum_{r=1}^{k}{\bf 1}_{E_r}\sum_{j=0}^{m_r-1}\left(\hbar(X_{t_r+j}, t_r+j)-\hbar(\tilde{X}_{0}, t_r)\right)$$

그 빅드랍이 일어나는 여러가지 경우는 차치하고,  $i^\star$에 떨어져서 바로 $j^\star$로 가는 경우만 잡아도 라운드별로 $p_0$ 이다. $k$의 라운드가 최소한 완전하게 있으니까 기대값의 센스에서 $kp_0$번정도 그러한 라운드가 포함될 것이라 기대할 수 있다. 따라서 
평균적인 센스에서 저값을 한번정리해주면 

$$\mathbb{E}\left[\sum_{r=1}^{k}{\bf 1}_{E_r}\sum_{j=0}^{m_r-1}\left(\hbar(\tilde{X}_{j}, t_r+j)-\hbar(\tilde{X}_{0}, t_r)\right) \right]\leq kp_0 \left[ -\frac{M(t)}{n-1}+\frac{Tb}{n-1}+C_{round}\right]$$

이제 $\sum_{r=1}^{k}\sum_{j=0}^{m_r-1} \hbar(\tilde{X}_0, t_r)=\sum_{r=1}^{k}m_r \hbar(\tilde{X}_0,t_r)$ 에 관심을 가지자. 기대값을 취하면 아래와 같다. 

$$\sum_{r=1}^km_r\mathbb{E}[\hbar(\tilde{X}_0,t_r)] \leq T\frac{\alpha M(t)}{2}+\frac{\alpha b T^2}{2}+ T C_{round}$$

