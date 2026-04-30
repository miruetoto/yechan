---
title: 연구 ▷ HST ▷ Appendix A 약한 에르고딕성 — Hajnal보다 reset coupling으로 쓰자
author: codex
date: 04/29/2026
draft: false
output-file: 260429_hst_appendixA_reset_coupling.html
---

# 0. 결론

HST 논문 Appendix A의 weak ergodicity 증명은 방향은 맞다. 핵심 직관은 다음이다.

> HST는 flow가 계속되더라도 $T_{\max}$번 후에는 강제로 reset되고, block이 생겨도 reset된다. reset 직후 다음 위치는 $\mu_0$에서 새로 뽑힌다. 따라서 어떤 시점들에서는 transition matrix의 모든 row가 같은 $\mu_0$가 된다.

이 아이디어는 좋다. 다만 현재 Appendix의 Hajnal/Dobrushin식 서술은 엄밀하게 쓰려면 손볼 부분이 있다. 더 깔끔한 방향은 **Hajnal theorem을 길게 쓰는 것보다 reset step의 Doeblin/coupling 구조를 직접 쓰는 것**이다.

---

# 1. Appendix A의 좋은 점

Appendix A가 잡은 핵심은 정확하다.

HST의 snow trace는 일반 random walk와 다르다. 일반 directed random walk에서는 graph가 reducible이면 한 component에 갇힐 수 있다. 하지만 HST는 막히거나 $T_{\max}$번 흐르면 다음 round가 시작되고, 그때 위치가

$$
X_{t+1}\sim\mu_0
$$

로 다시 뽑힌다. 즉 path가 어디에 있었는지와 무관한 global refresh가 반복적으로 들어온다.

따라서 $X_t$의 초기 위치 의존성은 reset step에서 끊긴다. 이건 weak ergodicity를 얻기에 매우 강한 구조다.

---

# 2. 현재 Appendix 증명의 약한 부분

## 2.1 $X_t$ 단독은 Markov chain이 아니다

HST의 transition은 $X_t$만으로 정해지지 않는다. 다음 이동은

$$
h(\cdot,t),\qquad X_t,
\qquad Z_t
$$

에 의존한다. 따라서 $X_t$ 단독은 homogeneous Markov chain이 아니다. 엄밀한 상태는 예컨대

$$
S_t=(\boldsymbol\delta_t,X_t,Z_t)
$$

또는

$$
S_t=(h_t,X_t,Z_t)
$$

이다.

Appendix는 $P(t)$를 $X_t$의 non-homogeneous transition matrix처럼 쓰는데, 이 matrix는 deterministic time sequence라기보다 history-dependent/random environment에 가깝다. 이 부분은 명시해야 한다.

## 2.2 Hajnal theorem 적용이 너무 빠르다

Hajnal/Dobrushin theorem은 보통 deterministic sequence of stochastic matrices

$$
P(0),P(1),P(2),\dots
$$

에 대해 서술된다. 그런데 HST에서 $P(t)$는 $h_t,Z_t$에 의존하고, 이는 sample path에 따라 달라진다.

따라서 이 노선을 유지하려면 다음처럼 써야 한다.

> Sample path 또는 filtration을 고정하고, conditional transition matrices의 product에 대해 Dobrushin coefficient를 평가한다.

이 문장이 없으면 theorem 적용이 약간 공중에 뜬다.

## 2.3 stopping time/indexing이 모호하다

Appendix의 Theorem A.5는 block/reset 시점 $M_k$를 잡고

$$
\delta(P(M_k,M_{k+1}))\le \delta(P(M_k))=0
$$

같은 식으로 닫는다. 이 결론은 $M_k$가 정확히 “다음 step이 $\mu_0$ fresh draw인 flag time”이면 맞다. 그때 transition matrix의 모든 row가 $\mu_0$라서 Dobrushin coefficient가 0이다.

하지만 현재 stopping time 정의와 indexing은 읽기에 모호하다. $M_k$가 block이 일어난 시간인지, block 직후 flag time인지, fresh draw 직전인지가 더 명확해야 한다.

---

# 3. 더 깔끔한 대체 방향 — reset coupling

Hajnal을 길게 쓰지 않고도 핵심을 바로 쓸 수 있다.

## 3.1 Flag time 정의

다음처럼 flag time을 정의한다.

$$
\tau_0:=0,
\qquad
\tau_{r+1}:=\inf\{t>\tau_r: Z_t=T_{\max}\}.
$$

또는 논문 notation에 맞춰 “다음 step이 fresh draw인 시점”으로 잡아도 된다. 핵심은

$$
Z_{\tau_r}=T_{\max}
$$

이면 다음 step에서

$$
X_{\tau_r+1}\sim\mu_0
$$

라는 점이다.

## 3.2 bounded return to flag

Algorithm 1의 구조상, 일단 어떤 시점에서든 다음 flag time까지의 시간은 deterministic하게 bounded된다.

$$
\tau_{r+1}-\tau_r\le T_{\max}+1.
$$

이유는 간단하다. flow가 계속되면 $Z$가 하나씩 증가하고, $T_{\max}$에 도달하면 flag가 된다. 중간에 block이 생기면 즉시 $Z=T_{\max}$가 된다.

## 3.3 Coupling

두 개의 HST trace를 생각하자. 서로 다른 초기 상태

$$
S_0,
\qquad S_0'
$$

에서 출발하더라도, 각각 flag time에 도달하면 다음 위치는 모두 $\mu_0$에서 뽑힌다. 이때 같은 uniform random variable을 써서 fresh draw를 coupling하면

$$
X_{\tau+1}=X'_{\tau'+1}
$$

이 되게 만들 수 있다.

더 중요한 점은, fresh draw의 law가 과거 상태와 무관하다는 것이다. 따라서 reset step은 initial $X_0$ 정보를 완전히 지운다. Dobrushin 관점에서는 reset transition matrix가

$$
Q(i,\cdot)=\mu_0(\cdot)
\qquad\text{for all }i
$$

이므로

$$
\delta(Q)=0.
$$

그런 reset kernel이 bounded gaps로 반복 등장하므로 weak ergodicity가 따라온다.

---

# 4. 논문용 정리 형태

아래는 Appendix A를 대체할 수 있는 더 짧은 정리 형태다.

::: {.callout-note title="Proposition: reset contraction of the snow trace"}
Let $\{S_t\}$ be the augmented HST process with state
$$
S_t=(\boldsymbol\delta_t,X_t,Z_t).
$$
Assume that $T_{\max}<\infty$. Let
$$
\mathcal F:=\{S:Z=T_{\max}\}
$$
be the block-flag set. If $S_t\in\mathcal F$, then the next trace location satisfies
$$
\mathcal L(X_{t+1}\mid S_t)=\mu_0.
$$
Equivalently, the conditional transition kernel of $X_{t+1}$ from any two flag states has Dobrushin coefficient zero. Moreover, every sample path returns to $\mathcal F$ within at most $T_{\max}+1$ steps. Hence the trace process has repeated reset contractions with uniformly bounded gaps.
:::

**Proof.** If $S_t\in\mathcal F$, then $Z_t=T_{\max}$. By the first branch of Algorithm 1, $X_{t+1}$ is drawn from $\mu_0$, independently of the previous trace location and height profile. Thus the transition kernel from any flag state has identical rows equal to $\mu_0$, and its Dobrushin coefficient is zero.

It remains to note that flag states recur with bounded gaps. If the process does not block, then each flow step increases $Z$ by one. After at most $T_{\max}$ consecutive flow steps, $Z$ reaches $T_{\max}$. If the process blocks earlier, then $Z$ is set to $T_{\max}$ immediately. Therefore the next visit to $\mathcal F$ occurs within at most $T_{\max}+1$ steps. $\square$

---

# 5. 해석

이 명제는 Appendix A의 약한 에르고딕성 논의를 더 HST답게 만든다.

- Hajnal theorem을 쓰려면 random/history-dependent transition matrix 문제를 정리해야 한다.
- Reset coupling을 쓰면 핵심 구조가 바로 보인다.
- HST는 directed graph에서 path가 막혀도 reset으로 전역 재시작한다.
- 따라서 weak forgetting은 graph connectivity보다 reset mechanism에서 나온다.

다만 이 명제만으로 $\overline{\mathrm{SD}}^2(t)$ 수렴이 끝나는 것은 아니다. 이것은 trace $X_t$의 초기 위치 망각에 관한 statement이고, height difference $\boldsymbol\delta_t$의 positive recurrence나 moment bound는 별도 Foster--Lyapunov 논증이 필요하다.


---

# 6. Hajnal은 꼭 필요한가?

결론부터 말하면, Appendix A의 현재 목적만 놓고 보면 Hajnal은 거의 필요 없다.

Hajnal/Dobrushin 이론이 유용한 전형적인 상황은 다음과 같다.

$$
P_1P_2P_3\cdots
$$

처럼 매 step의 transition matrix가 다르고, 각각의 block이 조금씩만 contraction을 줄 때, 그 product 전체가 초기분포를 잊는다는 것을 보이고 싶을 때다.

그런데 HST는 더 강한 구조를 가진다. reset step에서는 transition kernel이 아예

$$
Q(i,\cdot)=\mu_0(\cdot)\qquad\text{for all }i
$$

가 된다. 즉 모든 row가 같다. 따라서 Dobrushin coefficient는

$$
\delta(Q)=0.
$$

이건 “조금씩 섞는다”가 아니라, **한 번에 trace 위치의 초기 정보를 지운다**는 뜻이다.

비유하면 다음과 같다.

| 방식 | 의미 |
|---|---|
| Hajnal | 매 step 또는 block마다 조금씩 섞어서 결국 잊는다 |
| HST reset | 중간에 memory를 아예 reset한다 |

따라서 Appendix A는 Hajnal theorem을 호출하지 않고도 다음 네 줄로 충분히 정리될 수 있다.

1. $Z=T_{\max}$이면 다음 $X$는 $\mu_0$에서 뽑힌다.
2. 이 reset kernel은 모든 row가 같으므로 Dobrushin coefficient가 0이다.
3. 이런 reset step은 최대 $T_{\max}+1$ 간격으로 반복된다.
4. 따라서 trace $X_t$의 초기 위치 의존성은 reset step에서 사라진다.

다만 이 결론은 어디까지나 **trace $X_t$의 weak forgetting**에 관한 것이다. 논문의 진짜 어려운 부분인

$$
\overline{\mathrm{SD}}^2(t)
$$

수렴, 즉 height difference $\boldsymbol\delta_t$의 positive recurrence와 moment bound는 여전히 별도 문제다. 거기에는 Hajnal보다 Foster--Lyapunov drift가 핵심이다.

---

# 7. 최종 메모

Appendix A는 방향을 잘 잡았다. 하지만 논문용으로는 “Hajnal을 적용한다”보다 “reset kernel의 Dobrushin coefficient가 0이고, 그런 reset이 bounded gaps로 반복된다”라고 쓰는 편이 더 짧고 강하다.
