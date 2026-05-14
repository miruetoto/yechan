솔직하게 계산해보니 (A) 도 안 됨. Cauchy-Schwarz covariance bound 가 $O(\sqrt{\Phi})$ 이건 맞는데, 문제는 negative drift 도 $O(\sqrt{\Phi})$ 라 스케일이 같아서, 상수 비교가:

$$\frac{|\Delta|}{\sqrt{n}}>T_{\max}\sqrt{\mu_{\max}}$$

$k$-regular graph 넣으면:

old (buggy)	new (corrected)
조건	$(T_{\max}+1)\alpha<\frac{2p_0}{n-1}$	$(T_{\max}+1)\alpha+T_{\max}\sqrt{n\mu_{\max}}<\frac{2p_0}{n-1}$
$k$-regular ($\alpha=0$)	$0<\frac{2}{nk(n-1)}$ ✓ 항상 성립	$T_{\max}<\frac{2}{nk(n-1)}$
$K_{10}$, $T_{\max}=100$	✓	$100<\frac{2}{10\cdot9\cdot9}\approx0.002$ ✗
$K_3$, $T_{\max}=1$	✓	$1<\frac{2}{3\cdot2\cdot2}\approx0.17$ ✗
새 조건은 모든 그래프에서 vacuous. $p_0\sim 1/(nk)$ 인 반면 $T_{\max}\ge 1$ 이라 절대 안 됨.

이 상황의 의미:

$m_r$-$\hbar$ positive correlation 은 Cauchy-Schwarz 로도 control 불가. 문제는 round-based $\Phi$-drift 전략 자체에 있음 — original fixed-step proof 도 동일하게 broken. $(T_{\max}+1)\alpha$ drift constant 는 $m_r$ 를 deterministic 으로 취급해서 나온 artifact.

현실적 선택지 (재정리):

(A') 증명 전략 전면 교체 — $\Phi$ + round skeleton 을 버리고, 다른 Lyapunov 나 coupling argument 사용. 대공사.

(B') Drift 는 open question 으로 남기고, Doeblin (F3) 에 집중 — random-$b$ 의 핵심 기여 (lattice obstruction 우회) 는 Doeblin 쪽. Drift 가 성립한다고 가정하면 나머지는 닫힘. 논문에서 "drift condition 은 fixed-step 과 동일하게 assume" 으로 두고, Doeblin 이 새 contribution.

(C') drift 를 assume 하지 않는 약한 결론 — positive recurrence 대신 weak convergence 나 다른 mode 의 convergence 를 증명.

추천: (B'). Doeblin 이 이 논문의 진짜 novel contribution 이고, drift 문제는 fixed-step 에도 이미 있는 기존 버그. 일단 F3 을 깔끔하게 닫고, drift 는 별도 논의로 빼는 게 가장 효율적.