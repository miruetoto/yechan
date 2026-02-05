---
title: (책공부) [SOGT] General Topology
author: 신록예찬
date: 08/31/2025
draft: false
---

## Chap 5. Topological Spaces: Definitions

### TOPOLOGICAL SPACES

$X$를 공집합이 아닌 집합이라 하자. $X$의 부분집합들의 모임 $\mathcal{T}$가 $X$ 위의 위상이란, $\mathcal{T}$가 다음 공리들을 만족하는 것이다.

- [$O_1$] $X$ and $\varnothing$ belong to $\mathcal{T}$.
- [$O_2$] The union of any number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.
- [$O_3$] The intersection of any two sets in $\mathcal{T}$ belongs to $\mathcal{T}$.

$\mathcal{T}$의 원소들을 $\mathcal{T}$-열린 집합, 또는 간단히 열린 집합이라 부르며, $X$와 $\mathcal{T}$를 함께 쌍 $(X, \mathcal{T})$를 위상 공간이라 한다.

> 규빈: 위상이란, $X$의 부분집합 중 "모든 멤버들이 여유있는 집합들의 모임"이라 볼 수 있다. 여기에서 여유있다라는 개념이 좀 추상적인 편인데, 내가 이해하는 방식으로 설명하보겠다. 아래와 같은 상상을 하자. 
> 1. $X$의 임의의 부분집합 $A$를 생각하자. 
> 2. $A$의 원소들은 각각 스스로 "내가 여유있나?" 라고 반문한다. 
> 3. 내가 여유있는지 판단할때, 외부원소 없이 스스로 독립적으로 여유있다고 느낄수도 있고, 아니면 어떠한 원소의 존재때문에 "쟤랑 같이 있음 내가 여유있지" 라고 느낄수도 있다. 
> 4. 아무튼 어떠한 판단과정을 거쳐서 $A$의 모든 원소들이 "난 여유있어" 라고 생각한다면? 그 집합 $A$는 여유있는 집합이다. 
> 예를들어서 $X=\{a,b,c\}$ 라고 하자. 원소 $a$는 $b$가 있다면 여유있다고 느끼고, 원소 $b$는 본인스스로 여유있다고 느낀다. $c$는 $b$가 있어야 여유있다고 느낀다. 그렇다면 
> 1. $\{b\}$는 여유있는 집합이다. 왜냐하면 원소 $b$는 스스로 여유있다고 생각하니까. 
> 2. $\{a,b\}$ 역시 여유있는 집합이다. 왜냐하면 원소 $a$는 $b$가 있다면 여유있다고 생각하고, $b$는 그 스스로 여유있다고 생각하니까. 
> 3. $\{b,c\}$ 역시 여유있는 집합이다. 왜냐하면 원소 $b$는 스스로 여유있다고 생각하고, $c$는 $b$가 있다면 여유있다고 생각하니까. 
> 4. $\{a,b,c\}$ 역시 여유있는 집합이다. 왜냐하면 $a$는 $b$때문에 여유있다고 생각하고 $b$는 스스로 여유있다고 생각하고 $c$는 $b$때문에 여유있다고 생각하니까. 
> 5. $\emptyset$역시 여유있는집합이다. 왜냐하면 여유있음을 따질수 있는 주체가 없으니까. (그냥 여유있다 생각한다고 치는거) 
> 6. $\{a\}$는 여유있는 집합이 아니다. 왜냐하면 $a$는 스스로 여유있다고 생각하지 못하니까. 그리고 $\{c\}$역시 여유있는집합이 아니다. 왜냐하면 $c$ 역시 스스로 여유있다 생각하지 못하기 때문이다. 또한 $\{a,c\}$역시 여유있는 집합이 아니다. 
> 엄청 중요한 직관중 하나가 나를 여유있게 만들러주는 존재이다. $a$의 입장에서는 $b$만 있으면 숨통이 트이므로 "$b$만 있으면 괜찮아, 난 $b$랑 가까우니까" 라고 생각할 수 있다. 또한 $c$역시 $b$를 가깝게 생각할 것이다. 이렇게 "여유"라는 개념은 "가까움"이라는 개념과 추상적으로 닿아있다. 토폴로지는 사실 거리없이 가까움의 개념을 정의할 수 있는 도구라 해석할 수 있는데, 그 주된 테크닉이 "여유있는 집합들의 모임 $\to$ 가까움" 을 정의하는 방식이다. 이제 공리를 하나씩 따져보자. 
> 1. O1을 따져보자. 공집합은 항상 여유있는 집합이라 치자고 했다. $X$ 역시 비슷한 맥락인데 전체집합은 모든 원소들을 모은 집합인데 그래도 모든 원소가 모여있으면 각각의 원소가 "여유있어" 라고 생각한다고 가정하자. (만약에 그렇지 않다면 논의할게 없으니까.. 다 모여있어도 여유있다고 생각하지 않는데, 거기서 뭘 여유있는 집합의 모임을 만들어) 
> 2. O2를 따져보자. 이 성질은 합집합에 대한 성질인데, 이건 생각보다 당연하다. $A$모임의 각 원소가 스스로 여유있다고 느끼고, $B$모임의 각 원소가 스스로 여유있다고 느낀다고 가정하자. 그러면 당연히 $A\cup B$의 각 원소도 모두 스스로 여유있다고 여길것이다. 왜냐하면 최소한 $A\cup B$ 에서 $A$원소 출신은 $A\cup B$에 있는 $A$때문에 스스로 여유있다 여길것이고, $B$원소 출신역시 $A\cup B$에 있는 $B$때문에 스스로 여유있다고 여길테니까. 
> 3. O3조건을 따져보자. 이 조건이 킥이다. $A$모임의 각 원소가 스스로 여유있다고 여기고, $B$모임의 각 원소가 스스로 여유있다고 여긴다고 가정하자. $A\cap B$의 원소는 과거에 $A$에도, $B$에도 속해있던 원소이다. 이 원소를 편의상 $p$라고 하자. $p$는 $A\cap B$에서 여유있다고 느낄까? 그렇지 않을 수 있다. 예를들어 그 $p$가 $A$에서 느낀 여유있음은 $A,B$ 모두에 있는 원소때문일 수도 있지만 (즉 $A\cap B$에 있는 원소때문일수도 있지만) $A$에만 있는 원소때문일 수도 있다 (즉 $A-B$에 있는 원소때문일 수도 있다). 마찬가지로 그 $p$가 $B$에서 느낀 여유있음은 $A,B$ 모두에 있는 원소떄문일 수도 있지만 (즉 $A \cap B$에 있는 원소때문일수도 있지만) $B$에만 있는 원소때문일수도 있다 (즉 $B-A$에 있는 원소 때문일수도 있다). 결과적으로 $p$가 만약 (1) $A$에만 있는 원소때문에 $A$에서 편안함을 느꼈고, (2) $B$에만 있는 원소 때문에 $B$에서 편안함을 느꼈다면, 그 $p$와 가까웠던 두 원소가 교집합에서는 동시에 날아가는 상황이 생기는 것이다. 따라서 이러한 경우에는 $p$가 $A\cap B$에서 여유있음을 느낄 이유가 없다. 따라서 O3이 만족하는건 직관적이지 않을 수 있다. 여기서 하나의 센스가 필요한데 바로 "교집합은 가까운 원소를 남긴다" 이다. 즉 교집합으로 날아가는 원소는 가까운 원소가 아니라는 의미이다. 이 센스를 적용하면 $p$와 가까웠던 원소들이 교집합에 모두 날아가는 상황이 불가능함을 이해할 수 있다. 
> 요약하면 O1-O3는 여유있음(내부)과 가까움(이웃)의 개념을 동시에 정의하는 매우 영리한 원칙인 셈이다. 

Example 1.1: $\mathcal{U}$를 4장에서 다룬 실수의 모든 열린 집합들의 모임이라 하자. 그러면 $\mathcal{U}$는 $\mathbb{R}$ 위의 위상이며, 이를 $\mathbb{R}$ 위의 보통위상이라 한다. 마찬가지로 평면 $\mathbb{R}^2$에서 모든 열린 집합들의 모임 $\mathcal{U}$도 위상이며, $\mathbb{R}^2$ 위의 보통위상이라 한다. 별도로 명시하지 않는 한 $\mathbb{R}$과 $\mathbb{R}^2$에서는 항상 보통위상을 가정한다.

Example 1.2: $X = \{a, b, c, d, e\}$의 부분집합들의 다음 모임들을 생각하자.
$$\mathcal{T}_1 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$$
$$\mathcal{T}_2 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d\}\}$$
$$\mathcal{T}_3 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{a, b, d, e\}\}$$

$\mathcal{T}_1$은 세 공리 [$O_1$], [$O_2$], [$O_3$]을 만족하므로 $X$ 위의 위상임을 확인하라. 그러나 $\mathcal{T}_2$는 $X$위의 위상이 아닌데, $\mathcal{T}_2$의 두 원소의 합집합
$$\{a, c, d\} \cup \{b, c, d\} = \{a, b, c, d\}$$
이 $\mathcal{T}_2$에 속하지 않기 때문이다. 즉 $\mathcal{T}_2$는 공리 [$O_2$]를 만족하지 않는다.

또한 $\mathcal{T}_3$도 $X$ 위의 위상이 아닌데, $\mathcal{T}_3$의 두 집합의 교집합
$$\{a, c, d\} \cap \{a, b, d, e\} = \{a, d\}$$
이 $\mathcal{T}_3$에 속하지 않기 때문이다. 즉 $\mathcal{T}_3$은 공리 [$O_3$]을 만족하지 않는다.

> 규빈: 이것도 읽어보자. ${\cal T}_1$은 이런상황인 것이다. 
> 1. $a$는 스스로 여유있다 여긴다. 
> 2. $\{c,d\}$는 서로가 같이 있을때 여유있다 느낀다.
> 3. $\{b,e\}$는 $\{c,d\}$와 같이 있을때 여유있다 느낀다. 
> 그리고 ${\cal T}_2$는 이런상황인 것이다. 
> 1. $a$는 스스로 여유있다 여긴다. 
> 2. $\{c,d\}$는 서로가 같이 있을때 여유있다 느낀다.
> 3. $b$는 $\{c,d\}$와 같이 있을때 여유있다고 느낀다. 
> 4. 어? 그런데 왜 $\{a,b,c,d\}$는 여유있는 집합이 아니야???? 
> 그리고 ${\cal T}_3$은 이런상황인 것이다. 
> 1. $\{a,c,d\}$는 여유있는 집합이다. 
> 2. $\{a,b,d,e\}$는 여유있는 집합이다. 
> 3. 어?? 그럼 $\{a,d\}$는 가까운사이 아니야? 그런데 왜 $\{a,d\}$는 ${\cal T}_3$에 없어??

Example 1.3: $\mathcal{D}$를 $X$의 모든 부분집합들의 모임이라 하자. $\mathcal{D}$가 $X$ 위의 위상에 대한 공리들을 만족함을 확인하라. 이 위상을 이산위상이라 하며, $X$와 그 이산위상을 함께, 즉 쌍 $(X, \mathcal{D})$를 이산위상공간 또는 간단히 이산공간이라 한다.

> 규빈: 이산공간 = 모두 스스로 여유있음. 

Example 1.4: 공리 [$O_1$]에서 보듯이, $X$ 위의 위상은 반드시 집합 $X$와 $\varnothing$를 포함해야 한다. $X$와 $\varnothing$만으로 이루어진 모임 $\mathcal{J} = \{X, \varnothing\}$은 그 자체로 $X$ 위의 위상이다. 이를 비이산위상이라 하며, $X$와 그 비이산위상을 함께, 즉 $(X, \mathcal{J})$를 비이산위상공간 또는 간단히 비이산공간이라 한다.

> 규빈: 비이산공간 = 전체가 다 모인것 아니면 모두 여유없다고 생각함. 

Example 1.5: $\mathcal{T}$를 여집합이 유한인 $X$의 모든 부분집합들의 모임에 공집합 $\varnothing$를 함께 넣은 것이라 하자. 이 모임 $\mathcal{T}$도 $X$ 위의 위상이다. 이를 여유한위상 또는 $X$ 위의 $T_1$-위상이라 한다. ($T_1$의 의미는 이후 장에서 다룬다.)

> 규빈: $X$는 당연히 포함된다. $X^c=\emptyset$ 이고 $|\emptyset|=0<\infty$ 니까. 공집합은 정의에 의하여 포함된다. $A,B$가 포함된다고 하자. 그러면 $|A^c|<\infty$ 이고 $|B^c|<\infty$ 이다. 이제 $(A \cap B)^c = (A^c \cup B^c)$ 를 관찰하자. 우변에서 $A^c \cup B^c$ 의 원소가 유한므로 $A \cap B$ 역시 ${\cal T}$의 원소가 된다. 이 논리는 두개 집합의 intersection이 아니라 finite intersection까지 유효하다. 합집합은 더 따지기쉽다. $A,B$가 포함된다고 하자. 즉 $|A^c|<\infty$ 이고 $|B^c|<\infty$ 이라고 하자.  $(A\cup B)^c=A^c \cap B^c$ 을 관찰하면 당연히 우변의 집합은 finite하므로 $A\cup B$는 ${\cal T}$의 멤버가 된다. 이 논의는 uncountable union에 대해서도 성립한다. 

Example 1.6: $X$ 위의 임의의 두 위상 $\mathcal{T}_1$과 $\mathcal{T}_2$의 교집합 $\mathcal{T}_1 \cap \mathcal{T}_2$도 $X$ 위의 위상이다. [$O_1$]에 의해 $X$와 $\varnothing$는 각각 $\mathcal{T}_1$과 $\mathcal{T}_2$ 모두에 속하므로, $X$와 $\varnothing$는 교집합 $\mathcal{T}_1 \cap \mathcal{T}_2$에 속한다. 즉 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_1$]을 만족한다. 또한 $G, H \in \mathcal{T}_1 \cap \mathcal{T}_2$이면, 특히 $G, H \in \mathcal{T}_1$이고 $G, H \in \mathcal{T}_2$이다. 그런데 $\mathcal{T}_1$과 $\mathcal{T}_2$가 위상이므로 $G \cap H \in \mathcal{T}_1$이고 $G \cap H \in \mathcal{T}_2$이다. 따라서
$$G \cap H \in \mathcal{T}_1 \cap \mathcal{T}_2$$
즉 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_3$]을 만족한다. 마찬가지로 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_2$]를 만족한다.

앞의 예제에서의 진술은 사실 위상들의 임의의 모임으로 일반화할 수 있다. 즉,

Theorem 5.1: Let $\{\mathcal{T}_i : i \in I\}$ be any collection of topologies on a set $X$. Then the intersection $\bigcap_i \mathcal{T}_i$ is also a topology on $X$.

마지막 예제에서, 위상들의 합집합은 위상이 되지 않을 수 있음을 보인다.

Example 1.7: 각 모임
$$\mathcal{T}_1 = \{X, \varnothing, \{a\}\} \quad \text{and} \quad \mathcal{T}_2 = \{X, \varnothing, \{b\}\}$$
은 $X = \{a, b, c\}$ 위의 위상이다. 그러나 합집합
$$\mathcal{T}_1 \cup \mathcal{T}_2 = \{X, \varnothing, \{a\}, \{b\}\}$$
은 $X$ 위의 위상이 아닌데, [$O_2$]를 위반하기 때문이다. 즉 $\{a\} \in \mathcal{T}_1 \cup \mathcal{T}_2$, $\{b\} \in \mathcal{T}_1 \cup \mathcal{T}_2$이지만 $\{a\} \cup \{b\} = \{a, b\}$는 $\mathcal{T}_1 \cup \mathcal{T}_2$에 속하지 않는다.

$G$가 점 $p \in X$를 포함하는 열린 집합이면, $G$를 $p$의 열린 근방이라 한다. 또한 $p$를 제외한 $G$, 즉 $G \setminus \{p\}$를 $p$의 제거된 열린 근방이라 한다.

Remark: 공리 [$O_1$], [$O_2$], [$O_3$]은 다음 두 공리와 동치이다:

- [$O_1^*$] The union of any number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.
- [$O_2^*$] The intersection of any finite number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.

[$O_1^*$]이 $\varnothing$가 $\mathcal{T}$에 속함을 함의하는데, 이는
$$\cup \{G \in \mathcal{T} : G \in \varnothing\} = \varnothing$$
즉 집합들의 공합집합은 공집합이기 때문이다. 또한 [$O_2^*$]는 $X$가 $\mathcal{T}$에 속함을 함의하는데, 이는
$$\cap \{G \in \mathcal{T} : G \in \varnothing\} = X$$
즉 $X$의 부분집합들의 공교집합은 $X$ 자체이기 때문이다.

### ACCUMULATION POINTS

$X$를 위상 공간이라 하자. 점 $p \in X$가 $X$의 부분집합 $A$의 집적점 또는 극한점(군점 또는 도래점이라고도 함)이란, $p$를 포함하는 모든 열린 집합 $G$가 $p$와 다른 $A$의 점을 포함하는 것이다. 즉,
$$G \text{ open}, \ p \in G \quad \text{implies} \quad (G \setminus \{p\}) \cap A \neq \varnothing$$
$A$의 집적점들의 집합을 $A'$로 나타내며, $A$의 도래집합이라 한다.

Example 2.1: 모임 $\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$은 $X = \{a, b, c, d, e\}$ 위의 위상을 정의한다. $X$의 부분집합 $A = \{a, b, c\}$를 생각하자. $b \in X$는 $A$의 극한점인데, $b$를 포함하는 열린 집합들이 $\{b, c, d, e\}$와 $X$이고, 각각 $b$와 다른 $A$의 점, 즉 $c$를 포함하기 때문이다. 반면에 점 $a \in X$는 $A$의 극한점이 아닌데, $a$를 포함하는 열린 집합 $\{a\}$가 $a$와 다른 $A$의 점을 포함하지 않기 때문이다. 마찬가지로 점 $d$와 $e$는 $A$의 극한점이고, 점 $c$는 $A$의 극한점이 아니다. 따라서 $A' = \{b, d, e\}$가 $A$의 도래집합이다.

Example 2.2: $X$를 비이산 공간이라 하자. 즉 $X$와 $\varnothing$만이 $X$의 유일한 열린 부분집합이다. 그러면 $X$는 임의의 점 $p \in X$를 포함하는 유일한 열린 집합이다. 따라서 $p$는 공집합 $\varnothing$와 $p$만으로 이루어진 집합 즉 한원소 집합 $\{p\}$를 제외한 $X$의 모든 부분집합의 집적점이다. 따라서 $X$의 임의의 부분집합 $A$의 도래집합 $A'$은 다음과 같다:
$$A' = \begin{cases} \varnothing & \text{if } A = \varnothing \\ \{p\}^c = X \setminus \{p\} & \text{if } A = \{p\} \\ X & \text{if } A \text{ contains two or more points} \end{cases}$$

직선 $\mathbb{R}$과 평면 $\mathbb{R}^2$의 보통위상에 대해, 위의 집적점의 정의는 4장에서 주어진 것과 동일함을 확인하라.

### CLOSED SETS

$X$를 위상 공간이라 하자. $X$의 부분집합 $A$가 닫힌 집합이란, 그 여집합 $A^c$가 열린 집합인 것이다.

Example 3.1: 모임 $\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$은 $X = \{a, b, c, d, e\}$ 위의 위상을 정의한다. $X$의 닫힌 부분집합들은
$$\varnothing, \ X, \ \{b, c, d, e\}, \ \{a, b, e\}, \ \{b, e\}, \ \{a\}$$
즉 $X$의 열린 부분집합들의 여집합이다. $\{b, c, d, e\}$와 같이 열리면서 동시에 닫힌 $X$의 부분집합이 있고, $\{a, b\}$와 같이 열리지도 않고 닫히지도 않은 $X$의 부분집합도 있음에 주목하라.

Example 3.2: $X$를 이산 위상 공간이라 하자. 즉 $X$의 모든 부분집합이 열려 있다. 그러면 $X$의 모든 부분집합은 그 여집합이 항상 열려 있으므로 닫혀 있기도 하다. 다시 말해 $X$의 모든 부분집합은 열리면서 동시에 닫혀 있다.

공간 $X$의 임의의 부분집합 $A$에 대해 $A^{cc} = A$임을 상기하라. 따라서

Proposition 5.2: In a topological space $X$, a subset $A$ of $X$ is open if and only if its complement is closed.

위상 공간의 공리 [$O_1$], [$O_2$], [$O_3$]과 드 모르간 법칙으로부터 다음을 얻는다

Theorem 5.3: Let $X$ be a topological space. Then the class of closed subsets of $X$ possesses the following properties:

1. $X$ and $\varnothing$ are closed sets.
2. The intersection of any number of closed sets is closed.
3. The union of any two closed sets is closed.

닫힌 집합은 그 극한점을 이용하여 다음과 같이 특성화할 수도 있다:

Theorem 5.4: A subset $A$ of a topological space $X$ is closed if and only if $A$ contains each of its accumulation points.

다시 말해, 집합 $A$가 닫혀 있을 필요충분조건은 $A$의 도래집합 $A'$이 $A$의 부분집합인 것, 즉 $A' \subset A$이다.

### CLOSURE OF A SET

$A$를 위상 공간 $X$의 부분집합이라 하자. $A$의 폐포란,
$$\bar{A} \quad \text{or} \quad A^-$$
로 나타내며, $A$를 포함하는 모든 닫힌 상위집합들의 교집합이다. 다시 말해 $\{F_i : i \in I\}$가 $A$를 포함하는 $X$의 모든 닫힌 부분집합들의 모임이면,
$$\bar{A} = \cap_i F_i.$$
$\bar{A}$는 닫힌 집합들의 교집합이므로 닫힌 집합임을 먼저 확인하라. 더 나아가 $\bar{A}$는 $A$를 포함하는 가장 작은 닫힌 상위집합이다. 즉 $F$가 $A$를 포함하는 닫힌 집합이면
$$A \subset \bar{A} \subset F$$
따라서 집합 $A$가 닫혀 있을 필요충분조건은 $A = \bar{A}$이다. 이 결과들을 형식적으로 진술하면:

Proposition 5.5: Let $\bar{A}$ be the closure of a set $A$. Then: (i) $\bar{A}$ is closed; (ii) if $F$ is a closed superset of $A$, then $A \subset \bar{A} \subset F$; and (iii) $A$ is closed iff $A = \bar{A}$.

Example 4.1: Example 3.1의 $X = \{a, b, c, d, e\}$ 위의 위상 $\mathcal{T}$를 생각하자. 닫힌 부분집합들은
$$\varnothing, \ X, \ \{b, c, d, e\}, \ \{a, b, e\}, \ \{b, e\}, \ \{a\}$$
따라서,
$$\overline{\{b\}} = \{b, e\}, \quad \overline{\{a, c\}} = X, \quad \overline{\{b, d\}} = \{b, c, d, e\}.$$
Example 4.2: $X$를 여유한 위상 공간이라 하자. 즉 유한 집합의 여집합과 $\varnothing$가 열린 집합이다. 그러면 닫힌 집합들은 정확히 $X$의 유한 부분집합들과 $X$ 자체이다. 따라서 $A \subset X$가 유한이면, $A$ 자체가 닫혀 있으므로 그 폐포 $\bar{A}$는 $A$이다. 반면에 $A \subset X$가 무한이면, $X$가 $A$를 포함하는 유일한 닫힌 상위집합이므로 $\bar{A}$는 $X$이다. 더 간결하게, 여유한 공간 $X$의 임의의 부분집합 $A$에 대해,
$$\bar{A} = \begin{cases} A & \text{if } A \text{ is finite} \\ X & \text{if } A \text{ is infinite} \end{cases}$$
집합의 폐포는 그 극한점을 이용하여 다음과 같이 완전히 기술할 수 있다:

Theorem 5.6: Let $A$ be a subset of a topological space $X$. Then the closure of $A$ is the union of $A$ and its set of accumulation points, i.e.,
$$\bar{A} = A \cup A'$$

점 $p \in X$가 $A \subset X$의 폐포점 또는 부착점이란, $p$가 $A$의 폐포에 속하는 것, 즉 $p \in \bar{A}$이다. 위 정리에 의하면, $p \in X$가 $A \subset X$의 폐포점일 필요충분조건은 $p \in A$이거나 $p$가 $A$의 극한점인 것이다.

Example 4.3: 유리수의 집합 $\mathbb{Q}$를 생각하자. 앞서 보았듯이 $\mathbb{R}$의 보통위상에서 모든 실수 $a \in \mathbb{R}$은 $\mathbb{Q}$의 극한점이다. 따라서 $\mathbb{Q}$의 폐포는 실수 전체의 집합, 즉 $\bar{\mathbb{Q}} = \mathbb{R}$이다.

위상 공간 $X$의 부분집합 $A$가 $B \subset X$에서 조밀하다란, $B$가 $A$의 폐포에 포함되는 것, 즉 $B \subset \bar{A}$이다. 특히 $A$가 $X$에서 조밀하다, 또는 $X$의 조밀한 부분집합이란, $\bar{A} = X$인 것이다.

Example 4.4: Example 4.1에서 다음을 확인하라.
$$\overline{\{a, c\}} = X \quad \text{and} \quad \overline{\{b, d\}} = \{b, c, d, e\}$$
여기서 $X = \{a, b, c, d, e\}$이다. 따라서 집합 $\{a, c\}$는 $X$의 조밀한 부분집합이지만 $\{b, d\}$는 그렇지 않다.

Example 4.5: Example 4.3에서 언급했듯이 $\bar{\mathbb{Q}} = \mathbb{R}$이다. 다시 말해 보통위상에서 유리수의 집합 $\mathbb{Q}$는 $\mathbb{R}$에서 조밀하다.

"폐포" 연산자는, $X$의 각 부분집합 $A$에 그 폐포 $\bar{A} \subset X$를 대응시키며, 아래 명제에 나오는 네 가지 성질을 만족하는데, 이를 쿠라토프스키 폐포 공리라 한다. 사실 이 공리들은 이후에 증명하겠지만 $X$ 위의 위상을 정의하는 데 사용될 수 있다.

Proposition 5.7: (i) $\bar{\varnothing} = \varnothing$; (ii) $A \subset \bar{A}$; (iii) $\overline{A \cup B} = \bar{A} \cup \bar{B}$; and (iv) $\bar{\bar{A}} = \bar{A}$.

### INTERIOR, EXTERIOR, BOUNDARY

$A$를 위상 공간 $X$의 부분집합이라 하자. 점 $p \in A$가 $A$의 내점이란, $p$가 $A$에 포함된 어떤 열린 집합 $G$에 속하는 것이다:
$$p \in G \subset A \quad \text{where } G \text{ is open}$$
$A$의 내점들의 집합을
$$\text{int}(A), \quad \mathring{A} \quad \text{or} \quad A^\circ$$
로 나타내며, $A$의 내부라 한다. $A$의 내부는 다음과 같이 특성화할 수도 있다:

Proposition 5.8: The interior of a set $A$ is the union of all open subsets of $A$. Furthermore: (i) $A^\circ$ is open; (ii) $A^\circ$ is the largest open subset of $A$, i.e. if $G$ is an open subset of $A$ then $G \subset A^\circ \subset A$; and (iii) $A$ is open iff $A = A^\circ$.

$A$의 외부란, $\text{ext}(A)$로 쓰며, $A$의 여집합의 내부, 즉 $\text{int}(A^c)$이다. $A$의 경계란, $b(A)$로 쓰며, $A$의 내부에도 외부에도 속하지 않는 점들의 집합이다. 다음으로 내부, 외부, 폐포 사이의 중요한 관계가 이어진다.

Theorem 5.9: Let $A$ be any subset of a topological space $X$. Then the closure of $A$ is the union of the interior and boundary of $A$, i.e. $\bar{A} = A^\circ \cup b(A)$.

Example 5.1: 끝점이 $a$와 $b$인 네 구간 $[a, b]$, $(a, b)$, $(a, b]$, $[a, b)$를 생각하자. 각각의 내부는 열린 구간 $(a, b)$이고, 각각의 경계는 끝점들의 집합, 즉 $\{a, b\}$이다.

Example 5.2: $X = \{a, b, c, d, e\}$ 위의 위상
$$\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$$
과 $X$의 부분집합 $A = \{b, c, d\}$를 생각하자. 점 $c$와 $d$는 각각 $A$의 내점인데,
$$c, d \in \{c, d\} \subset A$$
이고 $\{c, d\}$는 열린 집합이기 때문이다. 점 $b \in A$는 $A$의 내점이 아니다. 따라서 $\text{int}(A) = \{c, d\}$이다. $X$의 점 중 오직 $a$만이 $A$에 대해 외점, 즉 $A$의 여집합 $A^c = \{a, e\}$의 내점이다. 따라서 $\text{int}(A^c) = \{a\}$이다. 따라서 $A$의 경계는 점 $b$와 $e$로 이루어진다. 즉 $b(A) = \{b, e\}$이다.

Example 5.3: 유리수의 집합 $\mathbb{Q}$를 생각하자. $\mathbb{R}$의 모든 열린 부분집합은 유리수와 무리수를 모두 포함하므로, $\mathbb{Q}$의 내점이나 외점은 존재하지 않는다. 따라서 $\text{int}(\mathbb{Q}) = \varnothing$이고 $\text{int}(\mathbb{Q}^c) = \varnothing$이다. 따라서 $\mathbb{Q}$의 경계는 실수 전체의 집합, 즉 $b(\mathbb{Q}) = \mathbb{R}$이다.

위상 공간 $X$의 부분집합 $A$가 $X$에서 어디에서도 조밀하지 않다란, $A$의 폐포의 내부가 공집합인 것, 즉 $\text{int}(\bar{A}) = \varnothing$인 것이다.

Example 5.4: $\mathbb{R}$의 부분집합 $A = \{1, \frac{1}{2}, \frac{1}{3}, \frac{1}{4}, ...\}$를 생각하자. 앞서 언급했듯이 $A$는 정확히 하나의 극한점 $0$을 가진다. 따라서 $\bar{A} = \{0, 1, \frac{1}{2}, \frac{1}{3}, \frac{1}{4}, ...\}$이다. $\bar{A}$는 내점을 갖지 않음을 확인하라. 따라서 $A$는 $\mathbb{R}$에서 어디에서도 조밀하지 않다.

Example 5.5: $A$를 0과 1 사이의 유리점들로 이루어진 집합, 즉 $A = \{x : x \in \mathbb{Q}, 0 < x < 1\}$이라 하자. $A$의 내부는 공집합, 즉 $\text{int}(A) = \varnothing$임을 확인하라. 그러나 $A$는 $\mathbb{R}$에서 어디에서도 조밀하지 않은 것은 아니다. $A$의 폐포가 $[0, 1]$이므로,
$$\text{int}(\bar{A}) = \text{int}([0, 1]) = (0, 1)$$
은 공집합이 아니기 때문이다.

### NEIGHBORHOODS AND NEIGHBORHOOD SYSTEMS

$p$를 위상 공간 $X$의 한 점이라 하자. $X$의 부분집합 $N$이 $p$의 근방이란, $N$이 $p$를 포함하는 어떤 열린 집합 $G$의 상위집합인 것이다:
$$p \in G \subset N \quad \text{where } G \text{ is an open set}$$
다시 말해 "$N$이 점 $p$의 근방이다"라는 관계는 "$p$가 $N$의 내점이다"라는 관계를 의미한다.  $p \in X$의 모든 근방들의 모임을 $\mathcal{N}_p$로 나타내며, $p$의 근방계라 한다.

Example 6.1: $a$를 임의의 실수, 즉 $a \in \mathbb{R}$이라 하자. 중심이 $a$인 각 닫힌 구간 $[a - \delta, a + \delta]$는 $a$를 포함하는 열린 구간 $(a - \delta, a + \delta)$를 포함하므로 $a$의 근방이다. 마찬가지로 $p$가 평면 $\mathbb{R}^2$의 한 점이면, 중심이 $p$인 모든 닫힌 원판 $\{q \in \mathbb{R}^2 : d(p, q) < \delta \neq 0\}$은 중심이 $p$인 열린 원판을 포함하므로 $p$의 근방이다.

임의의 점 $p \in X$의 근방계 $\mathcal{N}_p$에 관한 핵심적인 사실들은 아래 명제에 나오는 네 가지 성질이며, 이를 근방 공리라 한다. 사실 이 공리들은 이후에 언급하겠지만 $X$ 위의 위상을 정의하는 데 사용될 수 있다.

Proposition 5.10: 

1. $\mathcal{N}_p$ is not empty and $p$ belongs to each member of $\mathcal{N}_p$.
2. The intersection of any two members of $\mathcal{N}_p$ belongs to $\mathcal{N}_p$.
3. Every superset of a member of $\mathcal{N}_p$ belongs to $\mathcal{N}_p$.
4. Each member $N \in \mathcal{N}_p$ is a superset of a member $G \in \mathcal{N}_p$ where $G$ is a neighborhood of each of its points, i.e. $G \in \mathcal{N}_g$ for every $g \in G$.

### CONVERGENT SEQUENCES

위상공간 $X$에서 점들의 수열 $\langle a_1, a_2, \ldots \rangle$가 점 $b \in X$로 수렴(converge)한다는 것은, 또는 $b$가 수열 $\langle a_n \rangle$의 극한(limit)이라는 것은,
$$\lim_{n \to \infty} a_n = b, \quad \lim a_n = b \quad \text{or} \quad a_n \to b$$
로 표기하며, $b$를 포함하는 각 열린집합 $G$에 대해 양의 정수 $n_0 \in N$이 존재하여
$$n > n_0 \quad \text{implies} \quad a_n \in G$$
인 것과 동치이다. 즉, $G$가 수열의 항들 중 거의 모든(almost all), 다시 말해 유한 개를 제외한 모든 항을 포함하는 경우이다.

Example 7.1: $\langle a_1, a_2, \ldots \rangle$를 비이산 위상공간 $(X, \mathcal{G})$에서 점들의 수열이라 하자. 다음을 관찰하라: (i) $X$는 임의의 점 $b \in X$를 포함하는 유일한 열린집합이다; 그리고 (ii) $X$는 수열 $\langle a_n \rangle$의 모든 항을 포함한다. 따라서 수열 $\langle a_1, a_2, \ldots \rangle$는 모든 점 $b \in X$로 수렴한다.

Example 7.2: $\langle a_1, a_2, \ldots \rangle$를 이산 위상공간 $(X, \mathcal{D})$에서 점들의 수열이라 하자. 이제 모든 점 $b \in X$에 대해 한원소집합 $\{b\}$는 $b$를 포함하는 열린집합이다. 따라서 $a_n \to b$이면, 집합 $\{b\}$는 수열의 거의 모든 항을 포함해야 한다. 다시 말해, 수열 $\langle a_n \rangle$이 점 $b \in X$로 수렴하는 것은 수열이 $\langle a_1, a_2, \ldots, a_{n_0}, b, b, b, \ldots \rangle$의 형태인 것과 동치이다.

Example 7.3: $\mathcal{T}$를 무한집합 $X$ 위의 위상으로서 $\varnothing$과 가산집합들의 여집합으로 구성된 것이라 하자(Problem 56 참조). $X$에서 수열 $\langle a_1, a_2, \ldots \rangle$가 $b \in X$로 수렴하는 것은 수열이 $\langle a_1, a_2, \ldots, a_{n_0}, b, b, b, \ldots \rangle$의 형태인 것과 동치임을 주장한다. 즉, $\langle a_n \rangle$의 항 중 $b$와 다른 것들로 이루어진 집합 $A$가 유한이다. 이제 $A$는 가산이므로 $A^c$는 $b$를 포함하는 열린집합이다. 따라서 $a_n \to b$이면 $A^c$는 수열의 유한 개를 제외한 모든 항을 포함하고, 그러므로 $A$는 유한이다.

### COARSER AND FINER TOPOLOGIES

$\mathcal{T}_1$과 $\mathcal{T}_2$를 공집합이 아닌 집합 $X$ 위의 위상이라 하자. 각 $\mathcal{T}_1$-열린 부분집합이 또한 $X$의 $\mathcal{T}_2$-열린 부분집합이라 가정하자. 즉, $\mathcal{T}_1$이 $\mathcal{T}_2$의 부분족(subclass)이라 가정하자, 다시 말해 $\mathcal{T}_1 \subset \mathcal{T}_2$이다. 그러면 $\mathcal{T}_1$이 $\mathcal{T}_2$보다 더 거친(coarser) 또는 더 작은(smaller) (때로는 더 약한(weaker)이라고도 하는) 위상이라 하고, 또는 $\mathcal{T}_2$가 $\mathcal{T}_1$보다 더 세밀한(finer) 또는 더 큰(larger) 위상이라 한다. $X$ 위의 모든 위상들의 모임 $\mathbf{T} = \{\mathcal{T}_i\}$는 족 포함관계에 의해 반순서(partially ordered)됨을 관찰하라; 따라서 다음과 같이도 쓴다

$$\mathcal{T}_1 \precsim \mathcal{T}_2 \quad \text{for} \quad \mathcal{T}_1 \subset \mathcal{T}_2$$

그리고 $X$ 위의 두 위상이 어느 쪽도 다른 쪽보다 거칠지 않으면 비교불가능(not comparable)하다고 말한다.

Example 8.1: 임의의 집합 $X$ 위의 이산 위상 $\mathcal{D}$, 비이산 위상 $\mathcal{G}$, 그리고 임의의 다른 위상 $\mathcal{T}$를 생각하자. 그러면 $\mathcal{T}$는 $\mathcal{D}$보다 거칠고 $\mathcal{T}$는 $\mathcal{G}$보다 세밀하다. 즉, $\mathcal{G} \precsim \mathcal{T} \precsim \mathcal{D}$이다.

Example 8.2: 평면 $\mathbf{R}^2$ 위의 여유한 위상 $\mathcal{T}$와 보통 위상 $\mathcal{U}$를 생각하자. $\mathbf{R}^2$의 모든 유한 부분집합은 $\mathcal{U}$-닫힌집합임을 상기하라; 따라서 $\mathbf{R}^2$의 유한 부분집합의 여집합, 즉 $\mathcal{T}$의 임의의 원소는 또한 $\mathcal{U}$-열린집합이다. 다시 말해, $\mathcal{T}$는 $\mathcal{U}$보다 거칠다, 즉 $\mathcal{T} \precsim \mathcal{U}$이다.

### SUBSPACES, RELATIVE TOPOLOGIES

$A$를 위상공간 $(X, \mathcal{T})$의 공집합이 아닌 부분집합이라 하자. $A$와 $X$의 $\mathcal{T}$-열린 부분집합들의 모든 교집합으로 이루어진 족 $\mathcal{T}_A$는 $A$ 위의 위상이다; 이것을 $A$ 위의 상대 위상(relative topology) 또는 $\mathcal{T}$의 $A$로의 상대화(relativization)라 부르며, 위상공간 $(A, \mathcal{T}_A)$를 $(X, \mathcal{T})$의 부분공간(subspace)이라 부른다. 다시 말해, $A$의 부분집합 $H$가 $\mathcal{T}_A$-열린집합, 즉 $A$에 상대적으로 열린집합인 것은 $X$의 $\mathcal{T}$-열린 부분집합 $G$가 존재하여

$$H = G \cap A$$

인 것과 동치이다.

Example 9.1: 다음 위상을 생각하자
$$\mathcal{T} = \{X, \varnothing, \{a\}, \{c,d\}, \{a,c,d\}, \{b,c,d,e\}\}$$
$X = \{a,b,c,d,e\}$ 위의 위상이며, $X$의 부분집합 $A = \{a,d,e\}$를 생각하자. 다음을 관찰하라

$$X \cap A = A, \quad \{a\} \cap A = \{a\}, \quad \{a,c,d\} \cap A = \{a,d\}$$
$$\varnothing \cap A = \varnothing, \quad \{c,d\} \cap A = \{d\}, \quad \{b,c,d,e\} \cap A = \{d,e\}$$
따라서 $\mathcal{T}$의 $A$로의 상대화는
$$\mathcal{T}_A = \{A, \varnothing, \{a\}, \{d\}, \{a,d\}, \{d,e\}\}$$

Example 9.2: $\mathbf{R}$ 위의 보통 위상 $\mathcal{U}$와 닫힌 구간 $A = [3, 8]$ 위의 상대 위상 $\mathcal{T}_A$를 생각하자. 반닫힌-반열린 구간 $[3, 5)$는 $A$ 위의 상대 위상에서 열린집합, 즉 $\mathcal{T}_A$-열린집합임에 주목하라, 왜냐하면

$$[3, 5) = (2, 5) \cap A$$

이고 여기서 $(2, 5)$는 $\mathbf{R}$의 $\mathcal{T}$-열린 부분집합이기 때문이다. 이로써 어떤 집합이 부분공간에 상대적으로는 열린집합이지만 전체 공간에서는 열린집합도 닫힌집합도 아닐 수 있음을 알 수 있다.

### EQUIVALENT DEFINITIONS OF TOPOLOGIES

위상공간에 대한 우리의 정의는 위상공간에서 열린집합에 대한 공리를 제시하였다, 즉 열린집합을 위상의 원시적 개념(primitive notion)으로 사용하였다. 이제 집합 위에 위상을 정의하는 대안적 방법을 보여주는 두 정리를 제시하는데, "점의 근방(neighborhood of a point)"과 "집합의 폐포(closure of a set)"의 개념을 원시적 개념으로 사용한다.


Theorem 5.11: Let $X$ be a non-empty set and let there be assigned to each point $p \in X$ a class $\mathcal{A}_p$ of subsets of $X$ satisfying the following axioms:

- [$\mathbf{A_1}$] $\mathcal{A}_p$ is not empty and $p$ belongs to each member of $\mathcal{A}_p$.
- [$\mathbf{A_2}$] The intersection of any two members of $\mathcal{A}_p$ belongs to $\mathcal{A}_p$.
- [$\mathbf{A_3}$] Every superset of a member of $\mathcal{A}_p$ belongs to $\mathcal{A}_p$.
- [$\mathbf{A_4}$] Each member $N \in \mathcal{A}_p$ is a superset of a member $G \in \mathcal{A}_p$ such that $G \in \mathcal{A}_g$ for every $g \in G$.

Then there exists one and only one topology $\mathcal{T}$ on $X$ such that $\mathcal{A}_p$ is the $\mathcal{T}$-neighborhood system of the point $p \in X$.

> 규빈: 이 정리는 적당한 조건 A1-A4를 만족하는 $X$의 모임 ${\cal A}_p$를 점 $p$에 대한 근방모임(neighborhood system of point $p$)을 "위상없이" 정의할 수 있으며, ${\cal A}_p$를 이용하여 위상 ${\cal T}$를 유일하게 결정할 수 있다는 의미이다. 

Theorem 5.12: Let $X$ be a non-empty set and let $k$ be an operation which assigns to each subset $A$ of $X$ the subset $A^k$ of $X$, satisfying the following axioms, called the Kuratowski Closure Axioms:

- [$\mathbf{K_1}$] $\varnothing^k = \varnothing$
- [$\mathbf{K_2}$] $A \subset A^k$
- [$\mathbf{K_3}$] $(A \cup B)^k = A^k \cup B^k$
- [$\mathbf{K_4}$] $(A^k)^k = A^k$

Then there exists one and only one topology $\mathcal{T}$ on $X$ such that $A^k$ will be the $\mathcal{T}$-closure of the subset $A$ of $X$.

> 규빈: 이 정리는 Theorem 5.11과 같은 맥락으로, $A$의 closure $A^k$를 위상없이 정의할 수 있으며, $A^k$를 이용하여 위상 ${\cal T}$를 유일하게 결정할 수 있다는 의미이다. 

## Chap 6. Bases and Subbases

### BASE FOR A TOPOLOGY

$(X, \mathcal{T})$를 위상공간이라 하자. $X$의 열린 부분집합들의 족 $\mathcal{B}$, 즉 $\mathcal{B} \subset \mathcal{T}$가 위상 $\mathcal{T}$에 대한 기저(base)라 함은

(i) 모든 열린집합 $G \in \mathcal{T}$가 $\mathcal{B}$의 원소들의 합집합인 것이다.

동치로, $\mathcal{B} \subset \mathcal{T}$가 $\mathcal{T}$에 대한 기저인 것은

(ii) 열린집합 $G$에 속하는 임의의 점 $p$에 대해, $p \in B \subset G$인 $B \in \mathcal{B}$가 존재하는 것이다.

Example 1.1: 열린구간들은 직선 $\mathbf{R}$ 위의 보통 위상에 대한 기저를 이룬다. 왜냐하면 $G \subset \mathbf{R}$이 열린집합이고 $p \in G$이면, 정의에 의해 $p \in (a, b) \subset G$인 열린구간 $(a, b)$가 존재하기 때문이다. 마찬가지로, 열린 원판들은 평면 $\mathbf{R}^2$ 위의 보통 위상에 대한 기저를 이룬다.

Example 1.2: 평면 $\mathbf{R}^2$에서 $x$-축과 $y$-축에 평행한 변으로 둘러싸인 열린 직사각형들도 $\mathbf{R}^2$ 위의 보통 위상에 대한 기저 $\mathcal{B}$를 이룬다. $G \subset \mathbf{R}^2$를 열린집합이라 하고 $p \in G$라 하자. 그러면 $p$를 중심으로 하고  $p \in D_p \subset G$를 만족하는 열린 원판 $D_p$가 존재한다. 그러면 꼭짓점이 $D_p$의 경계 위에 놓이는 임의의 직사각형 $B \in \mathcal{B}$는

$p \in B \subset D_p \subset G$ or $p \in B \subset G$

를 만족한다. 즉 그림에서 보이는 바와 같다. 다시 말해, $\mathcal{B}$는 위의 (ii)를 만족한다.

![](attachments/250831_%EC%B1%85%EA%B3%B5%EB%B6%80_%5BSOGT%5D%20General%20Topology_01.png)
Example 1.3: 임의의 이산공간 $(X, \mathcal{D})$를 생각하자. 그러면 $X$의 모든 한원소 부분집합들의 족 $\mathcal{B} = \{\{p\} : p \in X\}$는 $X$ 위의 이산 위상 $\mathcal{D}$에 대한 기저이다. 각 한원소집합 $\{p\}$는 $\mathcal{D}$-열린이고, 모든 $A \subset X$가 $\mathcal{D}$-열린이므로, 모든 집합은 한원소집합들의 합집합이다. 사실 $X$의 부분집합들의 임의의 다른 족 $\mathcal{B}^*$가 $\mathcal{D}$에 대한 기저인 것은 $\mathcal{B}$의 상위족(superclass)인 것, 즉 $\mathcal{B}^* \supset \mathcal{B}$인 것과 동치이다.

> 규빈: 결국에는 ${\cal B}$가 최소기저의 느낌을 가진다는것. 

이제 다음 질문을 던진다: 집합 $X$의 부분집합들의 족 $\mathcal{B}$가 주어졌을 때, $\mathcal{B}$가 $X$ 위의 어떤 위상에 대한 기저가 되려면 언제인가? $X$는 $X$ 위의 모든 위상에서 열린집합이므로 $X = \bigcup\{B : B \in \mathcal{B}\}$가 필요조건임은 분명하다. 다음 예시는 다른 조건들도 필요함을 보여준다.